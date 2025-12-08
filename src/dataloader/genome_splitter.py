import hail as hl
import torch
from torch.utils.data import Dataset, DataLoader
from typing import Optional, Tuple, List, Dict, Any
import numpy as np
import pandas as pd
from dataclasses import dataclass
import os
import logging 
import time

from vcf_dataset import KGDataLoader

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

proxy_vars = [
    'http_proxy', 'https_proxy', 'ftp_proxy',
    'HTTP_PROXY', 'HTTPS_PROXY', 'FTP_PROXY',
    'no_proxy', 'NO_PROXY'
]

for var in proxy_vars:
    os.environ.pop(var, None)  # 使用pop避免KeyError

@dataclass
class ChromosomeWindow:
    """Data class to store chromosome window information."""
    vcz_path: str
    human_id: str
    chromosome: str
    snp_start: int
    snp_end: int
    variant_count: int = 0


class KGGenomeSplitter:
    """
    Split 1KG data by chromosome using sliding windows.
    Divides each chromosome into 1Mb windows with 1% sliding overlap.
    """
    
    def __init__(self, 
                 mt_path: str,
                 window_size: int = 1000000,  # 1Mb
                 slide_step: float = 0.01):   # 1% sliding
        """
        Initialize genome splitter.
        
        Args:
            mt_path: Path to Hail MatrixTable
            window_size: Window size in base pairs (default: 1Mb)
            slide_step: Sliding step as fraction of window size (default: 1%)
        """
        hl.stop()
        hl.init()
            
        self.mt_path = mt_path
        self.window_size = window_size
        self.slide_step = slide_step
        self.slide_bp = int(window_size * slide_step)  # Actual bp to slide
        
        # Load MatrixTable
        self.mt = hl.read_matrix_table(mt_path)
        
    def get_chromosome_windows(self, 
                             chromosome: str,
                             min_variants: int = 10) -> List[ChromosomeWindow]:
        """
        Generate sliding windows for a specific chromosome.
        
        Args:
            chromosome: Chromosome name (e.g., 'chr1'), GRCh37:1 ; GRCh38:chr1
            min_variants: Minimum number of variants required in a window
            
        Returns:
            List of ChromosomeWindow objects
        """
        # Filter to specific chromosome
        chr_mt = self.mt.filter_rows(self.mt.locus.contig == chromosome)
        
        if chr_mt.count_rows() == 0:
            print(f"No variants found for {chromosome}")
            return []
        
        # Get variant positions
        variants = chr_mt.rows()
        positions = variants.aggregate(hl.agg.collect(variants.locus.position))
        positions = sorted(positions)

        if not positions:
            return []
        
        windows = []
        
        # Generate sliding windows
        start_idx = 0
        min_vars = min(self.window_size*0.01, 100)
        while positions:
            window_variants = positions[start_idx:start_idx+self.window_size]
            # Count variants in current window
            variant_count = len(window_variants)

            if variant_count < min_vars:
                break

            start_pos, end_pos = min(window_variants), max(window_variants)
            samples = chr_mt.aggregate_cols(hl.agg.collect(chr_mt.s))
            for sample in samples:
                window = ChromosomeWindow(
                    vcz_path=self.mt_path
                    human_id=sample,
                    chromosome=chromosome,
                    snp_start=start_pos,
                    snp_end=end_pos,
                    variant_count=variant_count
                )
                windows.append(window)
            
            # Slide window
            start_idx += self.window_size
        
        print(f"Generated {len(windows)} windows for {chromosome} "
              f"({len(windows)//len(samples) if samples else 0} windows per sample)")
        return windows
    
    def get_all_chromosome_windows(self, 
                                 chromosomes: Optional[List[str]] = None,
                                 min_variants: int = 10) -> List[ChromosomeWindow]:
        """
        Generate windows for all chromosomes.
        
        Args:
            chromosomes: List of chromosomes to process (default: autosomes 1-22)
            min_variants: Minimum variants per window
            
        Returns:
            List of all ChromosomeWindow objects
        """
        if chromosomes is None:
            chromosomes = [f'chr{i}' for i in range(1, 23)]
        
        all_windows = []
        for chrom in chromosomes:
            chrom_windows = self.get_chromosome_windows(chrom, min_variants)
            all_windows.extend(chrom_windows)
            
        print(f"Total windows generated: {len(all_windows)}")
        return all_windows
    
    def save_windows_to_tsv(self, 
                          windows: List[ChromosomeWindow],
                          output_path: str) -> None:
        """
        Save window information to TSV file.
        
        Args:
            windows: List of ChromosomeWindow objects
            output_path: Output file path
        """
        data = []
        for window in windows:
            data.append({
                'human_id': window.human_id,
                'chromosome': window.chromosome,
                'snp_start': window.snp_start,
                'snp_end': window.snp_end,
                'variant_count': window.variant_count
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output_path, sep='\t', index=False)
        print(f"Saved {len(windows)} windows to {output_path}")
    
    def load_windows_from_tsv(self, tsv_path: str) -> List[ChromosomeWindow]:
        """
        Load windows from TSV file.
        
        Args:
            tsv_path: Path to TSV file
            
        Returns:
            List of ChromosomeWindow objects
        """
        df = pd.read_csv(tsv_path, sep='\t')
        windows = []
        
        for _, row in df.iterrows():
            window = ChromosomeWindow(
                human_id=row['human_id'],
                chromosome=row['chromosome'],
                snp_start=row['snp_start'],
                snp_end=row['snp_end'],
                variant_count=row['variant_count']
            )
            windows.append(window)
        
        print(f"Loaded {len(windows)} windows from {tsv_path}")
        return windows


class KGWindowDataLoader(KGDataLoader):
    """
    PyTorch DataLoader for 1KG data with chromosome window support.
    Extends KGDataLoader to handle window-based data loading.
    """
    
    def __init__(self,
                 mt_path: str,
                 windows: List[ChromosomeWindow],
                 window_index: Optional[List[int]] = None,
                 **kwargs):
        """
        Initialize window-based DataLoader.
        
        Args:
            mt_path: Path to Hail MatrixTable
            windows: List of ChromosomeWindow objects
            window_index: Optional list of window indices to use
            **kwargs: Additional KGDataLoader arguments
        """
        # Initialize base class
        super().__init__(mt_path, **kwargs)
        
        self.windows = windows
        self.window_index = window_index if window_index is not None else list(range(len(windows)))

    @staticmethod
    def encode_genotype(gt_call, missing_value=-1):
        """编码单个基因型"""
        if gt_call is None:
            return missing_value
        return gt_call.n_alt_alleles()


    def __len__(self) -> int:
        """Return number of windows in the dataset."""
        return len(self.windows)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """
        Get genotype data for a specific window.
        
        Args:
            idx: Window index
            
        Returns:
            Tuple of (genotypes, window_info) where:
            - genotypes: Tensor of shape [n_variants_in_window]
            - window_info: Dictionary with window metadata
        """
        # window_idx = self.window_index[idx]
        window = self.windows[idx]
        # variant_indices = self.window_variants[window_idx]
        
        # Get sample and filter to window
        sample_q_stime = time.time()
        sample_mt = self.mt.filter_cols(self.mt.s == window.human_id)
        sample_q_etiem = time.time()
        print(f'sample query time: {sample_q_etiem - sample_q_stime}')
        window_mt = sample_mt.filter_rows(
            (sample_mt.locus.contig == str(window.chromosome)) &
            (sample_mt.locus.position >= window.snp_start) &
            (sample_mt.locus.position <= window.snp_end)
        )
        sample_q_etime2 = time.time()
        print(f"snp quer time: {sample_q_etime2-sample_q_etiem}")
        # Convert genotypes to tensor
        genotypes = window_mt.GT.collect()
        print(f"collect time: {time.time() - sample_q_etime2}")
        genotype_tensor = torch.tensor([self.encode_genotype(gt, self.missing_value)
                                        for gt in genotypes], dtype=torch.int32)
        
        # Create window info dictionary
        window_info = {
            'human_id': window.human_id,
            'chromosome': window.chromosome,
            'snp_start': window.snp_start,
            'snp_end': window.snp_end,
            'variant_count': window.variant_count,
            'window_index': idx
        }
        
        return genotype_tensor, window_info
    
    def get_window_stats(self) -> Dict[str, Any]:
        """Get statistics about the windows."""
        variant_counts = [w.variant_count for w in self.windows]
        
        return {
            'total_windows': len(self.windows),
            'mean_variants_per_window': np.mean(variant_counts),
            'std_variants_per_window': np.std(variant_counts),
            'min_variants_per_window': np.min(variant_counts),
            'max_variants_per_window': np.max(variant_counts),
            'unique_samples': len(set(w.human_id for w in self.windows)),
            'unique_chromosomes': len(set(w.chromosome for w in self.windows))
        }


# Test cases for new functionality
def test_genome_splitter():
    """Test genome splitter functionality."""
    print("Testing KGGenomeSplitter...")
    from tqdm import tqdm
    
    # Create mock data
    # test_path = "/home/sukui/03.projects/01.pgt/snparray_analysis/work_dir/hail_data/1kg.mt"
    test_path = "/home/sukui/01.data/03.raw_data/1kg/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.mt"
    try:
        # Test splitter initialization
        splitter = KGGenomeSplitter(test_path, window_size=1_000_000, slide_step=0.1)
        print("✓ Splitter initialized successfully")
        
        # Test chromosome window generation
        windows = splitter.get_chromosome_windows('22', min_variants=5)
        assert len(windows) > 0, "Should generate windows"
        print(f"✓ Generated {len(windows)} windows for chr1")
        
        # Test window properties
        window = windows[0]
        assert hasattr(window, 'human_id'), "Window should have human_id"
        assert hasattr(window, 'chromosome'), "Window should have chromosome"
        assert hasattr(window, 'snp_start'), "Window should have snp_start"
        assert hasattr(window, 'snp_end'), "Window should have snp_end"
        assert window.snp_end > window.snp_start, "Invalid window coordinates"
        print("✓ Window properties correct")
        
        # Test TSV save/load
        tsv_path = "/home/sukui/03.projects/01.pgt/snparray_analysis/work_dir/hail_data/test_windows.tsv"
        splitter.save_windows_to_tsv(windows, tsv_path)
        loaded_windows = splitter.load_windows_from_tsv(tsv_path)
        assert len(loaded_windows) == len(windows), "Should load same number of windows"
        print("✓ TSV save/load works correctly")
        
        # Test window data loader
        window_loader = KGWindowDataLoader(test_path, loaded_windows[:10])  # Use first 10 windows
        assert len(window_loader) == 10, "Loader should have 10 windows"
        
        # Test data loading
        genotypes, window_info = window_loader[0]
        for i, data in tqdm(enumerate(window_loader)):
            print(f"{i} gt: {data[0].shape}")
        assert isinstance(genotypes, torch.Tensor), "Should return tensor"
        assert isinstance(window_info, dict), "Should return window info dict"
        assert 'human_id' in window_info, "Window info should contain human_id"
        print("✓ Window data loading works correctly")
        
        # Test statistics
        stats = window_loader.get_window_stats()
        assert 'total_windows' in stats, "Stats should contain total windows"
        print("✓ Window statistics work correctly")
        
        print("All genome splitter tests passed! ✓")
        
    except Exception as e:
        print(f"Test failed: {e}")
        raise
    finally:
        # Cleanup
        pass
            # os.remove(tsv_path)


def create_mock_matrix_table_with_chromosomes():
    """Create mock MatrixTable with multiple chromosomes for testing."""
    # Create samples
    samples = [f"sample_{i}" for i in range(5)]
    
    # Create variants across multiple chromosomes
    entries = []
    for chrom in ['chr1', 'chr2']:
        for pos in range(1, 1001, 10):  # Variants every 10bp
            for sample_idx, sample in enumerate(samples):
                # Create variant with random genotype
                gt = hl.Call([sample_idx % 2, (sample_idx + 1) % 2])  # Simple pattern
                entries.append(hl.Struct(
                    s=sample,
                    locus=hl.Locus(contig=chrom, position=pos),
                    alleles=['A', 'G'],
                    GT=gt
                ))
    
    # Convert to MatrixTable (simplified approach)
    return hl.Table.from_pandas(
        pd.DataFrame(entries)
    )._to_matrix(
        row_key=['locus', 'alleles'],
        col_key=['s'],
        row_fields=[],
        col_fields=[]
    )


# Example usage
def example_usage():
    """Example of how to use the genome splitting functionality."""
    # Initialize splitter
    splitter = KGGenomeSplitter("path/to/1kg.mt")
    
    # Generate windows for all autosomes
    windows = splitter.get_all_chromosome_windows()
    
    # Save windows to file
    splitter.save_windows_to_tsv(windows, "1kg_windows.tsv")
    
    # Create data loader for training
    train_loader = KGWindowDataLoader(
        mt_path="path/to/1kg.mt",
        windows=windows,
        split='train'
    )
    
    # Create batched loader
    batch_loader = BatchKGDataLoader(train_loader, batch_size=32, shuffle=True)
    
    return batch_loader


if __name__ == "__main__":
    # Run tests
    test_genome_splitter()