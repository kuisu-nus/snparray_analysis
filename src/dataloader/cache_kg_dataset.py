import hail as hl
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from dataclasses import dataclass
from typing import List, Dict, Any, Optional, Tuple
import time
import os
import logging
from collections import OrderedDict, defaultdict
import numpy as np


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
    human_id: str
    chromosome: str
    snp_start: int
    snp_end: int
    window_size: str
    variant_count: int = 0


class KGGenomeSplitter:
    """
    Split 1KG data by chromosome using configurable window sizes.
    Generates windows for specified sizes and saves as CSV.
    """
    
    def __init__(self, mt_path: str):
        """
        Initialize genome splitter.
        
        Args:
            mt_path: Path to Hail MatrixTable
        """
        hl.stop()
        hl.init()
        self.mt_path = mt_path
        self.mt = hl.read_matrix_table(mt_path)
    
    def generate_windows_for_size(self,
                                chromosome: str,
                                window_size: int,
                                min_variants: int = 10) -> List[ChromosomeWindow]:
        """
        Generate windows for a specific chromosome and window size.
        
        Args:
            chromosome: Chromosome name
            window_size: Window size in base pairs
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
        samples = chr_mt.aggregate_cols(hl.agg.collect(chr_mt.s))
        
        # Generate non-overlapping windows
        start_pos = positions[0]
        end_pos = positions[-1]
        
        current_start = start_pos
        while current_start <= end_pos:
            current_end = current_start + window_size - 1
            
            # Count variants in current window
            window_variants = [p for p in positions if current_start <= p <= current_end]
            variant_count = len(window_variants)
            
            if variant_count >= min_variants:
                for sample in samples:
                    window = ChromosomeWindow(
                        human_id=sample,
                        chromosome=chromosome,
                        snp_start=current_start,
                        snp_end=current_end,
                        window_size=self._format_window_size(window_size),
                        variant_count=variant_count
                    )
                    windows.append(window)
            
            # Move to next window
            current_start = current_end + 1
        
        print(f"Generated {len(windows)} windows for {chromosome} with size {window_size} "
              f"({len(windows)//len(samples) if samples else 0} windows per sample)")
        return windows
    
    def _format_window_size(self, size: int) -> str:
        """Format window size for human-readable output."""
        if size >= 1_000_000:
            return f"{size//1_000_000}M"
        elif size >= 1_000:
            return f"{size//1_000}K"
        else:
            return f"{size}"
    
    def generate_all_windows(self,
                           window_sizes: List[int] = [1000, 100000, 1000000],
                           chromosomes: Optional[List[str]] = None,
                           min_variants: int = 10) -> List[ChromosomeWindow]:
        """
        Generate windows for all specified sizes and chromosomes.
        
        Args:
            window_sizes: List of window sizes in base pairs
            chromosomes: List of chromosomes to process
            min_variants: Minimum variants per window
            
        Returns:
            List of all ChromosomeWindow objects
        """
        if chromosomes is None:
            chromosomes = [f'chr{i}' for i in range(1, 23)]
        
        all_windows = []
        
        for window_size in window_sizes:
            print(f"Generating windows with size {window_size}...")
            for chrom in chromosomes:
                chrom_windows = self.generate_windows_for_size(
                    chrom, window_size, min_variants
                )
                all_windows.extend(chrom_windows)
        
        print(f"Total windows generated: {len(all_windows)}")
        return all_windows
    
    def save_windows_to_csv(self,
                          windows: List[ChromosomeWindow],
                          output_path: str) -> None:
        """
        Save window information to CSV file.
        
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
                'window_size': window.window_size,
                'variant_count': window.variant_count
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
        print(f"Saved {len(windows)} windows to {output_path}")


class CachedKGDataset(Dataset):
    """
    PyTorch Dataset for 1KG data with sample-level caching.
    Efficiently loads and caches sample data to minimize Hail queries.
    """
    
    def __init__(self,
                 mt_path: str,
                 windows: List[ChromosomeWindow],
                 missing_value: int = -1,
                 cache_size:int=3):
        """
        Initialize cached dataset.
        
        Args:
            mt_path: Path to Hail MatrixTable
            windows: List of ChromosomeWindow objects
            missing_value: Value to represent missing genotypes
        """
        self.mt_path = mt_path
        self.windows = windows
        self.missing_value = missing_value
        
        # Initialize Hail\
        hl.stop()
        hl.init()
        self.mt = hl.read_matrix_table(mt_path)
        
        # Cache for loaded samples
        self.sample_cache = OrderedDict()
        self.cache_size = cache_size
        
        # Precompute window indices per sample for efficiency
        self._organize_windows_by_sample()
    
    def _organize_windows_by_sample(self) -> None:
        """Organize windows by sample for efficient data loading."""
        self.sample_windows = defaultdict(list)
        self.window_to_index = {}
        
        for idx, window in enumerate(self.windows):
            self.sample_windows[window.human_id].append((idx, window))
            self.window_to_index[idx] = window
    
    def _load_sample_data(self, human_id: str) -> Dict[str, Any]:
        """
        Load all genotype data for a specific sample into memory.
        
        Args:
            human_id: Sample identifier
            
        Returns:
            Dictionary with sample's genotype data organized by chromosome
        """
        if human_id in self.sample_cache:
            data = self.sample_cache.pop(human_id)
            self.sample_cache[human_id] = data
            return data
        
        print(f"Loading data for sample: {human_id}")
        load_start = time.time()
        
        # Filter to specific sample
        sample_mt = self.mt.filter_cols(self.mt.s == human_id) # NA20752, HG00096
        # human_ids = ['NA20752', 'HG00096']
        # sample_mt = sample_mt.filter_cols(hl.literal(human_ids).contains(sample_mt.s))

        # entries_table = sample_mt.entries()
        # all_data = entries_table.select(
        #                 'locus', 'alleles', 's', 'GT'
        #             ).collect()

        # Collect all variant positions and genotypes
        col_stime = time.time()
        entries_table = sample_mt.entries().key_by()
        sample_data = entries_table.select('locus', 'GT').order_by('locus').collect() #要减少collect()的调用，因为每次都是相当于把数据全部扫描一次。
        # TODO: 每次读取5个人的数据，
        col_etime = time.time()
        print(f"collect cost time: {col_etime - col_stime}")
        # Organize by chromosome for efficient window queries
        organized_data = defaultdict(dict)
        for row in sample_data:
            chrom = row.locus.contig
            pos = row.locus.position
            organized_data[chrom][pos] = self._encode_genotype(row.GT)
        
        load_time = time.time() - load_start
        print(f"Loaded {len(sample_data)} variants for {human_id} in {load_time:.2f}s")
        
        # Cache the data
        self.sample_cache[human_id] = organized_data
        self._add_to_cache(sample_id=human_id, data=organized_data)
        return organized_data
    
    def _add_to_cache(self, sample_id, data):
        # 如果缓存已满，移除最老的项
        if len(self.sample_cache) >= self.cache_size:
            oldest_key = next(iter(self.sample_cache))
            del self.sample_cache[oldest_key]
            print(f"从缓存中移除: {oldest_key}")
        
        # 添加新数据
        self.sample_cache[sample_id] = data
        print(f"添加到缓存: {sample_id}, 缓存大小: {len(self.sample_cache)}")

    def get_cache_info(self):
        """获取缓存信息"""
        return {
            'cache_size': len(self.sample_cache),
            'cached_samples': list(self.sample_cache.keys()),
            'max_cache_size': self.cache_size
        }
    
    def _encode_genotype(self, gt_call) -> int:
        """Encode single genotype call."""
        if gt_call is None:
            return self.missing_value
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
            Tuple of (genotypes, window_info)
        """
        window = self.window_to_index[idx]
        
        # Load sample data if not cached
        sample_data = self._load_sample_data(window.human_id)
        
        # Extract variants for the specific window
        chrom_data = sample_data.get(window.chromosome, {})
        window_variants = []
        
        # Collect genotypes for positions in the window
        for pos in sorted(chrom_data.keys()):
            if window.snp_start <= pos <= window.snp_end:
                window_variants.append(chrom_data[pos])
        
        # Convert to tensor
        if window_variants:
            genotype_tensor = torch.tensor(window_variants, dtype=torch.int32)
        else:
            # Empty window - return tensor with zeros
            genotype_tensor = torch.tensor([], dtype=torch.int32)
        
        # Create window info dictionary
        window_info = {
            'human_id': window.human_id,
            'chromosome': window.chromosome,
            'snp_start': window.snp_start,
            'snp_end': window.snp_end,
            'window_size': window.window_size,
            'variant_count': len(window_variants),
            'window_index': idx
        }
        
        return genotype_tensor, window_info
    
    def clear_cache(self) -> None:
        """Clear the sample cache to free memory."""
        self.sample_cache.clear()
    
    def preload_samples(self, sample_ids: List[str]) -> None:
        """
        Preload specific samples into cache.
        
        Args:
            sample_ids: List of sample IDs to preload
        """
        for sample_id in sample_ids:
            if sample_id not in self.sample_cache:
                self._load_sample_data(sample_id)
    
    def get_dataset_stats(self) -> Dict[str, Any]:
        """Get comprehensive statistics about the dataset."""
        variant_counts = [w.variant_count for w in self.windows]
        sample_ids = list(set(w.human_id for w in self.windows))
        chromosomes = list(set(w.chromosome for w in self.windows))
        
        # Count windows by size
        window_sizes = {}
        for window in self.windows:
            window_sizes[window.window_size] = window_sizes.get(window.window_size, 0) + 1
        
        return {
            'total_windows': len(self.windows),
            'unique_samples': len(sample_ids),
            'unique_chromosomes': len(chromosomes),
            'windows_by_size': window_sizes,
            'mean_variants_per_window': np.mean(variant_counts),
            'std_variants_per_window': np.std(variant_counts),
            'min_variants_per_window': np.min(variant_counts),
            'max_variants_per_window': np.max(variant_counts),
            'cached_samples': len(self.sample_cache)
        }


def create_kg_dataloader(mt_path: str,
                        windows_csv: str,
                        batch_size: int = 32,
                        shuffle: bool = True,
                        num_workers: int = 0,
                        **kwargs) -> DataLoader:
    """
    Convenience function to create a DataLoader from windows CSV.
    
    Args:
        mt_path: Path to Hail MatrixTable
        windows_csv: Path to windows CSV file
        batch_size: Batch size for DataLoader
        shuffle: Whether to shuffle the data
        num_workers: Number of worker processes
        **kwargs: Additional arguments for CachedKGDataset
        
    Returns:
        PyTorch DataLoader
    """
    # Load windows from CSV
    df = pd.read_csv(windows_csv)
    df = df.sort_values(['human_id', 'chromosome'])
    print(f"df data:\{df.head()}")
    windows = []
    
    for _, row in df.iterrows():
        window = ChromosomeWindow(
            human_id=row['human_id'],
            chromosome=row['chromosome'],
            snp_start=row['snp_start'],
            snp_end=row['snp_end'],
            window_size=row['window_size'],
            variant_count=row.get('variant_count', 0)
        )
        windows.append(window)
    
    # Create dataset
    dataset = CachedKGDataset(mt_path, windows, **kwargs)
    
    # Create DataLoader
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        collate_fn=collate_variable_length  # Custom collate for variable length sequences
    )
    
    return dataloader


def collate_variable_length(batch):
    """
    Custom collate function for variable length genotype sequences.
    
    Args:
        batch: List of (genotypes, window_info) tuples
        
    Returns:
        Batched tensors and metadata
    """
    genotypes, window_info = zip(*batch)
    
    # For variable length sequences, we can:
    # 1. Pad to max length in batch
    max_length = max(genotype.shape[0] for genotype in genotypes)
    padded_genotypes = []
    
    for genotype in genotypes:
        if genotype.shape[0] < max_length:
            # Pad with missing values
            pad_length = max_length - genotype.shape[0]
            padded = torch.cat([
                genotype, 
                torch.full((pad_length,), -1, dtype=genotype.dtype)
            ])
            padded_genotypes.append(padded)
        else:
            padded_genotypes.append(genotype)
    
    batched_genotypes = torch.stack(padded_genotypes)
    
    # Batch window info
    batched_info = []
    for info in window_info:
        batched_info.append(info)
    
    return batched_genotypes, batched_info


# Example usage
if __name__ == "__main__":
    # Step 1: Generate windows
    mt_path="/home/sukui/01.data/03.raw_data/1kg/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.mt"
    csv_path = "/home/sukui/03.projects/01.pgt/snparray_analysis/work_dir/hail_data/1kg_windows.csv"
    # splitter = KGGenomeSplitter(mt_path)
    # windows = splitter.generate_all_windows(
    #     window_sizes=[1_000_000, 10_000_000],
    #     chromosomes=['22']  # Process specific chromosomes
    # )
    # splitter.save_windows_to_csv(windows, csv_path)
    
    # Step 2: Create DataLoader
    dataloader = create_kg_dataloader(
        mt_path=mt_path,
        windows_csv=csv_path,
        batch_size=16,
        shuffle=False
    )
    
    # Step 3: Use in training loop
    for batch_idx, (genotypes, window_info) in enumerate(dataloader):
        print(f"Batch {batch_idx}: {genotypes.shape}")
        if batch_idx >= 5:  # Just show first 5 batches
            break