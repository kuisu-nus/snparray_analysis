from dataclasses import dataclass
import torch
from torch.utils.data import Dataset, DataLoader
from collections import OrderedDict, defaultdict
from typing import List, Dict, Any, Tuple
import xarray as xr
import time
import os
import sgkit as sg 
# from sgkit.io.vcf import vcf_to_zarr
import numpy as np
import logging
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ChromosomeWindow:
    """Data class to store chromosome window information."""
    human_id: str
    chromosome: str
    snp_start: int
    snp_end: int
    window_size: str
    variant_count: int = 0

class CachedKGSGDataset(Dataset):
    """
    PyTorch Dataset for 1KG data using sgkit with efficient caching.
    Supports both VCF and Zarr formats with optimized data loading.
    """
    
    def __init__(self,
                 data_path: str,
                 windows: List[Any],
                 missing_value: int = -1,
                 cache_size: int = 3,
                 is_zarr: bool = True):
        """
        Initialize cached dataset for genotype data.
        
        Args:
            data_path: Path to VCF or Zarr file
            windows: List of genomic windows to process
            missing_value: Value to represent missing genotypes
            cache_size: Maximum number of samples to cache in memory
            is_zarr: Whether the input is Zarr format (False for VCF)
        """
        self.data_path = data_path
        self.windows = windows
        self.missing_value = missing_value
        self.cache_size = cache_size
        self.is_zarr = is_zarr
        
        # Load dataset using sgkit
        self.ds = self._load_dataset(data_path, is_zarr)
        logger.info(f"Dataset loaded with {self.ds.samples.size} samples and {self.ds.variants.size} variants")
        
        # Cache for loaded samples
        self.sample_cache = OrderedDict()
        
        # Precompute window indices for efficient access
        self._organize_windows_by_sample()

    def convert_gt(self, ds_data)->xr.DataArray:
        genotype_data = ds_data['call_genotype'].values
        alt_count = np.sum(genotype_data, axis=2)
        missing_mask = np.any(genotype_data < 0, axis=2)
        alt_count[missing_mask] = -1
    
        # 创建新的 DataArray
        converted_data = xr.DataArray(
            alt_count,
            dims=['variants', 'samples'],
            coords={
                'variants': ds_data.variants,
                'samples': ds_data.samples
            },
            attrs={'description': 'Genotype encoded as -1(missing), 0(HomRef), 1(Het), 2(HomAlt)'}
        )
        return converted_data
    
    def _load_dataset(self, data_path: str, is_zarr: bool = True) -> xr.Dataset:
        """Load genotype dataset from VCF or Zarr format using sgkit."""
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"Data path not found: {data_path}")
        
        if is_zarr:
            # Load from Zarr format (optimized for sgkit)
            ds = sg.load_dataset(data_path)
        else:
            # Convert VCF to Zarr and load (one-time conversion)
            zarr_path = data_path.replace('.vcf', '.zarr').replace('.gz', '')
            # if not os.path.exists(zarr_path):
            #     logger.info(f"Converting VCF to Zarr format: {data_path} -> {zarr_path}")
            #     vcf_to_zarr(data_path, zarr_path, chunk_length=10000, chunk_width=100)
            # ds = sg.load_dataset(zarr_path)
            raise FileExistsError(f"not support vcf format:{data_path}")
        start_time = time.time()
        # ds = self.convert_gt(ds)
        logger.info(f"cost time for convert gt: {time.time()-start_time}")
        return ds
    
    def _organize_windows_by_sample(self) -> None:
        """Organize windows by sample for efficient data retrieval."""
        self.window_to_index = {}
        self.sample_windows = defaultdict(list)
        
        for idx, window in enumerate(self.windows):
            self.window_to_index[idx] = window
            self.sample_windows[window.human_id].append((idx, window))
    
    def _load_sample_data(self, human_id: str) -> Dict[str, Any]:
        """
        Load all genotype data for a specific sample with caching.
        
        Args:
            human_id: Sample identifier
            
        Returns:
            Dictionary with sample's genotype data organized by chromosome
        """
        # Check cache first (with LRU update)
        if human_id in self.sample_cache:
            data = self.sample_cache.pop(human_id)
            self.sample_cache[human_id] = data
            return data
        
        logger.info(f"Loading data for sample: {human_id}")
        load_start = time.time()
        
        # Find sample index
        sample_ids = self.ds['sample_id'].values
        sample_idx = np.where(sample_ids == human_id)[0]
        if len(sample_idx) == 0:
            raise ValueError(f"Sample {human_id} not found")
        sample_idx = sample_idx[0]
        genotype_data = self.ds["call_genotype"].values[:, sample_idx, :]

        alt_count = np.sum(genotype_data, axis=-1)
        missing_mask = np.any(genotype_data < 0, axis=-1)
        alt_count[missing_mask] = -1
        
        # Cache the data with LRU management
        self._add_to_cache(human_id, alt_count)
        return alt_count
    
    def _add_to_cache(self, sample_id: str, data: Dict[str, Any]) -> None:
        """Add sample data to cache with LRU eviction policy."""
        if len(self.sample_cache) >= self.cache_size:
            oldest_key = next(iter(self.sample_cache))
            del self.sample_cache[oldest_key]
            logger.debug(f"Evicted {oldest_key} from cache")
        
        self.sample_cache[sample_id] = data
        logger.debug(f"Cached {sample_id}, cache size: {len(self.sample_cache)}")
    
    def _encode_genotype(self, gt_call: Any) -> int:
        """Encode genotype call to alternate allele count."""
        if gt_call is None or np.any(gt_call < 0):
            return self.missing_value
        return np.sum(gt_call)
    
    def __len__(self) -> int:
        """Return number of windows in the dataset."""
        return len(self.windows)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """
        Get genotype data for a specific genomic window.
        
        Args:
            idx: Window index
            
        Returns:
            Tuple of (genotype_tensor, window_info_dict)
        """
        window = self.window_to_index[idx]
        
        # Load sample data (uses cache if available)
        sample_data = self._load_sample_data(window.human_id)
        
        # Extract variants for the specific genomic window
        # chrom_data = sample_data.get(window.chromosome, {})
        # window_variants = []
        
        # Collect genotypes for positions within window boundaries
        # for pos in sorted(chrom_data.keys()):
        #     if window.snp_start <= pos <= window.snp_end:
        #         window_variants.append(chrom_data[pos])
        window_variants = sample_data[window.snp_start:window.snp_end]
        # Convert to PyTorch tensor
        if window_variants:
            genotype_tensor = torch.tensor(window_variants, dtype=torch.int32)
        else:
            # Empty window - return empty tensor
            genotype_tensor = torch.tensor([], dtype=torch.int32)
        
        # Create comprehensive window metadata
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
        """Clear sample cache to free memory."""
        self.sample_cache.clear()
        logger.info("Sample cache cleared")
    
    def preload_samples(self, sample_ids: List[str]) -> None:
        """
        Preload specific samples into cache for faster access.
        
        Args:
            sample_ids: List of sample IDs to preload
        """
        for sample_id in sample_ids:
            if sample_id not in self.sample_cache:
                self._load_sample_data(sample_id)
        logger.info(f"Preloaded {len(sample_ids)} samples into cache")
    
    def get_cache_info(self) -> Dict[str, Any]:
        """Get current cache statistics and information."""
        return {
            'cache_size': len(self.sample_cache),
            'cached_samples': list(self.sample_cache.keys()),
            'max_cache_size': self.cache_size
        }
    
    def get_dataset_stats(self) -> Dict[str, Any]:
        """Get comprehensive statistics about the dataset."""
        variant_counts = []
        sample_ids = set()
        chromosomes = set()
        window_sizes = {}
        
        for window in self.windows:
            variant_counts.append(window.variant_count)
            sample_ids.add(window.human_id)
            chromosomes.add(window.chromosome)
            window_sizes[window.window_size] = window_sizes.get(window.window_size, 0) + 1
        
        variant_counts = np.array(variant_counts)
        
        return {
            'total_windows': len(self.windows),
            'unique_samples': len(sample_ids),
            'unique_chromosomes': len(chromosomes),
            'windows_by_size': window_sizes,
            'mean_variants_per_window': np.mean(variant_counts),
            'std_variants_per_window': np.std(variant_counts),
            'min_variants_per_window': np.min(variant_counts),
            'max_variants_per_window': np.max(variant_counts),
            'cached_samples': len(self.sample_cache),
            'total_samples': self.ds.samples.size,
            'total_variants': self.ds.variants.size
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
    dataset = CachedKGSGDataset(mt_path, windows, **kwargs)
    
    # Create DataLoader
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
    )
    
    return dataloader

# Example usage
if __name__ == "__main__":
    # Step 1: Generate windows
    mt_path="/home/sukui/01.data/03.raw_data/1kg/ALL.chr22.vcz"
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