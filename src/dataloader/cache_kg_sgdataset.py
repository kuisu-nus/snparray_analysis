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

def gt_data_collate_fn(batch):
    """
    自定义 collate 函数用于处理 GTData 对象
    """
    if isinstance(batch[0], tuple):
        # 如果 batch 是 (GTData, dict) 的元组
        gt_data_list = [item[0] for item in batch]
        info_list = [item[1] for item in batch]
        
        # 合并 GTData
        batched_gt_data = GTData.collate(gt_data_list)
        
        # 合并 info dicts (保持为列表)
        return batched_gt_data, info_list
    elif isinstance(batch[0], GTData):
        # 如果 batch 直接是 GTData 对象
        return GTData.collate(batch)
    else:
        # 回退到默认的 collate
        return torch.utils.data.default_collate(batch)

# 在 GTData 类中添加 collate 方法
@dataclass
class GTData:
    gt: torch.Tensor
    chr: torch.Tensor
    pos: torch.Tensor
    pos_encode: torch.Tensor = None
    window_info:dict = None

    @classmethod
    def from_numpy(cls, gt: np.ndarray, chr: np.ndarray, pos: np.ndarray):
        return cls(
            gt=torch.from_numpy(gt).to(torch.int32),
            chr=torch.from_numpy(chr).to(torch.int32),
            pos=torch.from_numpy(pos).to(torch.int32)
        )
    
    def encode_positions(self, resolution: int = 20_000):
        """编码位置信息"""
        lowres_pos = self.pos // resolution
        highres_pos = self.pos % resolution
        pos_encode = torch.stack([lowres_pos, highres_pos])
        
        return GTData(
            gt=self.gt,
            chr=self.chr,
            pos=self.pos,
            pos_encode=pos_encode
        )
    
    def get_window_by_position_exact(self, start_pos: int, end_pos: int, max_len:int=1024):
        # 找到在位置范围内的索引
        mask = (self.pos >= start_pos) & (self.pos <= end_pos)
        indices = torch.where(mask)[0]
        
        if len(indices) == 0:
            logger.warning(f"not find valid vars")
            return GTData(
                gt=torch.tensor([], dtype=torch.int32),
                chr=torch.tensor([], dtype=torch.int32),
                pos=torch.tensor([], dtype=torch.int32)
            )
        if len(indices) > max_len:
            indices = torch.randperm(len(indices))[:max_len].sort()[0]
        return GTData(
            gt=self.gt[indices],
            chr=self.chr[indices],
            pos=self.pos[indices],
            pos_encode=self.pos_encode[:, indices] if self.pos_encode is not None else None
        )
    
    @classmethod
    def collate(cls, batch: List['GTData']) -> 'GTData':
        """
        将多个 GTData 对象合并成一个
        
        Args:
            batch: GTData 对象列表
            
        Returns:
            合并后的 GTData 对象
        """
        if not batch:
            return cls(
                gt=torch.tensor([], dtype=torch.int32),
                chr=torch.tensor([], dtype=torch.int32),
                pos=torch.tensor([], dtype=torch.int32)
            )
        
        # 合并所有张量
        gt_tensors = [item.gt for item in batch]
        chr_tensors = [item.chr for item in batch]
        pos_tensors = [item.pos for item in batch]
        
        batched_gt = torch.stack(gt_tensors, dim=0)
        batched_chr = torch.stack(chr_tensors, dim=0)
        batched_pos = torch.stack(pos_tensors, dim=0)
        
        # 处理 pos_encode（如果有的话）
        batched_pos_encode = None
        if all(item.pos_encode is not None for item in batch):
            pos_encode_tensors = [item.pos_encode for item in batch]
            batched_pos_encode = torch.stack(pos_encode_tensors, dim=0)
        
        return cls(
            gt=batched_gt,
            chr=batched_chr,
            pos=batched_pos,
            pos_encode=batched_pos_encode
        )
    
    def __len__(self):
        return len(self.gt)

class CachedKGSGDataset(Dataset):
    """
    PyTorch Dataset for 1KG data using sgkit with efficient caching.
    Supports both VCF and Zarr formats with optimized data loading.
    """
    
    def __init__(self,
                 data_path: str,
                 windows: List[ChromosomeWindow],
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
    
    def _load_sample_data(self, human_id: str) -> GTData:
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
        
        # Find sample index
        sample_ids = self.ds['sample_id'].values
        sample_idx = np.where(sample_ids == human_id)[0]
        if len(sample_idx) == 0:
            raise ValueError(f"Sample {human_id} not found")
        sample_idx = sample_idx[0]
        genotype_data = self.ds["call_genotype"].values[:, sample_idx, :]
        gt_pos = self.ds["variant_position"].values
        chrom_indices = self.ds['variant_contig'].values

        alt_count = np.sum(genotype_data, axis=-1)
        missing_mask = np.any(genotype_data < 0, axis=-1)
        alt_count[missing_mask] = -1

        gt_data = GTData.from_numpy(gt=alt_count, chr=chrom_indices, pos=gt_pos)
        
        # Cache the data with LRU management
        self._add_to_cache(human_id, gt_data)
        return gt_data
    
    def _add_to_cache(self, sample_id: str, data: GTData) -> None:
        """Add sample data to cache with LRU eviction policy."""
        if len(self.sample_cache) >= self.cache_size:
            oldest_key = next(iter(self.sample_cache))
            del self.sample_cache[oldest_key]
            logger.debug(f"Evicted {oldest_key} from cache")
        
        self.sample_cache[sample_id] = data
        logger.debug(f"Cached {sample_id}, cache size: {len(self.sample_cache)}")
    
    def __len__(self) -> int:
        """Return number of windows in the dataset."""
        return len(self.windows)
    
    def __getitem__(self, idx: int) -> Tuple[GTData, Dict[str, Any]]:
        """
        Get genotype data for a specific genomic window.
        
        Args:
            idx: Window index
            
        Returns:
            Tuple of (genotype_tensor, window_info_dict)
        """
        window = self.windows[idx]
        
        # Load sample data (uses cache if available)
        gt_data = self._load_sample_data(window.human_id)
        
        window_gt = gt_data.get_window_by_position_exact(window.snp_start, window.snp_end, max_len=1024)
        window_gt = window_gt.encode_positions(resolution=20_000)
        
        # Create comprehensive window metadata
        window_info = {
            'human_id': window.human_id,
            'chromosome': window.chromosome,
            'snp_start': window.snp_start,
            'snp_end': window.snp_end,
            'window_size': window.window_size,
            'variant_count': len(window_gt.gt),
            'window_index': idx
        }
        window_gt.window_info = window_info
        return window_gt
    
    def clear_cache(self) -> None:
        """Clear sample cache to free memory."""
        self.sample_cache.clear()
        logger.info("Sample cache cleared")
    
    
    def get_cache_info(self) -> Dict[str, Any]:
        """Get current cache statistics and information."""
        return {
            'cache_size': len(self.sample_cache),
            'cached_samples': list(self.sample_cache.keys()),
            'max_cache_size': self.cache_size
        }
    
def gt_data_collate_fn(batch):
    """
    自定义 collate 函数用于处理 GTData 对象
    """
    if isinstance(batch[0], tuple):
        # 如果 batch 是 (GTData, dict) 的元组
        gt_data_list = [item[0] for item in batch]
        info_list = [item[1] for item in batch]
        
        # 合并 GTData
        batched_gt_data = GTData.collate(gt_data_list)
        
        # 合并 info dicts (保持为列表)
        return batched_gt_data, info_list
    elif isinstance(batch[0], GTData):
        # 如果 batch 直接是 GTData 对象
        return GTData.collate(batch)
    else:
        # 回退到默认的 collate
        return torch.utils.data.default_collate(batch)


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
        collate_fn=gt_data_collate_fn
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
        shuffle=False,
        num_workers=8
    )
    
    # Step 3: Use in training loop
    for batch_idx, genotypes in enumerate(dataloader):
        print(f"Batch {batch_idx}: shape:{genotypes.gt.shape} ")
        if batch_idx >= 64:  # Just show first 5 batches
            break