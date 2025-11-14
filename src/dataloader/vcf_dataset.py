import hail as hl
import torch
from torch.utils.data import Dataset, DataLoader
from typing import Optional, Tuple, List
import numpy as np
import os

proxy_vars = [
    'http_proxy', 'https_proxy', 'ftp_proxy',
    'HTTP_PROXY', 'HTTPS_PROXY', 'FTP_PROXY',
    'no_proxy', 'NO_PROXY'
]

for var in proxy_vars:
    os.environ.pop(var, None)  # 使用pop避免KeyError



class KGDataLoader(Dataset):
    """
    A PyTorch DataLoader for 1000 Genomes (1KG) data using Hail.
    Handles genotype data loading and preprocessing for ML models.
    """
    
    def __init__(self, 
                 mt_path: str,
                 sample_indices: Optional[List[int]] = None,
                 variant_indices: Optional[List[int]] = None,
                 missing_value=-1,
                 split: str = 'train',
                 train_ratio: float = 0.8,
                 val_ratio: float = 0.1,
                 random_seed: int = 42):
        """
        Initialize the 1KG DataLoader.
        
        Args:
            mt_path: Path to the Hail MatrixTable
            sample_indices: Optional list of sample indices to subset
            variant_indices: Optional list of variant indices to subset  
            split: Data split - 'train', 'val', or 'test'
            train_ratio: Proportion of data for training
            val_ratio: Proportion of data for validation
            random_seed: Random seed for reproducibility
        """
        # Initialize Hail if not already done
        # if not hl._initialized:
        #     hl.init()
        self.missing_value = -1
        self.mt_path = mt_path
        self.split = split
        self.random_seed = random_seed
        
        # Load MatrixTable
        self.mt = hl.read_matrix_table(mt_path)
        
        # Apply sample and variant filtering if provided
        if sample_indices is not None:
            samples = self.mt.cols().collect()
            selected_samples = [samples[i].s for i in sample_indices]
            self.mt = self.mt.filter_cols(hl.literal(selected_samples).contains(self.mt.s))
            
        if variant_indices is not None:
            variants = self.mt.rows().collect()
            selected_variants = [variants[i].locus for i in variant_indices]
            self.mt = self.mt.filter_rows(hl.literal(selected_variants).contains(self.mt.locus))
        
        # Split data
        self._setup_data_split(train_ratio, val_ratio)
        
        # Cache genotype matrix dimensions
        self.n_samples = self.mt.count_cols()
        self.n_variants = self.mt.count_rows()
        
    def _setup_data_split(self, train_ratio: float, val_ratio: float) -> None:
        """Split data into train/val/test sets."""
        # Add random split annotations
        self.mt = self.mt.annotate_cols(
            split=hl.rand_unif(0, 1, seed=self.random_seed)
        )
        
        # Filter based on split
        if self.split == 'train':
            self.mt = self.mt.filter_cols(self.mt.split < train_ratio)
        elif self.split == 'val':
            self.mt = self.mt.filter_cols(
                (self.mt.split >= train_ratio) & 
                (self.mt.split < train_ratio + val_ratio)
            )
        else:  # test
            self.mt = self.mt.filter_cols(self.mt.split >= train_ratio + val_ratio)
    
    def __len__(self) -> int:
        """Return number of samples in the dataset."""
        return self.n_samples
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Get genotype data for a specific sample.
        
        Args:
            idx: Sample index
            
        Returns:
            Tuple of (genotypes, sample_info) where:
            - genotypes: Tensor of shape [n_variants]
            - sample_info: Tensor containing sample metadata
        """
        # Get specific sample
        samples = self.mt.cols().collect()
        sample_id = samples[idx].s
        
        # Filter to specific sample and get genotypes
        sample_mt = self.mt.filter_cols(self.mt.s == sample_id)
        
        # Convert genotypes to tensor (GT: 0=hom_ref, 1=het, 2=hom_var)
        genotypes = sample_mt.GT.collect()
        genotype_tensor = torch.tensor([g.alleles for g in genotypes], dtype=torch.float32)
        
        # Get sample info (example: using population code if available)
        sample_info = torch.tensor([idx], dtype=torch.long)  # Placeholder
        
        return genotype_tensor, sample_info
    
    def get_variant_info(self) -> hl.Table:
        """Get variant information table."""
        return self.mt.rows()
    
    def get_sample_info(self) -> hl.Table:
        """Get sample information table."""
        return self.mt.cols()


class BatchKGDataLoader(DataLoader):
    """
    Batched DataLoader for 1KG data with optimized performance.
    """
    
    def __init__(self, 
                 dataset: KGDataLoader,
                 batch_size: int = 32,
                 shuffle: bool = True,
                 **kwargs):
        """
        Initialize batched data loader.
        
        Args:
            dataset: KGDataLoader instance
            batch_size: Number of samples per batch
            shuffle: Whether to shuffle data
            **kwargs: Additional DataLoader arguments
        """
        super().__init__(
            dataset=dataset,
            batch_size=batch_size,
            shuffle=shuffle,
            **kwargs
        )


# Test cases
def test_kg_dataloader():
    """Test cases for KGDataLoader functionality."""
    
    # Mock a small MatrixTable for testing
    print("Testing KGDataLoader...")
    
    # Create mock data
    # mock_mt = create_mock_matrix_table()
    test_path = "/home/sukui/03.projects/01.pgt/snparray_analysis/work_dir/hail_data/1kg_test.mt"
    # mock_mt.write(test_path)
    
    try:
        # Test initialization
        dataloader = KGDataLoader(test_path, split='train')
        print(f"✓ DataLoader initialized successfully")
        print(f"✓ Samples: {dataloader.n_samples}, Variants: {dataloader.n_variants}")
        
        # Test length
        assert len(dataloader) == dataloader.n_samples, "Length mismatch"
        print("✓ __len__ method works correctly")
        
        # Test item access
        genotype, sample_info = dataloader[0]
        assert isinstance(genotype, torch.Tensor), "Genotype should be tensor"
        assert isinstance(sample_info, torch.Tensor), "Sample info should be tensor"
        assert genotype.shape[0] == dataloader.n_variants, "Genotype shape mismatch"
        print("✓ __getitem__ method works correctly")
        
        # Test batch loader
        batch_loader = BatchKGDataLoader(dataloader, batch_size=2, shuffle=False)
        batch = next(iter(batch_loader))
        assert len(batch) == 2, "Batch should contain two elements"
        assert batch[0].shape[0] == 2, "Batch dimension incorrect"
        print("✓ Batch loading works correctly")
        
        # Test variant and sample info methods
        variant_info = dataloader.get_variant_info()
        sample_info = dataloader.get_sample_info()
        assert variant_info.count() > 0, "Variant info should not be empty"
        assert sample_info.count() > 0, "Sample info should not be empty"
        print("✓ Info methods work correctly")
        
        print("All tests passed! ✓")
        
    except Exception as e:
        print(f"Test failed: {e}")
        raise
    finally:
        # Cleanup
        pass
        # import shutil
        # shutil.rmtree(test_path)


def create_mock_matrix_table():
    """Create a mock MatrixTable for testing."""
    # Define sample data
    samples = [f"sample_{i}" for i in range(10)]
    
    # Define variant data  
    loci = [hl.locus('1', i) for i in range(1, 101)]  # 注意：应该是 hl.locus 不是 hl.Locus
    alleles = [['A', 'G'] for _ in range(100)]
    
    # Create mock entries
    entries = []
    for sample in samples:
        for j, (locus, allele) in enumerate(zip(loci, alleles)):
            # Random genotypes: 0=hom_ref, 1=het, 2=hom_var
            gt = hl.Call([0, 0]) if j % 3 == 0 else hl.Call([0, 1]) if j % 3 == 1 else hl.Call([1, 1])
            entries.append(hl.struct(s=sample, locus=locus, alleles=allele, GT=gt))
    
    # 直接使用 entries 创建 Table
    table = hl.Table.parallelize(entries)
    
    # 转换为 MatrixTable
    mt = table.to_matrix_table(
        row_key=['locus', 'alleles'],
        col_key=['s']
    )
    
    return mt


if __name__ == "__main__":
    # Run tests
    test_kg_dataloader()