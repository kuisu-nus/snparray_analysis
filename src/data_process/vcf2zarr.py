import allel
import zarr

# 方法2.1: 直接转换为Zarr格式（基于HDF5的现代替代）
def vcf_to_zarr(vcf_path, zarr_path):
    """使用scikit-allel将VCF转换为Zarr格式"""
    allel.vcf_to_zarr(vcf_path, zarr_path, group='/', fields='*', overwrite=True)

# 使用示例
vcf_to_zarr('input.vcf.gz', 'output.zarr')

# 方法2.2: 读取VCF并保存为HDF5
def vcf_to_h5_allel(vcf_path, h5_path):
    """使用scikit-allel读取VCF并保存为HDF5"""
    
    # 读取VCF文件
    callset = allel.read_vcf(vcf_path, fields=['variants/CHROM', 'variants/POS', 
                                              'variants/REF', 'variants/ALT',
                                              'calldata/GT'])
    
    # 转换为HDF5
    with h5py.File(h5_path, 'w') as h5f:
        for key, value in callset.items():
            if value is not None:
                h5f.create_dataset(key, data=value)
        
        h5f.attrs['source_vcf'] = vcf_path

# 使用示例
vcf_to_h5_allel('input.vcf.gz', 'output.h5')