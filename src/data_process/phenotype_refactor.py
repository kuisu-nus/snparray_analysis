"""
Refactor phenotype data for genetic analysis
Handles PED files and phenotype data integration with Chinese character conversion
"""

import pandas as pd
import logging
import argparse
import sys
from pathlib import Path
from typing import Dict, Optional, List
from pypinyin import lazy_pinyin

# Configure logging
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)


class Configuration:
    """Configuration management for the refactor process"""
    
    def __init__(
        self,
        ped_file: str,
        phenotype_file: str,
        full_phenotype_file: str,
        output_ped: str,
        output_phenotype: str,
        output_full_phenotype: str,
        column_rename: Optional[Dict[str, str]] = None,
        required_columns: Optional[List[str]] = None,
        log_level: str = "INFO"
    ):
        self.ped_file = ped_file
        self.phenotype_file = phenotype_file
        self.full_phenotype_file = full_phenotype_file
        self.output_ped = output_ped
        self.output_phenotype = output_phenotype
        self.output_full_phenotype = output_full_phenotype
        self.column_rename = column_rename or {}
        self.required_columns = required_columns or ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
        self.log_level = log_level
        
        # Validate paths
        self._validate_paths()
    
    def _validate_paths(self):
        """Validate that input files exist"""
        input_files = [
            (self.ped_file, "PED file"),
            (self.phenotype_file, "Phenotype file"), 
            (self.full_phenotype_file, "Full phenotype file")
        ]
        
        for file_path, file_type in input_files:
            if not Path(file_path).exists():
                raise FileNotFoundError(f"{file_type} not found: {file_path}")
        
        logger.info("All input files validated successfully")


class PedFileProcessor:
    """Process PED format files for genetic analysis"""
    
    @staticmethod
    def extract_sample_info(ped_file: str) -> pd.DataFrame:
        """
        Extract sample information from PED file
        
        Args:
            ped_file: Path to PED file
            
        Returns:
            DataFrame with columns: FID, IID, PID, MID, SEX, PHENO
        """
        try:
            with open(ped_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            sample_data = [line[:100].split()[:6] for line in lines]
            df_sample = pd.DataFrame(
                sample_data, 
                columns=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO']
            )
            
            logger.info(f"Extracted {len(df_sample)} samples from PED file")
            logger.debug(f"Sample data preview: {df_sample.head()}")
            
            return df_sample
            
        except Exception as e:
            logger.error(f"Error reading PED file {ped_file}: {e}")
            raise


class PhenotypeDataProcessor:
    """Process and transform phenotype data"""
    
    def __init__(self, config: Configuration):
        self.config = config
    
    @staticmethod
    def contains_chinese(text: str) -> bool:
        """
        Check if text contains Chinese characters
        
        Args:
            text: Input text to check
            
        Returns:
            True if Chinese characters are found
        """
        if pd.isna(text):
            return False
        return any('\u4e00' <= char <= '\u9fff' for char in str(text))
    
    @staticmethod
    def convert_chinese_to_pinyin(text: str) -> str:
        """
        Convert Chinese characters to Pinyin
        
        Args:
            text: Text containing Chinese characters
            
        Returns:
            Text with Chinese characters converted to Pinyin
        """
        if pd.isna(text):
            return text
        return '_'.join(lazy_pinyin(str(text)))
    
    def process_phenotype_data(self, sample_ids: pd.Series) -> pd.DataFrame:
        """
        Process phenotype data and merge with sample IDs
        
        Args:
            sample_ids: Series of sample IDs to filter by
            
        Returns:
            Processed phenotype DataFrame
        """
        try:
            # Read phenotype data
            df_phenotype = pd.read_csv(self.config.full_phenotype_file)
            logger.info(f"Loaded phenotype data: {df_phenotype.shape}")
            
            # Create IID column if not exists
            if 'IID' not in df_phenotype.columns:
                df_phenotype['IID'] = df_phenotype.apply(
                    lambda row: f"{row['idx']}_{row['chip_sub_idx']}", 
                    axis=1
                )
            
            # Merge with sample IDs
            df_merged = pd.merge(sample_ids, df_phenotype, on='IID', how='left')
            logger.info(f"Merged phenotype data: {df_merged.shape}")
            
            # Apply column renaming
            if self.config.column_rename:
                actual_rename = {
                    k: v for k, v in self.config.column_rename.items() 
                    if k in df_merged.columns
                }
                logging.info(f"df merge head: \n{df_merged.head()}")
                if actual_rename:
                    df_merged.rename(columns=actual_rename, inplace=True)
                    logger.info(f"Renamed columns: {actual_rename}")
            
            # Ensure required columns exist
            df_merged = self._ensure_required_columns(df_merged)
            
            # Convert Chinese text to Pinyin
            df_merged = self._convert_chinese_columns(df_merged)
            
            logger.debug(f"Processed phenotype data preview:\n{df_merged.head()}")
            return df_merged
            
        except Exception as e:
            logger.error(f"Error processing phenotype data {self.config.full_phenotype_file}: {e}")
            raise
    
    def _ensure_required_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure required columns exist in phenotype data"""
        # Add missing required columns with default values
        column_defaults = {
            "SEX": 0,
            "PHENO": -9,  # Standard missing value code in PLINK
            "FID": "FAM001",  # Default family ID
            "PID": "0",  # Default parent ID (0 = missing)
            "MID": "0"   # Default mother ID (0 = missing)
        }
        
        for col in self.config.required_columns:
            if col not in df.columns and col in column_defaults:
                df[col] = column_defaults[col]
                logger.info(f"Added missing {col} column with default value: {column_defaults[col]}")
        
        # Convert phenotype values if needed
        if "PHENO" in df.columns:
            df = self._standardize_phenotype_values(df)
        
        return df
    
    def _standardize_phenotype_values(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize phenotype values to PLINK format"""
        unique_values = df['PHENO'].unique()
        
        # Convert Chinese phenotype values
        if "是" in unique_values or "否" in unique_values:
            df['PHENO'] = df['PHENO'].map({'否': 1, '是': 2, 'unknown': -9})
            logger.info("Converted Chinese PHENO values to PLINK format (1=control, 2=case)")
        
        # Convert other common formats
        elif set(unique_values) & {"0", "1", "2"}:
            # Already in PLINK format
            pass
        else:
            logger.warning(f"Unrecognized PHENO values: {unique_values}. Please verify the mapping.")
        
        return df
    
    def _convert_chinese_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert Chinese characters to Pinyin in all columns"""
        chinese_columns = []
        for col in df.columns:
            if df[col].apply(self.contains_chinese).any():
                chinese_columns.append(col)
                df[col] = df[col].apply(
                    lambda x: self.convert_chinese_to_pinyin(x) if pd.notna(x) else x
                )
        
        if chinese_columns:
            logger.info(f"Converted Chinese characters to Pinyin in columns: {chinese_columns}")
        
        return df


class PhenotypeRefactor:
    """Main class for refactoring phenotype data"""
    
    def __init__(self, config: Configuration):
        """
        Initialize refactor with configuration
        
        Args:
            config: Configuration object with all settings
        """
        self.config = config
        self.ped_processor = PedFileProcessor()
        self.phenotype_processor = PhenotypeDataProcessor(config)
        
        # Initialize data
        self._load_data()
    
    def _load_data(self):
        """Load and process initial data"""
        logger.info("Loading sample information from PED file")
        df_sample = self.ped_processor.extract_sample_info(self.config.ped_file)
        self.sample_ids = df_sample["IID"]
        
        logger.info("Processing full phenotype data")
        self.full_phenotype_df = self.phenotype_processor.process_phenotype_data(
            self.sample_ids
        )
    
    def refactor_ped_file(self) -> None:
        """Refactor PED file with updated phenotype information"""
        try:
            logger.info(f"Refactoring PED file: {self.config.ped_file} -> {self.config.output_ped}")
            
            # Read original PED file
            with open(self.config.ped_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Extract phenotype information
            ped_phenotype_data = self.full_phenotype_df[
                ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
            ]
            
            # Ensure output directory exists
            Path(self.config.output_ped).parent.mkdir(parents=True, exist_ok=True)
            
            # Write refactored data
            with open(self.config.output_ped, 'w', encoding='utf-8') as f:
                for line, sample_info in zip(lines, ped_phenotype_data.to_numpy()):
                    # Extract SNP data
                    snp_data = "\t".join(line[:100].split()[6:]) + line[100:]
                    
                    # Create new line with updated sample info
                    sample_info_str = "\t".join(str(item) for item in sample_info)
                    new_line = f"{sample_info_str}\t{snp_data}"
                    
                    f.write(new_line)
            
            logger.info(f"Successfully refactored PED file: {self.config.output_ped}")
            
        except Exception as e:
            logger.error(f"Error refactoring PED file: {e}")
            raise
    
    def save_phenotype_data(self) -> None:
        """Save all phenotype data files"""
        try:
            # Save full phenotype data
            self.full_phenotype_df.to_csv(self.config.output_full_phenotype, index=False)
            logger.info(f"Saved full phenotype data to: {self.config.output_full_phenotype}")
            
            # Save simplified phenotype file (PLINK format)
            simple_phenotype = self.full_phenotype_df[["FID", "IID", "PHENO"]]
            simple_phenotype.to_csv(self.config.output_phenotype, index=False, sep="\t")
            logger.info(f"Saved simplified phenotype data to: {self.config.output_phenotype}")
            
        except Exception as e:
            logger.error(f"Error saving phenotype data: {e}")
            raise
    
    def run(self) -> None:
        """Execute complete refactoring pipeline"""
        logger.info("Starting phenotype data refactoring pipeline")
        
        self.refactor_ped_file()
        self.save_phenotype_data()
        
        logger.info("Phenotype data refactoring completed successfully")


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Refactor phenotype data for genetic analysis"
    )
    
    # Input files
    parser.add_argument("--ped-file", required=True, help="Input PED file path")
    parser.add_argument("--phenotype-file", required=True, help="Input phenotype file path")
    parser.add_argument("--full-phenotype-file", required=True, help="Input full phenotype file path")
    
    # Output files
    parser.add_argument("--output-ped", required=True, help="Output refactored PED file path")
    parser.add_argument("--output-phenotype", required=True, help="Output phenotype file path")
    parser.add_argument("--output-full-phenotype", required=True, help="Output full phenotype file path")
    
    # Optional parameters
    parser.add_argument("--rename", nargs='+', help="Column rename mappings (format: old:new)")
    parser.add_argument("--log-level", default="INFO", 
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       help="Logging level")
    
    return parser.parse_args()


def create_config_from_args(args) -> Configuration:
    """Create Configuration object from command line arguments"""
    
    # Parse rename mappings
    column_rename = {}
    if args.rename:
        for mapping in args.rename:
            if ':' in mapping:
                old, new = mapping.split(':', 1)
                column_rename[old.strip()] = new.strip()
    
    config = Configuration(
        ped_file=args.ped_file,
        phenotype_file=args.phenotype_file,
        full_phenotype_file=args.full_phenotype_file,
        output_ped=args.output_ped,
        output_phenotype=args.output_phenotype,
        output_full_phenotype=args.output_full_phenotype,
        column_rename=column_rename,
        log_level=args.log_level
    )
    
    return config


def get_default_config() -> Configuration:
    """Get default configuration for direct execution"""
    return Configuration(
        ped_file=r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS_ALL\PLINK_291025_0948\PGT_TLS_ALL.ped",
        phenotype_file=r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS_ALL\PLINK_291025_0948\PGT_TLS_ALL.phenotype",
        full_phenotype_file=r"D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype.csv",
        output_ped=r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS_ALL\PLINK_291025_0948\PGT_TLS_ALL.refactor.ped",
        output_phenotype=r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS_ALL\PLINK_291025_0948\PGT_TLS_ALL.refactor.phenotype",
        output_full_phenotype=r"D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype.refactor.csv",
        column_rename={
            "男方姓名": "PID",
            "女方姓名_clinical": "MID", 
            "诊断可移植": "PHENO"
        }
    )


def main():
    """Main execution function"""
    try:
        # Check if running with command line arguments
        if len(sys.argv) > 1:
            args = parse_arguments()
            config = create_config_from_args(args)
            logging.info(f"configs: \n{config}")
        else:
            # Use default configuration
            config = get_default_config()
            logger.info("Using default configuration")
        
        # Set log level
        logging.getLogger().setLevel(config.log_level)
        
        # Initialize and run refactor
        refactor = PhenotypeRefactor(config)
        refactor.run()
        
        logger.info("Phenotype data refactoring completed successfully")
        
    except Exception as e:
        logger.error(f"Phenotype refactoring failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()