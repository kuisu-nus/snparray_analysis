import logging
from pathlib import Path
import pandas as pd
import numpy as np
import os
import sys
import argparse
import glob
import re

def setup_logging() -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

class MergeData:
    def __init__(self, clinical_path:str, snparray_path, output_dir, output_name:str):
        self.clinical_path = Path(clinical_path)
        assert self.clinical_path.exists(), f"File {self.clinical_path} does not exist"
        self.snparray_path = Path(snparray_path)
        assert self.snparray_path.exists(), f"File {self.snparray_path} does not exist"
        self.output_dir = output_dir
        self.output_file = output_name
    
    def load_data(self, file_path, header=0):
        df_data = pd.read_csv(file_path, sep=',', header=header)
        if "女方姓名" and "PGD/PGS编号" in df_data.columns:
            df_data["女方姓名-PGD/PGS编号"] = df_data.apply(lambda row: f"{row['女方姓名']}-{row['PGD/PGS编号']}", axis=1)
            logging.info("Added new column '女方姓名-PGD/PGS编号' to the data")
        logging.info(f"Loaded {df_data.shape} data from {file_path}")
        logging.info(f"Data columns: {df_data.columns}")
        logging.info(f"Data head: \n{df_data.head()}")
        return df_data

    def merge_data(self, df_clinical, df_snparray, on='女方姓名-PGD/PGS编号', how='inner'):
        # Get all the CSV files in the input directory
        df_merge = pd.merge(df_clinical, df_snparray, on=on, how=how)
        logging.info(f"Merged {df_merge.shape} data")
        return df_merge
    
    def save_data(self, df_data, output_dir, output_file):
        df_data.to_csv(os.path.join(output_dir, output_file+".csv"), index=False)

    def run(self):
        df_clinical = self.load_data(self.clinical_path, header=1)
        df_snparray = self.load_data(self.snparray_path, header=0)
        df_merge = self.merge_data(df_clinical, df_snparray)
        self.save_data(df_merge, self.output_dir, self.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--clinical_path', type=str, 
                        default=r'D:\03.projects\AI.PGT\snparray_analysis\data\clinical_data_pgt.csv')
    parser.add_argument('--snparray_path', type=str, 
                        default=r'D:\03.projects\AI.PGT\snparray_analysis\data\5.芯片实验记录表2020.12.29_persons.csv')
    parser.add_argument('--output_dir', type=str, 
                        default=r'D:\03.projects\AI.PGT\snparray_analysis\work_dir')
    parser.add_argument('--output_name', type=str, 
                        default=r'merge_clinical_snparray')
    args = parser.parse_args()

    setup_logging()

    merge_data = MergeData(args.clinical_path, args.snparray_path, args.output_dir, args.output_name)
    merge_data.run()