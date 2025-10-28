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
        # df_merge = pd.merge(df_clinical, df_snparray, on=on, how=how)
        time_tolerance = "30D"
        df_clinical['活检日期'] = pd.to_datetime(df_clinical['活检日期'])
        df_snparray['活检日期'] = pd.to_datetime(df_snparray['实验时间'], format='%Y年%m月%d日', errors='coerce')
        
        df_clinical = df_clinical.dropna(subset=['活检日期', '女方姓名-PGD/PGS编号'])
        df_snparray = df_snparray.dropna(subset=['活检日期', '女方姓名-PGD/PGS编号'])

        # 先按患者编号合并所有可能的组合
        merged = pd.merge(
            df_clinical, 
            df_snparray, 
            on='女方姓名-PGD/PGS编号', 
            suffixes=('_clinical', '_snparray')
        )
    
        # 计算时间差
        merged['时间差'] = abs(merged['活检日期_clinical'] - merged['活检日期_snparray'])
        
        # 过滤在时间容差范围内的记录
        tolerance_td = pd.Timedelta(time_tolerance)
        merged = merged[merged['时间差'] <= tolerance_td]
        
        # 为每个临床记录找到时间最接近的实验记录
        merged = merged.loc[merged.groupby(['女方姓名-PGD/PGS编号', '活检日期_clinical'])['时间差'].idxmin()]
        logging.info(f"Merged {merged.shape} data from clinical data {df_clinical.shape} and snparray data {df_snparray.shape}")
        return merged
    
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
                        default=r'D:\03.projects\AI.PGT\snparray_analysis\data\multi_芯片实验记录表_all_persons.csv')
    parser.add_argument('--output_dir', type=str, 
                        default=r'D:\03.projects\AI.PGT\snparray_analysis\work_dir')
    parser.add_argument('--output_name', type=str, 
                        default=r'merge_clinical_snparray')
    args = parser.parse_args()

    setup_logging()

    merge_data = MergeData(args.clinical_path, args.snparray_path, args.output_dir, args.output_name)
    merge_data.run()