import json
import pandas as pd
import os
import numpy as np
from pathlib import Path
import argparse
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SnpArrayData:
    def __init__(self, data_root, file_name, output_root, sep='\t'):
        self.file_name = file_name
        self.data_root = Path(data_root)
        assert self.data_root.exists(), f"Data root path {data_root} does not exist."
        self.output_root = Path(output_root)
        assert self.output_root.exists(), f"Output root path {output_root} does not exist."
        self.sep = sep

        self.static_info = {}
        self.static_info["data_root"] = str(self.data_root)
        self.static_info["file_name"] = self.file_name

        logger.info(f"Initialized SnpArrayData with data_root: {data_root}, file_name: {file_name}, output_root: {output_root}, sep: {sep}")

    def load_data(self, file_type=".ped"):
        """
        Load SNP array data from specified file type. Currently supports .ped and .fam file types.
        file_path = self.data_root / f"{self.file_name}{file_type}"""
        file_path = self.data_root / f"{self.file_name}{file_type}"
        assert file_path.exists(), f"File {file_path} does not exist."

        if file_type == ".ped":
            df = pd.read_csv(file_path, sep=self.sep, header=None, dtype=str)
            self.static_info["ped_shape"] = str(df.shape)
            # self.static_info["ped_head"] = df.iloc[:5, :10]
            logger.info(f"Loaded .ped file with shape: {df.shape}")
        elif file_type == ".phenotype":
            df = pd.read_csv(file_path, sep=self.sep, header=None)
            self.static_info["fam_shape"] = str(df.shape)
            self.static_info["fam_head"] = str(df.head())
            logger.info(f"Loaded .fam file with shape: {df.shape}")
        elif file_type == ".map":
            df = pd.read_csv(file_path, sep=self.sep, header=None)
            self.static_info["map_shape"] = str(df.shape)
            self.static_info["map_head"] = str(df.head())
            logger.info(f"Loaded .map file with shape: {df.shape}")
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

        return df
    
    def process_data(self):
        """"""
        ped_df = self.load_data(file_type=".ped")
        fam_df = self.load_data(file_type=".phenotype")
        map_df = self.load_data(file_type=".map")
        self.dump_json()

    def dump_json(self, output_path=None):
        if output_path is None:
            output_path = f"{self.file_name}_static_info.json"
        output_path = self.output_root / output_path

        # transformer DataFrame to serializable format
        serializable_info = {}

        for key, value in self.static_info.items():
            if isinstance(value, pd.DataFrame):
                serializable_info[key] = {
                    'data': value.to_dict('records'),  # 'list', 'dict'
                    'columns': value.columns.tolist(),
                    'shape': value.shape
                }
            elif isinstance(value, (pd.Series, np.ndarray)):
                serializable_info[key] = value.tolist()
            else:
                serializable_info[key] = value

        print(f"static_info: {json.dumps(serializable_info, indent=2)}")
        with open(output_path, 'w') as f:
            json.dump(serializable_info, f, indent=2)
        logger.info(f"Dumped static info to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Process SNP array data.')
    parser.add_argument('--data_root', type=str, required=True, help='Path to the root directory of the data.')
    parser.add_argument('--file_name', type=str, required=True, help='Name of the file to process.')
    parser.add_argument('--output_root', type=str, required=True, help='Path to the root directory of the output.')
    parser.add_argument('--sep', type=str, default='\t', help='Separator used in the input file.')
    args = parser.parse_args()

    snp_array_data = SnpArrayData(args.data_root, args.file_name, args.output_root, args.sep)
    snp_array_data.process_data()


if __name__ == "__main__":
    main()