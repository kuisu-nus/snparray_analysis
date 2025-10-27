"""
SNP Experiment Data Parser

This module provides functionality to parse SNP experiment data from text files
and convert it to structured JSON format.
"""

import argparse
import json
import logging
import re
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any


class ChipInfoParser:
    """Parser for chip information in SNP experiment data."""
    
    def __init__(self) -> None:
        """Initialize the ChipInfoParser with required mappings."""
        self.experiments: Optional[Dict] = None
        self.chip_idx_dict = {
            "chip_num_1": "X1 number",
            "chip_num_2": "X2 number", 
            "chip_num_4": "X3 number",
            "chip_num_5": "X4 number"
        }
        self.chip_info = {}
        self.name_dict = {}

    def parse_name(self, text: str) -> Optional[Tuple[str, str, str]]:
        """
        Parse a name string containing Chinese characters, initials and numbers.
        
        Args:
            text: Input text to parse
            
        Returns:
            Tuple of (chinese_name, initial, number) if parsing successful, None otherwise
        """
        pattern = r'^([\u4e00-\u9fa5A-Z]+)\s*([A-Z])\s*(\d{1,2})$'
        match = re.match(pattern, text)
        
        if match:
            chinese = match.group(1)
            initial = match.group(2) if match.group(2) else ''
            number = match.group(3)
            
            # Validate number range (1-20)
            if 1 <= int(number) <= 20:
                return chinese, initial, number
        
        return None

    def parse_number(self, text: str) -> Optional[Tuple[str, str]]:
        """
        Parse a string containing initials and numbers.
        
        Args:
            text: Input text to parse
            
        Returns:
            Tuple of (initial, number) if parsing successful, None otherwise
        """
        pattern = r'^([A-Z])(\d{1,2})$'
        match = re.match(pattern, text)
        
        if match:
            initial = match.group(1)
            number = match.group(2)
            
            # Validate number range (2-20)
            if 2 <= int(number) <= 20:
                return initial, number
        
        return None

    def parse_sample_name(self, name: str, name_dict: Dict) -> Tuple[str, str]:
        """
        Parse a single sample name and extract person name and ID prefix.
        
        Args:
            name: Sample name to parse
            name_dict: Dictionary mapping initials to Chinese names
            
        Returns:
            Tuple of (person_name, id_prefix)
        """
        person_name = name
        id_prefix = name
        
        if "-" in name:
            id_prefix = name.split("-")[0]
            person_name = name.split("-")[1]
        elif self.parse_name(name):
            chinese, initial, number = self.parse_name(name)
            name_dict[initial] = chinese
            person_name = chinese
            id_prefix = number
        elif self.parse_number(name):
            initial, number = self.parse_number(name)
            person_name = name_dict.get(initial, name)
            id_prefix = number
        elif name == "RB":
            person_name = "RB"
            id_prefix = "1"
        else:
            logging.warning(f"Unable to parse name: {name}")
        
        return person_name, id_prefix

    def process_chip_data(self, chip_data: List[str], name_dict: Dict) -> Tuple[List[str], List[str]]:
        """
        Process chip data to extract names and name numbers.
        
        Args:
            chip_data: List of chip data entries
            name_dict: Dictionary mapping initials to Chinese names
            
        Returns:
            Tuple of (names_list, name_numbers_list)
        """
        names = []
        name_numbers = []
        
        for name in chip_data:
            person_name, id_prefix = self.parse_sample_name(name, name_dict)
            names.append(person_name)
            name_numbers.append(id_prefix)
            
            # Update name dictionary
            if person_name not in name_dict:
                name_dict[person_name] = []
            name_dict[person_name].append(id_prefix)
        
        return names, name_numbers

    def parse(self, experiments: Dict) -> List[Dict]:
        """
        Parse experiment data and extract chip information.
        
        Args:
            experiments: Dictionary containing experiment data
            
        Returns:
            List of dictionaries containing parsed chip information
        """
        chip_infos = []
        name_dict = {}
        experiment_info = {}
        
        for key, value in experiments.items():
            if key in self.chip_idx_dict:
                chip_info = self._process_chip_entry(key, value, name_dict, experiments)
                if chip_info:
                    chip_infos.append(chip_info)
            elif key in self.chip_idx_dict.values():
                continue
            else:
                experiment_info[key] = value
        
        # Add experiment info to all chip entries
        for chip_info in chip_infos:
            chip_info.update(experiment_info)
        
        return chip_infos

    def _process_chip_entry(self, key: str, value: List[str], name_dict: Dict, 
                           experiments: Dict) -> Optional[Dict]:
        """
        Process a single chip entry.
        
        Args:
            key: Chip identifier key
            value: Chip data values
            name_dict: Dictionary for name mappings
            experiments: Full experiments dictionary
            
        Returns:
            Dictionary with processed chip info or None if invalid
        """
        chip_info = {
            "type": value[0],
            "sample": value[1:-1]
        }
        
        # Process names and numbers
        names, name_numbers = self.process_chip_data(chip_info["sample"], name_dict)
        chip_info["names"] = names
        chip_info["name_numbers"] = name_numbers
        
        # Get chip index
        chip_idx_key = self.chip_idx_dict[key]
        if chip_idx_key in experiments:
            chip_info["idx"] = experiments[chip_idx_key]
        else:
            raise ValueError(f"Missing chip index information: {chip_idx_key}")
        
        # Validate sample count
        if not (len(names) == len(name_numbers) == 12):
            raise ValueError(
                f"Sample count mismatch (expected 12): "
                f"names({len(names)}), name_numbers({len(name_numbers)})"
            )
        
        return chip_info


class TextFileParser:
    """Parser for SNP experiment text files."""
    
    def __init__(self, chip_parser: ChipInfoParser) -> None:
        """
        Initialize the TextFileParser.
        
        Args:
            chip_parser: Instance of ChipInfoParser for processing chip data
        """
        self.chip_parser = chip_parser
        self.experiments = []

    def load_data(self, data_path: str) -> List[str]:
        """
        Load data from text file.
        
        Args:
            data_path: Path to the data file
            
        Returns:
            List of lines from the file
        """
        data_path_obj = Path(data_path)
        if not data_path_obj.exists():
            raise FileNotFoundError(f"Data path {data_path} does not exist.")
        
        with open(data_path_obj, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        return lines

    def parse_experiment_time(self, line_data: List[str], experiment_dict: Dict) -> None:
        """
        Parse experiment time and experimenter information.
        
        Args:
            line_data: Split line data
            experiment_dict: Dictionary to store experiment data
        """
        time_key, time_value = line_data[1].split("：")
        experiment_dict[time_key.strip()] = time_value.strip().replace(" ", "")
        
        experimenter_key, experimenter_value = line_data[7].split("：")
        experiment_dict[experimenter_key.strip()] = experimenter_value.strip()

    def parse_comment(self, line_data: List[str], experiment_dict: Dict) -> None:
        """
        Parse comment information.
        
        Args:
            line_data: Split line data
            experiment_dict: Dictionary to store experiment data
        """
        comment_key = line_data[0][:2].strip()
        comment_value = line_data[0][3:].strip()
        experiment_dict[comment_key] = comment_value

    def parse_chip_numbers(self, line_data: List[str], experiment_dict: Dict, 
                          chip_line: int) -> None:
        """
        Parse chip number information.
        
        Args:
            line_data: Split line data
            experiment_dict: Dictionary to store experiment data
            chip_line: Current chip line number
        """
        if chip_line in [1, 2, 4, 5]:  # chip_num_1 to chip_num_5
            experiment_dict[f"chip_num_{chip_line}"] = line_data

    def parse_chip_indices(self, line_data: List[str], experiment_dict: Dict) -> None:
        """
        Parse chip index information (X1-X4 numbers).
        
        Args:
            line_data: Split line data
            experiment_dict: Dictionary to store experiment data
        """
        # Parse X1 and X2
        X1_key, X1_value = line_data[1].split("：")
        experiment_dict[X1_key.strip()] = X1_value.strip().replace(" ", "")
        
        X2_key, X2_value = line_data[7].split("：")
        experiment_dict[X2_key.strip()] = X2_value.strip()
        
        # Parse X3 and X4 (assuming similar structure in next line)
        # This would need to be handled in the calling function

    def parse_txt(self, lines: List[str]) -> None:
        """
        Parse text data and extract experiment information.
        
        Args:
            lines: List of lines from the input file
        """
        experiment_dict = {}
        chip_line = -1
        
        for line in lines:
            line_data = [i.strip() for i in line.strip().split(",")]
            
            if "实验时间" in line:
                chip_line = 0
                if experiment_dict:
                    # Process completed experiment
                    parser_experiment = self.chip_parser.parse(experiment_dict)
                    self.experiments.extend(parser_experiment)
                    experiment_dict = {}
                
                self.parse_experiment_time(line_data, experiment_dict)
                chip_line += 1
            
            elif "备注" in line:
                self.parse_comment(line_data, experiment_dict)
                chip_line += 1
            
            elif chip_line in [0, 3]:
                chip_line += 1
                continue
            
            elif 1 <= chip_line <= 5:
                self.parse_chip_numbers(line_data, experiment_dict, chip_line)
                chip_line += 1
            
            elif chip_line in [6,7]:
                self.parse_chip_indices(line_data, experiment_dict)
                chip_line += 1
        
        # Process the last experiment
        if experiment_dict:
            parser_experiment = self.chip_parser.parse(experiment_dict)
            self.experiments.extend(parser_experiment)
        
        logging.info(f"Parsed total {len(self.experiments)} experiments")


class ParseSNPExperimentData:
    """Main class for parsing SNP experiment data."""
    
    def __init__(self, data_path: str) -> None:
        """
        Initialize the SNP experiment data parser.
        
        Args:
            data_path: Path to the data file
        """
        self.data_path = Path(data_path)
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data path {data_path} does not exist.")
        
        self.experiments = []
        self.chip_parser = ChipInfoParser()
        self.text_parser = TextFileParser(self.chip_parser)

    def load_data(self, data_path: str) -> List[str]:
        """
        Load data from the specified path.
        
        Args:
            data_path: Path to the data file
            
        Returns:
            List of lines from the file
        """
        return self.text_parser.load_data(data_path)

    def parse_txt(self, lines: List[str]) -> None:
        """
        Parse text data and extract experiment information.
        
        Args:
            lines: List of lines from the input file
        """
        self.text_parser.parse_txt(lines)
        self.experiments = self.text_parser.experiments

    def save_json(self, output_dir: str, file_name:str) -> None:
        """
        Save parsed experiments to JSON file.
        
        Args:
            output_path: Path where JSON file will be saved
        """
        output_path = f"{output_dir}/{file_name}.json"
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path_obj, 'w', encoding='utf-8') as file:
            json.dump(self.experiments, file, ensure_ascii=False, indent=4)
        
        logging.info(f"Saved experiments to {output_path}")
    
    def save_csv(self, output_dir: str, file_name:str) -> None:
        """
        Save parsed experiments to CSV file.
        
        Args:
            output_path: Path where CSV file will be saved
        """
        output_path = f"{output_dir}/{file_name}.csv"
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)

        chip_sub_idx = ["R01C01","R02C01","R03C01","R04C01","R05C01","R06C01","R01C02","R02C02","R03C02","R04C02","R05C02","R06C02"]
        
        results = []
        for experiment in self.experiments:
            result = {}
            info_result = {k: v for k, v in experiment.items() if k not in ["names", "name_numbers", "sample", "备注"]}
            for idx,sample_i, name, name_number, in zip(chip_sub_idx, experiment["sample"], experiment["names"], experiment["name_numbers"]):
                result.update(info_result)
                result["name"] = name
                result["name_number"] = name_number
                result["chip_sub_idx"] = idx
                result["sample"] = sample_i
                results.append(result.copy())
        
        df = pd.DataFrame(results)
        df.to_csv(output_path_obj, index=False, encoding='utf-8')
        
        logging.info(f"Saved experiments to {output_path}")


def setup_logging() -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(description='Process SNP experiment data.')
    parser.add_argument(
        '--data_path', 
        type=str, 
        default=r"D:\03.projects\AI.PGT\snparray_analysis\data\experiment_snparray.csv",
        required=False, 
        help='Path to the SNP experiment data text file.'
    )
    parser.add_argument(
        '--output_dir', 
        type=str, 
        required=False, 
        default=r"D:\03.projects\AI.PGT\snparray_analysis\data",
        help='Path to save the parsed JSON file.'
    )
    parser.add_argument(
        '--file_name', 
        type=str, 
        required=False, 
        default=r"experiment_snparray",
        help='file name for the output files (without extension).'
    )
    return parser.parse_args()


def main() -> None:
    """Main execution function."""
    setup_logging()
    args = parse_arguments()
    
    logging.info(f"Arguments: {args}")
    
    try:
        parser_snp = ParseSNPExperimentData(data_path=args.data_path)
        lines = parser_snp.load_data(data_path=args.data_path)
        parser_snp.parse_txt(lines)
        parser_snp.save_json(output_dir=args.output_dir, file_name=args.file_name)
        parser_snp.save_csv(output_dir=args.output_dir, file_name=args.file_name+"_persons")
        logging.info("SNP experiment data parsing completed successfully.")
    except Exception as e:
        logging.error(f"Error processing SNP experiment data: {e}")
        raise


if __name__ == "__main__":
    main()