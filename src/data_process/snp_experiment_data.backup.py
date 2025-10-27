import argparse
import json
from pathlib import Path
import re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ChipInfoParser:
    def __init__(self):
        self.experiments:dict = None
        self.chip_idx_dict = {"chip_num_1":"X1 number", "chip_num_2":"X2 number", "chip_num_4":"X3 number", "chip_num_5":"X4 number"}
        self.chip_info = {}
        self.name_dict ={}

    def parse_name(self, text):
        pattern = r'^([\u4e00-\u9fa5A-Z]+)\s*([A-Z])\s*(\d{1,2})$'
        match = re.match(pattern, text)
        
        if match:
            chinese = match.group(1)
            initial = match.group(2) if match.group(2) else ''
            number = match.group(3)
            
            # 验证数字范围在1-20之间
            if 1 <= int(number) <= 20:
                return chinese, initial, number
    
        return None

    def parse_number(self, text):
        pattern = r'^([A-Z])(\d{1,2})$'
        match = re.match(pattern, text)
        
        if match:
            initial = match.group(1)
            number = match.group(2) if match.group(2) else ''
            
            # 验证数字范围在1-20之间
            if 2 <= int(number) <= 20:
                return initial, number
        
        return None

    def parse(self, experiments:dict):
        chip_infos = []
        name_dict ={}
        experiment_info = {}
        for key, value in experiments.items():
            chip_info = {}
            if key  in self.chip_idx_dict.keys():
                chip_info["type"] = value[0]
                chip_info["sample"] = value[1:-1]

                names = []
                name_numbers = []
                for name in chip_info["sample"]:
                    if "-" in name:
                        person_name = name.split("-")[1]
                        id_prefix = name.split("-")[0]
                        names.append(person_name)
                        name_numbers.append(id_prefix)
                    elif self.parse_name(name):
                        chinese, initial, number = self.parse_name(name)
                        name_dict[initial] = chinese
                        names.append(chinese)
                        name_numbers.append(number)
                    elif self.parse_number(name):
                        initial, number = self.parse_number(name)
                        names.append(name_dict[initial])
                        name_numbers.append(number)
                    elif "RB" == name:
                        person_name = "RB"
                        id_prefix = "1"
                        names.append(person_name)
                        name_numbers.append(id_prefix)
                    else:
                        names.append(name)
                        name_numbers.append(name)
                        logger.error(f"无法解析的名字: {name}")

                    if person_name not in name_dict:
                        name_dict[person_name] = []

                    name_dict[person_name].append(id_prefix)
                
                chip_info["names"] = names
                chip_info["name_numbers"] = name_numbers

                if self.chip_idx_dict[key] in experiments:
                    chip_info["idx"] = experiments[self.chip_idx_dict[key]]
                else:
                    raise ValueError(f"缺少芯片编号信息: {self.chip_idx_dict[key]}")
                
                if len(chip_info["names"]) == len(chip_info["name_numbers"]) == 12:
                    chip_infos.append(chip_info)
                else:
                    raise ValueError(f"样本数量不匹配 or !=12: names({len(chip_info['names'])}), name_numbers({len(chip_info['name_numbers'])})")
            
            elif key in self.chip_idx_dict.values():
                continue
            else:
                experiment_info[key] = value

        for key, value in experiment_info.items():
            for chip_info in chip_infos:
                chip_info[key] = value

        return chip_infos
    

class ParseSNPExperimentData:
    def __init__(self, data_path:str):
        self.data_path = Path(data_path)
        assert self.data_path.exists(), f"Data path {data_path} does not exist."

        self.experiments = []
        self.chip_parser = ChipInfoParser()

    def load_data(self, data_path:str):
        with open(data_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        return lines
    
    def parse_txt(self, lines:str):
        experiment_dict = {}
        chip_line = -1
        for line in lines:
            # get experiment time
            line_data = [i.strip() for i in line.strip().split(",")]
            if "实验时间" in line:
                chip_line = 0
                if experiment_dict:
                    parser_experiment = self.chip_parser.parse(experiment_dict)
                    self.experiments.extend(parser_experiment)
                    experiment_dict = {}

                time_key, time_value = line_data[1].split("：")
                experiment_dict[time_key.strip()] = time_value.strip().replace(" ", "")
                experimenter_key, experimenter_value = line_data[7].split("：")
                experiment_dict[experimenter_key.strip()] = experimenter_value.strip()
            elif "备注" in line:
                comment_key, comment_value = line_data[0][:2].strip(), line_data[0][3:].strip()
                experiment_dict[comment_key.strip()] = comment_value.strip()
                chip_line += 1

            if chip_line in [0, 3]:
                chip_line += 1
                continue

            if chip_line >= 1 and chip_line <= 5:
                experiment_dict[f"chip_num_{chip_line}"] = line_data
                chip_line += 1
            elif chip_line == 6:
                X1_key, X1_value = line_data[1].split("：")
                experiment_dict[X1_key.strip()] = X1_value.strip().replace(" ", "")
                X2_key, X2_value = line_data[7].split("：")
                experiment_dict[X2_key.strip()] = X2_value.strip()
                chip_line += 1
            elif chip_line == 7:
                X3_key, X3_value = line_data[1].split("：")
                experiment_dict[X3_key.strip()] = X3_value.strip().replace(" ", "")
                X4_key, X4_value = line_data[7].split("：")
                experiment_dict[X4_key.strip()] = X4_value.strip()
                chip_line += 1


        logger.info(f"Parsed total {len(self.experiments)} experiments from {self.data_path}")
        

    def save_json(self, output_path:str):
        with open(output_path, 'w', encoding='utf-8') as file:
            json.dump(self.experiments, file, ensure_ascii=False, indent=4)
        logger.info(f"Saved experiments to {output_path}")


def parse():
    parser = argparse.ArgumentParser(description='Process SNP experiment data.')
    parser.add_argument('--data_path', type=str, required=True, help='Path to the SNP experiment data text file.')
    parser.add_argument('--output_path', type=str, required=True, help='Path to save the parsed JSON file.')
    args = parser.parse_args()

    logger.info(f"Arguments: {args}")
    return args

def main():
    args = parse()
    parser_snp = ParseSNPExperimentData(data_path=args.data_path)
    lines = parser_snp.load_data(data_path=args.data_path)
    parser_snp.parse_txt(lines)
    parser_snp.save_json(output_path=args.output_path)

if __name__ == "__main__":
    main()
