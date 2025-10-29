import pandas as pd
import logging
from pathlib import Path
from pypinyin import lazy_pinyin

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

ped_file=r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.ped"
phenotype_file=r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.phenotype"
full_phenotype_file=r"D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype.csv"

output_ped = r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.refactor.ped"
output_phenotype=r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.refactor.phenotype"
output_full_phenotype=r"D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype.refactor.csv"

class RefactorPhenotypeData:
    def __init__(self,ped_file,phenotype_file,full_phenotype_file,rename={}):
        self.ped_file=ped_file
        self.phenotype_file=phenotype_file
        self.full_phenotype_file=full_phenotype_file
        
        self.df_sample_idx = self.get_sample_idx(self.ped_file)
        self.df_full_phenotype = self.get_full_phenotype_by_sample(self.full_phenotype_file, rename=rename)
    
    def get_sample_idx(self, ped_file:str):
        """将ped文件分为两个文件，一个包含样本ID，一个包含样本信息
        家系ID 个体ID 父ID 母ID 性别 表型值 SNP1_allele1 SNP1_allele2 SNP2_allele1 SNP2_allele2 ...
        ——>
        家系ID 个体ID 父ID 母ID 性别 表型值 |
        SNP1_allele1 SNP1_allele2 SNP2_allele1 SNP2_allele2 ...
        """
        with open(ped_file, 'r') as f:
            lines = f.readlines()
            sample_ids = [line[:100].split()[:6] for line in lines]
        df_sample = pd.DataFrame(sample_ids, columns=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])
        logging.info(f'sample_ids: {df_sample.head()}')
        df_sample_idx = df_sample["IID"]
        return df_sample_idx
    
    def get_full_phenotype_by_sample(self, full_phenotype_file: str, rename={}):
        """根据样本ID获取全量phenotype数据"""
        df_full_phenotype = pd.read_csv(full_phenotype_file)
        df_full_phenotype['IID'] = df_full_phenotype.apply(
            lambda row: f"{row['idx']}_{row['chip_sub_idx']}", axis=1
        )
        # full_phenotype = full_phenotype[full_phenotype['sample_id'].isin(sample_ids)]
        logging.info(f'full_phenotype: num: {df_full_phenotype.shape}, head:\n{df_full_phenotype.head()}')

        # 根据IID进行向左合并
        df_full_phenotype = pd.merge(self.df_sample_idx, df_full_phenotype, on='IID', how='left')

        check_rename = {key:value for key, value in rename.items() if value in df_full_phenotype.columns}
        if check_rename:
            df_full_phenotype.rename(columns=rename, inplace=True)
        if "SEX" not in df_full_phenotype.columns:
            df_full_phenotype["SEX"] = 0
            logger.info(f"add column SEX to full_phenotype")

        if "是" in df_full_phenotype['PHENO'].unique():
            df_full_phenotype['PHENO'] = df_full_phenotype['PHENO'].map({'否': 0, '是': 1})
            logger.info(f"map column PHENO to full_phenotype")

        for col in df_full_phenotype.columns:
            if df_full_phenotype[col].apply(self.contains_chinese).any():
                df_full_phenotype[col] = df_full_phenotype[col].apply(
                    lambda x: '_'.join(lazy_pinyin(str(x))) if pd.notna(x) else x
                )
                logger.info(f"rename column {col} to full_phenotype to pinyin")

        logging.info(f'full_phenotype: num: {df_full_phenotype.shape}, head:\n{df_full_phenotype.head()}')
        return df_full_phenotype
            
    def refoctor_ped_file(self, ped_file:str, output_file:str):
        """根据df_full_phenotype进行重构ped文件"""
        with open(ped_file, 'r') as f:
            lines = f.readlines()
        ped_phenotype = self.df_full_phenotype[["MID", "IID", "PID", "MID", "SEX", "PHENO"]]
        
        ped_phenotype = ped_phenotype

        with open(output_file, 'w') as f:
                pass # delete the file content
        for line, sample in zip(lines, ped_phenotype.to_numpy()):
            snp_txt =  "\t".join(line[:100].split()[6:]) + line[100:]
            full_sample_txt = "\t".join(str(i) for i in sample)
            new_line = f"{full_sample_txt}\t{snp_txt}"
            with open(output_file, 'a') as f:
                f.write(new_line)
        
    @staticmethod
    def contains_chinese(text):
        if pd.isna(text):
            return False
        return any('\u4e00' <= char <= '\u9fff' for char in str(text))  

if __name__ == "__main__":
    rename = {
        "男方姓名":"PID",
        "女方姓名_clinical":"MID",
        "诊断可移植": "PHENO"
    }
    rpd = RefactorPhenotypeData(ped_file,phenotype_file,full_phenotype_file, rename)
    rpd.refoctor_ped_file(ped_file,output_ped)