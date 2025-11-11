"""
1. 根据p value选择SNPs features
2. 数据处理，one_hot编码, 标准化, 缺失值处理
3. 训练模型
4. 评估模型
"""
import pandas as pd
import numpy as np
import os
import torch
import logging

from sklearn.metrics import accuracy_score, roc_auc_score, average_precision_score, confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score, matthews_corrcoef
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier


# from ..data_process import logger

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SNPML:
    def __init__(self, output_dir:str):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            logger.info(f"create output dir: {output_dir}")
        self.output_dir = output_dir
        
        pass

    def get_snp_features(self, logistic_file:str):
        """根据p value选择SNPs features"""
        assert os.path.exists(logistic_file), f"{logistic_file} not exists"
        logistic_df = pd.read_csv(logistic_file, sep="\s+")
        assert "P" in logistic_df.columns, f"{logistic_file} not contains P column"
        # sort by p value
        logistic_df = logistic_df.sort_values(by="P")
        logistic_df = logistic_df[logistic_df["P"] < 0.05]
        logger.info(f"select {len(logistic_df)} snp features\n{logistic_df.head()}")

        # save into output_dir
        logistic_df.to_csv(os.path.join(self.output_dir, "snp_features_sort.csv"), index=False)
        return logistic_df
    
    def load_feature_phenotype_data(self, feature_file:str, phenotype_file:str, sep="\t"):
        """将vcf的数据进行处理, 行代表样本，列代表SNP，且将缺失值填充为-1，00为0, 01为1， 11为2"""
        assert os.path.exists(feature_file), f"{feature_file} not exists"
        assert os.path.exists(phenotype_file), f"{phenotype_file} not exists"
        feature_df = pd.read_csv(feature_file, sep=sep)
        # 将其转置，并这只main 列
        feature_df.columns = feature_df.columns.str.replace(r'^\d+_\d+_', '', regex=True)
        feature_df = feature_df.iloc[:, 9:].T
        genotype_map = {
            '0/0': 0,
            '0/1': 1,
            '1/0': 1,  # 杂合子的另一种表示
            '1/1': 2,
            './.': -1
        }
        feature_df = feature_df.applymap(lambda x: genotype_map.get(x, -1))

        phenotype_df = pd.read_csv(phenotype_file, sep=sep).iloc[:,1:]
        phenotype_df = phenotype_df.set_index(["IID"])
        phenotype_df = phenotype_df.add_prefix("label_")
        logger.info(f"feature_df shape: {feature_df.shape}\n{feature_df.head()}, \
                    phenotype_df shape: {phenotype_df.shape}\n{phenotype_df.head()}")
        
        merge_X_Y = pd.merge(feature_df, phenotype_df, left_index=True, right_index=True)

        logger.info(f"merge_X_Y shape: {merge_X_Y.shape}\n{merge_X_Y.head()}")
        merge_X_Y.to_csv(os.path.join(self.output_dir, "merge_X_Y.csv"), index=False, sep=sep)
        X, Y = merge_X_Y.iloc[:, :-1].to_numpy(), merge_X_Y.iloc[:, -1].to_numpy() - 1
        # Y = np.expand_dims(Y, axis=1)
        return X, Y
        

    def train_xgboost_classifier(self, X_train, y_train, X_test, random_state=42):
        """训练XGBoost分类器"""
        X_train_np = X_train.cpu().numpy() if torch.is_tensor(X_train) else X_train
        y_train_np = y_train.cpu().numpy() if torch.is_tensor(y_train) else y_train
        X_test_np = X_test.cpu().numpy() if torch.is_tensor(X_test) else X_test
        
        print("Training XGBoost classifier...")
        print("XGboost parameters: n_estimators=100, learning_rate=0.1, max_depth=6, random_state={}".format(random_state))
        xgb = XGBClassifier(
            n_estimators=100,
            random_state=random_state,
            # use_label_encoder=False, 
            eval_metric='mlogloss',
            learning_rate=0.1,
            max_depth=6        
        )
        xgb.fit(X_train_np, y_train_np)
        print("XGBoost training completed")
        
        probs = xgb.predict_proba(X_test_np)
        preds = xgb.predict(X_test_np)
        return preds, probs

    def evaluate(self, X_train, y_train, X_test, y_test):
        # 因为Xgboost不使用验证集，因此仅使用原始训练集和测试集
        print(f"\n[INFO] Since XGboost does not use validation set, using train and test sets without validation set.")
        print(f"Train set size: {len(X_train)}")
        print(f"Test set size: {len(X_test)}")

        # 使用训练集训练，在处理后的测试集上评估
        y_pred, y_probs = self.train_xgboost_classifier(
                X_train, y_train, X_test,
                random_state=42
            )
        y_true = y_test.cpu().numpy() if torch.is_tensor(y_test) else y_test


        # confusion maxtrix
        cm = confusion_matrix(y_true, y_pred)
        logger.info(f"Confusion Matrix:\n{cm}")
        # 计算评估指标
        acc = accuracy_score(y_true, y_pred)
        precision = precision_score(y_true, y_pred, average='macro', zero_division=0)
        recall = recall_score(y_true, y_pred, average='macro', zero_division=0)
        f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)

        try:
            mcc = matthews_corrcoef(y_true, y_pred)
        except:
            mcc = 0.0

        try:
            if y_probs.shape[1] == 2:
                auc = roc_auc_score(y_true, y_probs[:, 1])
                auprc = average_precision_score(y_true, y_probs[:, 1])
            else:
                auc = roc_auc_score(y_true, y_probs, multi_class='ovr', average='macro')
                auprc = average_precision_score(y_true, y_probs, average='macro')
        except Exception as e:
            print(f"Error calculating AUC/AUPRC: {e}")
            auc = 0.0
            auprc = 0.0

        logger.info(f"acc: {acc:.4f}, auc: {auc:.4f}, auprc: {auprc:.4f}, f1: {f1:.4f}, mcc: {mcc:.4f}, precision: {precision:.4f}, recall: {recall:.4f}")
        return acc, auc, auprc, f1, mcc, precision, recall


if __name__ == "__main__":
    logistic_file = r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS_ALL1\PLINK_311025_1032\1_QC\logistic_results.assoc_2.logistic"
    output_dir = r"D:\03.projects\AI.PGT\snparray_analysis\work_dir\snp_ml"
    feature_fiel = r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS_ALL1\PLINK_311025_1032\4_vcf\PGT_TLS_ALL1.feature500.table"
    phenotype_file = r"D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS_ALL1\PLINK_311025_1032\4_vcf\PGT_TLS_ALL1.refactor.phenotype"

    snpml = SNPML(output_dir)
    snpml.get_snp_features(logistic_file)
    X, Y = snpml.load_feature_phenotype_data(feature_fiel, phenotype_file)
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    snpml.train_xgboost_classifier(X_train, y_train, X_test, random_state=42)
    snpml.evaluate(X_train, y_train, X_test, y_test)
    pass