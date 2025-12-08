#!/usr/bin/env python3
"""
SNP-based Phenotype Prediction using XGBoost
Author: ML Assistant
Date: 2024
Description: Predict phenotype from VCF genotype data using XGBoost classifier
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import gzip
import warnings
from typing import Tuple, Dict, List, Optional, Any
from pathlib import Path
from collections import defaultdict

# Machine learning libraries
from tqdm import tqdm
import xgboost as xgb
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.feature_selection import (
    SelectKBest,
    mutual_info_classif,
    VarianceThreshold,
)
from sklearn.decomposition import PCA
import joblib

# For VCF parsing
import allel

warnings.filterwarnings("ignore")


class VCFProcessor:
    """Process VCF files to extract genotype data"""

    def __init__(self, vcf_path: str, n_workers: int = 4):
        """
        Initialize VCF processor

        Args:
            vcf_path: Path to VCF file (.vcf or .vcf.gz)
            n_workers: Number of workers for parallel processing
        """
        self.vcf_path = vcf_path
        self.n_workers = n_workers
        self.vcf_data = None
        self.sample_ids = None

    def load_vcf(
        self, max_snps: Optional[int] = None, quality_threshold: float = 20.0
    ) -> pd.DataFrame:
        """
        Load and filter VCF file

        Args:
            max_snps: Maximum number of SNPs to load (for testing)
            quality_threshold: Minimum QUAL score to keep variant

        Returns:
            DataFrame with genotype matrix
        """
        print(f"Loading VCF file: {self.vcf_path}")

        # try:
        # Use scikit-allel for efficient VCF reading
        vcf = allel.read_vcf(
            self.vcf_path,
            fields=[
                "samples",
                "variants/CHROM",
                "variants/POS",
                "variants/QUAL",
                "variants/ID",
                "calldata/GT",
            ],
            alt_number=16,
        )

        # Extract sample IDs
        print(vcf.keys())
        self.sample_ids = vcf["samples"]
        print(f"Found {len(self.sample_ids)} samples")

        # Extract genotype data
        gt = vcf["calldata/GT"]
        n_snps, n_samples, ploidy = gt.shape

        print(f"Loaded {n_snps} SNPs for {n_samples} samples (ploidy={ploidy})")

        # Convert genotype to numeric format
        # For diploid: 0=0/0, 1=0/1 or 1/0, 2=1/1, -1=./.
        genotype_matrix = np.sum(gt, axis=2)  # Sum alleles
        genotype_matrix[gt[:, :, 0] == -1] = -1  # Handle missing

        # Create SNP info DataFrame
        snp_info = pd.DataFrame(
            {
                "CHROM": vcf["variants/CHROM"],
                "POS": vcf["variants/POS"],
                "ID": vcf["variants/ID"],
                "QUAL": vcf["variants/QUAL"],
            }
        )

        # Filter by quality if provided
        if quality_threshold > 0:
            qual_mask = snp_info["QUAL"] >= quality_threshold
            genotype_matrix = genotype_matrix[qual_mask]
            snp_info = snp_info[qual_mask]
            print(f"Filtered to {len(snp_info)} SNPs with QUAL >= {quality_threshold}")

        # Limit SNPs if specified
        if max_snps and len(snp_info) > max_snps:
            genotype_matrix = genotype_matrix[:max_snps]
            snp_info = snp_info.iloc[:max_snps]
            print(f"Limited to {max_snps} SNPs for testing")

        # Create final DataFrame
        genotype_df = pd.DataFrame(
            genotype_matrix.T,
            columns=[f"SNP_{i}" for i in range(len(snp_info))],
            index=self.sample_ids,
        )

        self.vcf_data = {
            "genotypes": genotype_df,
            "snp_info": snp_info,
            "sample_ids": self.sample_ids,
        }

        return genotype_df

        # except Exception as e:
        #     print(f"Error loading VCF: {e}")
        #     sys.exit(1)

    def filter_snps(
        self,
        genotype_df: pd.DataFrame,
        missing_threshold: float = 0.1,
        maf_threshold: float = 0.01,
    ) -> pd.DataFrame:
        """
        Filter SNPs based on missing rate and minor allele frequency

        Args:
            genotype_df: Input genotype DataFrame
            missing_threshold: Maximum allowed missing rate
            maf_threshold: Minimum minor allele frequency

        Returns:
            Filtered genotype DataFrame
        """
        print("\nFiltering SNPs...")
        n_original = genotype_df.shape[1]

        # Calculate missing rate
        missing_rate = (genotype_df == -1).sum() / len(genotype_df)
        mask_missing = missing_rate <= missing_threshold
        df_filtered = genotype_df.loc[:, mask_missing]
        print(
            f"Removed {n_original - df_filtered.shape[1]} SNPs with missing rate > {missing_threshold}"
        )

        # Calculate MAF
        mafs = []
        for col in tqdm(df_filtered.columns, desc="Calculating MAF"):
            # Convert -1 to NaN for MAF calculation
            col_data = df_filtered[col].replace(-1, np.nan)
            # Calculate allele frequencies (0, 1, or 2)
            allele_counts = col_data.value_counts(normalize=True)
            if len(allele_counts) > 0:
                # For biallelic SNP, MAF = min(p, 1-p) where p is frequency of minor allele
                # Since we have genotype counts, we need allele frequency
                # Simplified: use minimum genotype frequency
                maf = allele_counts.min()
                mafs.append(maf)
            else:
                mafs.append(0)

        maf_series = pd.Series(mafs, index=df_filtered.columns)
        mask_maf = maf_series >= maf_threshold
        df_filtered = df_filtered.loc[:, mask_maf]

        print(
            f"Removed {mask_missing.sum() - df_filtered.shape[1]} SNPs with MAF < {maf_threshold}"
        )
        print(f"Final: {df_filtered.shape[1]} SNPs after filtering")

        return df_filtered


class PhenoProcessor:
    """Process PED phenotype files"""

    def __init__(self, ped_path: str):
        """
        Initialize PED processor

        Args:
            ped_path: Path to PED file
        """
        self.ped_path = ped_path
        self.pheno_data = None

    def load_ped(self) -> pd.DataFrame:
        """
        Load PED file

        Returns:
            DataFrame with phenotype data
        """
        print(f"\nLoading PED file: {self.ped_path}")

        # PED format: Family ID, Individual ID, Paternal ID, Maternal ID, Sex, Phenotype
        try:
            ped_df = pd.read_csv(
                self.ped_path,
                sep="\t",
                header=1,
                names=[
                    "Family ID",
                    "Individual ID",
                    "Paternal ID",
                    "Maternal ID",
                    "Gender",
                    "Phenotype",
                    "Population",
                    "Relationship",
                    "Siblings",
                    "Second Order",
                    "Third Order",
                    "Children",
                    "Other Comments",
                    "phase 3 genotypes",
                    "related genotypes",
                    "omni genotypes",
                    "affy_genotypes",
                ],
            )
            name_replace = {
                "Family ID": "FID",
                "Individual ID": "IID",
                "Paternal ID": "PID",
                "Maternal ID": "MID",
                "Gender": "SEX",
                "Population": "PHENO",
            }

            super_pop_mapping = {
                # 非洲 (AFR)
                "YRI": 0,  # "AFR",
                "LWK": 0,  # "AFR",
                "GWD": 0,  # "AFR",
                "MSL": 0,  # "AFR",
                "ESN": 0,  # "AFR",
                "ACB": 0,  # "AFR",
                # 欧洲 (EUR)
                "GBR": 1,  # "EUR",
                "FIN": 1,  # "EUR",
                "IBS": 1,  # "EUR",
                "CEU": 1,  # "EUR",
                "TSI": 1,  # "EUR",
                # 东亚 (EAS)
                "CHS": 2,  # "EAS",
                "CHB": 2,  # "EAS",
                "JPT": 2,  # "EAS",
                "CDX": 2,  # "EAS",
                "CHD": 2,  # "EAS",
                "KHV": 2,  # "EAS",
                # 南亚 (SAS)
                "GIH": 3,  # "SAS",
                "PJL": 3,  # "SAS",
                "BEB": 3,  # "SAS",
                "ITU": 3,  # "SAS",
                "STU": 3,  # "SAS",
                # 美洲 (AMR)
                "MXL": 4,  # "AMR",
                "PUR": 4,  # "AMR",
                "CLM": 4,  # "AMR",
                "PEL": 4,  # "AMR",
                # 混合群体 - 根据遗传背景分类
                "ASW": 5,  # "UNK",  # 非裔美国人，归为AFR
            }
            ped_df["Population"] = ped_df["Population"].map(super_pop_mapping)

            ped_df = ped_df.rename(columns=name_replace)
            required_columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
            ped_df = ped_df[required_columns]

            print(f"Loaded {len(ped_df)} individuals")
            print(f"head:\n{ped_df.head()}")
            print(
                f"Phenotype distribution:\n{ped_df['PHENO'].value_counts().sort_index()}"
            )

            self.pheno_data = ped_df
            return ped_df

        except Exception as e:
            print(f"Error loading PED file: {e}")
            sys.exit(1)

    def match_samples(
        self, ped_df: pd.DataFrame, vcf_samples: List[str]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Match samples between PED and VCF data

        Args:
            ped_df: PED DataFrame
            vcf_samples: List of sample IDs from VCF

        Returns:
            Tuple of (matched phenotype series, matched genotype DataFrame indices)
        """
        # Try different matching strategies
        matched_samples = []
        missing_samples = []

        for vcf_sample in vcf_samples:
            # Check if sample exists in PED
            # PED typically uses IID or FID_IID format
            mask = ped_df["IID"] == vcf_sample
            if mask.any():
                matched_samples.append(vcf_sample)
            else:
                missing_samples.append(vcf_sample)

        print(f"\nMatched {len(matched_samples)}/{len(vcf_samples)} samples")
        if missing_samples and len(missing_samples) < 10:
            print(f"Missing samples: {missing_samples}")

        # Filter phenotypes for matched samples
        pheno_matched = ped_df.set_index("IID").loc[matched_samples, "PHENO"]

        return pheno_matched, matched_samples


class FeatureSelector:
    """Perform feature selection on SNP data"""

    def __init__(self, n_features: int = 10000, method: str = "mutual_info"):
        """
        Initialize feature selector

        Args:
            n_features: Number of features to select
            method: Selection method ('mutual_info', 'variance', 'pca')
        """
        self.n_features = n_features
        self.method = method
        self.selector = None
        self.selected_features = None

    def fit_transform(self, X: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
        """
        Fit and transform features

        Args:
            X: Input features
            y: Target labels

        Returns:
            Selected features DataFrame
        """
        print(f"\nPerforming feature selection ({self.method})...")
        print(f"Input shape: {X.shape}")

        # Handle missing values
        X_filled = X.replace(-1, np.nan).fillna(X.replace(-1, np.nan).median())

        if self.method == "variance":
            # Remove low variance features
            selector = VarianceThreshold(threshold=0.01)
            X_selected = selector.fit_transform(X_filled)
            self.selected_features = X.columns[selector.get_support()]
            print(f"Selected {X_selected.shape[1]} features by variance")

        elif self.method == "mutual_info":
            # Use mutual information for feature selection
            if X_filled.shape[1] > self.n_features:
                # For large datasets, use iterative selection
                selector = SelectKBest(
                    mutual_info_classif, k=min(self.n_features, X_filled.shape[1])
                )
                X_selected = selector.fit_transform(X_filled, y)
                self.selected_features = X.columns[selector.get_support()]
                print(f"Selected {X_selected.shape[1]} features by mutual information")
            else:
                X_selected = X_filled
                self.selected_features = X.columns

        elif self.method == "pca":
            # Dimensionality reduction with PCA
            pca = PCA(n_components=min(100, X_filled.shape[1]))
            X_selected = pca.fit_transform(X_filled)
            print(
                f"PCA: {X_selected.shape[1]} components explain "
                f"{pca.explained_variance_ratio_.sum():.2%} variance"
            )
            self.selected_features = [f"PC{i}" for i in range(X_selected.shape[1])]

        else:
            X_selected = X_filled
            self.selected_features = X.columns

        # Convert back to DataFrame
        X_selected_df = pd.DataFrame(
            X_selected, index=X.index, columns=self.selected_features
        )

        return X_selected_df


class SNPClassifier:
    """Main classifier for SNP-based phenotype prediction"""

    def __init__(self, params: Optional[Dict] = None):
        """
        Initialize classifier

        Args:
            params: XGBoost parameters
        """
        self.params = params or {
            "objective": "multi:softprob",
            "num_class": 5,
            "max_depth": 6,
            "learning_rate": 0.1,
            "n_estimators": 100,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "random_state": 42,
            "n_jobs": -1,
            "verbosity": 0,
        }
        self.model = None
        self.feature_selector = None
        self.scaler = None
        self.label_encoder = None

    def prepare_data(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        test_size: float = 0.2,
        random_state: int = 42,
    ) -> Tuple:
        """
        Prepare data for training

        Args:
            X: Features
            y: Labels
            test_size: Test set proportion
            random_state: Random seed

        Returns:
            Tuple of (X_train, X_test, y_train, y_test)
        """
        # Encode labels if not numeric
        if y.dtype == "object":
            self.label_encoder = LabelEncoder()
            y_encoded = self.label_encoder.fit_transform(y)
        else:
            y_encoded = y.values

        # Scale features
        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(X.replace(-1, np.nan).fillna(0))

        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled,
            y_encoded,
            test_size=test_size,
            random_state=random_state,
            stratify=y_encoded,
        )

        print(f"Training set: {X_train.shape}, Test set: {X_test.shape}")

        return X_train, X_test, y_train, y_test

    def train(
        self,
        X_train: np.ndarray,
        y_train: np.ndarray,
        X_val: Optional[np.ndarray] = None,
        y_val: Optional[np.ndarray] = None,
        cv_folds: int = 5,
    ) -> Dict:
        """
        Train XGBoost classifier

        Args:
            X_train: Training features
            y_train: Training labels
            X_val: Validation features (optional)
            y_val: Validation labels (optional)
            cv_folds: Number of CV folds

        Returns:
            Dictionary with training results
        """
        print("\nTraining XGBoost classifier...")

        # Update num_class based on data
        n_classes = len(np.unique(y_train))
        if "num_class" in self.params and self.params["num_class"] != n_classes:
            print(f"Updating num_class from {self.params['num_class']} to {n_classes}")
            self.params["num_class"] = n_classes

        # Initialize model
        self.model = xgb.XGBClassifier(**self.params)

        # Train with or without validation set
        if X_val is not None and y_val is not None:
            eval_set = [(X_train, y_train), (X_val, y_val)]
            self.model.fit(
                X_train,
                y_train,
                eval_set=eval_set,
                verbose=False,
                early_stopping_rounds=10,
            )
        else:
            # Use cross-validation
            cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
            cv_scores = cross_val_score(
                self.model, X_train, y_train, cv=cv, scoring="accuracy", n_jobs=-1
            )
            print(f"CV Accuracy: {cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")

            # Train on full training set
            self.model.fit(X_train, y_train)

        # Get feature importance
        if hasattr(self.model, "feature_importances_"):
            importance = self.model.feature_importances_
            top_features = np.argsort(importance)[-10:]  # Top 10
            print(f"\nTop 10 important features:")
            for i, idx in enumerate(top_features[::-1]):
                print(f"  {i+1}. Feature {idx}: {importance[idx]:.4f}")

        return {"model": self.model}

    def evaluate(self, X_test: np.ndarray, y_test: np.ndarray) -> Dict:
        """
        Evaluate model performance

        Args:
            X_test: Test features
            y_test: Test labels

        Returns:
            Dictionary with evaluation metrics
        """
        print("\nEvaluating model...")

        # Make predictions
        y_pred = self.model.predict(X_test)
        y_pred_proba = self.model.predict_proba(X_test)

        # Calculate metrics
        accuracy = accuracy_score(y_test, y_pred)

        # Decode labels if encoder was used
        if self.label_encoder is not None:
            y_test_decoded = self.label_encoder.inverse_transform(y_test)
            y_pred_decoded = self.label_encoder.inverse_transform(y_pred)
        else:
            y_test_decoded = y_test
            y_pred_decoded = y_pred

        print(f"Test Accuracy: {accuracy:.4f}")
        print("\nClassification Report:")
        print(classification_report(y_test_decoded, y_pred_decoded))

        print("\nConfusion Matrix:")
        cm = confusion_matrix(y_test_decoded, y_pred_decoded)
        print(cm)

        return {
            "accuracy": accuracy,
            "predictions": y_pred,
            "probabilities": y_pred_proba,
            "confusion_matrix": cm,
        }

    def save_model(self, path: str):
        """Save model to disk"""
        joblib.dump(
            {
                "model": self.model,
                "scaler": self.scaler,
                "label_encoder": self.label_encoder,
                "params": self.params,
            },
            path,
        )
        print(f"Model saved to {path}")

    def load_model(self, path: str):
        """Load model from disk"""
        data = joblib.load(path)
        self.model = data["model"]
        self.scaler = data["scaler"]
        self.label_encoder = data["label_encoder"]
        self.params = data["params"]
        print(f"Model loaded from {path}")


class SNPPipeline:
    """Main pipeline for SNP phenotype prediction"""

    def __init__(self, args):
        """
        Initialize pipeline

        Args:
            args: Command line arguments
        """
        self.args = args
        self.vcf_processor = None
        self.pheno_processor = None
        self.feature_selector = None
        self.classifier = None
        self.results = {}

    def run(self):
        """Run complete pipeline"""
        print("=" * 60)
        print("SNP Phenotype Prediction Pipeline")
        print("=" * 60)

        # Step 1: Process VCF
        self.vcf_processor = VCFProcessor(self.args.vcf, self.args.n_workers)
        genotype_df = self.vcf_processor.load_vcf(
            max_snps=self.args.max_snps, quality_threshold=self.args.quality_threshold
        )

        # Filter SNPs
        genotype_filtered = self.vcf_processor.filter_snps(
            genotype_df,
            missing_threshold=self.args.missing_threshold,
            maf_threshold=self.args.maf_threshold,
        )

        # Step 2: Process phenotype
        self.pheno_processor = PhenoProcessor(self.args.pheno)
        pheno_df = self.pheno_processor.load_ped()

        # Match samples
        pheno_matched, matched_samples = self.pheno_processor.match_samples(
            pheno_df, genotype_filtered.index.tolist()
        )

        # Filter genotypes to matched samples
        X = genotype_filtered.loc[matched_samples]
        y = pheno_matched

        # Step 3: Feature selection
        self.feature_selector = FeatureSelector(
            n_features=self.args.n_features, method=self.args.selection_method
        )
        X_selected = self.feature_selector.fit_transform(X, y)

        # Step 4: Classification
        classifier_params = {
            "max_depth": self.args.max_depth,
            "learning_rate": self.args.learning_rate,
            "n_estimators": self.args.n_estimators,
            "objective": "multi:softprob",
            "num_class": 5,
            "subsample": self.args.subsample,
            "colsample_bytree": self.args.colsample_bytree,
            "random_state": self.args.random_seed,
            "n_jobs": self.args.n_workers,
            "verbosity": 0,
        }

        self.classifier = SNPClassifier(classifier_params)

        # Prepare data
        X_train, X_test, y_train, y_test = self.classifier.prepare_data(
            X_selected,
            y,
            test_size=self.args.test_size,
            random_state=self.args.random_seed,
        )

        # Train model
        train_results = self.classifier.train(
            X_train, y_train, cv_folds=self.args.cv_folds
        )

        # Evaluate model
        eval_results = self.classifier.evaluate(X_test, y_test)

        # Save results
        self.results = {
            "accuracy": eval_results["accuracy"],
            "n_samples": len(X),
            "n_features": X_selected.shape[1],
            "feature_names": self.feature_selector.selected_features.tolist(),
        }

        # Save model if requested
        if self.args.save_model:
            model_path = f"xgboost_snp_model_{self.args.random_seed}.joblib"
            self.classifier.save_model(model_path)

        # Save feature importance
        if self.args.save_features:
            features_df = pd.DataFrame(
                {
                    "feature": self.feature_selector.selected_features,
                    "importance": self.classifier.model.feature_importances_,
                }
            )
            features_df = features_df.sort_values("importance", ascending=False)
            features_df.to_csv(
                f"feature_importance_{self.args.random_seed}.csv", index=False
            )
            print(
                f"Feature importance saved to feature_importance_{self.args.random_seed}.csv"
            )

        print("\n" + "=" * 60)
        print(f"FINAL ACCURACY: {eval_results['accuracy']:.4f}")
        print("=" * 60)

        return eval_results["accuracy"]


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description="SNP-based Phenotype Prediction using XGBoost",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input/Output arguments
    parser.add_argument("--vcf", required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument(
        "--pheno", required=True, help="Input phenotype file (.ped format)"
    )

    # VCF processing arguments
    parser.add_argument(
        "--max-snps",
        type=int,
        default=None,
        help="Maximum number of SNPs to load (for testing)",
    )
    parser.add_argument(
        "--quality-threshold",
        type=float,
        default=20.0,
        help="Minimum QUAL score for SNP filtering",
    )
    parser.add_argument(
        "--missing-threshold",
        type=float,
        default=0.1,
        help="Maximum missing rate for SNP filtering",
    )
    parser.add_argument(
        "--maf-threshold",
        type=float,
        default=0.01,
        help="Minimum minor allele frequency for SNP filtering",
    )

    # Feature selection arguments
    parser.add_argument(
        "--n-features", type=int, default=10000, help="Number of features to select"
    )
    parser.add_argument(
        "--selection-method",
        choices=["mutual_info", "variance", "pca"],
        default="mutual_info",
        help="Feature selection method",
    )

    # Model arguments
    parser.add_argument("--max-depth", type=int, default=6, help="XGBoost max depth")
    parser.add_argument(
        "--learning-rate", type=float, default=0.1, help="XGBoost learning rate"
    )
    parser.add_argument(
        "--n-estimators", type=int, default=100, help="XGBoost number of estimators"
    )
    parser.add_argument(
        "--subsample", type=float, default=0.8, help="XGBoost subsample ratio"
    )
    parser.add_argument(
        "--colsample-bytree",
        type=float,
        default=0.8,
        help="XGBoost column sample ratio",
    )

    # Training arguments
    parser.add_argument(
        "--test-size", type=float, default=0.2, help="Test set proportion"
    )
    parser.add_argument(
        "--cv-folds", type=int, default=5, help="Number of cross-validation folds"
    )
    parser.add_argument(
        "--random-seed", type=int, default=42, help="Random seed for reproducibility"
    )

    # Output arguments
    parser.add_argument(
        "--save-model", action="store_true", help="Save trained model to disk"
    )
    parser.add_argument(
        "--save-features", action="store_true", help="Save feature importance to CSV"
    )

    # Performance arguments
    parser.add_argument(
        "--n-workers", type=int, default=4, help="Number of parallel workers"
    )

    args = parser.parse_args()

    # Check input files exist
    for file_path in [args.vcf, args.pheno]:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}")
            sys.exit(1)

    # Run pipeline
    pipeline = SNPPipeline(args)
    accuracy = pipeline.run()

    # Exit with code based on accuracy (for scripting purposes)
    if accuracy > 0.5:
        sys.exit(0)  # Success
    else:
        sys.exit(1)  # Poor performance


if __name__ == "__main__":
    main()
