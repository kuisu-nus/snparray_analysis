#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A professional, class-based Python script to import CSV data into a MySQL database.

This script uses pandas for efficient CSV reading and SQLAlchemy for robust
database writing, with special handling for long text fields and UTF-8 encoding.
"""

import argparse
import logging
import sys
from typing import Optional, Dict

import pandas as pd
from sqlalchemy import create_engine, engine, exc
from sqlalchemy.types import TEXT



def setup_logging() -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
# --- Configuration ---

# TODO: For production, load these from environment variables or a secure
# configuration file (e.g., using os.getenv, configparser, or python-dotenv).
DB_CONFIG = {
    "host": "localhost",       # Your MySQL server host
    "port": 3306,              # Your MySQL server port
    "user": "root",   # Your MySQL username
    "password": "sukui1016", # Your MySQL password
    "database": "ivf_data"  # The database to write to
}

# CSV file path and the target table name in MySQL
def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(description='add clinical data to mysql')
    parser.add_argument(
        '--csv_path', 
        type=str, 
        default=r"D:\03.projects\AI.PGT\snparray_analysis\data\5.芯片实验记录表2020.12.29.xlsx",
        required=False, 
        help='Path to the SNP experiment data text file.'
    )


CSV_FILE_PATH = r"D:\03.projects\AI.PGT\snparray_analysis\data\clinical_data_pgt.csv"
TABLE_NAME = "clinical_data"

# --- Logger Setup ---

# Configure a in-memory logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)  # Output logs to the console
        # You could also add a logging.FileHandler("import.log") here
    ]
)
logger = logging.getLogger(__name__)


class MySQLDataImporter:
    """
    Encapsulates the logic for importing a CSV file into a MySQL database table.
    """

    def __init__(self, db_config: Dict[str, any]):
        """
        Initializes the importer and creates the database engine.

        Args:
            db_config (Dict[str, any]): Database connection parameters.
        """
        self.db_config = db_config
        self.engine = self._create_db_engine()

    def _create_db_engine(self) -> Optional[engine.Engine]:
        """
        Creates a SQLAlchemy engine using the provided configuration.

        Uses 'utf8mb4' charset to ensure full Unicode support for
        Chinese characters and emojis.
        """
        # try:
        # Construct the database connection string
        # We use 'mysql+mysqlconnector' as the dialect
        connection_str = (
            f"mysql+mysqlconnector://{self.db_config['user']}:{self.db_config['password']}@"
            f"{self.db_config['host']}:{self.db_config['port']}/{self.db_config['database']}"
            f"?charset=utf8mb4"
        )
        
        engine = create_engine(connection_str)
        
        # Test the connection immediately
        with engine.connect() as conn:
            logger.info(
                f"Successfully connected to database: "
                f"{self.db_config['database']}@{self.db_config['host']}"
            )
        
        return engine
        
        # except exc.OperationalError as e:
        #     logger.error(f"Database connection failed: {e}")
        #     logger.error(
        #         "Please check your host, port, username, and password."
        #     )
        #     return None
        # except ImportError:
        #     logger.error(
        #         "The 'mysql-connector-python' driver was not found."
        #     )
        #     logger.error("Please install it by running: pip install mysql-connector-python")
        #     return None
        # except Exception as e:
        #     logger.error(f"An unexpected error occurred while creating DB engine: {e}")
        #     return None

    def read_csv_to_dataframe(self, file_path: str, header:int=0) -> Optional[pd.DataFrame]:
        """
        Reads the specified CSV file into a pandas DataFrame.

        Uses 'utf-8-sig' encoding to handle files with or without a
        Byte Order Mark (BOM), which is common for CSVs exported from Excel.

        Args:
            file_path (str): The path to the CSV file.

        Returns:
            Optional[pd.DataFrame]: A DataFrame if successful, else None.
        """
        try:
            df = pd.read_csv(file_path, encoding='utf-8-sig', header=header)
            logger.info(f"Successfully loaded {df.shape[0]} rows and {df.shape[1]} columns from '{file_path}'.")
            
            # Professional practice: Clean column names
            # Remove leading/trailing whitespace
            df.columns = [col.strip() for col in df.columns]
            logger.info("Cleaned column names (removed leading/trailing whitespace).")
            
            return df
        
        except FileNotFoundError:
            logger.error(f"Error: The file was not found at '{file_path}'.")
            return None
        except pd.errors.EmptyDataError:
            logger.error(f"Error: The file '{file_path}' is empty.")
            return None
        except UnicodeDecodeError:
            logger.error(f"Error: Could not decode '{file_path}' with 'utf-8-sig'.")
            logger.error("Please check the file's encoding. It might be GBK or another format.")
            return None
        except Exception as e:
            logger.error(f"An unexpected error occurred while reading the CSV: {e}")
            return None

    def write_dataframe_to_db(self, df: pd.DataFrame, table_name: str, if_exists: str = 'replace'):
        """
        Writes the DataFrame to the specified MySQL table.

        Args:
            df (pd.DataFrame): The data to write.
            table_name (str): The name of the target table.
            if_exists (str): How to behave if the table already exists.
                             'replace': Drop, re-create, and insert data.
                             'append': Insert new values into the existing table.
                             'fail': Raise a ValueError.
        """
        if self.engine is None:
            logger.error("Database engine is not initialized. Cannot write data.")
            return

        # CRITICAL STEP:
        # Map all 'object' (string) columns to sqlalchemy.types.TEXT.
        # This prevents 'DataError: (1406, "Data too long for column ...")'
        # which is very likely with your 'clinical_data.csv' file.
        dtype_mapping = {
            col: TEXT for col in df.select_dtypes(include=['object']).columns
        }
        logger.info(f"Mapped {len(dtype_mapping)} object columns to TEXT to prevent data truncation.")

        try:
            logger.info(f"Writing data to table '{table_name}' (mode: {if_exists})...")
            
            # Use pandas.to_sql with the SQLAlchemy engine
            # 'index=False' prevents pandas from writing the DataFrame index as a column
            # 'chunksize' writes data in batches, which is more memory-efficient
            df.to_sql(
                name=table_name,
                con=self.engine,
                if_exists=if_exists,
                index=False,
                dtype=dtype_mapping,
                chunksize=1000  # Insert 1000 rows at a time
            )
            
            logger.info(f"Successfully wrote {len(df)} rows to table '{table_name}'.")

        except exc.DataError as e:
            logger.error(f"Data error during insertion: {e}")
            logger.error("This might happen if data types still do not match.")
        except exc.SQLAlchemyError as e:
            logger.error(f"A database error occurred: {e}")
        except Exception as e:
            logger.error(f"An unexpected error occurred during database write: {e}")

    def run(self, csv_path: str, table_name: str, if_exists: str = 'replace', header:int=0):
        """
        Executes the full import process: read CSV, then write to DB.

        Args:
            csv_path (str): The path to the CSV file.
            table_name (str): The name of the target table.
            if_exists (str): Behavior if the table exists ('replace', 'append', 'fail').
        """
        logger.info(f"Starting import process for '{csv_path}'...")
        
        # Step 1: Read CSV
        df = self.read_csv_to_dataframe(csv_path, header)
        
        # Step 2: Write to DB, only if reading was successful and engine is valid
        if df is not None and self.engine is not None:
            self.write_dataframe_to_db(df, table_name, if_exists)
        elif self.engine is None:
            logger.error("Import process failed because database connection was not established.")
        else:
            logger.error("Import process failed because the CSV file could not be read.")
            
        logger.info("Import process finished.")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed arguments namespace

        python script.py --csv_path "/path/to/data.csv" --table_name "clinical_data" --if_exists replace
    """
    parser = argparse.ArgumentParser(description='Import clinical data to MySQL database')
    
    parser.add_argument(
        '--csv_path', 
        type=str,
        required=True,
        help='Path to the clinical data CSV file'
    )
    parser.add_argument(
        '--header', 
        type=int,
        required=True,
        help='Path to the clinical data CSV file'
    )
    
    parser.add_argument(
        '--table_name', 
        type=str,
        default='clinical_data',
        help='Target table name in MySQL (default: clinical_data)'
    )
    
    parser.add_argument(
        '--if_exists', 
        type=str,
        choices=['replace', 'append'],
        default='replace',
        help='Behavior if table exists: replace or append (default: replace)'
    )

    logging.info(f"Parsed arguments: {parser.parse_args()}")
    
    return parser.parse_args()

# --- Main Execution ---

if __name__ == "__main__":
    setup_logging()
    args = parse_arguments()

    try:
        # Create importer instance and run import process
        importer = MySQLDataImporter(DB_CONFIG)
        importer.run(args.csv_path, args.table_name, if_exists=args.if_exists, header=args.header)
        
        logging.info(f"Successfully imported data from {args.csv_path} to table {args.table_name}")
        
    except Exception as e:
        logging.info(f"Error during import: {str(e)}", file=sys.stderr)
        sys.exit(1)