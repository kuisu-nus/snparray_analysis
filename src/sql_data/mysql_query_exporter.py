#!/usr/bin/env python3
"""
MySQL Query Export Tool
A comprehensive tool for querying MySQL database and exporting results to CSV/XLSX format.
"""

import pandas as pd
import pymysql
import logging
import argparse
import sys
import os
from typing import Optional, Dict, Any
import json


class MySQLQueryExporter:
    """
    A class to handle MySQL database queries and export results to various formats.
    """
    
    def __init__(self, config_file: Optional[str] = None, **kwargs):
        """
        Initialize the MySQL Query Exporter.
        
        Args:
            config_file (str, optional): Path to JSON configuration file
            **kwargs: Database connection parameters
        """
        self.setup_logging()
        self.db_config = self.load_config(config_file, kwargs)
        self.connection = None
        self.logger = logging.getLogger(__name__)
        
    def setup_logging(self, log_level: str = "INFO") -> None:
        """
        Set up logging configuration.
        
        Args:
            log_level (str): Logging level (DEBUG, INFO, WARNING, ERROR)
        """
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler('mysql_query_export.log')
            ]
        )
        
    def load_config(self, config_file: Optional[str], kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Load configuration from file and/or keyword arguments.
        
        Args:
            config_file (str, optional): Path to JSON config file
            kwargs: Additional configuration parameters
            
        Returns:
            Dict[str, Any]: Combined configuration dictionary
        """
        config = {
            'host': 'localhost',
            'port': 3306,
            'user': 'root',
            'password': '',
            'database': 'ivf_data',
            'charset': 'utf8mb4'
        }
        
        # Load from config file if provided
        if config_file and os.path.exists(config_file):
            try:
                with open(config_file, 'r') as f:
                    file_config = json.load(f)
                config.update(file_config)
                logging.info(f"Configuration loaded from {config_file}")
            except Exception as e:
                logging.error(f"Error loading config file: {e}")
                raise
        
        # Update with keyword arguments
        config.update(kwargs)
        
        return config
    
    def connect(self) -> bool:
        """
        Establish connection to MySQL database.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        # try:``
        self.connection = pymysql.connect(
                host=self.db_config.get('host'),
                port=int(self.db_config.get('port', 3306)),
                user=self.db_config.get('user'),
                password=self.db_config.get('password'),
                database=self.db_config.get('database'),
                charset=self.db_config.get('charset', 'utf8mb4'),
                cursorclass=pymysql.cursors.DictCursor
            )
        self.logger.info(f"Successfully connected to database: {self.db_config['database']}")
        return True
        # except Exception as e:
        #     self.logger.error(f"Database connection failed: {e}")
        #     return False
    
    def execute_query(self, query: str, params: Optional[tuple] = None) -> Optional[pd.DataFrame]:
        """
        Execute SQL query and return results as DataFrame.
        
        Args:
            query (str): SQL query to execute
            params (tuple, optional): Query parameters
            
        Returns:
            pd.DataFrame: Query results as DataFrame, None if error occurs
        """
        if not self.connection:
            if not self.connect():
                return None
        
        try:
            with self.connection.cursor() as cursor:
                self.logger.info(f"Executing query: {query}")
                if params:
                    cursor.execute(query, params)
                    self.logger.debug(f"Query parameters: {params}")
                else:
                    cursor.execute(query)
                
                results = cursor.fetchall()
                self.logger.info(f"Query returned {len(results)} rows")
                
                # Convert to DataFrame
                df = pd.DataFrame(results)
                if "idx" in df.columns:
                    df['idx'] = df['idx'].astype(int)
                return df
                
        except Exception as e:
            self.logger.error(f"Query execution failed: {e}")
            return None
    
    def export_to_file(self, df: pd.DataFrame, output_file: str, 
                      file_format: str = 'auto') -> bool:
        """
        Export DataFrame to file in specified format.
        
        Args:
            df (pd.DataFrame): DataFrame to export
            output_file (str): Output file path
            file_format (str): File format ('csv', 'xlsx', or 'auto')
            
        Returns:
            bool: True if export successful, False otherwise
        """
        if df is None or df.empty:
            self.logger.warning("No data to export")
            return False
        
        try:
            # Determine file format
            if file_format == 'auto':
                if output_file.lower().endswith('.xlsx'):
                    file_format = 'xlsx'
                else:
                    file_format = 'csv'
            
            # Export based on format
            if file_format.lower() == 'csv':
                df.to_csv(output_file, index=False)
                self.logger.info(f"Data exported to CSV: {output_file}")
            elif file_format.lower() == 'xlsx':
                df.to_excel(output_file, index=False)
                self.logger.info(f"Data exported to Excel: {output_file}")
            else:
                self.logger.error(f"Unsupported file format: {file_format}")
                return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"Export failed: {e}")
            return False
    
    def run_export(self, query: str, output_file: str, 
                  params: Optional[tuple] = None,
                  file_format: str = 'auto') -> bool:
        """
        Complete workflow: query database and export results.
        
        Args:
            query (str): SQL query to execute
            output_file (str): Output file path
            params (tuple, optional): Query parameters
            file_format (str): Output file format
            
        Returns:
            bool: True if entire process successful, False otherwise
        """
        self.logger.info("Starting database query and export process")
        
        # Execute query
        df = self.execute_query(query, params)
        if df is None:
            return False
        
        # Export results
        success = self.export_to_file(df, output_file, file_format)
        
        if success:
            self.logger.info("Query and export completed successfully")
        else:
            self.logger.error("Query and export process failed")
        
        return success
    
    def close(self) -> None:
        """Close database connection."""
        if self.connection:
            self.connection.close()
            self.logger.info("Database connection closed")
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


def main():
    """Main function to handle command line arguments and execute export."""
    parser = argparse.ArgumentParser(description='MySQL Query Export Tool')
    
    # Database connection arguments
    parser.add_argument('--host', default='localhost', help='MySQL host')
    parser.add_argument('--port', type=int, default=3306, help='MySQL port')
    parser.add_argument('--user', required=True, help='MySQL username')
    parser.add_argument('--password', required=True, help='MySQL password')
    parser.add_argument('--database', default='ivf_data', help='MySQL database name')
    
    # Query and export arguments
    parser.add_argument('--query', required=True, help='SQL query to execute')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--format', choices=['csv', 'xlsx', 'auto'], 
                       default='auto', help='Output file format')
    parser.add_argument('--config', help='JSON configuration file path')
    
    # Additional options
    parser.add_argument('--log-level', default='INFO', 
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')
    
    args = parser.parse_args()

    logging.info(f"Starting MySQL Query Export Tool with arguments: {args}")
    
    # Create exporter instance
    try:
        with MySQLQueryExporter(
            config_file=args.config,
            host=args.host,
            port=args.port,
            user=args.user,
            password=args.password,
            database=args.database
        ) as exporter:
            
            # Set logging level
            exporter.setup_logging(args.log_level)
            
            # Execute query and export
            success = exporter.run_export(
                query=args.query,
                output_file=args.output,
                file_format=args.format
            )
            
            sys.exit(0 if success else 1)
            
    except Exception as e:
        logging.error(f"Application error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # Example usage
    if len(sys.argv) == 1:
        # Demonstrate class usage without command line arguments
        print("MySQL Query Export Tool")
        print("Usage examples:")
        print("1. Command line: python script.py --user root --password pass --query 'SELECT * FROM merge_clinical_snparray' --output results.csv")
        print("2. Programmatic usage:")
        
        # Example programmatic usage
        exporter = MySQLQueryExporter(
            host='localhost',
            user='root',
            password='sukui1016',
            database='ivf_data'
        )
        
        # Example query
        query = "SELECT * FROM merge_clinical_snparray LIMIT 10"
        success = exporter.run_export(
            query=query,
            output_file='sample_output.csv',
            file_format='csv'
        )
        
        exporter.close()
        
        if success:
            print("Example execution completed successfully")
        else:
            print("Example execution failed")
    else:
        main()