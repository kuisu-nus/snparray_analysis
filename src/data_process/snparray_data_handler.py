import pandas as pd
import pymysql
import logging
import os
import subprocess
import paramiko
from typing import Optional, List, Dict, Any, Callable
import configparser
from pathlib import Path

from rsync_downloader import RsyncDownloader


class SNPArrayDataHandler:
    """
    A comprehensive class to handle SNP array data operations including 
    MySQL connection, data processing, and secure file downloading from Linux systems.
    """
    
    def __init__(self, config_file: str = 'config.ini', log_level: str = 'INFO'):
        """
        Initialize the SNPArrayDataHandler.
        
        Args:
            config_file (str): Path to configuration file
            log_level (str): Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.setup_logging(log_level)
        self.logger = logging.getLogger(__name__)
        self.config = self.load_config(config_file)
        self.connection = None
        self.dataframe = None
        self.ssh_client = None
        self.sftp_client = None
        
    def setup_logging(self, log_level: str) -> None:
        """Setup logging configuration."""
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('snp_array_handler.log'),
                logging.StreamHandler()
            ]
        )
    
    def load_config(self, config_file: str) -> configparser.ConfigParser:
        """
        Load configuration from INI file.
        
        Args:
            config_file (str): Path to configuration file
            
        Returns:
            configparser.ConfigParser: Configuration object
        """
        config = configparser.ConfigParser()
        try:
            if not os.path.exists(config_file):
                self._create_default_config(config_file)
            config.read(config_file, encoding='utf-8')
            self.logger.info(f"Configuration loaded from {config_file}")
            return config
        except Exception as e:
            self.logger.error(f"Error loading configuration: {e}")
            raise
    
    def _create_default_config(self, config_file: str) -> None:
        """Create default configuration file with comprehensive settings."""
        config = configparser.ConfigParser()
        
        # Database configuration
        config['DATABASE'] = {
            'host': 'localhost',
            'port': '3306',
            'user': 'your_username',
            'password': 'your_password',
            'database': 'ivf_data',
            'charset': 'utf8mb4',
            'table_name': 'merge_clinical_snparray'
        }
        
        # Linux/SSH connection configuration
        config['LINUX_SERVER'] = {
            'host': 'your_linux_server',
            'port': '22',
            'username': 'your_username',
            'password': 'your_password',
            'private_key_path': '~/.ssh/id_rsa',
            'timeout': '30'
        }
        
        # Path configuration
        config['PATHS'] = {
            'remote_base_path': '/remote/path/to/snp_array_data/',
            'local_download_dir': './downloaded_files/',
            'output_csv': 'snp_array_data.csv',
            'temp_dir': './temp/'
        }
        
        # Download configuration
        config['DOWNLOAD'] = {
            'max_retries': '3',
            'timeout': '300',
            'chunk_size': '8192',
            'parallel_downloads': '5'
        }
        
        # File patterns and path generation
        config['FILE_PATTERNS'] = {
            'path_template': '{base_path}/{idx}/{sub_chip_idx}.cel',
            'file_extension': '.cel',
            'supported_extensions': '.cel,.txt,.csv'
        }
        
        # Query configuration
        config['QUERY'] = {
            'custom_where_clause': '',
            'order_by': 'idx ASC',
            'limit': '0'
        }
        
        with open(config_file, 'w') as f:
            config.write(f)
        self.logger.info(f"Default configuration created at {config_file}")
    
    def connect_to_mysql(self) -> bool:
        """
        Establish connection to MySQL database.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            db_config = self.config['DATABASE']
            self.connection = pymysql.connect(
                host=db_config.get('host'),
                port=int(db_config.get('port', 3306)),
                user=db_config.get('user'),
                password=db_config.get('password'),
                database=db_config.get('database'),
                charset=db_config.get('charset', 'utf8mb4'),
                cursorclass=pymysql.cursors.DictCursor
            )
            self.logger.info("Successfully connected to MySQL database")
            return True
        except Exception as e:
            self.logger.error(f"Failed to connect to MySQL: {e}")
            return False
    
    def connect_to_linux_server(self) -> bool:
        """
        Establish SSH connection to Linux server for file downloads.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            server_config = self.config['LINUX_SERVER']
            self.ssh_client = paramiko.SSHClient()
            self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            
            # Connection parameters
            connect_params = {
                'hostname': server_config.get('host'),
                'port': int(server_config.get('port', 22)),
                'username': server_config.get('username'),
                'timeout': int(server_config.get('timeout', 30))
            }
            
            # Use password or private key authentication
            password = server_config.get('password')
            private_key_path = server_config.get('private_key_path')
            
            if password and password != 'your_password':
                connect_params['password'] = password
            elif private_key_path and private_key_path != '~/.ssh/id_rsa':
                private_key = paramiko.RSAKey.from_private_key_file(
                    os.path.expanduser(private_key_path)
                )
                connect_params['pkey'] = private_key
            
            self.ssh_client.connect(**connect_params)
            self.sftp_client = self.ssh_client.open_sftp()
            
            self.logger.info("Successfully connected to Linux server via SSH")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to connect to Linux server: {e}")
            return False
    
    def fetch_snparray_data(self) -> Optional[pd.DataFrame]:
        """
        Fetch data from the configured MySQL table.
        
        Returns:
            Optional[pd.DataFrame]: DataFrame containing SNP array data or None if failed
        """
        if not self.connection:
            self.logger.error("No database connection established")
            return None
        
        # Get table name from config
        table_name = self.config['DATABASE'].get('table_name', 'merge_clinical_snparray')
        #  检查表记录数
        
        
        # Build query with configurable options
        query = f"SELECT * FROM {table_name}"
        
        # Add WHERE clause if specified
        where_clause = self.config['QUERY'].get('custom_where_clause')
        if where_clause:
            query += f" WHERE {where_clause}"
        
        # Add ORDER BY if specified
        order_by = self.config['QUERY'].get('order_by')
        if order_by:
            query += f" ORDER BY {order_by}"
        
        # Add LIMIT if specified
        limit = self.config['QUERY'].get('limit', '0')
        if limit and limit != '0':
            query += f" LIMIT {limit}"
        
        try:
            with self.connection.cursor() as cursor:
                cursor.execute(query)
                query_data = cursor.fetchall()
                self.logger.info(f"获取到 {len(query_data)} 行数据")
                
            self.dataframe = pd.DataFrame(query_data)
            self.logger.info(f"Successfully fetched {len(self.dataframe)} records from {table_name}\n{self.dataframe.head()}")
            return self.dataframe
        except Exception as e:
            self.logger.error(f"Error fetching data from MySQL: {e}")
            return None
    
    def generate_download_paths(self, custom_path_function: Callable = None) -> None:
        """
        Generate download paths based on configuration or custom function.
        
        Args:
            custom_path_function (Callable, optional): Custom function to generate paths
        """
        if self.dataframe is None:
            self.logger.error("No data available. Please fetch data first.")
            return
        
        # try:
        if custom_path_function:
            # Use custom path generation function
            self.dataframe['remote_path'] = self.dataframe.apply(
                lambda row: custom_path_function(row), 
                axis=1
            )
            self.logger.info("Custom download paths generated")
        else:
            # Use configured path template
            path_template = self.config['FILE_PATTERNS'].get(
                'path_template', 
                '{base_path}/{idx}/{sub_chip_idx}.cel'
            )
            remote_base_path = self.config['PATHS'].get('remote_base_path', '')
            if "idx" in self.dataframe.columns:
                self.dataframe['idx'] = self.dataframe['idx'].astype(int)
            self.dataframe['remote_path'] = self.dataframe.apply(
                lambda row: self._apply_path_template(row, path_template, remote_base_path), 
                axis=1
            )
            self.logger.info("Config-based download paths generated")
        
        # Generate local paths
        local_download_dir = Path(self.config['PATHS'].get('local_download_dir', './downloaded_files/'))
        os.makedirs(local_download_dir, exist_ok=True)
        
        self.dataframe['local_path'] = self.dataframe['remote_path'].apply(
            lambda remote_path: self._generate_local_path(remote_path, local_download_dir)
        )
            
        # except Exception as e:
        #     self.logger.error(f"Error generating download paths: {e}")
    
    def _apply_path_template(self, row: pd.Series, template: str, base_path: str) -> str:
        """
        Apply path template to generate remote file path.
        
        Args:
            row (pd.Series): DataFrame row
            template (str): Path template string
            base_path (str): Base remote path
            
        Returns:
            str: Generated remote path
        """
        try:
            # Replace base_path placeholder
            path = template.replace('{base_path}', base_path)
            
            # Replace other column placeholders
            for column in row.index:
                placeholder = f'{{{column}}}'
                if placeholder in path:
                    path = path.replace(placeholder, str(row[column]))
            
            return path
        except Exception as e:
            self.logger.error(f"Error applying path template: {e}")
            return ""
    
    def _generate_local_path(self, remote_path: str, local_download_dir: str) -> str:
        """
        Generate local path from remote path.
        
        Args:
            remote_path (str): Remote file path
            local_download_dir (str): Local download directory
            
        Returns:
            str: Local file path
        """
        try:
            # Extract filename from remote path
            filename = os.path.basename(os.path.dirname(remote_path))
            
            # Create local path
            local_path = os.path.join(local_download_dir, filename)
            
            return local_path
        except Exception as e:
            self.logger.error(f"Error generating local path: {e}")
            return ""
    
    def save_to_csv(self, output_file: str = None) -> bool:
        """
        Save dataframe to CSV file.
        
        Args:
            output_file (str, optional): Output CSV file path
            
        Returns:
            bool: True if save successful, False otherwise
        """
        if self.dataframe is None:
            self.logger.error("No data available to save")
            return False
        
        try:
            if output_file is None:
                output_file = self.config['PATHS'].get('output_csv', 'snp_array_data.csv')
            
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', 
                       exist_ok=True)
            
            self.dataframe.to_csv(output_file, index=False)
            self.logger.info(f"Data successfully saved to {output_file}")
            return True
        except Exception as e:
            self.logger.error(f"Error saving to CSV: {e}")
            return False
    
    def download_files(self, max_retries: int = None, timeout: int = None) -> Dict[str, Any]:
        """
        Download SNP array files from Linux server to local machine.
        
        Args:
            max_retries (int, optional): Maximum number of download retries
            timeout (int, optional): Timeout for download commands
            
        Returns:
            Dict[str, Any]: Download statistics
        """
        if self.dataframe is None or 'remote_path' not in self.dataframe.columns:
            self.logger.error("No download paths available. Please generate paths first.")
            return {'success': False, 'message': 'No download paths available'}
        
        if not self.sftp_client:
            if not self.connect_to_linux_server():
                return {'success': False, 'message': 'Failed to connect to Linux server'}
        
        if max_retries is None:
            max_retries = int(self.config['DOWNLOAD'].get('max_retries', 3))
        
        if timeout is None:
            timeout = int(self.config['DOWNLOAD'].get('timeout', 300))
        
        download_stats = {
            'total_files': 0,
            'successful_downloads': 0,
            'failed_downloads': 0,
            'failed_paths': [],
            'downloaded_files': []
        }
        
        # Filter valid paths
        valid_records = self.dataframe[
            self.dataframe['remote_path'].notna() & 
            (self.dataframe['remote_path'] != '')
        ]
        download_stats['total_files'] = len(valid_records)
        
        rsync_config = self.config["RSYNC"]
        rsync_downloader = RsyncDownloader(rsync_config)
        for _, record in valid_records.iterrows():
            remote_path = record['remote_path']
            local_path = record['local_path']
            
            rsync_status = rsync_downloader.sync_directory(
                remote_source=f"{self.config['LINUX_SERVER']['username']}@{self.config['LINUX_SERVER']['host']}:{remote_path}",
                local_destination=local_path,
                sync_type=rsync_config.get('sync_type', 'mirror')
            )
            if rsync_status["success"]:
                download_stats['successful_downloads'] += 1
                download_stats['downloaded_files'].append({
                    'remote_path': remote_path,
                    'local_path': local_path
                })
            else:
                download_stats['failed_downloads'] += 1
                download_stats['failed_paths'].append(remote_path)
        
        self.logger.info(
            f"Download completed: {download_stats['successful_downloads']}/"
            f"{download_stats['total_files']} files downloaded successfully"
        )
        
        return download_stats
    
    def _download_single_file_sftp(self, remote_path: str, local_path: str, max_retries: int) -> bool:
        """
        Download a single file using SFTP.
        
        Args:
            remote_path (str): Remote file path on Linux server
            local_path (str): Local file path to save to
            max_retries (int): Maximum retry attempts
            
        Returns:
            bool: True if download successful, False otherwise
        """
        for attempt in range(max_retries):
            try:
                # Create local directory if it doesn't exist
                local_dir = os.path.dirname(local_path)
                os.makedirs(local_dir, exist_ok=True)
                
                # Download file via SFTP
                self.sftp_client.get(remote_path, local_path)
                
                # Verify file was downloaded
                if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
                    self.logger.info(f"Successfully downloaded: {remote_path} -> {local_path}")
                    return True
                else:
                    self.logger.warning(f"Downloaded file verification failed: {local_path}")
                    
            except Exception as e:
                self.logger.warning(
                    f"Download attempt {attempt + 1} failed for {remote_path}: {e}"
                )
        
        self.logger.error(f"All download attempts failed for {remote_path}")
        return False
    
    def download_files_rsync_batch(self, download_configs: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        批量rsync下载。
        
        Args:
            download_configs: 下载配置列表
                Example:
                [{
                    'remote_host': 'user@server',
                    'remote_path': '/path/to/snp_data',
                    'file_pattern': 'AAA*.jpg',
                    'local_destination': './downloads/snp_images/',
                    'sync_type': 'mirror'
                }]
        """
        rsync_downloader = RsyncDownloader(self.config)
        
        results = {
            'total_tasks': len(download_configs),
            'successful_tasks': 0,
            'failed_tasks': 0,
            'task_results': []
        }
        
        for config in download_configs:
            try:
                if 'file_pattern' in config:
                    # 模式同步
                    result = rsync_downloader.sync_pattern(
                        remote_host=config['remote_host'],
                        remote_path=config['remote_path'],
                        patterns=[config['file_pattern']],
                        local_destination=config['local_destination']
                    )
                else:
                    # 目录同步
                    result = rsync_downloader.sync_directory(
                        remote_source=f"{config['remote_host']}:{config['remote_path']}/",
                        local_destination=config['local_destination'],
                        sync_type=config.get('sync_type', 'mirror')
                    )
                
                task_result = {
                    'config': config,
                    'success': result['success'],
                    'details': result
                }
                
                if result['success']:
                    results['successful_tasks'] += 1
                    self.logger.info(f"rsync任务成功: {config['remote_path']}")
                else:
                    results['failed_tasks'] += 1
                    self.logger.error(f"rsync任务失败: {config['remote_path']} - {result.get('stderr', '')}")
                
                results['task_results'].append(task_result)
                
            except Exception as e:
                self.logger.error(f"rsync任务异常: {config['remote_path']} - {e}")
                results['failed_tasks'] += 1
                results['task_results'].append({
                    'config': config,
                    'success': False,
                    'error': str(e)
                })
        
        return results

    def close_connections(self) -> None:
        """Close all connections (MySQL and SSH)."""
        if self.connection:
            self.connection.close()
            self.logger.info("MySQL connection closed")
        
        if self.sftp_client:
            self.sftp_client.close()
        
        if self.ssh_client:
            self.ssh_client.close()
            self.logger.info("SSH connection closed")
    
    def run_pipeline(self, output_file: str = None, custom_path_function: Callable = None) -> Dict[str, Any]:
        """
        Run the complete SNP array data processing pipeline.
        
        Args:
            output_file (str, optional): Output CSV file path
            custom_path_function (Callable, optional): Custom path generation function
            
        Returns:
            Dict[str, Any]: Pipeline execution results
        """
        results = {
            'database_connection': False,
            'data_fetched': False,
            'paths_generated': False,
            'csv_saved': False,
            'linux_connection': False,
            'download_stats': None
        }
        
        # try:
        # Step 1: Connect to MySQL
        if self.connect_to_mysql():
            results['database_connection'] = True
            
            # Step 2: Fetch data
            if self.fetch_snparray_data() is not None:
                results['data_fetched'] = True
                
                # Step 3: Generate download paths
                self.generate_download_paths(custom_path_function)
                results['paths_generated'] = True
                
                # Step 4: Save to CSV
                if self.save_to_csv(output_file):
                    results['csv_saved'] = True
                    
                    # Step 5: Connect to Linux and download files
                    if self.connect_to_linux_server():
                        results['linux_connection'] = True
                        results['download_stats'] = self.download_files()
        
        return results
            
        # except Exception as e:
        #     self.logger.error(f"Pipeline execution failed: {e}")
        #     return results
        # finally:
        #     self.close_connections()


# Example custom path generation function
def custom_path_generator(row: pd.Series) -> str:
    """
    Example custom function for generating download paths.
    
    Args:
        row (pd.Series): DataFrame row containing data
        
    Returns:
        str: Generated remote file path
    """
    # Custom logic for path generation
    base_path = "/custom/base/path"
    idx = row.get('idx', 'unknown')
    sub_chip_idx = row.get('sub_chip_idx', 'unknown')
    
    # Example: /custom/base/path/IDX_123/SUBCHIP_456.cel
    return f"{base_path}/IDX_{idx}/SUBCHIP_{sub_chip_idx}.cel"


# Example usage and main execution
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='SNP Array Data Handler')
    parser.add_argument('--config', 
                        default=r'D:\03.projects\AI.PGT\snparray_analysis\src\configs\snparray_data_handler_config.ini', 
                        help='Configuration file path')
    parser.add_argument('--output', help='Output CSV file path')
    parser.add_argument('--log-level', default='INFO', help='Logging level')
    parser.add_argument('--use-custom-paths', action='store_true', 
                       help='Use custom path generation function')
    
    args = parser.parse_args()
    
    # Initialize handler
    handler = SNPArrayDataHandler(config_file=args.config, log_level=args.log_level)
    
    # Run complete pipeline
    custom_function = custom_path_generator if args.use_custom_paths else None
    results = handler.run_pipeline(output_file=args.output, custom_path_function=custom_function)
    
    # Print results
    print("\nPipeline Execution Results:")
    for key, value in results.items():
        print(f"{key}: {value}")