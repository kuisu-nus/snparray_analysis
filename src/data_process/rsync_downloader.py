import subprocess
import logging
from typing import List, Dict, Any, Optional

class RsyncDownloader:
    """
    专门的rsync下载器类
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.logger = logging.getLogger(__name__)
    
    def sync_directory(self, remote_source: str, local_destination: str, 
                     sync_type: str = 'mirror') -> Dict[str, Any]:
        """
        同步整个目录。
        
        Args:
            remote_source (str): 远程源
            local_destination (str): 本地目标
            sync_type (str): 同步类型 - 'mirror', 'update', 'backup'
        """
        # base_options = ['-avz', '--progress']
        base_options = ['-r']
        
        # if sync_type == 'mirror':
        #     base_options.extend(['--delete', '--delete-excluded'])
        # elif sync_type == 'update':
        #     base_options.append('--update')
        # elif sync_type == 'backup':
        #     base_options.extend(['--backup', '--backup-dir=../backup'])
        
        # # 添加SSH配置
        # ssh_options = self._get_ssh_options()
        # if ssh_options:
        #     base_options.extend(['-e', ssh_options])
        
        return self._execute_rsync(base_options, remote_source, local_destination)
    
    def sync_pattern(self, remote_host: str, remote_path: str, patterns: List[str],
                   local_destination: str, recursive: bool = True) -> Dict[str, Any]:
        """
        同步匹配特定模式的文件。
        """
        options = ['-avz', '--progress']
        
        # 添加包含/排除模式
        for pattern in patterns:
            options.extend(['--include', pattern])
        options.extend(['--exclude', '*'])
        
        if recursive:
            options.append('--prune-empty-dirs')
        else:
            options.append('--no-recursive')
        
        # SSH配置
        ssh_options = self._get_ssh_options()
        if ssh_options:
            options.extend(['-e', ssh_options])
        
        remote_source = f"{remote_host}:{remote_path}/"
        
        return self._execute_rsync(options, remote_source, local_destination)
    
    def _get_ssh_options(self) -> str:
        """获取SSH配置选项"""
        ssh_config = self.config.get('LINUX_SERVER', {})
        options = []
        
        if ssh_config.get('port'):
            options.append(f"-p {ssh_config['port']}")
        
        if ssh_config.get('private_key_path'):
            options.append(f"-i {ssh_config['private_key_path']}")
        
        return ' '.join(options) if options else None
    
    def _execute_rsync(self, options: List[str], source: str, destination: str) -> Dict[str, Any]:
        """执行rsync命令"""
        # cmd = ['rsync'] + options + [source, destination]
        cmd = ['scp'] + options + [source, destination]
        
        self.logger.info(f"执行rsync: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2小时超时
            )
            
            return {
                'success': result.returncode == 0,
                'returncode': result.returncode,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'command': ' '.join(cmd)
            }
            
        except Exception as e:
            self.logger.error(f"rsync执行失败: {e}")
            return {
                'success': False,
                'error': str(e),
                'command': ' '.join(cmd)
            }