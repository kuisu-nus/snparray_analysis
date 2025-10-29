
# SNP Array Analysis

一个用于处理和分析SNP芯片数据的Python工具包，支持数据导入、处理、合并和存储等功能。

## 功能特性

- SNP芯片数据处理（PLINK格式）
- 实验数据解析和转换
- 临床数据与SNP数据合并
- MySQL数据库集成
- 文件同步和下载

## 项目结构

```
snparray_analysis/
├── data/                    # 数据目录
├── docs/                    # 文档目录
├── scripts/                 # 脚本目录
├── src/                     # 源代码目录
│   ├── data_process/        # 数据处理模块
│   ├── sql_data/           # SQL数据模块
│   └── notebook/           # Jupyter notebooks
└── work_dir/               # 工作目录
```

## 核心模块

### 1. 数据处理模块 (src/data_process/)

- `snparray_data.py`: SNP芯片数据处理
  - 支持PLINK格式数据加载
  - 数据格式转换和导出

- `snp_experiment_data.py`: 实验数据解析
  - 芯片实验记录解析
  - 支持多种数据格式输出

- `merger_data.py`: 数据合并
  - 临床数据与SNP数据合并
  - 支持自定义合并策略

- `snparray_data_handler.py`: 数据处理器
  - MySQL数据库连接
  - 文件下载和同步
  - 数据处理流水线

- `rsync_downloader.py`: 文件同步工具
  - 支持rsync协议
  - 批量文件下载

### 2. SQL数据模块 (src/sql_data/)

- `clinical_sql.py`: MySQL数据导入工具
  - CSV数据导入
  - 支持UTF-8编码
  - 数据库表管理

- `pymysql.py`: MySQL操作封装
  - 数据库连接管理
  - CRUD操作封装
  - 事务处理

## 使用说明

### 1. 环境要求

- Python 3.6+
- MySQL 5.7+
- 必要的Python包（见requirements.txt）

### 2. 数据处理流程

1. SNP芯片数据处理
```bash
python src/data_process/snparray_data.py \
    --data_root "path/to/data" \
    --file_name "filename" \
    --output_root "path/to/output"
```

2. 实验数据解析
```bash
python src/data_process/snp_experiment_data.py \
    --data_path "path/to/experiment_data" \
    --file_name "output_name" \
    --output_dir "path/to/output"
```

3. 数据合并
```bash
python src/data_process/merger_data.py \
    --clinical_path "path/to/clinical_data" \
    --snparray_path "path/to/snparray_data" \
    --output_dir "path/to/output" \
    --output_name "merge_result"
```

### 3. 数据库操作

1. 创建表
```bash
python src/sql_data/clinical_sql.py \
    --csv_path "path/to/data.csv" \
    --table_name "table_name" \
    --if_exists "replace"
```

2. 数据导入
```bash
python src/sql_data/pymysql.py
```

## 注意事项

1. 数据格式要求
   - SNP数据：PLINK格式
   - 实验数据：Excel或CSV格式
   - 临床数据：CSV格式

2. 数据清洗规则
   - 同一病人可能存在多个周期
   - 使用"姓名-PGD编号-年份-遗传材料接收者"作为唯一标识

3. 异常处理
   - 详细的错误日志记录
   - 异常数据标记和处理

## 开发说明

1. 代码规范
   - 遵循PEP 8规范
   - 完整的文档字符串
   - 类型提示支持

2. 测试
   - 单元测试覆盖
   - 集成测试验证

## 贡献指南

1. Fork项目
2. 创建特性分支
3. 提交变更
4. 推送到分支
5. 创建Pull Request

## 许可证

MIT License

## 联系方式

如有问题或建议，请提交Issue或联系项目维护者。