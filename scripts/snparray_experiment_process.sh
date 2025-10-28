
ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}

python data_process/snp_experiment_data.py \
    --data_path "D:\03.projects\AI.PGT\snparray_analysis\data\experiment_record\*芯片实验记录表*.xlsx" \
    --file_name multi_芯片实验记录表_all \
    --output_dir "D:\03.projects\AI.PGT\snparray_analysis\data"