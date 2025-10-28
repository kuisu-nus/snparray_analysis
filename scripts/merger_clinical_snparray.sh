
ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}

python data_process/merger_data.py \
    --clinical_path "D:\03.projects\AI.PGT\snparray_analysis\data\clinical_data_pgt.csv" \
    --snparray_path "D:\03.projects\AI.PGT\snparray_analysis\data\multi_芯片实验记录表_all_persons.csv" \
    --output_dir "D:\03.projects\AI.PGT\snparray_analysis\work_dir" \
    --output_name merge_clinical_snparray