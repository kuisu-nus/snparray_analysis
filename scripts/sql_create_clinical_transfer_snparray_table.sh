ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}


CSV_PATH="D:\03.projects\AI.PGT\snparray_analysis\work_dir\merge_clinical_transfer_snparray.csv"
TABLE_NAME="merge_clinical_transfer_snparray"
IF_EXISTS="replace"

python sql_data/clinical_sql.py \
    --csv_path "${CSV_PATH}" \
    --table_name "${TABLE_NAME}" \
    --if_exists "${IF_EXISTS}" \
    --header 0