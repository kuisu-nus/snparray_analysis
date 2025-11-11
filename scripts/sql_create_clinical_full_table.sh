ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}


CSV_PATH="D:\01.data\01.gene\clinical_transfer_data.csv"
TABLE_NAME="clinical_transfer_data"
IF_EXISTS="replace"

python sql_data/clinical_sql.py \
    --csv_path "${CSV_PATH}" \
    --table_name "${TABLE_NAME}" \
    --if_exists "${IF_EXISTS}" \
    --header 0