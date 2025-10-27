ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}


CSV_PATH="D:\03.projects\AI.PGT\snparray_analysis\data\5.芯片实验记录表2020.12.29_persons.csv"
TABLE_NAME="snparray_data"
IF_EXISTS="replace"

python sql_data/clinical_sql.py \
    --csv_path "${CSV_PATH}" \
    --table_name "${TABLE_NAME}" \
    --if_exists "${IF_EXISTS}" \
    --header 0