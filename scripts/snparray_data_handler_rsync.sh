
ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}

# python data_process/snparray_data_handler.py \
#     --config "D:\03.projects\AI.PGT\snparray_analysis\src\configs\snparray_data_handler_config.ini" \
#     --log-level "INFO"

python data_process/snparray_data_handler.py \
    --config "D:\03.projects\AI.PGT\snparray_analysis\src\configs\snparray_data_handler_config.clinical_trainsfer_snparray.ini" \
    --log-level "INFO"