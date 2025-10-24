
ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}

python data_process/snparray_data.py \
    --data_root "D:\03.projects\AI.PGT\data\pengjingjing\SNP_result\pengjingjign_203368710006\PLINK_241025_0935" \
    --file_name pengjingjign_203368710006 \
    --output_root "D:\03.projects\AI.PGT\snparray_analysis\work_dir"
