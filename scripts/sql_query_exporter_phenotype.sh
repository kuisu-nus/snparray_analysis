ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}

QUERY="
SELECT 
    女方姓名_clinical,
    idx,
    chip_sub_idx, 
    男方姓名,
    女方姓名_clinical,
    诊断可移植
FROM merge_clinical_snparray
"

python sql_data/mysql_query_exporter.py \
    --user root \
    --password sukui1016 \
    --query "${QUERY}" \
    --output "D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype.csv" \
    --config "D:\03.projects\AI.PGT\snparray_analysis\src\configs\mysql_query_exporter.json"