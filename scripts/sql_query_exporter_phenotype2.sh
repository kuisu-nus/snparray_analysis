ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}

    QUERY="
    SELECT 
    m.女方姓名_clinical,
    m.idx,
    m.chip_sub_idx, 
    m.男方姓名,
    m.女方姓名_clinical,
    c.\`是否妊娠\`
FROM merge_clinical_snparray AS m
LEFT JOIN clinical_transfer_data AS c ON m.E_IDNo = c.E_IDNo
    "

    python sql_data/mysql_query_exporter.py \
        --user root \
        --password sukui1016 \
        --query "${QUERY}" \
        --output "D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype2.csv" \
        --config "D:\03.projects\AI.PGT\snparray_analysis\src\configs\mysql_query_exporter.json"