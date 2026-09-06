process CHECKDE {

    input:
    path diff_dir

    output:
    tuple path("checkDE_result.txt"), path(diff_dir), emit: checked

    script:
    def de_number = params.de_number ?: 100
    """
    set +e
    bash /workplace/pipeline/WTSS/scripts/checkDE.sh ${diff_dir}/ ${de_number} checkDE.txt > checkDE.log 2>&1
    rc=\$?
    if [ \$rc -eq 0 ]; then
        echo "PASS" > checkDE_result.txt
    else
        echo "FAIL" > checkDE_result.txt
        cat checkDE.txt.error >> checkDE_result.txt 2>/dev/null || true
    fi
    """
}
