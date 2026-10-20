process MERGE_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path typed
    path secondary
    path zscore
    path tertiary
    path annotation

    output:
    path 'final_dataset.csv', emit: csv

    script:
    def tert = tertiary.name.startsWith('NO_FILE') ? '' : "--tertiary ${tertiary}"
    def annot = annotation.name.startsWith('NO_FILE') ? '' : "--annotation ${annotation}"
    """
    python3 ${projectDir}/bin/merge_results.py --typed ${typed} --secondary ${secondary} --zscore ${zscore} ${tert} ${annot} \\
        --output final_dataset.csv
    """
}

process EXPORT_EXCEL {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_csv

    output:
    path 'final_dataset.xlsx'

    script:
    """
    python3 ${projectDir}/bin/export_excel.py --input ${final_csv} --output final_dataset.xlsx
    """
}
