process VERIFY_STRUCTURES {
    publishDir "${params.outdir}/06_qc", mode: 'copy'

    input:
    path typed
    path secondary
    path zscore

    output:
    path 'qc_checks.csv', emit: csv
    path 'qc_report.txt', emit: report

    script:
    def strict = params.strict ? '--strict' : ''
    """
    python3 ${projectDir}/bin/verify_structures.py --typed ${typed} --secondary ${secondary} --zscore ${zscore} \\
        --shuffles ${params.shuffles} --output qc_checks.csv ${strict} | tee qc_report.txt
    """
}
