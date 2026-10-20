process FOLD_SECONDARY {
    publishDir "${params.outdir}/03_secondary_structure", mode: 'copy'

    input:
    path typed

    output:
    path 'secondary.csv', emit: csv

    script:
    """
    python3 ${projectDir}/bin/fold_secondary.py --input ${typed} --output secondary.csv --temperature ${params.temperature}
    """
}
