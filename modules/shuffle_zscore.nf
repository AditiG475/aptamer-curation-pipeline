process SHUFFLE_ZSCORE {
    publishDir "${params.outdir}/04_zscore", mode: 'copy'

    input:
    path typed

    output:
    path 'zscore.csv', emit: csv

    script:
    """
    python3 ${projectDir}/bin/shuffle_zscore.py --input ${typed} --output zscore.csv \\
        --shuffles ${params.shuffles} --seed ${params.seed} --temperature ${params.temperature}
    """
}
