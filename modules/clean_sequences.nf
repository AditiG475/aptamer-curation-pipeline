process CLEAN_SEQUENCES {
    tag "${input.name}"
    publishDir "${params.outdir}/01_cleaned", mode: 'copy'

    input:
    path input

    output:
    path 'cleaned.csv', emit: csv

    script:
    """
    python3 ${projectDir}/bin/clean_sequences.py --input ${input} --output cleaned.csv \\
        --id-column ${params.id_column} --seq-column ${params.seq_column} \\
        --min-length ${params.min_length} --max-length ${params.max_length}
    """
}
