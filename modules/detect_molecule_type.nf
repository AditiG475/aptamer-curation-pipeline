process DETECT_MOLECULE_TYPE {
    publishDir "${params.outdir}/02_molecule_type", mode: 'copy'

    input:
    path cleaned

    output:
    path 'typed.csv', emit: csv

    script:
    """
    python3 ${projectDir}/bin/detect_molecule_type.py --input ${cleaned} --output typed.csv --default-type ${params.default_type}
    """
}
