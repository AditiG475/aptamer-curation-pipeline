// Literature / compound annotation branch (needs internet access and a CSV input
// with target, paper_title, journal and binding_text columns).

process EXTRACT_BINDING {
    input:
    path entries

    output:
    path 'binding.csv', emit: csv

    script:
    """
    python3 ${projectDir}/bin/extract_binding.py --input ${entries} --output binding.csv --text-column ${params.binding_text_column}
    """
}

process DOI_LOOKUP {
    publishDir "${params.outdir}/07_annotation", mode: 'copy', pattern: '*needs_review.csv'

    input:
    path binding

    output:
    path 'doi.csv', emit: csv
    path '*needs_review.csv', optional: true, emit: review

    script:
    """
    python3 ${projectDir}/bin/doi_lookup.py --input ${binding} --output doi.csv --min-similarity ${params.doi_min_similarity}
    """
}

process PUBCHEM_ANNOTATE {
    input:
    path doi

    output:
    path 'pubchem.csv', emit: csv

    script:
    """
    python3 ${projectDir}/bin/pubchem_annotate.py --input ${doi} --output pubchem.csv --name-column ${params.target_column}
    """
}

process CLASSIFY {
    publishDir "${params.outdir}/07_annotation", mode: 'copy'

    input:
    path pubchem
    path rules

    output:
    path 'annotated.csv', emit: csv

    script:
    """
    python3 ${projectDir}/bin/classify.py --input ${pubchem} --output annotated.csv --rules ${rules}
    """
}
