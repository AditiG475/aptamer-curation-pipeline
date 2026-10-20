process PREDICT_TERTIARY {
    tag "${entry_id}"
    label 'tertiary'
    publishDir "${params.outdir}/05_tertiary_structure", mode: 'copy'

    input:
    tuple val(entry_id), val(sequence), val(molecule_type)

    output:
    path "${entry_id}_3d/${entry_id}_3d_summary.csv", emit: summary
    path "${entry_id}_3d/*", emit: files

    script:
    """
    python3 ${projectDir}/bin/predict_tertiary.py --entry-id ${entry_id} --sequence ${sequence} --molecule-type ${molecule_type} \\
        --rhofold-dir ${params.rhofold_dir} --ckpt ${params.rhofold_ckpt} \\
        --outdir ${entry_id}_3d --device ${params.tertiary_device} \\
        --relax-steps ${params.relax_steps} --max-length ${params.tertiary_max_length} \\
        --python ${params.rhofold_python}
    """

    stub:
    def flag = molecule_type == 'DNA' ? 'rna_analogue_of_dna' : 'ok'
    """
    mkdir -p ${entry_id}_3d
    touch ${entry_id}_3d/${entry_id}_3d.pdb
    printf 'entry_id,model_sequence,tertiary_tool,tertiary_flag,pdb_file,atom_count,tertiary_status\\n${entry_id},${sequence},RhoFold+,${flag},${entry_id}_3d.pdb,,stub\\n' \\
        > ${entry_id}_3d/${entry_id}_3d_summary.csv
    """
}
