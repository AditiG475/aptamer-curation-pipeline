#!/usr/bin/env nextflow
/*
 * Aptamer curation & structure workflow
 *
 *   input (FASTA or CSV)
 *     -> CLEAN_SEQUENCES        whitespace / case / IUPAC clean-up, duplicates
 *     -> DETECT_MOLECULE_TYPE   DNA vs RNA from sequence evidence (+ metadata)
 *     -> FOLD_SECONDARY  ┐      ViennaRNA, molecule-specific parameters
 *     -> SHUFFLE_ZSCORE  ┘      fold stability vs shuffled sequences (parallel)
 *     -> PREDICT_TERTIARY       RhoFold+ 3D models, one task per sequence (optional)
 *     -> VERIFY_STRUCTURES      consistency checks + warnings
 *     -> annotation branch      Kd, DOI, PubChem, classification (CSV input, optional)
 *     -> MERGE_RESULTS -> EXPORT_EXCEL
 */
nextflow.enable.dsl = 2

include { CLEAN_SEQUENCES      } from './modules/clean_sequences'
include { DETECT_MOLECULE_TYPE } from './modules/detect_molecule_type'
include { FOLD_SECONDARY       } from './modules/fold_secondary'
include { SHUFFLE_ZSCORE       } from './modules/shuffle_zscore'
include { PREDICT_TERTIARY     } from './modules/predict_tertiary'
include { VERIFY_STRUCTURES    } from './modules/verify_structures'
include { EXTRACT_BINDING; DOI_LOOKUP; PUBCHEM_ANNOTATE; CLASSIFY } from './modules/annotate'
include { MERGE_RESULTS; EXPORT_EXCEL } from './modules/report'

def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf --input <aptamers.fasta|entries.csv> [options]

    Main options:
      --outdir              results folder                       [${params.outdir}]
      --temperature         folding temperature, Celsius         [${params.temperature}]
      --shuffles            shuffles per sequence for Z-score     [${params.shuffles}]
      --default_type        DNA/RNA when a sequence has no T or U [${params.default_type}]
      --run_3d              predict 3D structures with RhoFold+   [${params.run_3d}]
      --rhofold_dir/--rhofold_ckpt  RhoFold+ install + checkpoint
      --skip_annotation     skip Kd / DOI / PubChem / classification [${params.skip_annotation}]
      --strict              stop if a QC check fails              [${params.strict}]
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); return }
    if (!params.input) { helpMessage(); error "Please provide --input" }

    input_file = file(params.input, checkIfExists: true)
    is_csv     = input_file.name.toLowerCase().endsWith('.csv')

    // 1-2: clean and type the sequences
    cleaned = CLEAN_SEQUENCES(input_file).csv
    typed   = DETECT_MOLECULE_TYPE(cleaned).csv

    // 3-4: secondary structure and Z-score run in parallel
    secondary = FOLD_SECONDARY(typed).csv
    zscore    = SHUFFLE_ZSCORE(typed).csv

    // 5: tertiary structure, one task per sequence
    if (params.run_3d) {
        if (!params.rhofold_dir || !params.rhofold_ckpt) error "--run_3d needs --rhofold_dir and --rhofold_ckpt"
        per_sequence = typed
            .splitCsv(header: true, quote: '"')
            .filter { row -> row.clean_status == 'ok' && !row.duplicate_of }
            .map    { row -> tuple(row.entry_id, row.fold_sequence, row.detected_type) }
        tertiary = PREDICT_TERTIARY(per_sequence).summary
            .collectFile(name: 'tertiary_summary.csv', keepHeader: true, skip: 1,
                         storeDir: "${params.outdir}/05_tertiary_structure")
    } else {
        tertiary = Channel.value(file("${projectDir}/assets/NO_FILE_tertiary"))
    }

    // 6: QC
    VERIFY_STRUCTURES(typed, secondary, zscore)

    // 7: literature / compound annotation (CSV input only)
    if (is_csv && !params.skip_annotation) {
        binding    = EXTRACT_BINDING(input_file).csv
        doi        = DOI_LOOKUP(binding).csv
        pubchem    = PUBCHEM_ANNOTATE(doi).csv
        annotation = CLASSIFY(pubchem, file(params.rules, checkIfExists: true)).csv
    } else {
        annotation = Channel.value(file("${projectDir}/assets/NO_FILE_annotation"))
    }

    // 8: merge and export
    final_csv = MERGE_RESULTS(typed, secondary, zscore, tertiary, annotation).csv
    EXPORT_EXCEL(final_csv)
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nDone. Results in: ${params.outdir}  (final_dataset.xlsx, 06_qc/qc_report.txt)"
        : "\nWorkflow failed - see .nextflow.log")
}
