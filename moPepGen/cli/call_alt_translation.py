""" `callAltTranslation` calls peptide sequences from coding transcripts that
harbor any alternative translation event. """
from __future__ import annotations
import argparse
import os
from pathlib import Path
from contextlib import ExitStack
from moPepGen import constant, params, aa, get_logger, svgraph
from moPepGen.cli import common
from moPepGen.pipeline.call_alt_translation_worker import call_alt_translation_for_transcript
from moPepGen.svgraph.VariantPeptideTable import (
    VariantPeptideTable,
    get_peptide_table_path,
    get_peptide_table_path_temp
)

OUTPUT_FILE_FORMATS = ['.fa', '.fasta']

# pylint: disable=W0212
def add_subparser_call_alt_translation(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callAltTranslation """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callAltTranslation',
        help='Call peptides with alternative translation from coding transcripts.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-o', '--output-path',
        type=Path,
        help='Output path to the alternative translation peptide FASTA.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--w2f-reassignment',
        action='store_true',
        help='Include peptides with W > F (Tryptophan to Phenylalanine) '
        'reassignment.'
    )
    p.add_argument(
        '--selenocysteine-termination',
        action='store_true',
        help='Include peptides of selenoprotiens that the UGA is treated as '
        'termination instead of Sec.'
    )
    p.add_argument(
        '--peptide-finding-mode',
        type=str,
        help='Peptide finding mode to use when calling peptides.',
        choices=[
            mode.value.replace('_', '-') for mode in constant.PeptideFindingMode
            if mode is not constant.PeptideFindingMode.ARCHIPEL
        ],
        default=constant.PeptideFindingMode.MISC.value,
    )

    common.add_args_reference(p)
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)

    p.set_defaults(func=call_alt_translation)
    common.print_help_if_missing_args(p)
    return p


def call_alt_translation(args: argparse.Namespace) -> None:
    """Main entrypoint for calling alternative translation peptides."""
    logger = get_logger()

    common.validate_file_format(
        args.output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )

    # Archipel mode is not supported for alt translation peptide calling
    supported_modes = [
        mode.value.replace('_', '-') for mode in constant.PeptideFindingMode
        if mode is not constant.PeptideFindingMode.ARCHIPEL
    ]
    if args.peptide_finding_mode not in supported_modes:
        raise ValueError(
            f"Peptide finding mode '{args.peptide_finding_mode}' is not supported "
            "for alt translation peptide calling."
        )

    cp = params.CleavageParams(
        enzyme=args.cleavage_rule,
        exception=args.cleavage_exception,
        miscleavage=args.miscleavage,
        min_mw=args.min_mw,
        min_length=args.min_length,
        max_length=args.max_length,
        peptide_finding_mode=args.peptide_finding_mode,
        flanking_size=args.flanking_size
    )

    if not (args.selenocysteine_termination or args.w2f_reassignment):
        raise ValueError(
            'At least one of --selenocysteine-termination and --w2f-reassignment'
            ' must be given.'
        )

    common.print_start_message(args)

    # Load references in CLI layer
    ref_data = common.load_references(
        args=args, load_proteome=True, cleavage_params=cp,
        load_codon_tables=True
    )

    peptide_table_temp_path = get_peptide_table_path_temp(args.output_path)
    peptide_table_output_path = get_peptide_table_path(args.output_path)

    with ExitStack() as stack:
        # Open temporary file for peptide table
        seq_anno_handle = stack.enter_context(open(peptide_table_temp_path, 'w+'))
        peptide_table = VariantPeptideTable(seq_anno_handle)
        peptide_table.write_header()

        # Process each coding transcript
        for tx_id in ref_data.anno.transcripts:
            tx_model = ref_data.anno.transcripts[tx_id]
            if not tx_model.is_protein_coding:
                continue

            codon_table = ref_data.codon_tables[tx_model.transcript.chrom]

            try:
                peptide_anno = call_alt_translation_for_transcript(
                    tx_id=tx_id,
                    tx_model=tx_model,
                    genome=ref_data.genome,
                    anno=ref_data.anno,
                    codon_table=codon_table,
                    cleavage_params=cp,
                    w2f_reassignment=args.w2f_reassignment,
                    sec_truncation=args.selenocysteine_termination
                )
            except:
                logger.error('Exception raised from %s', tx_id)
                raise

            # Write peptides to table immediately
            for peptide in peptide_anno:
                is_valid = peptide_table.is_valid(
                    seq=peptide,
                    canonical_peptides=ref_data.canonical_peptides,
                    cleavage_params=cp
                )
                if is_valid:
                    for seq_anno in peptide_anno[peptide]:
                        peptide_table.add_peptide(peptide, seq_anno)

        # Write FASTA from table
        logger.info('Writing peptides to FASTA...')
        peptide_table.write_fasta(args.output_path)
        logger.info('Alternative translation peptide FASTA written.')
        peptide_table.sort_table(peptide_table_output_path)
        logger.info('Alternative translation peptide table sorted.')

    os.remove(peptide_table_temp_path)

    logger.info('Alternative translation peptide FASTA file written to disk.')

