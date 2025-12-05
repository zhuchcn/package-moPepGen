"""
Brute force variant peptide caller module.

This module provides a brute force algorithm for calling variant peptides from GVF files.
It has been refactored from a monolithic script into a modular structure for better
maintainability and to support multiple peptide finding modes.
"""

# Import all public APIs for backward compatibility
from .orf_tracker import ORFTracker
from .cli import parse_args, main
from .utils import create_mnvs, fix_indel_after_start_codon
from .caller import BruteForceVariantPeptideCaller
