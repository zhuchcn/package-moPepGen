# Brute Force Variant Peptide Caller Module

## Overview

This module provides a brute force algorithm for calling variant peptides from GVF files. It was refactored from a monolithic 1403-line script into a modular structure.

## Module Structure

### `__init__.py`
Main module interface that exports all public APIs for backward compatibility.

**Exports:**
- `ORFTracker` - ORF tracking across reading frames
- `parse_args` - CLI argument parsing
- `main` - Main entry point
- `create_mnvs` - Create MNVs from adjacent variants
- `fix_indel_after_start_codon` - Fix indel positioning
- `BruteForceVariantPeptideCaller` - Main caller class

### `orf_tracker.py`
Contains the `ORFTracker` class for managing ORF positions.

**Class:** `ORFTracker`
- Tracks multiple ORF start positions in a transcript
- Organizes ORFs by reading frame (0, 1, 2)
- Provides methods to navigate between in-frame ORFs

### `cli.py`
CLI interface for the brute force algorithm.

**Functions:**
- `parse_args(subparsers)` - Parse command-line arguments
- `main(args)` - Main execution function

**Command-line arguments:**
- `-i, --input-gvf` - Input GVF file(s)
- `-r, --reference-dir` - Reference directory
- Cleavage parameters (enzyme, miscleavage, min/max length, etc.)
- Codon table options
- Translation options (selenocysteine, W2F reassignment)

### `utils.py`
Utility functions for variant processing.

**Functions:**
- `create_mnvs(pool, max_adjacent_as_mnv)` - Merge adjacent SNVs/INDELs into MNVs
- `fix_indel_after_start_codon(pool, ref_data)` - Adjust INDEL positioning for consistency

### `caller.py`
Main `BruteForceVariantPeptideCaller` class with ~1130 lines.

**Responsibilities:**
- Generate variant sequences (fusion, circRNA, standard)
- Analyze variant effects (stop loss/gain, silent mutations)
- Validate variant combinations
- Call variant peptides through combinatorial enumeration
- Handle translational modifications (selenocysteine, W2F)

**Key Methods:**
- `call_peptides()` - Main entry point
- `call_peptides_main()` - Core peptide calling logic
- `get_variant_sequence()` - Apply variants to sequence
- `get_variant_sequence_fusion()` - Handle fusion variants
- `get_variant_sequence_circ_rna()` - Handle circRNA
- `check_variant_effect()` - Analyze variant translation effects
- `has_incompatible_variants()` - Validate variant combinations
- `generate_variant_comb()` - Generate all variant combinations

## Usage

### As a Module

```python
from moPepGen.util import brute_force

# Access classes and functions
caller = brute_force.BruteForceVariantPeptideCaller()
tracker = brute_force.ORFTracker(orf_starts=[0, 3, 6])

# Use utility functions
pool = brute_force.create_mnvs(variant_pool, max_adjacent_as_mnv=2)
```

### As a CLI Tool

```bash
moPepGen-util bruteForce \
    -i variants.gvf \
    -r reference_dir/ \
    --cleavage-rule trypsin \
    --min-length 7 \
    --max-length 25
```

## Future Enhancements

### Planned: peptide_finding_mode Support

Add support for different peptide finding strategies:
- `misc` - Enzymatic cleavage (current behavior)
- `sliding_window` - 8-11mer enumeration for neoantigens
- `archipel` - Variant islands with flanking regions

### Planned: Further Modularization

Extract from `caller.py`:
- `sequence_builder.py` - Variant sequence generation (~300 lines)
- `effect_analyzer.py` - Variant effect analysis (~150 lines)
- `validator.py` - Variant compatibility validation (~200 lines)
- `peptide_generator.py` - Mode-aware peptide generation (~300 lines)

This would reduce `caller.py` to ~300 lines of orchestration logic.

## Testing

```bash
# Test imports
python3 -c "from moPepGen.util import brute_force; print('✓ Import successful')"

# Run module
moPepGen-util bruteForce --help

# Run tests (if available)
pytest test/util/test_brute_force.py
```

## Backward Compatibility

All existing imports continue to work:

```python
# These all work as before
from moPepGen.util import brute_force
from moPepGen.util.brute_force import BruteForceVariantPeptideCaller
from moPepGen.util.brute_force import ORFTracker, create_mnvs
```

## Migration from Legacy

The old monolithic `brute_force.py` was renamed to `brute_force_legacy.py` and can be safely deleted after confirming the new module works correctly.

## See Also

- `REFACTORING_SUMMARY.md` - Detailed refactoring implementation summary
- `../../../REFACTORING_PLAN.md` - Overall refactoring plan for brute_force and fuzz_test
