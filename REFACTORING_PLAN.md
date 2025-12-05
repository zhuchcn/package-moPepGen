# Refactoring Plan: brute_force.py and fuzz_test.py

## Implementation Status

✅ **brute_force.py** → **brute_force/** module (IN PROGRESS)
- ✅ Renamed old file to `brute_force_legacy.py`
- ✅ Created `brute_force/` directory structure
- 🔄 Extracting classes into separate modules

## Executive Summary

Both modules need refactoring to:
1. Add support for the new `peptide_finding_mode` (misc, sliding_window, archipel)
2. Reduce complexity by extracting focused classes
3. Improve maintainability and testability
4. Ensure consistency between the two modules

---

## Part 1: brute_force.py Refactoring

### Current State
- **Total lines**: 1403
- **Main class**: `BruteForceVariantPeptideCaller` (1125 lines - god class)
- **Problem**: No peptide_finding_mode support, monolithic design, complex nested logic

### Proposed Structure

```
brute_force.py (main module)
├── ORFTracker (keep as-is, 43 lines) ✓
├── VariantSequenceBuilder (NEW, ~250 lines)
├── VariantEffectAnalyzer (NEW, ~150 lines)
├── VariantCompatibilityValidator (NEW, ~200 lines)
├── PeptideGenerator (NEW, ~300 lines) - MODE AWARE
└── BruteForceVariantPeptideCaller (refactored, ~200 lines)
```

### 1.1 Extract VariantSequenceBuilder

**Purpose**: Handle all variant sequence generation logic

**Methods to extract**:
- `get_variant_sequence()` - Apply variants to sequence
- `get_variant_sequence_fusion()` - Handle fusion variants
- `get_variant_sequence_circ_rna()` - Handle circRNA variants
- `get_variant_ref_seq()` - Get reference sequence
- `get_sec_positions()` - Track selenocysteine positions
- `get_gene_seq()` - Cache gene sequence

**Constructor**:
```python
def __init__(self, reference_data, tx_model, tx_seq, variant_pool):
    self.reference_data = reference_data
    self.tx_model = tx_model
    self.tx_seq = tx_seq
    self.variant_pool = variant_pool
    self.gene_seq = None  # cached
```

### 1.2 Extract VariantEffectAnalyzer

**Purpose**: Analyze variant effects on translation

**Methods to extract**:
- `check_variant_effect()` - Check stop_lost, stop_gain, silent_mutation
- `is_stop_lost()` - Check if variant causes stop loss
- `has_any_stop_codon_between()` - Check for stop codons in range

**Constructor**:
```python
def __init__(self, tx_seq, tx_model):
    self.tx_seq = tx_seq
    self.tx_model = tx_model
```

### 1.3 Extract VariantCompatibilityValidator

**Purpose**: Validate variant combinations

**Methods to extract**:
- `has_incompatible_variants()` - Check variant compatibility
- `has_overlapping_variants()` - Check for overlaps (static)
- `has_any_invalid_variants_on_inserted_sequences()` - Validate fusion/splicing
- `has_any_invalid_variant_on_circ()` - Validate circRNA (static)
- `should_clip_trailing_nodes()` - Check if nodes should be clipped

**Constructor**:
```python
def __init__(self, reference_data, tx_model, tx_seq, variant_pool, tx_id):
    self.reference_data = reference_data
    self.tx_model = tx_model
    self.tx_seq = tx_seq
    self.variant_pool = variant_pool
    self.tx_id = tx_id
```

### 1.4 Create PeptideGenerator (NEW - MODE AWARE)

**Purpose**: Generate peptides based on peptide_finding_mode

**Methods**:
```python
class PeptideGenerator:
    def __init__(self, cleavage_params, mode='misc'):
        self.cleavage_params = cleavage_params
        self.mode = mode

    def generate_peptide_candidates(self, aa_seq, cds_start,
                                    variant_coordinates, stop_lost,
                                    stop_gain, silent_mutation):
        """Route to mode-specific generator"""
        if self.mode == 'misc':
            return self._generate_misc(...)
        elif self.mode == 'sliding_window':
            return self._generate_sliding_window(...)
        elif self.mode == 'archipel':
            return self._generate_archipel(...)

    def _generate_misc(self, aa_seq, cds_start, ...):
        """Current enzymatic cleavage logic"""
        # Existing logic from call_peptides_main
        sites = aa_seq.find_all_enzymatic_cleave_sites_with_ranges(
            rule=self.cleavage_params.enzyme,
            exception=self.cleavage_params.exception
        )
        # ... rest of existing logic

    def _generate_sliding_window(self, aa_seq, cds_start, ...):
        """NEW: Generate all 8-11mers for neoantigens"""
        for length in range(self.cleavage_params.min_length,
                           self.cleavage_params.max_length + 1):
            for i in range(len(aa_seq) - length + 1):
                peptide = aa_seq[i:i+length]
                tx_lhs = cds_start + i * 3
                tx_rhs = cds_start + (i + length) * 3
                # Check if contains variants
                effective_variants = self._get_effective_variants(
                    tx_lhs, tx_rhs, ...
                )
                if effective_variants:
                    yield (peptide, tx_lhs, tx_rhs, effective_variants)

    def _generate_archipel(self, aa_seq, cds_start, ...):
        """NEW: Generate variant islands with flanking regions"""
        # Find variant positions
        variant_positions = [v.location.start for v in variant_coordinates]

        # For each variant, create island with flanking
        for var_pos in variant_positions:
            var_pos_aa = (var_pos - cds_start) // 3
            start = max(0, var_pos_aa - self.cleavage_params.flanking_size)
            end = min(len(aa_seq), var_pos_aa + self.cleavage_params.flanking_size)

            peptide = aa_seq[start:end]
            tx_lhs = cds_start + start * 3
            tx_rhs = cds_start + end * 3

            if self._validate_length(peptide):
                yield (peptide, tx_lhs, tx_rhs, ...)
```

### 1.5 Refactor BruteForceVariantPeptideCaller

**New role**: Coordinator/orchestrator

**Composition**:
```python
class BruteForceVariantPeptideCaller:
    def __init__(self, reference_data=None, cleavage_params=None, ...):
        # ... existing init

        # New: inject dependencies
        self.sequence_builder = None  # created when needed
        self.effect_analyzer = None
        self.compatibility_validator = None
        self.peptide_generator = None

    def _initialize_helpers(self):
        """Lazy initialization of helper classes"""
        if not self.sequence_builder:
            self.sequence_builder = VariantSequenceBuilder(
                self.reference_data, self.tx_model,
                self.tx_seq, self.variant_pool
            )
        if not self.effect_analyzer:
            self.effect_analyzer = VariantEffectAnalyzer(
                self.tx_seq, self.tx_model
            )
        if not self.compatibility_validator:
            self.compatibility_validator = VariantCompatibilityValidator(
                self.reference_data, self.tx_model,
                self.tx_seq, self.variant_pool, self.tx_id
            )
        if not self.peptide_generator:
            mode = self.cleavage_params.peptide_finding_mode
            self.peptide_generator = PeptideGenerator(
                self.cleavage_params, mode
            )

    def call_peptides_main(self, variants, denylist, ...):
        """Simplified orchestration"""
        self._initialize_helpers()

        # Delegate to sequence builder
        if variants[self.tx_id].fusion:
            seq, variant_coordinates = self.sequence_builder.get_variant_sequence_fusion(...)
        elif variants[self.tx_id].circ_rna:
            seq, variant_coordinates = self.sequence_builder.get_variant_sequence_circ_rna(...)
        else:
            seq, variant_coordinates = self.sequence_builder.get_variant_sequence(...)

        # Delegate to effect analyzer
        stop_lost, stop_gain, silent_mutation = self.effect_analyzer.check_variant_effect(
            seq, variant_coordinates
        )

        # Delegate to peptide generator
        for peptide, tx_lhs, tx_rhs, effective_variants in \
                self.peptide_generator.generate_peptide_candidates(
                    aa_seq, cds_start, variant_coordinates,
                    stop_lost, stop_gain, silent_mutation
                ):
            # Apply translational modifications
            peptide_seqs = self.translational_modification(
                peptide, lhs, tx_lhs, effective_variants,
                denylist, check_variants, check_canonical,
                selenocysteine_termination
            )
            for peptide_seq in peptide_seqs:
                yield peptide_seq
```

### 1.6 Add peptide_finding_mode to parse_args and main

**parse_args modifications**:
```python
def parse_args(subparsers):
    # ... existing args

    parser.add_argument(
        '--peptide-finding-mode',
        type=str,
        default='misc',
        choices=['misc', 'sliding-window', 'archipel'],
        help='Peptide finding mode. "misc" for enzyme-based cleavage, '
             '"sliding-window" for 8-11mer enumeration (neoantigens), '
             '"archipel" for variant islands with flanking regions.'
    )

    parser.add_argument(
        '--flanking-size',
        type=int,
        default=9,
        help='Flanking size for archipel mode (default: 9).'
    )
```

**main() modifications**:
```python
def main(args):
    # ... existing reference loading

    caller.cleavage_params = params.CleavageParams(
        enzyme=args.cleavage_rule,
        exception=args.cleavage_exception,
        miscleavage=int(args.miscleavage),
        min_mw=float(args.min_mw) if args.min_mw else None,
        min_length=args.min_length if args.min_length else None,
        max_length=args.max_length if args.max_length else None,
        peptide_finding_mode=args.peptide_finding_mode,  # NEW
        flanking_size=args.flanking_size  # NEW
    )
```

---

## Part 2: fuzz_test.py Refactoring

### Current State
- **Total lines**: 879
- **Main classes**: FuzzRecord (200), FuzzTestCase (320), FuzzTestConfig (70), Fuzzer (80)
- **Problem**: No peptide_finding_mode support, hard-coded parameters, monolithic FuzzTestCase

### Proposed Structure

```
fuzz_test.py
├── FuzzRecordStatus (keep as-is) ✓
├── FuzzRecord (keep as-is) ✓
├── PeptideCallerArgsBuilder (NEW, ~100 lines) - MODE AWARE
├── FuzzTestCase (refactored, ~250 lines)
├── FuzzTestConfig (extended, ~100 lines)
├── Fuzzer (keep as-is) ✓
└── parse_args (extended)
```

### 2.1 Create PeptideCallerArgsBuilder (NEW)

**Purpose**: Centralize argument construction for call_variants and brute_force

```python
class PeptideCallerArgsBuilder:
    """Builder for constructing peptide caller arguments"""

    def __init__(self, config: FuzzTestConfig, record: FuzzRecord):
        self.config = config
        self.record = record

    def build_call_variant_args(self) -> argparse.Namespace:
        """Build args for callVariantPeptide based on peptide_finding_mode"""
        args = argparse.Namespace()
        args.index_dir = None
        args.command = 'callPeptides'
        args.input_path = [self.record.var_gvf_file]
        if self.record.circ_gvf_file.exists():
            args.input_path.append(self.record.circ_gvf_file)

        # Reference files
        args.genome_fasta = self.config.path_genome_fasta
        args.annotation_gtf = self.config.path_annotation_gtf
        args.proteome_fasta = self.config.path_proteome_fasta
        args.reference_source = None

        # Codon tables
        args.codon_table = 'Standard'
        args.chr_codon_table = []
        args.start_codons = ['ATG']
        args.chr_start_codons = []

        # Output
        args.output_path = self.record.call_variant_fasta
        if self.config.save_graph:
            args.graph_output_dir = self.record.work_dir/'graph'
            args.graph_output_dir.mkdir(exist_ok=True)
        else:
            args.graph_output_dir = None

        # Variant processing
        args.quiet = True
        args.backsplicing_only = False
        args.max_adjacent_as_mnv = 2
        args.selenocysteine_termination = True
        args.w2f_reassignment = True

        # Cleavage parameters - MODE AWARE
        args.cleavage_rule = self.config.cleavage_rule
        args.cleavage_exception = self.config.cleavage_exception
        args.miscleavage = self.config.miscleavage
        args.min_mw = self.config.min_mw
        args.min_length = self.config.min_length
        args.max_length = self.config.max_length
        args.peptide_finding_mode = self.config.peptide_finding_mode  # NEW
        args.flanking_size = self.config.flanking_size  # NEW

        # Graph parameters
        args.threads = 1
        args.max_variants_per_node = (9, )
        args.additional_variants_per_misc = (2, )
        args.min_nodes_to_collapse = 30
        args.naa_to_collapse = 5

        # Other settings
        args.noncanonical_transcripts = False
        args.invalid_protein_as_noncoding = False
        args.debug_level = 1
        args.timeout_seconds = 300
        args.coding_novel_orf = False
        args.skip_failed = False

        return args

    def build_brute_force_args(self) -> argparse.Namespace:
        """Build args for bruteForce"""
        args = argparse.Namespace()
        args.input_gvf = [self.record.var_gvf_file]
        if self.record.circ_gvf_file.exists():
            args.input_gvf.append(self.record.circ_gvf_file)

        args.reference_dir = self.config.ref_dir

        # Codon tables
        args.codon_table = 'Standard'
        args.chr_codon_table = []
        args.start_codons = ['ATG']
        args.chr_start_codons = []

        # Variant processing
        args.force = True
        args.variant_ids = []
        args.max_adjacent_as_mnv = 2
        args.selenocysteine_termination = True
        args.w2f_reassignment = True

        # Cleavage parameters - MODE AWARE (MUST MATCH call_variant_args)
        args.cleavage_rule = self.config.cleavage_rule
        args.cleavage_exception = self.config.cleavage_exception
        args.miscleavage = self.config.miscleavage
        args.min_mw = self.config.min_mw
        args.min_length = self.config.min_length
        args.max_length = self.config.max_length
        args.peptide_finding_mode = self.config.peptide_finding_mode  # NEW
        args.flanking_size = self.config.flanking_size  # NEW

        return args

    def validate_mode_compatibility(self):
        """Validate that mode and enzyme are compatible"""
        if self.config.peptide_finding_mode == 'sliding_window':
            if self.config.cleavage_rule is not None \
                    and self.config.cleavage_rule.lower() != 'none':
                raise ValueError(
                    "sliding_window mode requires cleavage_rule to be None or 'none'"
                )
```

### 2.2 Extend FuzzTestConfig

**Add new parameters**:
```python
class FuzzTestConfig:
    def __init__(self, tx_id, n_iter, max_size, max_variants, min_variants,
                 exonic_only, snv_frac, snv_only_frac, fusion_frac,
                 circ_rna_frac, ci_ratio, alt_splice_frac, keep_succeeded,
                 save_graph, cleavage_rule, cleavage_exception, miscleavage,
                 min_mw, min_length, max_length, temp_dir, output_dir, ref_dir,
                 peptide_finding_mode='misc',  # NEW
                 flanking_size=9,  # NEW
                 fuzz_start=None, fuzz_end=None, seed=None, nthreads=1):
        # ... existing assignments

        # NEW: peptide finding mode
        self.peptide_finding_mode = peptide_finding_mode
        self.flanking_size = flanking_size

        # Apply mode-aware defaults (similar to params.CleavageParams)
        self._apply_mode_defaults()

    def _apply_mode_defaults(self):
        """Apply default min/max length based on mode if not specified"""
        if self.min_length is None or self.max_length is None or self.min_mw is None:
            defaults = self._get_mode_defaults()
            if self.min_length is None:
                self.min_length = defaults['min_length']
            if self.max_length is None:
                self.max_length = defaults['max_length']
            if self.min_mw is None:
                self.min_mw = defaults['min_mw']

    def _get_mode_defaults(self):
        """Get default parameters based on peptide_finding_mode"""
        if self.peptide_finding_mode == 'misc':
            return {
                'min_mw': 500,
                'min_length': 7,
                'max_length': 25
            }
        elif self.peptide_finding_mode == 'sliding_window':
            return {
                'min_mw': 0,
                'min_length': 8,
                'max_length': 11
            }
        elif self.peptide_finding_mode == 'archipel':
            return {
                'min_mw': 0,
                'min_length': self.flanking_size * 2 + 1,
                'max_length': self.flanking_size * 2 + 10
            }
        else:
            raise ValueError(f"Unknown peptide_finding_mode: {self.peptide_finding_mode}")
```

### 2.3 Refactor FuzzTestCase

**Simplify call_variants() and brute_force()**:
```python
class FuzzTestCase:
    def __init__(self, config, record=None):
        self.config = config
        self.record = record or FuzzRecord.generate_record()
        if not record:
            self.record.set_work_dir(self.config.temp_dir)

        # NEW: create args builder
        self.args_builder = PeptideCallerArgsBuilder(self.config, self.record)

    def call_variants(self):
        """call variants using moPepGen's graph algorithm"""
        # Validate mode compatibility
        self.args_builder.validate_mode_compatibility()

        # Build args using the builder
        args = self.args_builder.build_call_variant_args()

        # Call the peptide caller
        call_variant_peptide(args=args)

    def brute_force(self):
        """call the brute force variant peptide caller"""
        # Build args using the builder
        args = self.args_builder.build_brute_force_args()

        # Call brute force
        with open(self.record.brute_force_fasta, 'wt') as handle:
            with redirect_stdout(handle):
                brute_force.main(args)
```

### 2.4 Extend parse_args

**Add new arguments**:
```python
def parse_args(subparsers):
    # ... existing args

    parser.add_argument(
        '--peptide-finding-mode',
        type=str,
        default='misc',
        choices=['misc', 'sliding-window', 'archipel'],
        help='Peptide finding mode. "misc" for enzyme-based cleavage '
             '(default), "sliding-window" for 8-11mer enumeration '
             '(neoantigens), "archipel" for variant islands with flanking.'
    )

    parser.add_argument(
        '--flanking-size',
        type=int,
        default=9,
        help='Flanking size for archipel mode. Default: 9.'
    )

    # ... existing cleavage args (add_args_cleavage)
```

### 2.5 Update main()

**Pass new parameters**:
```python
def main(args):
    config = FuzzTestConfig(
        # ... existing params
        cleavage_rule=args.cleavage_rule,
        cleavage_exception=args.cleavage_exception,
        miscleavage=args.miscleavage,
        min_mw=args.min_mw,
        min_length=args.min_length,
        max_length=args.max_length,
        peptide_finding_mode=args.peptide_finding_mode,  # NEW
        flanking_size=args.flanking_size,  # NEW
        # ... rest of params
    )
    with Fuzzer(config) as fuzzer:
        fuzzer.fuzz()
```

---

## Part 3: Implementation Order

### Phase 1: Add peptide_finding_mode support (no major refactoring)
**Goal**: Get moving-window working with minimal changes

1. **brute_force.py**:
   - Add `peptide_finding_mode` and `flanking_size` to parse_args
   - Add PeptideGenerator class with mode routing
   - Modify `call_peptides_main()` to use PeptideGenerator
   - Update main() to pass peptide_finding_mode to CleavageParams

2. **fuzz_test.py**:
   - Add `peptide_finding_mode` and `flanking_size` to parse_args
   - Add these to FuzzTestConfig
   - Update `call_variants()` to pass peptide_finding_mode
   - Update `brute_force()` to pass peptide_finding_mode
   - Validate mode in main()

**Estimated changes**: ~200 lines added/modified
**Risk**: Low - additive changes, backward compatible

### Phase 2: Extract helper classes from brute_force.py
**Goal**: Improve maintainability and testability

3. **brute_force.py**:
   - Extract VariantSequenceBuilder
   - Extract VariantEffectAnalyzer
   - Extract VariantCompatibilityValidator
   - Refactor BruteForceVariantPeptideCaller to use helpers

**Estimated changes**: ~400 lines moved/refactored
**Risk**: Medium - structural changes, needs thorough testing

### Phase 3: Extract PeptideCallerArgsBuilder from fuzz_test.py
**Goal**: Centralize parameter construction

4. **fuzz_test.py**:
   - Create PeptideCallerArgsBuilder
   - Update FuzzTestCase to use builder
   - Add mode validation

**Estimated changes**: ~150 lines moved/refactored
**Risk**: Low - isolated changes

---

## Part 4: Benefits

### Immediate Benefits (Phase 1)
- ✅ Support for sliding_window mode (neoantigen calling)
- ✅ Support for archipel mode (variant island detection)
- ✅ Consistent peptide_finding_mode handling across both modules
- ✅ Mode-aware parameter defaults

### Long-term Benefits (Phase 2 & 3)
- ✅ Reduced complexity (god class broken down)
- ✅ Better testability (isolated components)
- ✅ Easier to maintain and extend
- ✅ Clear separation of concerns
- ✅ Reusable components

---

## Part 5: Testing Strategy

### Test Cases for peptide_finding_mode

1. **misc mode (enzyme-based)**:
   - Test with trypsin (existing behavior)
   - Test with other enzymes
   - Validate miscleavage handling

2. **sliding_window mode**:
   - Test 8-11mer enumeration
   - Validate all kmers contain variants
   - Check that enzyme is ignored

3. **archipel mode**:
   - Test flanking region generation
   - Validate variant island detection
   - Test with different flanking_size values

### Fuzz Test Scenarios

1. Run fuzz test with all three modes
2. Mix different variant types (SNV, INDEL, fusion, circRNA, alt-splice)
3. Compare call_variant and brute_force outputs for each mode
4. Validate mode-specific edge cases

---

## Part 6: Backward Compatibility

### Ensuring Compatibility

1. **Default mode**: `misc` (existing behavior)
2. **Optional parameters**: All new parameters have defaults
3. **Existing tests**: Should pass without modification
4. **CLI**: New options are optional

### Migration Path

Users can gradually adopt new modes:
- Continue using `misc` mode (default)
- Try `sliding_window` for neoantigen discovery
- Experiment with `archipel` for variant-focused analysis

---

## Summary

This refactoring plan provides:

1. **Minimal initial changes** to add peptide_finding_mode support (Phase 1)
2. **Structural improvements** for long-term maintainability (Phase 2 & 3)
3. **Clear separation** between different peptide finding strategies
4. **Consistency** between fuzz_test and brute_force modules
5. **Backward compatibility** with existing functionality

The phased approach allows you to:
- Get moving-window support quickly (Phase 1)
- Optionally improve structure later (Phase 2 & 3)
- Test thoroughly at each phase
