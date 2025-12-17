""" This module defined classes in order to group certain parameters together. """
from __future__ import annotations
from typing import TYPE_CHECKING
import dataclasses
from moPepGen import aa, constant, get_logger


if TYPE_CHECKING:
    from typing import Set, Dict, List
    from moPepGen import dna, gtf


@dataclasses.dataclass
class DefaultPeptideParams:
    """ Default peptide parameters for different modes. """
    min_mw: int
    min_length: int
    max_length: int

class CleavageParams():
    """ Cleavage related parameters.

    ## Attributes:
        - enzyme (str): Enzyme name.
        - exception (str): Enzymatic cleavage exception.
        - miscleavage (int): Number of cleavages to allow per non-canonical peptide.
        - min_mw (int): The minimal molecular weight of the non-canonical peptides.
        - min_length (int): The minimal length of non-canonical peptides, inclusive.
        - max_length (int): The maximum length of non-canonical peptides, inclusive.
        - max_variants_per_node (int): Maximal number of variants per node. This
            can be useful when there are local regions that are heavily mutated.
            When creating the cleavage graph, nodes containing variants larger
            than this value are skipped. Setting to -1 will avoid checking for
            this.
        - additional_variants_per_misc (int): Additional variants allowed for
            every miscleavage. This argument is used together with
            max_variants_per_node to handle hypermutated regions. Setting to -1
            will avoid checking for this.
        - min_nodes_to_collapse (int): When making the cleavage graph, the minimal
            number of nodes to trigger pop collapse.
        - naa_to_collapse (int): The number of bases used for pop collapse.
    """
    def __init__(self, enzyme:str=None, exception:str=None, miscleavage:int=2,
            min_mw:int=None, min_length:int=None, max_length:int=None,
            max_variants_per_node:int=7, additional_variants_per_misc:int=2,
            in_bubble_cap_step_down:int=0, min_nodes_to_collapse:int=30,
            naa_to_collapse:int=5, flanking_size:int=10, peptide_finding_mode:str='misc'):
        """ constructor """
        self.enzyme = enzyme
        if self.enzyme and self.enzyme.lower() == 'none':
            self.enzyme = None
        peptide_finding_mode = peptide_finding_mode.lower().replace('-', '_')
        logger = get_logger()
        if self.enzyme is None and peptide_finding_mode == constant.PeptideFindingMode.MISC.value:
            logger.warning(
                "Cleavage enzyme is not specified, but peptide finding mode is 'misc'. "
                "Setting peptide finding mode to 'sliding_window'."
            )
            peptide_finding_mode = constant.PeptideFindingMode.SLIDING_WINDOW.value
        elif self.enzyme is not None \
                and peptide_finding_mode != constant.PeptideFindingMode.MISC.value:
            logger.warning(
                "Cleavage enzyme is specified as '%s', but peptide finding "
                "mode is '%s'. Setting peptide finding mode to 'misc'.",
                self.enzyme, peptide_finding_mode
            )
            peptide_finding_mode = constant.PeptideFindingMode.MISC.value

        self.peptide_finding_mode = peptide_finding_mode
        self.exception = exception
        self.miscleavage = int(miscleavage)

        default_params = self.get_default_peptide_params(
            peptide_finding_mode=peptide_finding_mode,
            flanking_size=flanking_size
        )
        if min_mw is None:
            logger.info(
                "Using default min_mw = %i for peptide finding mode '%s'.",
                default_params.min_mw, peptide_finding_mode
            )
            min_mw = default_params.min_mw
        else:
            min_mw = float(min_mw)

        if min_length is None:
            logger.info(
                "Using default min_length = %i for peptide finding mode '%s'.",
                default_params.min_length, peptide_finding_mode
            )
            min_length = default_params.min_length
        else:
            min_length = int(min_length)

        if max_length is None:
            logger.info(
                "Using default max_length = %i for peptide finding mode '%s'.",
                default_params.max_length, peptide_finding_mode
            )
            max_length = default_params.max_length
        else:
            max_length = int(max_length)

        self.min_mw = min_mw
        self.min_length = min_length
        self.max_length = max_length
        self.max_variants_per_node = max_variants_per_node
        self.additional_variants_per_misc = additional_variants_per_misc
        self.in_bubble_cap_step_down = in_bubble_cap_step_down
        self.min_nodes_to_collapse = min_nodes_to_collapse
        self.naa_to_collapse = naa_to_collapse
        self.flanking_size = flanking_size
        if self.exception == 'auto':
            if enzyme == 'trypsin':
                self.exception = 'trypsin_exception'
            else:
                self.exception = None

    def get_default_peptide_params(self, peptide_finding_mode:str, flanking_size:int
            ) -> DefaultPeptideParams:
        """ Get default peptide lengths based on mode """
        if peptide_finding_mode == constant.PeptideFindingMode.MISC.value:
            min_mw = 500
            min_length = 7
            max_length = 25
        elif peptide_finding_mode == constant.PeptideFindingMode.SLIDING_WINDOW.value:
            min_mw = 0
            min_length = 8
            max_length = 11
        else:
            min_mw = 0
            min_length = flanking_size
            max_length = float("inf")

        return DefaultPeptideParams(
            min_mw=min_mw,
            min_length=min_length,
            max_length=max_length
        )

    def jsonfy(self, graph_params:bool=False):
        """ jsonfy """
        data = {
            'peptide_finding_mode': self.peptide_finding_mode,
            'enzyme': self.enzyme,
            'exception': self.exception,
            'miscleavage': self.miscleavage,
            'min_mw': self.min_mw,
            'min_length': self.min_length,
            'max_length': self.max_length,
            'flanking_size': self.flanking_size
        }
        if graph_params:
            data.update({
                'max_variants_per_node': self.max_variants_per_node,
                'additional_variants_per_misc': self.additional_variants_per_misc,
                'min_nodes_to_collapse': self.min_nodes_to_collapse,
                'naa_to_collapse': self.naa_to_collapse
            })
        return data

@dataclasses.dataclass
class ReferenceData:
    """ Reference related parameters

    ## Attributes
        - genome (dna.DNASeqDict)
        - anno (GenomicAnnotation)
        - canonical_peptides (Set[str])
        - proteome (aa.AminoAcidSeqDict)
        - codon_tables (Dict[str, CodonTableInfo])
    """
    genome: dna.DNASeqDict
    anno: gtf.GenomicAnnotation
    canonical_peptides: Set[str] = dataclasses.field(default_factory=set)
    proteome: aa.AminoAcidSeqDict = dataclasses.field(default_factory=aa.AminoAcidSeqDict)
    codon_tables: Dict[str, CodonTableInfo] = dataclasses.field(default_factory=dict)

@dataclasses.dataclass
class CodonTableInfo:
    """ Codon table info """
    codon_table: str
    start_codons: List[str] = dataclasses.field(default_factory=list)
