""" This module defined classes in order to group certain parameters together. """
from __future__ import annotations
from typing import TYPE_CHECKING
import dataclasses
from moPepGen import aa, constant


if TYPE_CHECKING:
    from typing import Set, Dict, List, Tuple
    from moPepGen import dna, gtf


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
            min_mw:int=500, min_length:int=None, max_length:int=None,
            max_variants_per_node:int=7, additional_variants_per_misc:int=2,
            in_bubble_cap_step_down:int=0, min_nodes_to_collapse:int=30,
            naa_to_collapse:int=5, flanking_size:int=9, peptide_finding_mode:str='misc'):
        """ constructor """
        self.enzyme = enzyme
        if self.enzyme and self.enzyme.lower() == 'none':
            self.enzyme = None
        if self.enzyme is None and peptide_finding_mode == constant.PeptideFindingMode.MISC.value:
            self.peptide_finding_mode = constant.PeptideFindingMode.SLIDING_WINDOW.value
        self.exception = exception
        self.miscleavage = miscleavage
        self.min_mw = min_mw

        min_length_default, max_length_default = self.get_default_peptide_lengths(
            peptide_finding_mode
        )
        if min_length is None:
            min_length = min_length_default
        if max_length is None:
            max_length = max_length_default
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

    def get_default_peptide_lengths(self, peptide_finding_mode) -> Tuple[int, int]:
        """ Get default peptide lengths based on mode """
        if peptide_finding_mode == constant.PeptideFindingMode.MISC.value:
            return (7, 25)
        else:
            return (8, 11)

    def jsonfy(self, graph_params:bool=False):
        """ jsonfy """
        data = {
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
