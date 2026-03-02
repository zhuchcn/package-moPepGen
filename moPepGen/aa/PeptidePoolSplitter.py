""" module for peptide pool splitter """
from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path
import itertools
import copy
from moPepGen.seqvar import GVFMetadata
from moPepGen import seqvar, circ, constant, VARIANT_PEPTIDE_SOURCE_DELIMITER, \
    SPLIT_DATABASE_KEY_SEPARATER
from .VariantPeptidePool import VariantPeptidePool
from .VariantPeptideLabel import VariantPeptideInfo, VariantSourceSet, \
    LabelSourceMapping

if TYPE_CHECKING:
    from typing import Dict, IO, List, Set, Union, FrozenSet
    from .AminoAcidSeqRecord import AminoAcidSeqRecord
    Databases = Dict[str,VariantPeptidePool]

class PeptidePoolSplitter():
    """ PeptidePoolSplitter to split peptide pool into separate databases.

    Attributes:
        peptides (VariantPeptidePool): The variant peptide pool.
        order (Dict[str,str]): Order of sources.
        label_map (Dict[str,Dict[str,str]]): This a nested dict of transcripts,
            variant lables, and their corresponding source type. This is to
            quickly retrieve the source type of a given variant label.
        group_map (Dict[str,str]): A map of certain variant sources being
            grouped together.
        databases (Dict[str,VariantPeptidePool]): This is the splitted
            databases. It is a dict that keys are the string representation of
            variant sources, and values are the splitted VariantPeptidePool.
        sources (Set[str]): Variant sources.
    """
    def __init__(self, peptides:VariantPeptidePool=None,
            order:Dict[Union[str,FrozenSet[str]],int]=None,
            label_map:LabelSourceMapping=None, group_map:Dict[str,str]=None,
            databases:Databases=None, sources:Set[str]=None):
        self.peptides = peptides
        self.databases = databases or {}
        self.label_map = label_map or LabelSourceMapping()
        self.group_map = group_map or {}
        self.order = order or {}
        if sources:
            self.sources = sources
        else:
            self.sources = set()
            for source_group in self.order:
                if isinstance(sources, str):
                    self.sources.add(source_group)
                else:
                    self.sources.update([s for s in source_group if s not in ['+', '*']])

    def get_reversed_group_map(self) -> Dict[str, List[str]]:
        """ Reverse group map """
        group_map: Dict[str, List[str]] = {}
        for k,v in self.group_map.items():
            if v in group_map:
                group_map[v].append(k)
            else:
                group_map[v] = [k]
        return group_map

    def append_order_internal_sources(self):
        """ Add internal sources that are not present in any GTFs, including
        novel ORF, sec termination, and codon reassignment. """
        sources = [
            constant.SOURCE_NOVEL_ORF,
            constant.SOURCE_SEC_TERMINATION,
            constant.SOURCE_CODON_REASSIGNMENT
        ]
        for source in sources:
            if source in self.group_map:
                source = self.group_map[source]
            if source not in self.order:
                self.append_order(source)

    def append_order(self, source:str):
        """ Add a group to the end of the order. """
        if source in self.order:
            return
        if source in self.group_map:
            source = self.group_map[source]
        if source in self.order:
            return
        self.order[source] = max(self.order.values()) + 1 if self.order else 0
        self.sources.add(source)

    def load_database(self, handle:IO) -> None:
        """ load peptide database """
        pool = VariantPeptidePool.load(handle)
        if not self.peptides:
            self.peptides = pool
            return
        for peptide in pool.peptides:
            self.peptides.add_peptide(peptide, None, skip_checking=True)

    def load_gvf(self, handle:IO) -> None:
        """ Load variant lables from GVF file. """
        metadata = GVFMetadata.parse(handle)
        source = metadata.source

        if source not in self.order:
            self.append_order(source)

        if metadata.parser == 'parseCIRCexplorer':
            for record in circ.io.parse(handle):
                self.label_map.add_circ_rna(record, source)
        else:
            for variant in seqvar.io.parse(handle):
                self.label_map.add_variant(variant, source)

    def add_peptide_to_database(self, database_key:str,
            peptide:AminoAcidSeqRecord) -> None:
        """ Add peptide to database """
        if database_key not in self.databases:
            self.databases[database_key] = VariantPeptidePool()
        self.databases[database_key].peptides.add(peptide)

    def get_all_peptide_sources(self, tx2gene:Dict[str,str], coding_tx:Set[str]):
        """ Get sources of all peptides """
        peptide_sources = []
        for peptide in self.peptides.peptides:
            peptide_infos = VariantPeptideInfo.from_variant_peptide(
                peptide=peptide, tx2gene=tx2gene, coding_tx=coding_tx,
                label_map=self.label_map
            )
            peptide_infos.sort()
            peptide_sources.append([x.sources for x in peptide_infos])
        return peptide_sources

    def create_wildcard_map(self):
        """ Create a wildcard map """
        wildcard_map = {}
        wildcard_chrs = ['+', '*']
        for sources in sorted(self.order, key=lambda x: self.order[x]):
            if isinstance(sources, str):
                sources = frozenset([sources])
            if not any(x in wildcard_chrs for x in sources):
                if sources not in wildcard_map:
                    wildcard_map[sources] = sources
                continue
            individual_sources = [x for x in self.sources if x not in sources]
            if all(x in sources for x in wildcard_chrs):
                raise ValueError(f"Invalid wildcard in souce order: {sources}")
            start = 0 if '*' in sources else 1
            for i in range(start, len(individual_sources)):
                for extra_sources in itertools.combinations(individual_sources, i):
                    expanded_sources = [x for x in sources if x not in wildcard_chrs] \
                        + list(extra_sources)
                    expanded_sources = frozenset(expanded_sources)
                    if expanded_sources not in wildcard_map:
                        wildcard_map[expanded_sources] = sources
        return wildcard_map

    def split(self, max_groups:int, additional_split:List[Set],
            tx2gene:Dict[str,str], coding_tx:Set[str],
            split_mode:str='peptide'):
        """ Split peptide pool into separate databases """
        if split_mode not in {'peptide', 'entry'}:
            raise ValueError(
                f"Unsupported split mode: {split_mode}. Use 'peptide' or 'entry'."
            )
        self.append_order_internal_sources()
        VariantSourceSet.set_levels(self.order)
        delimiter = VARIANT_PEPTIDE_SOURCE_DELIMITER
        additional_split = [VariantSourceSet(x) for x in additional_split]
        wildcard_map = self.create_wildcard_map()
        for peptide in self.peptides.peptides:
            peptide_infos = VariantPeptideInfo.from_variant_peptide(
                peptide=peptide,
                tx2gene=tx2gene,
                coding_tx=coding_tx,
                label_map=self.label_map,
                group_map=self.group_map,
                wildcard_map=wildcard_map
            )
            peptide_infos.sort()
            if split_mode == 'entry':
                self.add_peptide_to_database_entry_mode(
                    peptide=peptide,
                    peptide_infos=peptide_infos,
                    max_groups=max_groups,
                    additional_split=additional_split
                )
            else:
                self.add_peptide_to_database_peptide_mode(
                    peptide=peptide,
                    peptide_infos=peptide_infos,
                    delimiter=delimiter,
                    max_groups=max_groups,
                    additional_split=additional_split
                )

    def add_peptide_to_database_peptide_mode(self, peptide:AminoAcidSeqRecord,
            peptide_infos:List[VariantPeptideInfo], delimiter:str, max_groups:int,
            additional_split:List[VariantSourceSet]) -> None:
        """Assign one peptide sequence to exactly one database tier."""
        peptide.description = delimiter.join([str(x) for x in peptide_infos])
        peptide.id = peptide.description
        peptide.name = peptide.description

        database_key = self.get_database_key(
            sources=peptide_infos[0].sources,
            max_groups=max_groups,
            additional_split=additional_split
        )
        self.add_peptide_to_database(database_key, peptide)

    def add_peptide_to_database_entry_mode(self, peptide:AminoAcidSeqRecord,
            peptide_infos:List[VariantPeptideInfo], max_groups:int,
            additional_split:List[VariantSourceSet]) -> None:
        """Assign each header entry of a peptide sequence to its own tier."""
        for peptide_info in peptide_infos:
            database_key = self.get_database_key(
                sources=peptide_info.sources,
                max_groups=max_groups,
                additional_split=additional_split
            )
            entry_peptide = copy.copy(peptide)
            entry_peptide.description = str(peptide_info)
            entry_peptide.id = entry_peptide.description
            entry_peptide.name = entry_peptide.description
            if database_key not in self.databases:
                self.databases[database_key] = VariantPeptidePool()
            # Keep all entry labels for identical peptide sequences.
            self.databases[database_key].add_peptide(
                entry_peptide, None, skip_checking=True
            )

    def get_database_key(self, sources:VariantSourceSet, max_groups:int,
            additional_split:List[VariantSourceSet]) -> str:
        """Resolve output database key for a source set."""
        if len(sources) <= max_groups:
            return str(sources)
        for additional_set in additional_split:
            if additional_set.issubset(sources):
                return self.get_additional_database_key(additional_set)
        return self.get_remaining_database_key()

    @staticmethod
    def get_additional_database_key(additional:str) -> str:
        """ Get the database key for priority """
        return f'{str(additional)}{SPLIT_DATABASE_KEY_SEPARATER}additional'

    @staticmethod
    def get_remaining_database_key() -> str:
        """ Get the database key for remanent peptides """
        return 'Remaining'

    @staticmethod
    def get_database_filename(database_key:str, prefix:str) -> str:
        """ Get the database filename """
        return f"{prefix}_{database_key}.fasta"

    def write(self, output_prefix:str) -> None:
        """ Write databases to disk """
        output_dir = Path(output_prefix).parent
        output_dir.mkdir(exist_ok=True)
        filename_prefix = str(Path(output_prefix).name)
        for database_key, database in self.databases.items():
            filename = self.get_database_filename(database_key, filename_prefix)
            database.write(output_dir/filename)


def get_index_contains_source_set(peptide_infos:List[VariantPeptideInfo],
        target_set:VariantSourceSet) -> List[int]:
    """ Find the indices of a list of source sets,  """
    for i, info in enumerate(peptide_infos):
        if target_set.issubset(info.sources):
            return i
    return None
