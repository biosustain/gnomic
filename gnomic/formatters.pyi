from abc import abstractmethod, ABCMeta
from typing import Dict, Tuple


class Formatter(metaclass=ABCMeta):
    format: str

    def format_genotype(self, genotype: 'gnomic.Genotype') -> str: ...

    @abstractmethod
    def format_change(self, change: 'gnomic.types.Change') -> str: ...

    def format_annotation(self, annotation: 'gnomic.types.Annotation') -> str: ...

    @staticmethod
    def format_accession(accession: 'gnomic.types.Accession') -> str: ...

    def format_variant(self, variant: Tuple[str]) -> str: ...

    def format_feature(self, feature: 'gnomic.types.Feature') -> str: ...

    def format_fusion(self, fusion: 'gnomic.types.Fusion') -> str: ...

    def format_plasmid(self, plasmid: 'gnomic.types.Plasmid') -> str: ...

    def format_at_locus(self, at_locus: 'gnomic.types.AtLocus') -> str: ...


class GnomicFormatter(Formatter):
    def format_change(self, change: 'gnomic.types.Change') -> str: ...


class TextFormatter(Formatter):
    def format_change(self, change: 'gnomic.types.Change') -> str: ...


class HTMLFormatter(TextFormatter):
    def format_variant(self, variant: Tuple[str]) -> str: ...

    def format_feature(self, feature: 'gnomic.types.Feature') -> str: ...

    def format_fusion(self, fusion: 'gnomic.types.Fusion') -> str: ...

    def format_plasmid(self, plasmid: 'gnomic.types.Plasmid') -> str: ...


BUILTIN_FORMATTERS: Dict[str, Formatter]