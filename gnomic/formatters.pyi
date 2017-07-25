from abc import ABC, abstractmethod
from typing import Dict


class Formatter(ABC):
    @abstractmethod
    def format_genotype(self, genotype: 'gnomic.Genotype') -> str: ...

    @abstractmethod
    def format_change(self, change: 'gnomic.types.Change') -> str: ...

    @abstractmethod
    def format_annotation(self, annotation: 'gnomic.types.Annotation') -> str: ...

    @abstractmethod
    def format_feature(self, feature: 'gnomic.types.Feature') -> str: ...

    # TODO fusion, plasmid, etc.


class GnomicFormatter(Formatter):
    def format_genotype(self, genotype: 'gnomic.Genotype') -> str: ...

    def format_change(self, change: 'gnomic.types.Change') -> str: ...

    def format_annotation(self, annotation: 'gnomic.types.Annotation') -> str: ...


class TextFormatter(Formatter):
    def format_genotype(self, genotype: 'gnomic.Genotype') -> str: ...

    def format_change(self, change: 'gnomic.types.Change') -> str: ...

    def format_annotation(self, annotation: 'gnomic.types.Annotation') -> str: ...


class HTMLFormatter(Formatter):
    def format_genotype(self, genotype: 'gnomic.Genotype') -> str: ...

    def format_change(self, change: 'gnomic.types.Change') -> str: ...

    def format_annotation(self, annotation: 'gnomic.types.Annotation') -> str: ...


BUILTIN_FORMATTERS: Dict[str, Formatter]