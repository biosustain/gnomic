from abc import abstractmethod, ABC

from gnomic.types import Feature


class Formatter(ABC):
    @abstractmethod
    def format_genotype(self, genotype):
        pass

    @abstractmethod
    def format_change(self, change):
        pass

    @abstractmethod
    def format_annotation(self, annotation):
        pass

    @abstractmethod
    def format_feature(self, feature):
        pass


class GnomicFormatter(Formatter):
    def format_genotype(self, genotype):
        pass

    def format_change(self, change):
        pass

    def format_annotation(self, annotation):
        if isinstance(annotation, Feature):
            return self.format_feature(annotation)
        raise NotImplementedError

    def format_feature(self, feature):
        s = ''
        if feature.organism:
            s += '{}/'.format(feature.organism)
        if feature.type:
            s += '{}.'.format(feature.type)
        if feature.name:
            s += feature.name
        if feature.accession:
            s += '{}#'.format(feature.accession)
        if feature.variant:
            s += '({})'.format('; '.join(feature.variant))
        return s


BUILTIN_FORMATTERS = {
    'gnomic': GnomicFormatter(),
    # 'text': TextFormatter(),
    # 'html': HTMLFormatter()
}
