from abc import abstractmethod, ABC

from gnomic.types import Feature, Fusion, Plasmid


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

    @abstractmethod
    def format_fusion(self, fusion):
        pass

    @abstractmethod
    def format_plasmid(self, plasmid):
        pass


class GnomicFormatter(Formatter):
    def format_genotype(self, genotype):
        return ' '.join(self.format_change(change) for change in genotype.changes())

    def format_change(self, change):
        after = self.format_annotation(change.after) if change.after is not None else None
        before = self.format_annotation(change.before) if change.before is not None else None
        if after is None:
            return '-{}'.format(before)
        elif before is None:
            return '+{}'.format(after)
        elif change.multiple:
            return '{}>>{}'.format(before, after)
        else:
            return '{}>{}'.format(before, after)

    def format_annotation(self, annotation):
        if isinstance(annotation, Feature):
            return self.format_feature(annotation)
        elif isinstance(annotation, Fusion):
            return self.format_fusion(annotation)
        elif isinstance(annotation, Plasmid):
            return self.format_plasmid(annotation)
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
            if feature.accession.database:
                s += '#{}:{}'.format(feature.accession.database, feature.accession.identifier)
            else:
                s += '#{}'.format(feature.accession.identifier)
        if feature.variant:
            s += '({})'.format('; '.join(feature.variant))
        return s

    def format_fusion(self, fusion):
        return ':'.join(map(self.format_annotation, fusion.annotations))

    def format_plasmid(self, plasmid):
        if plasmid.annotations:
            return '({} {})'.format(plasmid.name, ' '.join(map(self.format_annotation, plasmid.annotations)))
        return '({})'.format(plasmid.name)


BUILTIN_FORMATTERS = {
    'gnomic': GnomicFormatter(),
    # 'text': TextFormatter(),
    # 'html': HTMLFormatter()
}
