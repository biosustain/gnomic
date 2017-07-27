from abc import abstractmethod, ABC

from gnomic.types import Feature, Fusion, Plasmid, AtLocus

VARIANT_MAP = {
    'wild-type':{
        'text': u'\u207A',
        'html': '+',
        'string': '+',
    },
    'mutant': {
        'text': u'\u207B',
        'html': '-',
        'string': '-',
    }

}

DELTA = u'\u0394'


class Formatter(ABC):
    def format_genotype(self, genotype):
        return ' '.join(self.format_change(change) for change in genotype.changes())

    @abstractmethod
    def format_change(self, change):
        pass

    def format_annotation(self, annotation):
        if isinstance(annotation, Feature):
            return self.format_feature(annotation)
        elif isinstance(annotation, Fusion):
            return self.format_fusion(annotation)
        elif isinstance(annotation, Plasmid):
            return self.format_plasmid(annotation)
        elif isinstance(annotation, AtLocus):
            return self.format_at_locus(annotation)
        raise NotImplementedError

    @staticmethod
    def format_accession(accession):
        if accession.database:
            return '#{}:{}'.format(accession.database, accession.identifier)
        else:
            return '#{}'.format(accession.identifier)

    def format_feature(self, feature):
        s = ''
        if feature.organism:
            s += '{}/'.format(feature.organism)
        if feature.type:
            s += '{}.'.format(feature.type)
        if feature.name:
            s += feature.name
        if feature.accession:
            s += self.format_accession(feature.accession)
        if feature.variant:
            if any(variant in VARIANT_MAP for variant in feature.variant):
                s += '{}'.format(''.join(VARIANT_MAP[v]['string'] for v in feature.variant if v in VARIANT_MAP))
            if any(variant not in VARIANT_MAP for variant in feature.variant):
                s += '({})'.format('; '.join(v for v in feature.variant if v not in VARIANT_MAP))
        return s

    def format_fusion(self, fusion):
        return ':'.join(map(self.format_annotation, fusion.annotations))

    def format_plasmid(self, plasmid):
        if plasmid.annotations:
            return '({} {})'.format(plasmid.name, ' '.join(map(self.format_annotation, plasmid.annotations)))
        return '({})'.format(plasmid.name)

    def format_at_locus(self, at_locus):
        return '{}@{}'.format(self.format_annotation(at_locus.annotation), self.format_annotation(at_locus.locus))


class GnomicFormatter(Formatter):
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


class TextFormatter(Formatter):
    def format_change(self, change):
        after = self.format_annotation(change.after) if change.after is not None else None
        before = self.format_annotation(change.before) if change.before is not None else None
        if before is not None:
            return '{}{}'.format(DELTA, before)
        else:
            return '{}'.format(after)


class HTMLFormatter(Formatter):
    def format_change(self, change):
        pass


BUILTIN_FORMATTERS = {
    'gnomic': GnomicFormatter(),
    'text': TextFormatter(),
    'html': HTMLFormatter()
}
