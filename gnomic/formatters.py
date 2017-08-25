# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from abc import abstractmethod, ABCMeta

import six

from gnomic.types import Feature, Fusion, Plasmid, AtLocus, CompositeAnnotation

DELTA = '\u0394'

RIGHTWARDS_ARROW = '\u2192'

RIGHTWARDS_PAIRED_ARROW = '\u21c9'


def escape_html(s, quote=False):
    s = s.replace('&', '&amp;')
    s = s.replace('<', '&lt;')
    s = s.replace('>', '&gt;')
    if quote:
        s = s.replace('"', '&quot;')
    return s


class Formatter(six.with_metaclass(ABCMeta)):
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
        elif isinstance(annotation, CompositeAnnotation):
            return self.format_composite_annotation(annotation)
        raise NotImplementedError

    @staticmethod
    def format_accession(accession):
        if accession.database:
            return '#{}:{}'.format(accession.database, accession.identifier)
        else:
            return '#{}'.format(accession.identifier)

    def format_variant(self, variant):
        return '({})'.format('; '.join(v for v in variant))

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
            s += self.format_variant(feature.variant)
        return s

    def format_fusion(self, fusion):
        return ':'.join(map(self.format_annotation, fusion.annotations))

    def format_plasmid(self, plasmid):
        if plasmid.annotations:
            return '({} {})'.format(plasmid.name, ' '.join(map(self.format_annotation, plasmid.annotations)))
        return '({})'.format(plasmid.name)

    def format_at_locus(self, at_locus):
        return '{}@{}'.format(self.format_annotation(at_locus.annotation), self.format_annotation(at_locus.locus))

    def format_composite_annotation(self, composite_annotation):
        return '{{{}}}'.format(', '.join(map(self.format_annotation, composite_annotation.annotations)))


class GnomicFormatter(Formatter):
    format = 'string'

    def format_change(self, change):
        after = self.format_annotation(change.after) if change.after is not None else None
        before = self.format_annotation(change.before) if change.before is not None else None

        if change.is_presence():
            return after
        elif after is None:
            return '-{}'.format(before)
        elif before is None:
            return '+{}'.format(after)
        elif change.multiple:
            return '{}>>{}'.format(before, after)
        else:
            return '{}>{}'.format(before, after)


class TextFormatter(Formatter):
    format = 'text'

    def format_change(self, change):
        after = self.format_annotation(change.after) if change.after is not None else None
        before = self.format_annotation(change.before) if change.before is not None else None
        if after is None:
            return '{}{}'.format(DELTA, before)
        elif before is None:
            return '{}'.format(after)
        elif change.multiple:
            return '{}{}{}{}'.format(DELTA, before, RIGHTWARDS_PAIRED_ARROW, after)
        else:
            return '{}{}{}{}'.format(DELTA, before, RIGHTWARDS_ARROW, after)


class HTMLFormatter(TextFormatter):
    format = 'html'

    def format_variant(self, variant):
        return '<sup>{}</sup>'.format('; '.join(escape_html(v) for v in variant))

    def format_feature(self, feature):
        s = ''
        if feature.organism:
            s += '{}/'.format(escape_html(feature.organism))
        if feature.type:
            s += '{}.'.format(escape_html(feature.type))
        if feature.name:
            s += escape_html(feature.name)
        if feature.accession:
            s += escape_html(self.format_accession(feature.accession))
        if feature.variant:
            s += self.format_variant(feature.variant)
        return '<span class="gnomic-feature">{}</span>'.format(s)

    def format_fusion(self, fusion):
        return '<span class="gnomic-fusion">{}</span>'.format(':'.join(map(self.format_annotation, fusion.annotations)))

    def format_plasmid(self, plasmid):
        if plasmid.annotations:
            s = '(<span class="gnomic-plasmid-name">{}</span> {})' \
                .format(escape_html(plasmid.name), ' '.join(map(self.format_annotation, plasmid.annotations)))
        else:
            s = '(<span class="gnomic-plasmid-name">{}</span>)'.format(escape_html(plasmid.name))
        return '<span class="gnomic-plasmid">{}</span>'.format(s)


BUILTIN_FORMATTERS = {
    'gnomic': GnomicFormatter(),
    'text': TextFormatter(),
    'html': HTMLFormatter()
}
