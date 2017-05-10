# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import collections

from gnomic.models import Mutation, Plasmid, Fusion, FeatureTree, FeatureSet


VARIANT_MAP = {
    'wild-type': {
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

FORMAT = {
    'substitution': {
        'text': '{delta}{old}::{new}',
        'html': '{delta}{old}::{new}',
        'string': '{old}>{multiple}{new}',
    },
    'indel': {
        'text': '{delta}{feature}',
        'html': '{delta}{feature}',
        'string': '{operator}{feature}',
    },
    'plasmid_contents_integrated': {
        'text': '{name}({content})',
        'html': '<span class="gnomic-plasmid"><span class="gnomic-plasmid-name">{name}</span>({content})</span>',
        'string': r'({name} {contents}',
    },
    'plasmid_contents': {
        'text': '({name} {content})',
        'html': '<span class="gnomic-plasmid">(<span class="gnomic-plasmid-name">{name}</span> {content})</span>',
        'string': r'({name} {content})',
    },
    'plasmid_integrated': {
        'text': '{name}',
        'html': '<span class="gnomic-plasmid"><span class="gnomic-plasmid-name">{name}</span></span>',
        'string': '({name})',
    },
    'plasmid': {
        'text': '({name})',
        'html': '<span class="gnomic-plasmid">(<span class="gnomic-plasmid-name">{name}</span>)</span>',
        'string': '({name})',
    },
    'variant': {
        'text': '({variant})',
        'html': '<sup>{variant}</sup>',
        'string': '({variant})',
    },
    'variant_in_map': {
        'text': '{variant}',
        'html': '<sup>{variant}</sup>',
        'string': '{variant}',
    },
    'fusion': {
        'text': '{fusion}',
        'html': '<span class="gnomic-fusion">{fusion}</span>',
        'string': '{fusion}',
    },
    'feature-set': {
        'text': '{{{contents}}}',
        'html': '<span class="gnomic-feature-set">{{{contents}}}</span>',
        'string': '{{{contents}}}',
    },
    'feature': {
        'text': '{is_marker}{feature}',
        'html': '{is_marker}<span class="gnomic-feature">{feature}</span>',
        'string': '{is_marker}{feature}'
    },
}


def html_escape(s, quote=False):
    s = s.replace('&', '&amp;')
    s = s.replace('<', '&lt;')
    s = s.replace('>', '&gt;')
    if quote:
        s = s.replace('"', '&quot;')
    return s


def namedtuple_with_defaults(typename, field_names, default_values=()):
    T = collections.namedtuple(typename, field_names)
    T.__new__.__defaults__ = (None,) * len(T._fields)
    if isinstance(default_values, collections.Mapping):
        prototype = T(**default_values)
    else:
        prototype = T(*default_values)
    T.__new__.__defaults__ = tuple(prototype)
    return T


def genotype_to_text(genotype, fusions=False):
    '''
    A method to create a more biologist friendly representation
    of the genotype string
    :param genotype: An instance of the Genotype class
    :param bool fusions: Whether to keep fusions together.
    :param delta_char: This symbol is used to display a deletion
    :return: str
    '''
    return ' '.join(change_to_text(c) for c in genotype.changes(fusions=fusions))


def genotype_to_string(genotype, fusions=True):
    return ' '.join(change_to_string(c) for c in genotype.changes(fusions=fusions))


def genotype_to_html(genotype, fusions=False):
    return ' '.join(change_to_html(c) for c in genotype.changes(fusions=fusions))


def change_to_text(change, delta_char=u'\u0394'):
    return change_to_format(change, output='text', delta_char=delta_char)


def change_to_string(change):
    return change_to_format(change, output='string')


def change_to_html(change, delta_char=u'\u0394'):
    return change_to_format(change, output='html', delta_char=delta_char)


def change_to_format(change, output, delta_char=u'\u0394'):
    if output not in ('text', 'html', 'string'):
        raise ValueError('output must be one of: "html", "text", "string"')

    feature_to_format = {'text': feature_to_text, 'string': feature_to_string, 'html': feature_to_html}[output]

    if isinstance(change, Mutation):
        if change.old and change.new:
            s = FORMAT['substitution'][output].format(old=feature_to_format(change.old),
                                                      new=feature_to_format(change.new),
                                                      delta=delta_char,
                                                      multiple='>' if change.multiple else '')
        elif change.old is None: # insertion
            s = FORMAT['indel'][output].format(delta='', feature=feature_to_format(change.new), operator='+')
        elif change.new is None: # deletion
            s = FORMAT['indel'][output].format(delta=delta_char,
                                               feature=feature_to_format(change.old, integrated=False),
                                               operator='-')
    elif isinstance(change, Plasmid):
        s = feature_to_format(change, integrated=False)
    else:
        raise TypeError()

    if change.markers:
        s += '::{}'.format(feature_to_format(change.markers, is_marker=True))
    return s


def feature_to_text(feature, integrated=True, is_marker=False):
    '''
    A method to transform a genotype feature into text
    :param feature: Genotype feature
    :param integrated: boolean
    :param is_marker: boolean
    :return: str
    '''
    return feature_to_output(feature, output='text', integrated=integrated, is_marker=is_marker)


def feature_to_string(feature, **kwargs):
    return feature_to_output(feature, output='string', integrated=False, is_marker=False)


def feature_to_html(feature, integrated=True, is_marker=False):
    return feature_to_output(feature, output='html', integrated=integrated, is_marker=is_marker)


def feature_to_output(feature, output, integrated=True, is_marker=False):
    if output not in ('text', 'html', 'string'):
        raise ValueError('output must be one of: "html", "text", "string"')

    escape_html = (output == 'html')

    if isinstance(feature, Plasmid):
        name = html_escape(feature.name) if escape_html else feature.name
        if feature.contents:
            content = ' '.join(feature_to_output(f, output) for f in feature.contents)
            if integrated:
                return FORMAT['plasmid_contents_integrated'][output].format(name=name, content=content)
            else:
                return FORMAT['plasmid_contents'][output].format(name=name, content=content)
        elif integrated:
            return FORMAT['plasmid_integrated'][output].format(name=name)
        else:
            return FORMAT['plasmid'][output].format(name=name)

    elif isinstance(feature, Fusion):
        return FORMAT['fusion'][output].format(fusion=':'.join(feature_to_output(f, output) for f in feature.contents))

    elif isinstance(feature, FeatureTree):
        contents = ' '.join(feature_to_output(f, output) for f in feature.contents)
        if len(feature) != 1:
            if isinstance(feature, FeatureSet):
                return  FORMAT['feature-set'][output].format(contents=contents)
            else:
                return '{{{contents}}}'.format(contents=contents)
        elif isinstance(feature, FeatureSet) and output == 'string':
            return FORMAT['feature-set'][output].format(contents=contents)
        else:
            return contents
    else:
        text = ''
        if feature.organism:
            text += '{}/'.format(feature.organism.name)
        if feature.type and feature.type.name != 'phene':
            text += '{}.'.format(feature.type.name)

        text += feature.name
        if escape_html:
            text = html_escape(text)

        if feature.variant in VARIANT_MAP:
            text += FORMAT['variant_in_map'][output].format(variant=VARIANT_MAP[feature.variant][output])
        elif feature.variant:
            text += FORMAT['variant'][output].format(variant=html_escape(feature.variant) if escape_html
                                                     else feature.variant)

        if feature.accession:
            accession = '#' + (feature.accession.database + ':' if feature.accession.database else '') \
                    + feature.accession.identifier
            text += html_escape(accession) if escape_html else accession

        return FORMAT['feature'][output].format(feature=text, is_marker='::' if is_marker else '')