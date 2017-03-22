import collections
from gnomic.models import Mutation, Plasmid, Fusion, FeatureTree, FeatureSet, Feature


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
    """
    A method to create a more biologist friendly representation
    of the genotype string
    :param genotype: An instance of the Genotype class
    :param bool fusions: Whether to keep fusions together.
    :param delta_char: This symbol is used to display a deletion
    :return: str
    """
    parts = []
    for change in genotype.changes(fusions=fusions):
        parts.append(change_to_text(change))

    return " ".join(parts)


def change_to_text(change, delta_char=u"\u0394"):
    s = None
    change_type = type(change)
    if change_type is Mutation:
        if change.old and change.new:  # substitution
            s = delta_char + feature_to_text(change.old) + '::' + feature_to_text(change.new)
        elif change.old is None:  # insertion
            s = feature_to_text(change.new)
        elif change.new is None:  # deletion
            s = delta_char + feature_to_text(change.old, integrated=False)
    elif change_type is Plasmid:
        s = feature_to_text(change, integrated=False)

    if change.markers:
        s += "::" + feature_to_text(change, is_maker=True)

    return s


def feature_to_text(feature, integrated=True, is_maker=False):
    """
    A method to transform a genotype feature into text
    :param feature: Genotype feature
    :param integrated: boolean
    :param is_maker: boolean
    :return: str
    """
    feature_type = type(feature)
    if feature_type is Plasmid:
        name = feature.name
        if feature.contents:
            content = ' '.join(map(feature_to_text, feature.contents))

            if integrated:
                return '%s(%s)' % (name, content)
            else:
                return '(%s %s)' % (name, content)

        elif integrated:
            return name
        else:
            return '(%s)' % name

    elif feature_type is Fusion:
        return ':'.join(map(feature_to_text, feature.contents))

    elif feature_type is FeatureTree:
        return ' '.join(map(feature_to_text, feature.contents))

    else:
        text = ''
        if is_maker:
            text += '::'

        if feature.organism:
            text += '%s/' % feature.organism.name

        text += feature.name

        variant_map = {'wild-type': u"\u207A",
                       'mutant': u"\u207B"}
        variant = feature.variant

        if variant:
            if variant in variant_map:
                text += variant_map[variant]
            else:
                text += "(%s)" % variant

        return text


def genotype_to_string(genotype, fusions=True):
    parts = []

    for change in genotype.changes(fusions=fusions):
        parts.append(change_to_string(change))
    return ' '.join(parts)


def change_to_string(change):
    s = None

    if isinstance(change, Mutation):
        if change.old and change.new:
            if change.multiple:
                s = '{}>>{}'.format(feature_to_string(change.old), feature_to_string(change.new))
            else:
                s = '{}>{}'.format(feature_to_string(change.old), feature_to_string(change.new))
        elif change.old is None:
            # FIXME phenotypes should not have a + or - prefix; only a postfix.
            s = '+{}'.format(feature_to_string(change.new))
        elif change.new is None:
            s = '-{}'.format(feature_to_string(change.old))
    elif isinstance(change, Plasmid):
        s = feature_to_string(change)
    else:
        raise TypeError()

    if change.markers:
        s += '::{}'.format(feature_to_string(change.markers))

    return s


def feature_to_string(feature):
    if isinstance(feature, Plasmid):
        name = feature.name
        if feature.contents:
            contents = ' '.join(feature_to_string(f) for f in feature.contents)
            return r'({name} {contents})'.format(name=name, contents=contents)
        else:
            return '({name})'.format(name=name)
    elif isinstance(feature, Fusion):
        return ':'.join(feature_to_string(f) for f in feature.contents)
    elif isinstance(feature, FeatureTree):
        contents = ' '.join(feature_to_string(f) for f in feature.contents)
        if isinstance(feature, FeatureSet):
            return '{{{contents}}}'.format(contents=contents)
        elif len(feature) != 1:
            return '{{{contents}}}'.format(contents=contents)
        else:
            return contents
    else:
        s = ''
        if feature.type and feature.type.name != 'phene':
            s += '{}.'.format(feature.type.name)

        if feature.organism:
            s += '{}/'.format(feature.organism.name)

        s += feature.name

        variant = feature.variant
        variant_map = {'wild-type': '+', 'mutant': '-'}

        if variant:
            try:
                s += variant_map[variant]
            except KeyError:
                s += '({})'.format(variant)
        return s
