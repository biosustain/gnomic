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


def genotype_to_text(genotype, fusions=False, delta_char=u"\u0394"):
    """
    A method to create a more biologist friendly representation
    of the genotype string
    :param genotype: An instance of the Genotype class
    :param bool fusions: Whether to keep fusions together.
    :param delta_char: This symbol is used to display a deletion
    :return: str
    """
    result = []
    for change in genotype.changes(fusions=fusions):
        change_type = type(change)
        if change_type is Mutation:
            if change.old and change.new:
                # Substitution
                result.append(delta_char + feature_to_text(change.old) + '::' + feature_to_text(change.new))
            elif not change.old:
                # Insertion
                result.append(feature_to_text(change.new))
            elif not change.new:
                # Deletion
                result.append(delta_char + feature_to_text(change.old))

        elif change_type is Plasmid:
            result.append(feature_to_text(change, integrated=False))

        if change.markers:
            result.append(feature_to_text(change, is_maker=True))

        result_string = ""

    for i, item in enumerate(result):
        if i + 1 < len(result) - 1:
            if "::" in result[i]:
                result_string += item
                continue
        result_string += " %s" % item

    return " ".join(result)


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
                # Added the caret char to show it should be superscript
                text += "^%s" % variant

        return text


def genotype_to_string(genotype, fusions=True):
    parts = []

    for change in genotype.changes(fusions=fusions):
        s = None

        if isinstance(change, Mutation):
            if change.old and change.new:
                if change.multiple:
                    s = '{}>>{}'.format(feature_to_string(change.old), feature_to_string(change.old))
                else:
                    s = '{}>{}'.format(feature_to_string(change.old), feature_to_string(change.old))
            elif change.old is None:
                # FIXME phenotypes should not have a + or - prefix; only a postfix.
                s = '+{}'.format(feature_to_string(change.new))
            elif change.new is None:
                s = '-{}'.format(feature_to_string(change.old))
        elif isinstance(change, Plasmid):
            s = feature_to_string(change)

        if change.markers:
            s += '::{}'.format(feature_to_string(change.markers))

        parts.append(s)
    return ' '.join(parts)


def feature_to_string(feature):
    if isinstance(feature, Plasmid):
        name = feature.name
        if feature.contents:
            contents = ' '.join(feature_to_string(f) for f in feature.contents)
            return r'{name}{{{contents}}}'.format(name=name, contents=contents)
        else:
            return '{name}{{}}'.format(name=name)
    elif isinstance(feature, Fusion):
        return ':'.join(feature_to_string(f) for f in feature.contents)
    elif isinstance(feature, FeatureTree):
        contents = ' '.join(feature_to_string(f) for f in feature.contents)
        if isinstance(feature, FeatureSet):
            return '{{{contents}}}'.format(contents=contents)
        elif len(feature) != 1:
            raise NotImplementedError
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
