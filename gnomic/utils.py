import collections
from gnomic.models import Mutation, Plasmid, Fusion, FeatureTree

def namedtuple_with_defaults(typename, field_names, default_values=()):
    T = collections.namedtuple(typename, field_names)
    T.__new__.__defaults__ = (None,) * len(T._fields)
    if isinstance(default_values, collections.Mapping):
        prototype = T(**default_values)
    else:
        prototype = T(*default_values)
    T.__new__.__defaults__ = tuple(prototype)
    return T


def genotype_to_text(genotype, delta_char=u"\u0394"):
    """
    A method to create a more biologist friendly representation
    of the genotype string
    :param genotype: An instance of the Genotype class
    :param delta_char: This symbol is used to display a deletion
    :return: str
    """
    result = []
    for change in genotype.changes():
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
