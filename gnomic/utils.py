from gnomic import Genotype
from gnomic.formatters import BUILTIN_FORMATTERS


def chain(*gnomic_strings, **kwargs):
    parent = kwargs.pop('parent', None)
    genotype = Genotype.parse(gnomic_strings[0], parent=parent, **kwargs)
    for gnomic_string in gnomic_strings[1:]:
        genotype = Genotype.parse(gnomic_string, parent=genotype, **kwargs)
    return genotype


def genotype_to_string(genotype):
    return BUILTIN_FORMATTERS['gnomic'].format_genotype(genotype)


def genotype_to_text(genotype):
    return BUILTIN_FORMATTERS['text'].format_genotype(genotype)


def genotype_to_html(genotype):
    return BUILTIN_FORMATTERS['html'].format_genotype(genotype)


def change_to_string(change):
    return BUILTIN_FORMATTERS['gnomic'].format_change(change)


def change_to_text(change):
    return BUILTIN_FORMATTERS['text'].format_change(change)


def change_to_html(change):
    return BUILTIN_FORMATTERS['html'].format_change(change)


def feature_to_string(feature):
    return BUILTIN_FORMATTERS['gnomic'].format_feature(feature)


def feature_to_text(feature):
    return BUILTIN_FORMATTERS['text'].format_feature(feature)


def feature_to_html(feature):
    return BUILTIN_FORMATTERS['html'].format_feature(feature)
