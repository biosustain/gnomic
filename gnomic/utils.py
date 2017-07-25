from gnomic import Genotype


def chain(*gnomic_strings, parent=None, **kwargs):
    genotype = Genotype.parse(gnomic_strings[0], parent=parent, **kwargs)
    for gnomic_string in gnomic_strings[1:]:
        genotype = Genotype.parse(gnomic_string, parent=genotype, **kwargs)
    return genotype
