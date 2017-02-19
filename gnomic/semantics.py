from gnomic.models import Mutation, Fusion, Plasmid, Feature, Organism, Accession, Type, FeatureTree, Range, \
    FeatureSet, Presence
from gnomic.grammar import GnomicSemantics


class DefaultSemantics(GnomicSemantics):
    def __init__(self,
                 organisms=None,
                 types=None):
        self._organisms = {} if organisms is None else Organism.map(organisms)
        self._types = {} if types is None else Type.map(types)

    def FUSION(self, ast):
        return Fusion(*ast)

    def FEATURE_FUSION(self, ast):
        return Fusion(*ast)

    def FEATURE_SET(self, ast):
        return FeatureSet(*ast)

    def ORGANISM(self, name):
        try:
            return self._organisms[name]
        except KeyError:
            self._organisms[name] = organism = Organism(name)
            return organism

    def BINARY_VARIANT(self, variant):
        if variant == '+':
            return 'wild-type'
        else:
            return 'mutant'

    def change(self, ast):
        if isinstance(ast, Mutation) or isinstance(ast, Presence):  # Mutation or deleted Plasmid
            return ast
        else:  # inserted Plasmid or Phene
            return Presence(ast, True)

    def insertion(self, ast):
        return Mutation(None, ast.after, locus=ast.locus, markers=ast.markers)

    def replacement(self, ast):
        return Mutation(ast.before,
                        ast.after.contents if isinstance(ast.after, Plasmid) else ast.after,
                        locus=ast.locus,
                        markers=ast.markers,
                        multiple=ast.op == '>>')

    def deletion(self, ast):
        if isinstance(ast.before, Plasmid):
            return Presence(ast.before, False, markers=ast.markers)
        else:
            return Mutation(ast.before, None, locus=ast.locus, markers=ast.markers)

    def RANGE(self, ast):
        level = {
            'c': 'coding',
            'r': 'RNA',
            'p': 'protein'
        }[ast.level]

        if ast.pos:
            return Range(level, ast.pos, ast.pos)
        return Range(level, ast.start, ast.end)

    def INTEGER(self, ast):
        return int(ast)

    def ACCESSION(self, ast):
        return Accession(ast['id'], ast['db'])

    def PLASMID(self, ast):
        return Plasmid(ast.name, ast.contents, markers=ast.markers)

    def PHENE(self, ast):
        return self.FEATURE(ast, default_type='phene')

    def FEATURE(self, ast, default_type=None):
        if ast.type or default_type:
            name = ast.type or default_type
            try:
                type = self._types[name]
            except KeyError:
                self._types[name] = type = Type(name)
        else:
            type = None

        return Feature(ast.name, type,
                       accession=ast.accession,
                       organism=ast.organism,
                       variant=', '.join(ast.variant) if isinstance(ast.variant, list) else ast.variant,
                       range=ast.range)
