from gnomic.grammar import GnomicSemantics
from gnomic.types import Fusion, Feature, Plasmid, Change, Accession, AtLocus, CompositeAnnotation


class DefaultSemantics(GnomicSemantics):
    def FUSION(self, ast):
        return Fusion(*ast)

    def FEATURE_FUSION(self, ast):
        return Fusion(*ast)

    def COMPOSITE_ANNOTATION(self, ast):
        return CompositeAnnotation(*ast)

    def DNA_SEQUENCE_VARIANT(self, ast):
        return ''.join(map(str, ast))

    def PROTEIN_SEQUENCE_VARIANT(self, ast):
        return ''.join(map(str, ast))

    def VARIABLE_VARIANT(self, ast):
        return ''.join(ast)

    def SEQUENCE_VARIANT(self, ast):
        return ''.join(ast)

    def INSERTION(self, ast):
        return Change(after=ast.after)

    def REPLACEMENT(self, ast):
        return Change(before=ast.before, after=ast.after, multiple=ast.op == '>>')

    def DELETION(self, ast):
        return Change(before=ast.before)

    def INTEGER(self, ast):
        return int(ast)

    def CHANGE(self, ast):
        if isinstance(ast, (Plasmid, Feature)):
            return Change(after=ast)
        return ast

    def PLASMID(self, ast):
        return Plasmid(ast.name, ast.annotations or ())

    def PHENE(self, ast):
        return self.FEATURE(ast)

    def FEATURE(self, ast):
        return Feature(ast.name,
                       ast.type,
                       accession=ast.accession,
                       organism=ast.organism,
                       variant=tuple(ast.variant) if ast.variant else None)

    def NUCLEOTIDE_SEQUENCE(self, ast):
        return ''.join(ast)

    def AMINO_ACID_SEQUENCE(self, ast):
        return ''.join(ast)

    def ACCESSION(self, ast):
        return Accession(ast.id, ast.db)

    def ANNOTATION_AT_LOCUS(self, ast):
        return AtLocus(ast.annotation, ast.locus)
