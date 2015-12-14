import json
from pprint import pprint

from gnomic import Fusion, Organism, Accession
from gnomic.grammar import GnomicParser
from semantics import DefaultSemantics


class GnomicJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, Fusion):
            return {
                "fusion": o.contents
            }
        elif isinstance(o, Organism):
            return o.default_alias
        elif isinstance(o, Accession):
            return {"db": o.database, "id": o.identifier}
        # TODO full JSON encoding
        try:
            iterable = iter(o)
        except TypeError:
            pass
        else:
            return list(iterable)
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, o)


parser = GnomicParser()
ast = parser.parse('a(r) -b +c+'
                   # '+{y:} '       # insertion
                   # 'p1{x:y }::m+ '   # plasmid
                   # 'p2{a}::m(R) '  # plasmid alt
                   # 'p3{}'          # plasmid no content
                   # 'site>p4{}'     # plasmid integrated
                   # 'a > b '        # replacement
                   # 'a>ecoli/b(variant)::m+ '
                   # '-del '          # deletion
                   # '-del::bX+ '
                   # '-p5{} '
                   # 'X(v)'          # phenotype
                   , whitespace='', semantics=DefaultSemantics(organisms=[Organism('E.coli', ('e.coli', 'ec', 'ecoli'))]), rule_name='start')
print()
print()
pprint(ast)
#print(json.dumps(ast, indent=2, cls=GnomicJSONEncoder))  # ASTs are JSON-friendy


# TODO type bindings, organism bindings
