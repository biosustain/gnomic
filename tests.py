import json

from gnomic.grammar import GnomicParser

parser = GnomicParser()
ast = parser.parse('+{y:z} '       # insertion
                   'p1{x:y }::m+ '   # plasmid
                   'p2{a}::m(R) '  # plasmid alt
                   'p3{}'          # plasmid no content
                   'site>p4{}'     # plasmid integrated
                   'a > b '        # replacement
                   'a>ecoli/b(variant)::m+ '
                   '-del '          # deletion
                   '-del::bX+ '
                   '-p5{} '
                   'X(v)'          # phenotype
                   '', rule_name='start')
print(ast)
print(json.dumps(ast, indent=2)) # ASTs are JSON-friendy