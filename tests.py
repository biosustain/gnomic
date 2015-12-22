from pprint import pprint
from gnomic import Genotype, Organism

ast = Genotype.parse('+{P.y:gene.dd} '       # insertion
                           # 'p1{x:y }::m+ '   # plasmid
                           # 'p2{a}::m(R) '  # plasmid alt
                           # 'p3{} '          # plasmid no content
                           # 'site>p4{} '     # plasmid integrated
                           # 'a>b '        # replacement
                           # 'a>E.coli/b(variant)::m+ '
                           # '-del '          # deletion
                           # '-del::bX+ '
                           # '-p5{} '
                           # 'X(v)'
                     )         # phenotype)


print()
print()
pprint(ast)
#print(json.dumps(ast, indent=2, cls=GnomicJSONEncoder))  # ASTs are JSON-friendy


# TODO type bindings, organism bindings
