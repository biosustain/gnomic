from pprint import pprint
from gnomic import Genotype, Organism

ast = Genotype.parse('site>plasmidA{abc} -{a:b} -{first second} -plasmidB{} +{P.y:gene.dd} +X Y(variant)'       # insertion
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

print(Genotype.is_valid('sdassfs df=fwer34'))
print()
pprint(ast)

print()
pprint(ast.added_features)
pprint(ast.removed_features)
print()
pprint(ast.added_plasmids)
pprint(ast.removed_plasmids)
print()
pprint(ast.added_fusion_features)
pprint(ast.removed_fusion_features)
print()
#print(json.dumps(ast, indent=2, cls=GnomicJSONEncoder))  # ASTs are JSON-friendy


# TODO type bindings, organism bindings
