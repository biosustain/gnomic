from typing import Tuple, Optional, Iterator, Sequence, Union, Iterable


class Change(object):
    before: Optional[Union['Annotation', 'AtLocus']]
    after: Optional['Annotation']
    multiple: bool

    def __init__(self,
                 before: Union['Annotation', 'AtLocus'] = None,
                 after: 'Annotation' = None,
                 multiple: bool = False) -> None: ...

    @classmethod
    def parse(self, gnomic_change_string: str) -> 'Change': ...

    def is_presence(self) -> bool:
        ...

    def __mod__(self, locus: 'Annotation') -> 'Change': ...

    def __matmul__(self, locus: 'Annotation') -> 'Change': ...


class Present(Change):
    def __init__(self, annotation: 'Annotation') -> None:
        ...


class Annotation(object):
    def match(self, other) -> bool: ...

    def __pos__(self) -> 'Mutation': ...

    def __neg__(self) -> 'Mutation': ...

    def __gt__(self, other: 'Annotation') -> 'Mutation': ...

    def __rshift__(self, other: 'Annotation') -> 'Mutation': ...

    def __mod__(self, locus: 'Annotation') -> 'AtLocus': ...

    def __matmul__(self, locus: 'Annotation') -> 'AtLocus': ...

    def __pow__(self, other: 'Annotation') -> 'Annotation': ...

    def __str__(self) -> str: ...


class AtLocus(object):
    annotation: Annotation
    locus: Annotation

    def __init__(self, annotation: 'Annotation', locus: 'Annotation') -> None: ...

    def __neg__(self) -> 'Mutation': ...

    def __gt__(self, other: 'Annotation') -> 'Change': ...

    def __rshift__(self, other: 'Annotation') -> 'Change': ...


class Feature(Annotation):
    name: Optional[str]
    type: Optional[str]
    accession: Optional['Accession']
    organism: Optional[str]
    variant: Optional[Tuple[str]]

    def __init__(self,
                 name: str = None,
                 type: str = None,
                 accession: Accession = None,
                 organism: str = None,
                 variant: Tuple[str] = None
                 ) -> None: ...

    @classmethod
    def parse(cls, gnomic_feature_string: str) -> 'Feature': ...

    def match(self, other, match_variant: bool = True) -> bool: ...


class CompositeAnnotationBase(Annotation):
    annotations: Tuple[Annotation]

    def __init__(self, *annotations: Annotation) -> None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, item: int) -> Annotation: ...

    def __iter__(self) -> Iterator[Annotation]: ...

    def contains(self, other: Annotation) -> bool: ...


class CompositeAnnotation(CompositeAnnotationBase):
    ...


class AnnotationSetCompositeAnnotation(CompositeAnnotationBase):
    ...


class Fusion(CompositeAnnotationBase):
    def index(self, other: Annotation) -> int: ...

    @classmethod
    def fuse(cls, annotations: Sequence[Annotation]) -> Optional[Annotation]: ...


class Plasmid(CompositeAnnotationBase):
    name: str

    def __init__(self, name: str, annotations: Iterable[Annotation] = ()) -> None: ...


class Accession(object):
    identifier: str
    database: Optional[str]

    def __init__(self, identifier: str, database: str = None) -> None: ...
