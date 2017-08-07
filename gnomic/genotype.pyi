from typing import Iterable, Tuple, Optional, Sequence, List, Union, Set

from gnomic.types import Change, Annotation, AtLocus, Plasmid, Feature, CompositeAnnotation, Fusion


def partial_match(search: Annotation, target: Annotation) -> bool: ...


def change_annotation(annotation: Annotation,
                      site: Annotation = None,
                      replacement: Annotation = None) -> Optional[Annotation]:
    pass


class GenotypeState(object):
    _changes: List[Tuple[Union[Annotation, AtLocus, None], Union[Annotation, None]]]

    def __init__(self, changes: Sequence[Change] = ()) -> None:
        ...

    def insert(self, annotation: Annotation, multiple: bool = False) -> None:
        ...

    def remove(self, site: Union[Annotation, AtLocus], multiple: bool = False) -> None:
        ...

    def replace(self, site: Union[Annotation, AtLocus], replacement: Annotation, multiple: bool = False) -> None:
        ...

    def change(self, change: Change) -> None:
        ...

    @property
    def changes(self) -> Tuple[Change]: ...


class Genotype(object):
    parent: Optional['Genotype']
    state: GenotypeState

    # TODO implement these (as @property or method)
    added_plasmids: Set[Plasmid]
    removed_plasmids: Set[Plasmid]
    added_features: Set[Feature]
    removed_features: Set[Feature]
    added_fusion_features: Set[Union[Feature, Fusion, CompositeAnnotation]]
    removed_fusion_features: Set[Feature, Fusion, CompositeAnnotation]

    def __init__(self,
                 changes: Iterable[Change],
                 parent: 'Genotype' = None) -> None:
        ...

    @classmethod
    def _parse_gnomic_string(cls, gnomic_string: str, *args, **kwargs) -> List[Change]: ...

    @classmethod
    def parse(cls, gnomic_string: str, parent: 'Genotype' = None) -> 'Genotype': ...

    @classmethod
    def is_valid(cls, gnomic_string: str) -> bool: ...

    def changes(self) -> Tuple[Change]: ...

    @property
    def added_features(self) -> Tuple[Feature]: ...

    @property
    def removed_features(self) -> Tuple[Feature]: ...

    @property
    def added_plasmids(self) -> Tuple[Plasmid]: ...

    @property
    def removed_plasmids(self) -> Tuple[Plasmid]: ...

    @property
    def added_fusions(self) -> Tuple[Fusion]: ...

    @property
    def removed_fusions(self) -> Tuple[Fusion]: ...


    def format(self, output: 'str' = 'text') -> str: ...
