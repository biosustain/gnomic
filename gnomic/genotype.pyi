from typing import Iterable, Tuple, Optional, Sequence, List, Union, Set

from gnomic.types import Change, Annotation, AtLocus, Plasmid, Feature, CompositeAnnotation, Fusion


def partial_match(search: Annotation, target: Annotation) -> bool: ...


def change_annotation(annotation: Annotation,
                      site: Annotation = None,
                      replacement: Annotation = None) -> Optional[Annotation]:
    ...


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

    def __init__(self, changes: Iterable[Change], parent: 'Genotype' = None) -> None: ...

    @classmethod
    def _parse_gnomic_string(cls, gnomic_string: str, *args, **kwargs) -> List[Change]: ...

    @classmethod
    def parse(cls, gnomic_string: str, parent: 'Genotype' = None) -> 'Genotype': ...

    @classmethod
    def is_valid(cls, gnomic_string: str) -> bool: ...

    def changes(self) -> Tuple[Change]: ...

    @property
    def added_features(self) -> Set[Feature]: ...

    @property
    def removed_features(self) -> Set[Feature]: ...

    @property
    def added_plasmids(self) -> Set[Plasmid]: ...

    @property
    def removed_plasmids(self) -> Set[Plasmid]: ...

    @property
    def added_fusions(self) -> Set[Fusion]: ...

    @property
    def removed_fusions(self) -> Set[Fusion]: ...

    @property
    def added_fusion_features(self) -> Set[Union[Feature, Fusion, CompositeAnnotation]]: ...

    @property
    def removed_fusion_features(self) -> Set[Union[Feature, Fusion, CompositeAnnotation]]: ...

    def format(self, output: 'str' = 'text') -> str: ...
