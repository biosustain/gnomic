from typing import Any

from gnomic import Genotype
from gnomic.types import Change, Feature


def chain(*gnomic_strings: str, parent: Genotype = None, **kwargs: Any) -> Genotype: ...


def genotype_to_string(genotype: Genotype) -> str: ...


def genotype_to_text(genotype: Genotype) -> str: ...


def genotype_to_html(genotype: Genotype) -> str: ...


def change_to_string(change: Change) -> str: ...


def change_to_text(change: Change) -> str: ...


def change_to_html(change: Change) -> str: ...


def feature_to_string(feature: Feature) -> str: ...


def feature_to_text(feature: Feature) -> str: ...


def feature_to_html(feature: Feature) -> str: ...
