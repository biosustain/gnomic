from typing import Any

from gnomic import Genotype


def chain(*gnomic_strings: str, parent: Genotype = None, **kwargs: Any) -> Genotype: ...
