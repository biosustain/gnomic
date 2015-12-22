import json

from gnomic.models import Fusion, Organism, Accession


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