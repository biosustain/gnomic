

class Mutation(object):
    def __init__(self, old, new, marker=None, multiple=False):
        self.old = old
        self.new = new
        self.marker = marker
        self.multiple = multiple

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.__dict__))


class FeatureTree(object):
    def __init__(self, *contents):
        self.contents = contents

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.contents))


class Fusion(FeatureTree):
    pass


class Plasmid(FeatureTree):
    def __init__(self, name, contents, site=None, marker=None):
        super(Plasmid, self).__init__(*contents if contents else ())
        self.name = name
        self.site = site
        self.marker = marker

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.__dict__))


class Feature(object):
    def __init__(self,
                 name=None,
                 type=None,
                 accession=None,
                 organism=None,
                 variant=None,
                 range=None):
        self.name = name
        self.type = type
        self.accession = accession
        self.organism = organism
        self.variant = variant
        self.range = range

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.__dict__))


class Organism(object):
    def __init__(self, name, aliases=None):
        self.name = name
        self.aliases = aliases or ()

    @property
    def default_alias(self):
        if len(self.aliases):
            return self.aliases[0]
        return self.name

    @classmethod
    def map(cls, organisms):
        dct = {}
        for organism in organisms:
            dct[organism.name] = organism
            if organism.aliases:
                for alias in organism.aliases:
                    dct[alias] = organism
        return dct

    def __eq__(self, other):
        return isinstance(other, Organism) and self.name == other.name

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.name))


class Accession(object):
    def __init__(self, identifier, database=None):
        self.identifier = identifier
        self.database = database

    def __eq__(self, other):
        return isinstance(other, Accession) and \
               self.database == other.database and \
               self.identifier == other.identifier

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.__dict__))


class Type(object):
    def __init__(self, name, aliases=None):
        self.name = name
        self.aliases = aliases or ()

    @property
    def default_alias(self):
        if len(self.aliases):
            return self.aliases[0]
        return self.name

    @classmethod
    def map(cls, types):
        dct = {}
        for type in types:
            dct[type.name] = type
            if type.aliases:
                for alias in type.aliases:
                    dct[alias] = type
        return dct

    def __eq__(self, other):
        return isinstance(other, Type) and self.name == other.name

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.name))