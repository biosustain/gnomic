

class Mutation(object):
    """

    .. attribute:: old

        Used either in a deletion or replacement. When used in a replacement, this must
        be a :class:`Feature`; when used in a deletion this can be any *insertable*.

        This may be a :class:`Plasmid` to signify a removed non-integrated plasmid. (Plasmid vectors that have been
        integrated cannot be removed directly --- their contents must be removed instead.)

    .. attribute:: new

        Used either in an insertion, or replacement.

    .. attribute:: markers

        A :class:`FeatureSet`, acting as a set of selection markers. The selection markers are generally inserted as part
        of the mutation.

    .. attribute:: multiple

        A :class:`bool`, specifying whether :attr:`old` is a multiple integration site.

        This attribute is only used with a replacement.

    """

    def __init__(self, old, new, markers=None, multiple=False):
        if isinstance(old, (list, tuple)):
            old = FeatureTree(*old)
        elif old and not isinstance(old, Plasmid):
            old = FeatureTree(old)
        if isinstance(new, (list, tuple)):
            new = FeatureTree(*new)
        elif new and not isinstance(new, Plasmid):
            new = FeatureTree(new)
        if markers is not None:
            markers = FeatureTree(*markers)
        self.old = old
        self.new = new
        self.markers = markers
        self.multiple = multiple

    def __eq__(self, other):
        return isinstance(other, Mutation) and \
            self.old == other.old and \
            self.new == other.new and \
            self.markers == other.markers and \
            self.multiple == other.multiple

    def __hash__(self):
        return hash(self.old) + \
               hash(self.new) + \
               hash(self.markers) + \
               hash(self.multiple)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join('{}={}'.format(key, repr(value))
                                         for key, value in self.__dict__.items() if value is not None))


def Ins(insert, **kwargs):
    return Mutation(None, insert, **kwargs)


def Sub(before, after, **kwargs):
    return Mutation(before, after, **kwargs)


def Del(delete, **kwargs):
    return Mutation(delete, None, **kwargs)


class MatchableMixin(object):
    def match(self, other):
        return self == other


class FeatureTree(object):
    """
    A grouping or collection of features and fusions.
    """
    def __init__(self, *contents):
        self.contents = contents

    def features(self):
        for item in self.contents:
            if isinstance(item, FeatureTree):
                # yield from item.features()
                for feature in item.features():
                    yield feature
            else:
                yield item

    def __len__(self):
        return len(self.contents)

    def __getitem__(self, item):
        return self.contents[item]

    def __iter__(self):
        return iter(self.contents)

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.contents == other.contents

    def __hash__(self):
        return hash(self.contents)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, ', '.join(map(repr, self.contents)))


class FeatureSet(FeatureTree, MatchableMixin):

    def match(self, other):
        if not isinstance(other, FeatureSet):
            return False

        # feature sets must have the same number of features to match
        if len(self) != len(other):
            return False

        return all(a.match(b) for a, b in zip(self.contents, other.contents))


class Fusion(FeatureTree, MatchableMixin):
    """
    A fusion must contain at least two features.
    """

    def __hash__(self):
        return hash(self.contents)

    def contains(self, other):
        if isinstance(other, Fusion):
            for i, feature in enumerate(self.contents):
                if feature == other[0]:
                    if self.contents[i:i + len(other)] == other.contents:
                        return True
            return False
        else:
            return other in self.contents

    def match(self, other):
        if not isinstance(other, Fusion):
            return False

        # fusions must have the same number of features to match
        if len(self) != len(other):
            return False

        return all(a.match(b) for a, b in zip(self.contents, other.contents))

    def __eq__(self, other):
        return isinstance(other, Fusion) and self.contents == other.contents


class Plasmid(FeatureTree, MatchableMixin):
    def __init__(self, name, contents, site=None, markers=None):
        super(Plasmid, self).__init__(*contents if contents else ())
        if markers is not None:
            markers = FeatureSet(*markers)
        self.name = name
        self.site = site
        self.markers = markers

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return isinstance(other, Plasmid) and self.name == other.name

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.__dict__))

    def match(self, other, match_contents=False):
        if not self == other:
            return False

        if match_contents and not super(Plasmid, self).match(other):
            return False

        return True


class Feature(MatchableMixin):
    def __init__(self, name=None, type=None, accession=None, organism=None, variant=None, range=None):
        self.name = name
        self.type = type
        self.accession = accession
        self.organism = organism
        self.variant = variant
        self.range = range

    @classmethod
    def parse(cls, string, *args, **kwargs):
        """

        ::

            >>> Feature.parse('MYO1#SGD:S000001065')
            Feature(name='MYO1', accession=Accession(database='SGD', identifier='S000001065'}))

        :param string:
        :param args:
        :param kwargs:
        :return:
        """
        from gnomic.grammar import GnomicParser
        from gnomic.semantics import DefaultSemantics

        parser = GnomicParser()
        semantics = DefaultSemantics(*args, **kwargs)
        return parser.parse(string,
                            whitespace='',
                            semantics=semantics,
                            rule_name='FEATURE')

    def match(self, other, match_variant=True):
        if not isinstance(other, Feature):
            return False

        if self.accession and other.accession:
            return self.accession == other.accession
        elif self.name:
            # only match features with the same name
            if self.name != other.name:
                return False

            # if an organism is specified, match only features with the same organism
            if self.organism and self.organism != other.organism:
                return False

            # TODO ranges
            if self.range and self.range != other.range:
                return False

            # if this feature has no variant, match any other feature; otherwise, match only features with the same
            # variant
            if not self.variant or match_variant is False:
                return True
            if self.variant == other.variant:
                return True
            return False
        else:
            # not enough information for any match
            return False

    def __eq__(self, other):
        if not isinstance(other, Feature):
            return False

        if self.accession and other.accession:
            return self.accession == other.accession and self.range == other.range
        elif self.name:
            return self.name == other.name and \
                   self.type == other.type and \
                   self.organism == other.organism and \
                   self.variant == other.variant and \
                   self.range == other.range
        # TODO range
        return False

    def __hash__(self):
        return hash(self.name) + \
               hash(self.type) + \
               hash(self.accession) + \
               hash(self.organism) + \
               hash(self.variant) + \
               hash(self.range)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join('{}={}'.format(key, repr(value))
                                         for key, value in self.__dict__.items() if value))


class Range(object):
    """
    An inclusive range at a coding (DNA), RNA or protein level.
    """
    def __init__(self, level, start, end):
        self.level = level
        self.start = start
        self.end = end

    def __hash__(self):
        return hash(self.level) + \
               hash(self.start) + \
               hash(self.end)

    def __len__(self):
        return self.end - self.start + 1

    def __eq__(self, other):
        return isinstance(other, Range) \
               and self.level == other.level \
               and self.start == other.start\
               and self.end == other.end

    def __repr__(self):
        if self.start == self.end:
            return '{}({}, {})'.format(self.__class__.__name__, repr(self.level), self.start)
        return '{}({}, {}, {})'.format(self.__class__.__name__, repr(self.level), self.start, self.end)


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

    def __hash__(self):
        return hash(self.name)

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

    def __hash__(self):
        return hash(self.identifier) + hash(self.database)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join('{}={}'.format(key, repr(value))
                                         for key, value in self.__dict__.items() if value))

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

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self.name))
