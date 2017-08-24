from __future__ import unicode_literals

from abc import ABCMeta
from itertools import chain

import six


class Change(object):
    def __init__(self, before=None, after=None, multiple=False):
        if before is None and after is None:
            raise ValueError()

        self.before = before
        self.after = after
        self.multiple = multiple

    @classmethod
    def parse(cls, gnomic_change_string):
        from gnomic.grammar import GnomicParser
        from gnomic.semantics import DefaultSemantics

        parser = GnomicParser()
        semantics = DefaultSemantics()
        return parser.parse(gnomic_change_string,
                            whitespace='',
                            semantics=semantics,
                            rule_name='CHANGE')

    def is_presence(self):
        if self.before is None and isinstance(self.after, Plasmid):
            return True
        return self.before and self.after and self.before.match(self.after)

    def __eq__(self, other):
        return isinstance(other, Change) and \
               self.before == other.before and \
               self.after == other.after and \
               self.multiple == other.multiple

    def __ne__(self, other):
        return (not isinstance(other, Change)) or \
               self.before != other.before or \
               self.after != other.after or \
               self.multiple != other.multiple

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join('{}={}'.format(key, repr(value))
                                         for key, value in self.__dict__.items() if value is not None))

    def __mod__(self, locus):
        return self.__matmul__(locus)

    def __matmul__(self, locus):
        if not isinstance(locus, Annotation):
            raise ValueError()
        if not self.before:
            raise ValueError('A change can only occur at a locus if it has a "before" part.')
        if isinstance(self.before, AtLocus):
            raise ValueError('Change "before" is already at a locus')
        return Change(before=AtLocus(self.before, locus),
                      after=self.after,
                      multiple=self.multiple)

    def __str__(self):
        if self.after is None:
            return '-{!s}'.format(self.before)
        elif self.before is None:
            return '+{!s}'.format(self.after)
        elif self.multiple:
            return '{!s}>>{!s}'.format(self.before, self.after)
        else:
            return '{!s}>{!s}'.format(self.before, self.after)

    def __hash__(self):
        return hash(self.before) + hash(self.after) + hash(self.multiple)


class Present(Change):
    def __init__(self, annotation):
        super(Present, self).__init__(after=annotation)


class Annotation(object):
    def match(self, other, match_variants=True):
        return False

    def __hash__(self):
        return hash(id(self))

    def __pow__(self, other):
        return Fusion(self, other)

    def __pos__(self):
        return Change(after=self)

    def __neg__(self):
        return Change(before=self)

    def __gt__(self, other):
        if not isinstance(other, Annotation):
            raise ValueError()
        return Change(self, other)

    def __rshift__(self, other):
        if not isinstance(other, Annotation):
            raise ValueError()
        return Change(self, other, multiple=True)

    def __mod__(self, locus):
        return self.__matmul__(locus)

    def __matmul__(self, locus):
        if not isinstance(locus, Annotation):
            raise ValueError()
        return AtLocus(self, locus)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join('{}={}'.format(key, repr(value))
                                         for key, value in self.__dict__.items() if value is not None))


class AtLocus(Annotation):
    def __init__(self, annotation, locus):
        if not isinstance(annotation, Annotation):
            raise ValueError()

        if not isinstance(locus, Annotation):
            raise ValueError()

        self.annotation = annotation
        self.locus = locus

    def __neg__(self):
        return Change(self)

    def __gt__(self, other):
        if not isinstance(other, Annotation):
            raise ValueError()
        return Change(self, other)

    def __rshift__(self, other):
        if not isinstance(other, Annotation):
            raise ValueError()
        return Change(self, other, multiple=True)

    def __hash__(self):
        return hash(self.annotation) + hash(self.locus)

    def __eq__(self, other):
        return isinstance(other, AtLocus) and \
               self.annotation == other.annotation and \
               self.locus == other.locus

    def __ne__(self, other):
        return (not isinstance(other, AtLocus)) or \
               self.annotation != other.annotation or \
               self.locus != other.locus

    def __str__(self):
        return '{!s}@{!s}'.format(self.annotation, self.locus)


class Feature(Annotation):
    """

    -Feature("foo")\

    Feature("site") >> Feature("insertion")

    """

    def __init__(self, name=None, type=None, accession=None, organism=None, variant=None):
        self.name = name
        self.type = type
        self.accession = accession
        self.organism = organism
        self.variant = variant

    @classmethod
    def parse(cls, gnomic_feature_string):
        from gnomic.grammar import GnomicParser
        from gnomic.semantics import DefaultSemantics

        if not isinstance(gnomic_feature_string, six.string_types):
            raise ValueError('"gnomic_feature_string" must a string, got {}'.format(repr(gnomic_feature_string)))

        parser = GnomicParser()
        semantics = DefaultSemantics()
        return parser.parse(gnomic_feature_string,
                            whitespace='',
                            semantics=semantics,
                            rule_name='FEATURE')

    def match(self, other, match_variants=True):
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

            # if this feature has no variant, match any other feature; otherwise, match only features with the same
            # variant
            if not self.variant or match_variants is False:
                return True
            if self.variant == other.variant:
                return True
            return False
        else:
            # not enough information for any match
            return False

    def __hash__(self):
        return hash(self.name) + \
               hash(self.type) + \
               hash(self.accession) + \
               hash(self.organism) + \
               hash(self.variant)

    def __eq__(self, other):
        if not isinstance(other, Feature):
            return False

        if self.accession and other.accession:
            return self.accession == other.accession
        elif self.name:
            return self.name == other.name and \
                   self.type == other.type and \
                   self.organism == other.organism and \
                   self.variant == other.variant
        return False

    def __ne__(self, other):
        if not isinstance(other, Feature):
            return True

        if self.accession and other.accession:
            return self.accession != other.accession
        elif self.name:
            return self.name != other.name or \
                   self.type != other.type or \
                   self.organism != other.organism or \
                   self.variant != other.variant
        return True

    def __str__(self):
        s = ''
        if self.organism:
            s += '{}/'.format(self.organism)
        if self.type:
            s += '{}.'.format(self.type)
        if self.name:
            s += self.name
        if self.accession:
            s += '#{}'.format(self.accession)
        if self.variant:
            s += '({})'.format('; '.join(self.variant))
        return s


class CompositeAnnotationBase(six.with_metaclass(ABCMeta, Annotation)):
    def __init__(self, *annotations):
        if not all(isinstance(annotation, Annotation) for annotation in annotations):
            raise ValueError()

        self.annotations = annotations

    def __hash__(self):
        return hash(self.annotations)

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.annotations == other.annotations

    def __ne__(self, other):
        return (not isinstance(other, self.__class__)) or (self.annotations != other.annotations)

    def __len__(self):
        return len(self.annotations)

    def __getitem__(self, item):
        return self.annotations[item]

    def __iter__(self):
        for annotation in self.annotations:
            yield annotation

    def contains(self, other):
        return other in self.annotations

    def features(self):
        return chain(*(a.features() if isinstance(a, CompositeAnnotationBase) else [a] for a in self.annotations))


class CompositeAnnotation(CompositeAnnotationBase):
    def __init__(self, *annotations):
        super(CompositeAnnotation, self).__init__(*chain(*(annotation.annotations
                                                           if isinstance(annotation, CompositeAnnotation)
                                                           else (annotation,)
                                                           for annotation in annotations)))

    def __str__(self):
        return '{{{}}}'.format(' '.join(map(str, self.annotations)))


class Fusion(CompositeAnnotationBase):
    def __init__(self, *annotations):
        annotations = tuple(chain(*(annotation.annotations
                                    if isinstance(annotation, Fusion)
                                    else (annotation,)
                                    for annotation in annotations)))

        if len(annotations) < 2:
            raise ValueError()

        super(Fusion, self).__init__(*annotations)

    @classmethod
    def fuse(cls, annotations):
        if len(annotations) == 0:
            return None
        if len(annotations) == 1:
            return annotations[0]
        return Fusion(*annotations)

    def match(self, other, match_variants=True):
        if not isinstance(other, Fusion):
            return False

        # fusions must have the same number of features to match
        if len(self) != len(other):
            return False

        return all(a.match(b, match_variants=match_variants) for a, b in zip(self.annotations, other.annotations))

    def index(self, other):
        if isinstance(other, Fusion):
            for i, annotation in enumerate(self.annotations):
                if annotation == other[0]:
                    if self.annotations[i:i + len(other)] == other.annotations:
                        return i
            raise ValueError('{} is not in Fusion'.format(other))
        else:
            return self.annotations.index(other)

    def contains(self, other):
        if isinstance(other, Fusion):
            for i, annotation in enumerate(self.annotations):
                if annotation == other[0]:
                    if self.annotations[i:i + len(other)] == other.annotations:
                        return True
            return False
        else:
            return other in self.annotations

    def __str__(self):
        return ':'.join(map(str, self.annotations))


class Plasmid(CompositeAnnotationBase):
    def __init__(self, name, annotations=()):
        if any(isinstance(annotation, Plasmid) for annotation in annotations):
            raise ValueError()

        super(Plasmid, self).__init__(*annotations)
        self.name = name

    def __eq__(self, other):
        return isinstance(other, Plasmid) and self.name == other.name

    def __ne__(self, other):
        return (not isinstance(other, Plasmid)) or (self.name != other.name)

    def __hash__(self):
        return hash(self.name)

    def match(self, other, match_variants=True):
        if not isinstance(other, Plasmid):
            return False

        return self.name == other.name

    def __bool__(self):
        if self.name:
            return True
        return False

    __nonzero__ = __bool__

    def __str__(self):
        if self.annotations:
            return '({} {})'.format(self.name, ' '.join(map(str, self.annotations)))
        return '({})'.format(self.name)


class Accession(object):
    def __init__(self, identifier, database=None):
        self.identifier = identifier
        self.database = database

    def __eq__(self, other):
        return isinstance(other, Accession) and \
               self.database == other.database and \
               self.identifier == other.identifier

    def __ne__(self, other):
        return (not isinstance(other, Accession)) or \
               self.database != other.database or \
               self.identifier != other.identifier

    def __hash__(self):
        return hash(self.identifier) + hash(self.database)

    def __repr__(self):
        if self.database:
            return '{}:{}'.format(self.database, self.identifier)
        else:
            return '{}'.format(self.identifier)
