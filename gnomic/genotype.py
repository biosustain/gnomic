import itertools

import six
from grako.exceptions import GrakoException

from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics
from gnomic.types import Plasmid, Change, Fusion, CompositeAnnotation, AtLocus, Feature, CompositeAnnotationBase
from gnomic.formatters import BUILTIN_FORMATTERS


def partial_match(search, target):
    """
    Check that something in ``target`` matches ``search``.

    :param search:
    :param target:
    :return:
    """
    if isinstance(target, CompositeAnnotationBase):
        if target.contains(search):
            return True
    return search.match(target)


def change_annotation(annotation, site=None, replacement=None):
    """
    Apply a change to an annotation.

    - If ``site`` is specified, looks for ``site`` and replaces with ``replacement``.
    - If ``site`` is not specified, create a :class:`CompositeAnnotation` of ``annotation`` and ``replacement``.
    """
    if not (site or replacement):
        raise ValueError()

    if site:
        if isinstance(annotation, (Fusion, CompositeAnnotation, Plasmid)):
            if not annotation.contains(site):
                if replacement:
                    return CompositeAnnotation(annotation, replacement)
                return annotation
            elif isinstance(annotation, (CompositeAnnotation, Plasmid)) \
                    or isinstance(site, (Feature, CompositeAnnotation)):
                if replacement is None:
                    annotations = tuple(b for b in annotation if b != site)
                else:
                    annotations = tuple(replacement if b == site else b for b in annotation if b != replacement)

                if isinstance(annotation, Fusion):
                    return Fusion.fuse(annotations)
                return CompositeAnnotation(*annotations)
            elif isinstance(site, Fusion):
                before_index = annotation.index(site)
                upstream = annotation.annotations[:before_index]
                downstream = annotation.annotations[before_index + len(site):]

                if replacement is None:
                    return Fusion.fuse(upstream + downstream)
                else:
                    return Fusion.fuse(upstream + (replacement,) + downstream)
            else:
                raise NotImplementedError()
        elif isinstance(annotation, Feature):
            if annotation == site:
                return replacement
            elif replacement:
                return CompositeAnnotation(annotation, replacement)
            return annotation
        else:
            raise NotImplementedError()

    elif isinstance(annotation, (Feature, Fusion)):
        return CompositeAnnotation(annotation, replacement)
    elif isinstance(annotation, (CompositeAnnotation, Plasmid)):
        if any(a == replacement for a in annotation):  # ignore duplicates
            return annotation
        return CompositeAnnotation(annotation, replacement)
    else:
        raise NotImplementedError()


class GenotypeState(object):
    def __init__(self, changes=()):
        self._changes = [(change.before, change.after) for change in changes]

    def insert(self, annotation, multiple=False):
        # skip repeated insertions
        # e.g. +annotation
        for before, after in self._changes:
            if before is None and after == annotation:
                return

        # e.g. -annotation
        matches = [(before, after) for before, after in self._changes
                   if after is None and annotation.match(before)]

        if len(matches) == 1 or len(matches) > 1 and multiple:
            for match in matches:
                self._changes.remove(match)
            return

        # e.g. +annotation(foo)
        matches = [(before, after) for before, after in self._changes
                   if before is None and annotation.match(after, match_variants=False)]

        # if len(matches) == 0:
        #     matches = [(before, after) for before, after in self._changes
        #                if before is None and annotation.match(after)]
        #
        # else:
        if len(matches) <= 1 or multiple:
            for match in matches:
                self._changes.remove(match)

        self._changes.append((None, annotation))

    def remove(self, site, multiple=False):
        # skip repeated deletions
        for before, after in self._changes:
            if before == site and after is None:
                return  # e.g. -gene.A, -gene.A

        if isinstance(site, AtLocus):
            # e.g. gene.A>gene.X -gene.X@gene.A
            matches = [(before, after) for before, after in self._changes
                       if before and after
                       and (isinstance(before, AtLocus) and site.locus == before.locus
                            or site.locus == before)]

            if len(matches) == 0:
                self._changes.append((site, None))
            elif len(matches) == 1:
                before, after = matches[0]

                after = change_annotation(after, site.annotation, None)

                self._changes.remove(matches[0])

                if isinstance(before, AtLocus):
                    if before.annotation != after:
                        self._changes.append((before, after))
                elif before != after:
                    self._changes.append((before, after))
            else:
                raise NotImplementedError
        else:
            # e.g. +site, foo>site
            matches = [(before, after) for before, after in self._changes
                       if site == after]

            if len(matches) == 1 or len(matches) > 1 and multiple:
                for (before, after) in matches:
                    self._changes.remove((before, after))
                    # recursive?
                    if before:
                        self._changes.append((before, None))
                return

            # e.g. -site
            matches = [(before, after) for before, after in self._changes
                       if after is None and site.match(before, match_variants=False)]

            if len(matches) == 1 or len(matches) > 1 and multiple:
                for match in matches:
                    self._changes.remove(match)
                return

            # e.g. +site
            matches = [(before, after) for before, after in self._changes
                       if before is None and site.match(after)]

            if len(matches) <= 1 or multiple:
                for match in matches:
                    self._changes.remove(match)

            if len(matches) == 0:
                self._changes.append((site, None))

    def replace(self, site, replacement, multiple=False):
        # e.g. gene.A>gene.B
        assert site and replacement

        # skip repeated replacements
        # XXX possibility of different behavior with multiple=True
        for before, after in self._changes:
            if site == before and replacement == after:
                return

        # skip changes without effect
        # e.g. gene.A>gene.A or gene.A@foo>gene.A
        if site == replacement \
                or isinstance(site, AtLocus) and site.annotation == replacement:
            return

        # e.g. gene.A>gene.X gene.A>gene.Y or
        #      gene.X@gene.A>gene.Y gene.X@gene.A>gene.Y
        matches = [(before, after) for before, after in self._changes
                   if before and site == before]

        if len(matches) == 0:
            if isinstance(site, AtLocus):
                # e.g. gene.A>gene.X gene.X@gene.A>gene.Y
                matches = [(before, after) for before, after in self._changes
                           if before and (isinstance(before, AtLocus) and site.locus == before.locus
                                          or site.locus == before)]
            else:
                # e.g. gene.A>gene.X gene.X>gene.Y
                matches = [(before, after) for before, after in self._changes
                           if after and partial_match(site, after)]

            if len(matches) == 0:
                self._changes.append((site, replacement))
            elif len(matches) == 1:

                before, after = matches[0]

                if isinstance(site, AtLocus):
                    after = change_annotation(after, site.annotation, replacement)
                else:
                    after = change_annotation(after, site, replacement)

                self._changes.remove(matches[0])

                if isinstance(before, AtLocus):
                    if before.annotation != after:
                        self._changes.append((before, after))
                elif before != after:
                    self._changes.append((before, after))
            else:
                # TODO
                raise NotImplementedError()

        elif len(matches) == 1:
            self._changes.remove(matches[0])
            self._changes.append((site, replacement))
        else:
            # TODO
            raise NotImplementedError()

    def change(self, change):
        # XXX consider not simplifying (foo@foo>... and -foo@foo) if there is a logical use
        # simplify change at locus where annotation and locus are the same (e.g. foo@foo>... or -foo@foo)
        if change.before and isinstance(change.before, AtLocus) and change.before.annotation == change.before.locus:
            change = Change(change.before.annotation, change.after, multiple=change.multiple)

        if change.before is None:
            # e.g. +gene.A
            self.insert(change.after, multiple=change.multiple)
        elif change.after is None:
            # e.g. -gene.A
            self.remove(change.before, multiple=change.multiple)
        else:
            # e.g. gene.A>gene.B
            assert change.before and change.after
            self.replace(change.before, change.after, multiple=change.multiple)

    @property
    def changes(self):
        return tuple(Change(before, after) for before, after in self._changes)


class Genotype(object):
    def __init__(self, changes, parent=None):
        if parent:
            state = GenotypeState(parent.state.changes)
        else:
            state = GenotypeState()

        for change in changes:
            state.change(change)

        self.parent = parent
        self.state = state

    @classmethod
    def _parse_gnomic_string(cls, gnomic_string, *args, **kwargs):
        parser = GnomicParser()
        semantics = DefaultSemantics(*args, **kwargs)
        return parser.parse(gnomic_string,
                            whitespace='',
                            semantics=semantics,
                            rule_name='start')

    @classmethod
    def parse(cls, gnomic_string, parent=None, **kwargs):
        if not isinstance(gnomic_string, six.string_types):
            raise ValueError('"gnomic_string" must a string, got {}'.format(repr(gnomic_string)))

        changes = Genotype._parse_gnomic_string(gnomic_string, **kwargs)
        return Genotype(changes, parent=parent, **kwargs)

    @classmethod
    def is_valid(cls, gnomic_string, **kwargs):
        """
        Tests whether a gnomic genotype definition can be parsed.
        """
        try:
            cls._parse_gnomic_string(gnomic_string, **kwargs)
            return True
        except GrakoException:
            return False

    @property
    def added_features(self):
        return {feature for feature in itertools.chain(*(
            iter(change.after.features()) if isinstance(change.after, CompositeAnnotationBase) else (change.after,)
            for change in self.state.changes if change.after is not None))}

    @property
    def removed_features(self):
        return {feature for feature in itertools.chain(*(
            iter(change.before.features()) if isinstance(change.before, CompositeAnnotationBase) else (change.before,)
            for change in self.state.changes if change.before is not None))}

    @property
    def added_plasmids(self):
        return {change.after
                for change in self.state.changes
                if change.before is None and isinstance(change.after, Plasmid)}

    @property
    def removed_plasmids(self):
        return {change.before
                for change in self.state.changes
                if change.after is None and isinstance(change.before, Plasmid)}

    @property
    def added_fusions(self):
        return {change.after
                for change in self.state.changes
                if change.before is None and isinstance(change.after, Fusion)}

    @property
    def removed_fusions(self):
        return {change.before
                for change in self.state.changes
                if change.after is None and isinstance(change.before, Fusion)}

    @property
    def added_fusion_features(self):
        return {change.after
                for change in self.state.changes
                if change.before is None and isinstance(change.after, (Feature, Fusion, CompositeAnnotation))}

    @property
    def removed_fusion_features(self):
        return {change.before
                for change in self.state.changes
                if change.after is None and isinstance(change.before, (Feature, Fusion, CompositeAnnotation))}

    def changes(self):
        return self.state.changes

    def format(self, output='text'):
        return BUILTIN_FORMATTERS[output].format_genotype(self)
