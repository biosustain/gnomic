#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CAVEAT UTILITOR
#
# This file was automatically generated by Grako.
#
#    https://pypi.python.org/pypi/grako/
#
# Any changes you make to it will be overwritten the next time
# the file is generated.


from __future__ import print_function, division, absolute_import, unicode_literals

from grako.buffering import Buffer
from grako.parsing import graken, Parser
from grako.util import re, RE_FLAGS, generic_main  # noqa


__all__ = [
    'GnomicParser',
    'GnomicSemantics',
    'main'
]

KEYWORDS = {}


class GnomicBuffer(Buffer):
    def __init__(
        self,
        text,
        whitespace=re.compile('[\\t ]+', RE_FLAGS | re.DOTALL),
        nameguard=None,
        comments_re=None,
        eol_comments_re=None,
        ignorecase=None,
        namechars='',
        **kwargs
    ):
        super(GnomicBuffer, self).__init__(
            text,
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            namechars=namechars,
            **kwargs
        )


class GnomicParser(Parser):
    def __init__(
        self,
        whitespace=re.compile('[\\t ]+', RE_FLAGS | re.DOTALL),
        nameguard=None,
        comments_re=None,
        eol_comments_re=None,
        ignorecase=None,
        left_recursion=False,
        parseinfo=True,
        keywords=None,
        namechars='',
        buffer_class=GnomicBuffer,
        **kwargs
    ):
        if keywords is None:
            keywords = KEYWORDS
        super(GnomicParser, self).__init__(
            whitespace=whitespace,
            nameguard=nameguard,
            comments_re=comments_re,
            eol_comments_re=eol_comments_re,
            ignorecase=ignorecase,
            left_recursion=left_recursion,
            parseinfo=parseinfo,
            keywords=keywords,
            namechars=namechars,
            buffer_class=buffer_class,
            **kwargs
        )

    @graken()
    def _SEQUENCE_VARIANT_(self):
        with self._choice():
            with self._option():
                self._DNA_SEQUENCE_VARIANT_()
            with self._option():
                self._PROTEIN_SEQUENCE_VARIANT_()
            self._error('no available options')

    @graken()
    def _DNA_SEQUENCE_VARIANT_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._token('g')
                with self._option():
                    self._token('c')
                with self._option():
                    self._token('n')
                self._error('expecting one of: c g n')
        self._token('.')
        with self._group():
            with self._choice():
                with self._option():
                    self._INTEGER_()
                    self._NUCLEOTIDE_()
                    self._token('>')
                    self._NUCLEOTIDE_()
                with self._option():
                    self._INTEGER_()
                    self._token('+')
                    self._INTEGER_()
                    self._NUCLEOTIDE_()
                    self._token('>')
                    self._NUCLEOTIDE_()
                with self._option():
                    self._token('[')
                    self._INTEGER_()
                    self._NUCLEOTIDE_()
                    self._token('>')
                    self._NUCLEOTIDE_()
                    self._token(';')
                    self._INTEGER_()
                    self._NUCLEOTIDE_()
                    self._token('>')
                    self._NUCLEOTIDE_()
                    self._token(']')
                with self._option():
                    self._INTEGER_()
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('=//')
                            with self._option():
                                self._token('=/')
                            self._error('expecting one of: =/ =//')
                    self._NUCLEOTIDE_()
                    self._token('>')
                    self._NUCLEOTIDE_()
                with self._option():
                    self._INTEGER_()
                    self._token('=')
                with self._option():
                    self._INTEGER_()
                    self._token('+')
                    self._INTEGER_()
                    self._token('del')
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('del=//del')
                            with self._option():
                                self._token('=/del')
                            self._error('expecting one of: =/del del=//del')
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('ins')
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._NUCLEOTIDE_SEQUENCE_()
                            with self._option():
                                self._token('(')
                                self._INTEGER_()
                                self._token(')')
                            self._error('no available options')
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('ins')
                    self._pattern(r'\w')
                    self._INTEGER_()
                    self._token('.')
                    self._INTEGER_()
                    self._token(':')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('ins')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('inv')
                    with self._optional():
                        self._INTEGER_()
                        self._token('_')
                        self._INTEGER_()
                        self._token('inv')
                with self._option():
                    self._token('(')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token(')ins')
                    self._NUCLEOTIDE_()
                    self._token('(')
                    self._PROTEIN_SEQUENCE_VARIANT_()
                    self._token(')')
                with self._option():
                    self._token('(')
                    self._INTEGER_()
                    self._token('+')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('-')
                    self._INTEGER_()
                    self._token(')_(')
                    self._INTEGER_()
                    self._token('+')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('-')
                    self._INTEGER_()
                    self._token(')')
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('dup')
                            with self._option():
                                self._token('[')
                                self._INTEGER_()
                                self._token(']')
                            self._error('expecting one of: dup')
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('con')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                with self._option():
                    self._INTEGER_()
                    self._token('delins')
                    self._NUCLEOTIDE_SEQUENCE_()
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('delins')
                    self._NUCLEOTIDE_SEQUENCE_()
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('[')
                    self._INTEGER_()
                    self._token(']')
                    with self._optional():
                        self._token(';[')
                        self._INTEGER_()
                        self._token(']')
                with self._option():
                    self._token('-')
                    self._INTEGER_()
                    self._token('_-')
                    self._INTEGER_()
                    self._token('[')
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('(')
                                self._INTEGER_()
                                self._token('_')
                                self._INTEGER_()
                                self._token(')')
                            with self._option():
                                self._INTEGER_()
                            self._error('no available options')
                    self._token(']')
                with self._option():
                    self._INTEGER_()
                    self._pattern(r'([ACGTBDHKMNRSVWY]{3}\[\d+\])+')
                with self._option():
                    self._INTEGER_()
                    self._token('-')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('-')
                    self._INTEGER_()
                    self._DEL_DUP_()
                with self._option():
                    self._INTEGER_()
                    self._DEL_DUP_()
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._DEL_DUP_()
                            with self._option():
                                self._token('inv')
                            self._error('expecting one of: inv')
                with self._option():
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('+')
                    self._INTEGER_()
                    self._DEL_DUP_()
                with self._option():
                    self._token('(?_-')
                    self._INTEGER_()
                    self._token(')_(*')
                    self._INTEGER_()
                    self._token('_?)')
                    self._DEL_DUP_()
                with self._option():
                    self._token('(?_-')
                    self._INTEGER_()
                    self._token(')_(')
                    self._INTEGER_()
                    self._token('+')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token('-')
                    self._INTEGER_()
                    self._token(')')
                    self._DEL_DUP_()
                self._error('no available options')

    @graken()
    def _DEL_DUP_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._token('del')
                with self._option():
                    self._token('dup')
                self._error('expecting one of: del dup')

    @graken()
    def _NUCLEOTIDE_(self):
        self._pattern(r'[ACGTBDHKMNRSVWY]')

    @graken()
    def _NUCLEOTIDE_SEQUENCE_(self):

        def block0():
            self._NUCLEOTIDE_()
        self._positive_closure(block0)

    @graken()
    def _PROTEIN_SEQUENCE_VARIANT_(self):
        self._token('p.')
        with self._group():
            with self._choice():
                with self._option():
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._AMINO_ACID_()
                            with self._option():
                                self._token('*')
                            self._error('expecting one of: *')
                    self._INTEGER_()
                    with self._optional():
                        self._AMINO_ACID_SEQUENCE_()
                    self._token('ext')
                    with self._group():
                        with self._choice():
                            with self._option():
                                with self._group():
                                    with self._choice():
                                        with self._option():
                                            self._token('-')
                                        with self._option():
                                            self._token('*')
                                        self._error('expecting one of: * -')
                                self._INTEGER_()
                            with self._option():
                                self._token('*?')
                            with self._option():
                                self._token('*')
                            self._error('expecting one of: * *?')
                with self._option():
                    self._token('(')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._AMINO_ACID_()
                    self._token('ext')
                    self._AMINO_ACID_SEQUENCE_()
                    self._token(')')
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._AMINO_ACID_()
                    self._token('fs')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._token('fs')
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._AMINO_ACID_()
                    self._token('fs*')
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._INTEGER_()
                            with self._option():
                                self._token('?')
                            self._error('expecting one of: ?')
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._token('delins')
                    self._AMINO_ACID_SEQUENCE_()
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._token('_')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('ins')
                            with self._option():
                                self._token('delins')
                            self._error('expecting one of: delins ins')
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._AMINO_ACID_SEQUENCE_()
                            with self._option():
                                self._INTEGER_()
                            self._error('no available options')
                with self._option():
                    self._token('(')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._token('_')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('ins')
                            with self._option():
                                self._token('delins')
                            self._error('expecting one of: delins ins')
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._AMINO_ACID_SEQUENCE_()
                            with self._option():
                                self._INTEGER_()
                            self._error('no available options')
                    self._token(')')
                with self._option():
                    self._token('[')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._AMINO_ACID_()
                    self._token(';')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._AMINO_ACID_()
                    self._token(']')
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._AMINO_ACID_()
                with self._option():
                    self._token('(')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('*')
                            with self._option():
                                self._token('Ter')
                            with self._option():
                                self._token('=')
                            with self._option():
                                self._token('?')
                            with self._option():
                                self._AMINO_ACID_()
                            self._error('expecting one of: * = ? Ter')
                    self._token(')')
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    with self._group():
                        with self._choice():
                            with self._option():
                                self._token('*')
                            with self._option():
                                self._token('Ter')
                            with self._option():
                                self._token('=')
                            with self._option():
                                self._token('?')
                            with self._option():
                                self._AMINO_ACID_()
                            self._error('expecting one of: * = ? Ter')
                with self._option():
                    self._token('0')
                with self._option():
                    self._token('?')
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._DEL_DUP_()
                with self._option():
                    self._token('(')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._DEL_DUP_()
                    self._token(')')
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._token('_')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._DEL_DUP_()
                with self._option():
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._token('[')
                    self._INTEGER_()
                    self._token(']')
                    with self._optional():
                        self._token(';[')
                        self._INTEGER_()
                        self._token(']')
                with self._option():
                    self._token('(')
                    self._AMINO_ACID_()
                    self._INTEGER_()
                    self._token(')[(')
                    self._INTEGER_()
                    self._token('_')
                    self._INTEGER_()
                    self._token(')]')
                self._error('expecting one of: 0 ?')

    @graken()
    def _AMINO_ACID_(self):
        self._pattern(r'[A-Z]([a-z]{2})?')

    @graken()
    def _AMINO_ACID_SEQUENCE_(self):

        def block0():
            self._AMINO_ACID_()
        self._positive_closure(block0)

    @graken()
    def _VARIABLE_VARIANT_(self):
        self._VARIABLE_VARIANT_IDENTIFIER_()
        self._token('=')
        self._VARIABLE_VARIANT_VALUE_()

    @graken()
    def _VARIABLE_VARIANT_IDENTIFIER_(self):
        self._pattern(r'[a-z0-9][a-zA-Z0-9]*(\.[a-z0-9][a-zA-Z0-9]*)*')

    @graken()
    def _VARIABLE_VARIANT_VALUE_(self):
        with self._choice():
            with self._option():
                self._NUMBER_()
            with self._option():
                self._QUOTED_STRING_()
            with self._option():
                self._UNQUOTED_STRING_()
            self._error('no available options')

    @graken()
    def _NUMBER_(self):
        self._pattern(r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?')

    @graken()
    def _QUOTED_STRING_(self):
        self._pattern(r'"(?:[^"\\]|\\.)*"')

    @graken()
    def _UNQUOTED_STRING_(self):
        self._pattern(r'[a-zA-Z0-9]+')

    @graken()
    def _start_(self):
        with self._optional():
            self._SEP_()
        with self._group():
            with self._choice():
                with self._option():
                    self._CHANGE_()
                    self.add_last_node_to_name('@')

                    def block1():
                        self._LIST_SEPARATOR_()
                        self._CHANGE_()
                        self.add_last_node_to_name('@')
                    self._closure(block1)
                with self._option():
                    self._empty_closure()
                self._error('no available options')
        with self._optional():
            self._SEP_()
        self._check_eof()

    @graken()
    def _CHANGE_(self):
        with self._choice():
            with self._option():
                self._INSERTION_()
            with self._option():
                self._REPLACEMENT_()
            with self._option():
                self._DELETION_()
            with self._option():
                self._PLASMID_()
            with self._option():
                self._PHENE_()
            self._error('no available options')

    @graken()
    def _INSERTION_(self):
        self._token('+')
        self._ANNOTATION_()
        self.name_last_node('after')
        self.ast._define(
            ['after'],
            []
        )

    @graken()
    def _REPLACEMENT_(self):
        with self._choice():
            with self._option():
                with self._group():
                    with self._choice():
                        with self._option():
                            self._ANNOTATION_AT_LOCUS_()
                        with self._option():
                            self._ANNOTATION_()
                        self._error('no available options')
                self.name_last_node('before')
                self._token('>')
                self.name_last_node('op')
                with self._group():
                    with self._choice():
                        with self._option():
                            self._PLASMID_()
                        with self._option():
                            self._ANNOTATION_()
                        self._error('no available options')
                self.name_last_node('after')
            with self._option():
                with self._group():
                    with self._choice():
                        with self._option():
                            self._ANNOTATION_AT_LOCUS_()
                        with self._option():
                            self._ANNOTATION_()
                        self._error('no available options')
                self.name_last_node('before')
                self._token('>>')
                self.name_last_node('op')
                with self._group():
                    with self._choice():
                        with self._option():
                            self._PLASMID_()
                        with self._option():
                            self._ANNOTATION_()
                        self._error('no available options')
                self.name_last_node('after')
            self._error('no available options')
        self.ast._define(
            ['after', 'before', 'op'],
            []
        )

    @graken()
    def _DELETION_(self):
        self._token('-')
        with self._group():
            with self._choice():
                with self._option():
                    self._PLASMID_()
                with self._option():
                    self._ANNOTATION_AT_LOCUS_()
                with self._option():
                    self._ANNOTATION_()
                self._error('no available options')
        self.name_last_node('before')
        self.ast._define(
            ['before'],
            []
        )

    @graken()
    def _PLASMID_(self):
        with self._choice():
            with self._option():
                self._token('(')
                self._IDENTIFIER_()
                self.name_last_node('name')
                self._SEP_()
                self._ANNOTATIONS_()
                self.name_last_node('annotations')
                self._token(')')
            with self._option():
                self._token('(')
                self._IDENTIFIER_()
                self.name_last_node('name')
                self._token(')')
            self._error('no available options')
        self.ast._define(
            ['annotations', 'name'],
            []
        )

    @graken()
    def _ANNOTATION_AT_LOCUS_(self):
        self._ANNOTATION_()
        self.name_last_node('annotation')
        self._token('@')
        self._FEATURE_()
        self.name_last_node('locus')
        self.ast._define(
            ['annotation', 'locus'],
            []
        )

    @graken()
    def _ANNOTATION_(self):
        with self._choice():
            with self._option():
                self._FUSION_()
            with self._option():
                self._FEATURE_()
            with self._option():
                self._FEATURE_SET_()
            with self._option():
                self._ANNOTATIONS_()
            self._error('no available options')

    @graken()
    def _ANNOTATIONS_(self):
        with self._optional():
            self._SEP_()
        with self._group():
            with self._choice():
                with self._option():
                    self._FEATURE_FUSION_()
                with self._option():
                    self._FEATURE_()
                self._error('no available options')
        self.add_last_node_to_name('@')

        def block2():
            self._LIST_SEPARATOR_()
            with self._group():
                with self._choice():
                    with self._option():
                        self._FEATURE_FUSION_()
                    with self._option():
                        self._FEATURE_()
                    self._error('no available options')
            self.add_last_node_to_name('@')
        self._closure(block2)
        with self._optional():
            self._SEP_()

    @graken()
    def _FEATURE_SET_(self):
        self._token('{')
        self._ANNOTATIONS_()
        self.name_last_node('@')
        self._token('}')

    @graken()
    def _FUSION_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._FEATURE_SET_()
                    self.add_last_node_to_name('@')
                with self._option():
                    self._FEATURE_()
                    self.add_last_node_to_name('@')
                self._error('no available options')

        def block3():
            self._token(':')
            with self._group():
                with self._choice():
                    with self._option():
                        self._FEATURE_SET_()
                        self.add_last_node_to_name('@')
                    with self._option():
                        self._FEATURE_()
                        self.add_last_node_to_name('@')
                    self._error('no available options')
        self._positive_closure(block3)

    @graken()
    def _FEATURE_FUSION_(self):
        self._FEATURE_()
        self.add_last_node_to_name('@')

        def block1():
            self._token(':')
            self._FEATURE_()
            self.add_last_node_to_name('@')
        self._positive_closure(block1)

    @graken()
    def _FEATURE_(self):
        with self._choice():
            with self._option():
                with self._optional():
                    self._FEATURE_ORGANISM_()
                    self.name_last_node('organism')
                with self._optional():
                    self._IDENTIFIER_()
                    self.name_last_node('type')
                    self._token('.')
                self._IDENTIFIER_()
                self.name_last_node('name')
                with self._optional():
                    self._ACCESSION_()
                    self.name_last_node('accession')
                with self._optional():
                    self._FEATURE_VARIANT_()
                    self.name_last_node('variant')
            with self._option():
                self._ACCESSION_()
                self.name_last_node('accession')
                with self._optional():
                    self._FEATURE_VARIANT_()
                    self.name_last_node('variant')
            self._error('no available options')
        self.ast._define(
            ['accession', 'name', 'organism', 'type', 'variant'],
            []
        )

    @graken()
    def _PHENE_(self):
        with self._choice():
            with self._option():
                with self._optional():
                    self._FEATURE_ORGANISM_()
                    self.name_last_node('organism')
                with self._optional():
                    self._IDENTIFIER_()
                    self.name_last_node('type')
                    self._token('.')
                self._IDENTIFIER_()
                self.name_last_node('name')
                with self._optional():
                    self._ACCESSION_()
                    self.name_last_node('accession')
                self._FEATURE_VARIANT_()
                self.name_last_node('variant')
            with self._option():
                self._ACCESSION_()
                self.name_last_node('accession')
                self._FEATURE_VARIANT_()
                self.name_last_node('variant')
            self._error('no available options')
        self.ast._define(
            ['accession', 'name', 'organism', 'type', 'variant'],
            []
        )

    @graken()
    def _FEATURE_ORGANISM_(self):
        self._ORGANISM_IDENTIFIER_()
        self.name_last_node('@')
        self._token('/')

    @graken()
    def _ORGANISM_IDENTIFIER_(self):
        self._pattern(r'[a-zA-Z0-9]+(\.[a-zA-Z0-9]+)?')

    @graken()
    def _ACCESSION_(self):
        with self._choice():
            with self._option():
                self._token('#')
                self._DATABASE_()
                self.name_last_node('db')
                self._token(':')
                with self._group():
                    with self._choice():
                        with self._option():
                            self._INTEGER_()
                        with self._option():
                            self._IDENTIFIER_()
                        self._error('no available options')
                self.name_last_node('id')
            with self._option():
                self._token('#')
                with self._group():
                    with self._choice():
                        with self._option():
                            self._INTEGER_()
                        with self._option():
                            self._IDENTIFIER_()
                        self._error('no available options')
                self.name_last_node('id')
            self._error('no available options')
        self.ast._define(
            ['db', 'id'],
            []
        )

    @graken()
    def _DATABASE_(self):
        self._pattern(r'[A-Za-z0-9-][A-Za-z0-9]+')

    @graken()
    def _INTEGER_(self):
        self._pattern(r'[0-9]+')
        self.name_last_node('@')

    @graken()
    def _FEATURE_VARIANT_(self):
        self._token('(')
        self._VARIANT_()
        self.add_last_node_to_name('@')

        def block1():
            with self._group():
                with self._choice():
                    with self._option():
                        self._token(',')
                    with self._option():
                        self._token(';')
                    self._error('expecting one of: , ;')
            with self._optional():
                self._SEP_()
            self._VARIANT_()
            self.add_last_node_to_name('@')
        self._closure(block1)
        self._token(')')

    @graken()
    def _VARIANT_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._VARIABLE_VARIANT_()
                with self._option():
                    self._SEQUENCE_VARIANT_()
                with self._option():
                    self._VARIANT_IDENTIFIER_()
                self._error('no available options')

    @graken()
    def _VARIANT_IDENTIFIER_(self):
        self._pattern(r'[A-Za-z0-9]+([A-Za-z0-9_\-]+[A-Za-z0-9])?')

    @graken()
    def _IDENTIFIER_(self):
        self._pattern(r'[a-zA-Z0-9]+([A-Za-z0-9_-]+[A-Za-z0-9])?')

    @graken()
    def _LIST_SEPARATOR_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._SEP_()
                with self._option():
                    self._token(',')
                    with self._optional():
                        self._SEP_()
                self._error('expecting one of: ,')

    @graken()
    def _SEP_(self):
        self._pattern(r'[\t ]+')


class GnomicSemantics(object):
    def SEQUENCE_VARIANT(self, ast):
        return ast

    def DNA_SEQUENCE_VARIANT(self, ast):
        return ast

    def DEL_DUP(self, ast):
        return ast

    def NUCLEOTIDE(self, ast):
        return ast

    def NUCLEOTIDE_SEQUENCE(self, ast):
        return ast

    def PROTEIN_SEQUENCE_VARIANT(self, ast):
        return ast

    def AMINO_ACID(self, ast):
        return ast

    def AMINO_ACID_SEQUENCE(self, ast):
        return ast

    def VARIABLE_VARIANT(self, ast):
        return ast

    def VARIABLE_VARIANT_IDENTIFIER(self, ast):
        return ast

    def VARIABLE_VARIANT_VALUE(self, ast):
        return ast

    def NUMBER(self, ast):
        return ast

    def QUOTED_STRING(self, ast):
        return ast

    def UNQUOTED_STRING(self, ast):
        return ast

    def start(self, ast):
        return ast

    def CHANGE(self, ast):
        return ast

    def INSERTION(self, ast):
        return ast

    def REPLACEMENT(self, ast):
        return ast

    def DELETION(self, ast):
        return ast

    def PLASMID(self, ast):
        return ast

    def ANNOTATION_AT_LOCUS(self, ast):
        return ast

    def ANNOTATION(self, ast):
        return ast

    def ANNOTATIONS(self, ast):
        return ast

    def FEATURE_SET(self, ast):
        return ast

    def FUSION(self, ast):
        return ast

    def FEATURE_FUSION(self, ast):
        return ast

    def FEATURE(self, ast):
        return ast

    def PHENE(self, ast):
        return ast

    def FEATURE_ORGANISM(self, ast):
        return ast

    def ORGANISM_IDENTIFIER(self, ast):
        return ast

    def ACCESSION(self, ast):
        return ast

    def DATABASE(self, ast):
        return ast

    def INTEGER(self, ast):
        return ast

    def FEATURE_VARIANT(self, ast):
        return ast

    def VARIANT(self, ast):
        return ast

    def VARIANT_IDENTIFIER(self, ast):
        return ast

    def IDENTIFIER(self, ast):
        return ast

    def LIST_SEPARATOR(self, ast):
        return ast

    def SEP(self, ast):
        return ast


def main(filename, startrule, **kwargs):
    with open(filename) as f:
        text = f.read()
    parser = GnomicParser(parseinfo=False)
    return parser.parse(text, startrule, filename=filename, **kwargs)


if __name__ == '__main__':
    import json
    ast = generic_main(main, GnomicParser, name='Gnomic')
    print('AST:')
    print(ast)
    print()
    print('JSON:')
    print(json.dumps(ast, indent=2))
    print()
