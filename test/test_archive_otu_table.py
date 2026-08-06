#!/usr/bin/env python3

#=======================================================================
# Authors: Ben Woodcroft
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import sys, os, unittest, json, random
from io import StringIO

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

from singlem.archive_otu_table import ArchiveOtuTable, InsufficientArchiveOtuTableVersionException

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')


def otu(gene='4.11.22seqs', sample='sample1', sequence='ACGT', num_hits=1, coverage=1.0,
        taxonomy='Root', read_names=None, nucleotides_aligned=None, taxonomy_by_known=False,
        read_unaligned_sequences=None, equal_best_hit_taxonomies=None,
        taxonomy_assignment_method='diamond', reads_with_repaired_deletions=None):
    return [gene, sample, sequence, num_hits, coverage, taxonomy,
            read_names if read_names is not None else ['read1'],
            nucleotides_aligned if nucleotides_aligned is not None else [60],
            taxonomy_by_known,
            read_unaligned_sequences if read_unaligned_sequences is not None else ['AAAACCCC'],
            equal_best_hit_taxonomies, taxonomy_assignment_method,
            reads_with_repaired_deletions]


class Tests(unittest.TestCase):
    maxDiff = None

    def write(self, data, fields=None):
        table = ArchiveOtuTable()
        table.alignment_hmm_sha256s = ['hmm_sha']
        table.singlem_package_sha256s = ['spkg_sha']
        if fields is not None:
            table.fields = fields
        table.data = data
        output = StringIO()
        table.write_to(output)
        return output.getvalue()

    def round_trip(self, data, fields=None):
        return ArchiveOtuTable.read(StringIO(self.write(data, fields)))

    def assertRoundTrips(self, data):
        observed = self.round_trip(data)
        self.assertEqual(data, observed.data)
        self.assertEqual(ArchiveOtuTable.version, observed.version)
        self.assertEqual(ArchiveOtuTable.FIELDS, observed.fields)

    def test_writes_the_current_version(self):
        j = json.loads(self.write([otu()]))
        self.assertEqual(5, j['version'])
        self.assertEqual(ArchiveOtuTable.FIELDS, j['fields'])

    def test_round_trip_single_otu(self):
        self.assertRoundTrips([otu()])

    def test_round_trip_multiple_otus(self):
        self.assertRoundTrips([
            otu(gene='gene1', sequence='AAAA', read_names=['r1'], read_unaligned_sequences=['AAAACCCC']),
            otu(gene='gene2', sequence='CCCC', read_names=['r2', 'r3'], nucleotides_aligned=[60, 57],
                read_unaligned_sequences=['GGGGTTTT', 'TTTTGGGG'], num_hits=2),
        ])

    def test_round_trip_empty_table(self):
        observed = self.round_trip([])
        self.assertEqual([], observed.data)

    def test_round_trip_none_valued_fields(self):
        self.assertRoundTrips([otu(taxonomy_assignment_method=None, read_unaligned_sequences=None)])

    def test_columns_are_stored_column_wise(self):
        j = json.loads(self.write([otu(gene='gene1'), otu(gene='gene2')]))
        self.assertEqual(['gene1', 'gene2'], j['otus']['gene'])
        self.assertEqual(2, j['num_otus'])

    def test_constant_fields_are_hoisted_out_of_the_columns(self):
        j = json.loads(self.write([otu(gene='gene1'), otu(gene='gene2')]))
        # sample is the same in both OTUs, gene is not.
        self.assertEqual('sample1', j['constant_fields']['sample'])
        self.assertNotIn('sample', j['otus'])
        self.assertNotIn('gene', j['constant_fields'])

    def test_constant_fields_do_not_collapse_bool_onto_int(self):
        # True == 1 in Python, so a column of [1, True] must not be hoisted.
        data = [otu(num_hits=1), otu(num_hits=True)]
        self.assertNotIn('num_hits', json.loads(self.write(data))['constant_fields'])
        self.assertEqual([1, True], [row[3] for row in self.round_trip(data).data])

    def test_list_valued_fields_are_never_hoisted(self):
        j = json.loads(self.write([otu(read_names=['r1']), otu(read_names=['r1'])]))
        self.assertNotIn('read_names', j['constant_fields'])

    def test_reads_are_dereplicated(self):
        data = [
            otu(gene='gene1', read_names=['r1', 'r2'], nucleotides_aligned=[60, 60],
                read_unaligned_sequences=['AAAA', 'CCCC']),
            # r1 also hits gene2, so its sequence is stored once, not twice.
            otu(gene='gene2', read_names=['r1'], read_unaligned_sequences=['AAAA']),
        ]
        j = json.loads(self.write(data))
        self.assertEqual(['AAAA', 'CCCC'], j['reads'])
        self.assertEqual([[ArchiveOtuTable.NEXT_NEW_READ, ArchiveOtuTable.NEXT_NEW_READ], [1]],
                         j['otus']['read_unaligned_sequences'])
        self.assertEqual(data, self.round_trip(data).data)

    def test_round_trip_read_repeated_within_one_otu(self):
        data = [otu(read_names=['r1', 'r2'], nucleotides_aligned=[60, 60],
                    read_unaligned_sequences=['AAAA', 'AAAA'])]
        self.assertEqual(['AAAA'], json.loads(self.write(data))['reads'])
        self.assertEqual(data, self.round_trip(data).data)

    def test_round_trip_otu_with_no_reads_between_otus_with_reads(self):
        data = [
            otu(gene='gene1', read_unaligned_sequences=['AAAA']),
            otu(gene='gene2', read_unaligned_sequences=None),
            otu(gene='gene3', read_unaligned_sequences=['CCCC']),
        ]
        self.assertEqual(data, self.round_trip(data).data)

    def test_short_rows_are_padded_to_the_field_width(self):
        # Entries parsed from an older archive can be added to a table declaring
        # the current field list; the fields they lack are unknown, not empty.
        version1_row = ['4.11.22seqs', 'sample1', 'ACGT', 1, 1.0, 'Root', ['read1'], [60], False]
        observed = self.round_trip([version1_row])
        self.assertEqual([version1_row + [None, None, None, None]], observed.data)

    def test_row_longer_than_fields_is_rejected(self):
        with self.assertRaises(Exception):
            self.write([otu() + ['extra']])

    def test_randomised_round_trips(self):
        random.seed(1)
        reads = ['AAAA', 'CCCC', 'GGGG', 'TTTT']
        for _ in range(200):
            data = []
            for _ in range(random.randint(0, 8)):
                n = random.randint(1, 4)
                data.append(otu(
                    gene=random.choice(['gene1', 'gene2']),
                    sample=random.choice(['s1', 's2']),
                    sequence=random.choice(['ACGT', 'TGCA']),
                    num_hits=n,
                    taxonomy_by_known=random.choice([True, False]),
                    read_names=['r%i' % random.randint(0, 20) for _ in range(n)],
                    nucleotides_aligned=[60] * n,
                    read_unaligned_sequences=random.choice(
                        [None, [random.choice(reads) for _ in range(n)]])))
            self.assertEqual(data, self.round_trip(data).data)

    def test_reads_version_4(self):
        with open(os.path.join(path_to_data, 'SRR8653040.json')) as f:
            observed = ArchiveOtuTable.read(f)
        self.assertEqual(4, observed.version)
        self.assertEqual(95, len(observed.data))
        self.assertEqual(ArchiveOtuTable.FIELDS_VERSION4, observed.fields)

    def test_version_4_and_version_5_hold_the_same_otus(self):
        with open(os.path.join(path_to_data, 'SRR8653040.json')) as f:
            version_4 = ArchiveOtuTable.read(f)
        version_5 = ArchiveOtuTable.read(StringIO(self.write(version_4.data)))
        # Version 5 adds reads_with_repaired_deletions, which a version 4 archive
        # does not record, so its OTUs gain a None for it.
        self.assertEqual([row + [None] for row in version_4.data], version_5.data)

    def test_version_4_tables_are_still_written_row_wise(self):
        table = ArchiveOtuTable()
        table.version = 4
        table.fields = ArchiveOtuTable.FIELDS_VERSION4
        table.alignment_hmm_sha256s = ['hmm_sha']
        table.singlem_package_sha256s = ['spkg_sha']
        table.data = [otu()[:len(ArchiveOtuTable.FIELDS_VERSION4)]]
        output = StringIO()
        table.write_to(output)
        j = json.loads(output.getvalue())
        self.assertEqual(4, j['version'])
        self.assertEqual([otu()[:len(ArchiveOtuTable.FIELDS_VERSION4)]], j['otus'])
        self.assertNotIn('constant_fields', j)

    def test_min_version(self):
        with self.assertRaises(InsufficientArchiveOtuTableVersionException):
            ArchiveOtuTable.read(StringIO(self.write([otu()])), min_version=6)
        ArchiveOtuTable.read(StringIO(self.write([otu()])), min_version=5)

    def test_unknown_version_rejected(self):
        j = json.loads(self.write([otu()]))
        j['version'] = 99
        with self.assertRaises(Exception):
            ArchiveOtuTable.read(StringIO(json.dumps(j)))

    def test_missing_field_rejected(self):
        j = json.loads(self.write([otu(gene='gene1'), otu(gene='gene2')]))
        del j['otus']['gene']
        with self.assertRaises(Exception):
            ArchiveOtuTable.read(StringIO(json.dumps(j)))

    def test_wrong_length_column_rejected(self):
        j = json.loads(self.write([otu(gene='gene1'), otu(gene='gene2')]))
        j['otus']['gene'] = ['gene1']
        with self.assertRaises(Exception):
            ArchiveOtuTable.read(StringIO(json.dumps(j)))

    def test_out_of_range_read_reference_rejected(self):
        j = json.loads(self.write([otu()]))
        j['otus']['read_unaligned_sequences'] = [[99]]
        with self.assertRaises(Exception):
            ArchiveOtuTable.read(StringIO(json.dumps(j)))


if __name__ == "__main__":
    unittest.main()
