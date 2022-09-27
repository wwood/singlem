import itertools
from .sequence_database import SequenceDatabase

class SparseResultFormatter:
    def write(self, query_results, output_io, sequence_type, streaming=False):
        ''' If streaming, then do not sort the results before writing.'''
        if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
            output_io.write("\t".join([
                'query_name','query_sequence',
                'divergence','num_hits','coverage','sample',
                'marker','hit_sequence','taxonomy'])+"\n")
        elif sequence_type == SequenceDatabase.PROTEIN_TYPE:
            output_io.write("\t".join([
                'query_name','query_sequence','query_protein_sequence',
                'divergence','num_hits','coverage','sample',
                'marker','hit_sequence','hit_protein_sequence','taxonomy'])+"\n")
        else:
            raise Exception("Unexpected sequence type: %s" % sequence_type)


        if streaming:
            all_results = query_results
        else:
            all_results = list(query_results)
            all_results.sort(
                key=lambda res: (res.query.name, res.divergence, res.subject.count))
        for res in all_results:
            if sequence_type == SequenceDatabase.NUCLEOTIDE_TYPE:
                output_io.write("\t".join([
                    res.query.name,
                    res.query.sequence,
                    str(res.divergence),
                    str(res.subject.count),
                    str(res.subject.coverage),
                    res.subject.sample_name,
                    res.subject.marker,
                    res.subject.sequence,
                    str(res.subject.taxonomy)
                ])+"\n")
            else: #protein
                output_io.write("\t".join([
                    res.query.name,
                    res.query.sequence,
                    res.query_protein_sequence,
                    str(res.divergence),
                    str(res.subject.count),
                    str(res.subject.coverage),
                    res.subject.sample_name,
                    res.subject.marker,
                    res.subject.sequence,
                    res.subject_protein_sequence,
                    str(res.subject.taxonomy)
                ])+"\n")

