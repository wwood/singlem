import itertools


class SparseResultFormatter:
    def write(self, query_results, output_io, streaming=False):
        ''' If streaming, then do not sort the results before writing.'''
        output_io.write("\t".join([
            'query_name','query_sequence',
            'divergence','num_hits','sample',
            'marker','hit_sequence','taxonomy'])+"\n")

        if streaming:
            all_results = query_results
        else:
            all_results = list(query_results)
            all_results.sort(
                key=lambda res: (res.query.name, res.divergence, res.subject.count))
        for res in all_results:
            output_io.write("\t".join([
                res.query.name,
                res.query.sequence,
                str(res.divergence),
                str(res.subject.count),
                res.subject.sample_name,
                res.subject.marker,
                res.subject.sequence,
                str(res.subject.taxonomy)
            ])+"\n")
