import itertools


class SparseResultFormatter:
    def write(self, query_results, output_io):
        output_io.write("\t".join([
            'query_name','query_sequence',
            'divergence','num_hits','sample',
            'marker','hit_sequence','taxonomy'])+"\n")

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
                res.subject.taxonomy
            ])+"\n")
