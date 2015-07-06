import itertools


class QueryResultFormatter:
    def __init__(self, divergences_and_subjects):
        self.divergences_and_subjects = divergences_and_subjects

class SparseResultFormatter(QueryResultFormatter):
    def write(self, output_io):
        # sort the array according to increasing divergence, then increasing number of hits
        sortable_list = [(x[0], x[1].count, x[1]) for x in self.divergences_and_subjects]
        sortable_list.sort(key=lambda x: [x[0],-x[1]])
        
        # print out array
        output_io.write("\t".join(['divergence','num_hits','sample',
                                   'marker','sequence','taxonomy'])+"\n")
        for hit in sortable_list:
            subject = hit[2]
            output_io.write("\t".join([str(hit[0]),
                             str(subject.count), 
                             subject.sample_name,
                             subject.marker, 
                             subject.sequence,
                             subject.taxonomy])+"\n")

class DenseResultFormatter(QueryResultFormatter):
    def write(self, output_io):
        # make array of sequence arrays, where index in the First
        # array is divergence
        max_divergence = max([x[0] for x in self.divergences_and_subjects])
        seq_sets = []
        for _ in range(max_divergence+1):
            seq_sets.append(set())
        seq_to_subjects = {}
        samples = set()
        for x in self.divergences_and_subjects:
            i = x[0]
            subject = x[1]
            samples.add(subject.sample_name)
            seq = subject.sequence
            
            seq_sets[i].add(seq)
            try:
                seq_to_subjects[seq].append(subject)
            except KeyError:
                seq_to_subjects[seq] = [subject]
            
        # Build each row as best as possible
        sample_name_order = []
        to_print = []
        for divergence, seq_set in enumerate(seq_sets):
            if len(seq_set)==0: continue
            
            # Determine order of sequences - put most common at top
            total_count_and_seq = []
            for seq in seq_set:
                total = 0
                for subject in seq_to_subjects[seq]:
                    total += subject.count
                total_count_and_seq.append([-total, seq])
            
            # For each row printed finally, now that all the sorting is done
            for arr in sorted(total_count_and_seq):
                seq = arr[1]
                
                taken_sample_names = set(sample_name_order)
                sample_count_and_subject = []
                sample_name_to_subject = {}
                max_count = 0
                max_count_subject = None
                for subject in seq_to_subjects[seq]:
                    sample_name_to_subject[subject.sample_name] = subject
                    if subject.count > max_count:
                        max_count = subject.count
                        max_count_subject = subject
                    if subject.sample_name in taken_sample_names: continue
                    sample_count_and_subject.append([-subject.count, subject])
                sample_count_and_subject.sort()
                
                # Cannot print directly here since order of the columns not 
                # entirely known (the ones that lack this sequence or any less divergent/common).
                row = [seq, divergence, max_count_subject.taxonomy]
                
                # add columns from previous column runs
                for sample_name in sample_name_order:
                    try:
                        count = str(sample_name_to_subject[sample_name].count)
                    except KeyError:
                        count = '0'
                    row.append(count)
                # add columns that are new this row
                for arr in sample_count_and_subject:
                    count, subject = arr
                    row.append(str(-count))
                    sample_name_order.append(subject.sample_name)
                to_print.append(row)
                
        # write headers
        output_io.write("\t".join(itertools.chain(['sequence','divergence'],
                                                  sample_name_order,
                                                  ['taxonomy'])))
        output_io.write("\n")
        for row in to_print:
            output_io.write("\t".join(itertools.chain([row[0], str(row[1])],
                                                      row[3:],
                                                      ['0']*(len(sample_name_order)-len(row[3:])),
                                                      [row[2]]))+"\n")
        