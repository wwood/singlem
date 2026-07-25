import csv

with open(snakemake.input.fscore_list, 'r') as f:
    with open(snakemake.output.resolved_fscores, 'w') as out:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(out, delimiter='\t')
        headers = next(reader)
        writer.writerow(headers)  # write header
        for row in reader:
            gene, first_fscore, sum_best_three, count_iterations = row
            first_fscore = float(first_fscore)
            sum_best_three = float(sum_best_three)
            if first_fscore >= 0.8:
                writer.writerow(row)
            elif first_fscore >= 0.6 and sum_best_three >= 1.8:
                writer.writerow(row)

with open(snakemake.output.done, 'w') as __:
    pass