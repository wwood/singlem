import csv

headers = ["gene","first_fscore","sum_best_three", "count_iterations"]
with open(snakemake.output.fscore_list, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(headers)
    for fscore in snakemake.input.fscores:
        with open(fscore, 'r') as fs:
            reader = csv.reader(fs, delimiter='\t')
            next(reader)  # skip header
            for row in reader:
                writer.writerow(row)

with open(snakemake.output.done, 'w') as __:
    pass