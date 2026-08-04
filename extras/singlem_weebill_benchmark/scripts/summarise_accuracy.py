import csv
from pathlib import Path


WANTED_RANKS = {"genus", "species"}
WANTED_METRICS = {
    "Bray-Curtis distance",
    "F1 score",
    "True positives",
    "False positives",
    "False negatives",
}

output_rows = []
for filename in snakemake.input:
    path = Path(filename)
    by_rank = {rank: {} for rank in WANTED_RANKS}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["tool"] == "Gold standard" or row["rank"] not in WANTED_RANKS:
                continue
            if row["metric"] in WANTED_METRICS:
                by_rank[row["rank"]][row["metric"]] = float(row["value"])

    for rank in sorted(WANTED_RANKS):
        metrics = by_rank[rank]
        missing = WANTED_METRICS - metrics.keys()
        if missing:
            raise ValueError(f"Missing OPAL metrics in {path} at {rank}: {sorted(missing)}")
        tp = metrics["True positives"]
        fp = metrics["False positives"]
        fn = metrics["False negatives"]
        output_rows.append({
            "method": path.parent.name,
            "dataset": path.name.removesuffix(".results.tsv"),
            "rank": rank,
            "bray_curtis": metrics["Bray-Curtis distance"],
            "f1": metrics["F1 score"],
            "false_positive_rate": fp / (tp + fp) if tp + fp else 0.0,
            "false_negative_rate": fn / (tp + fn) if tp + fn else 0.0,
            "true_positives": tp,
            "false_positives": fp,
            "false_negatives": fn,
        })

with open(snakemake.output[0], "w") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(output_rows[0]))
    writer.writeheader()
    writer.writerows(output_rows)
