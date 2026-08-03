import csv
from pathlib import Path


rows = []
for filename in snakemake.input:
    path = Path(filename)
    with path.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    rows.append({"method": path.parent.name, "dataset": path.stem, **row})

fieldnames = ["method", "dataset"] + [
    field for field in rows[0] if field not in {"method", "dataset"}
]
with open(snakemake.output[0], "w") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
