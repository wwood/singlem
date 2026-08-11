#!/usr/bin/env python3
"""Add static legacy redirects to an already-built Zensical site."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"
REDIRECTS = {
    "Installation": "../",
    "Glossary": "../explanation/otu-tables/",
    "FAQ": "../how-to/process-multiple-samples/",
    "Lyrebird": "../how-to/profile-phages/",
    "tools/data": "../../reference/singlem-data/",
    "tools/pipe": "../../reference/singlem-pipe/",
    "tools/appraise": "../../reference/singlem-appraise/",
    "tools/summarise": "../../reference/singlem-summarise/",
    "tools/renew": "../../reference/singlem-renew/",
    "tools/supplement": "../../reference/singlem-supplement/",
    "tools/prokaryotic_fraction": "../../reference/singlem-prokaryotic_fraction/",
    "tools/lyrebird_data": "../../reference/lyrebird-data/",
    "tools/lyrebird_pipe": "../../reference/lyrebird-pipe/",
    "advanced/makedb": "../../reference/singlem-makedb/",
    "advanced/query": "../../reference/singlem-query/",
    "advanced/condense": "../../reference/singlem-condense/",
    "advanced/seqs": "../../reference/singlem-seqs/",
    "advanced/create": "../../reference/singlem-create/",
    "advanced/metapackage": "../../reference/singlem-metapackage/",
    "advanced/regenerate": "../../reference/singlem-regenerate/",
    "advanced/lyrebird_condense": "../../reference/lyrebird-condense/",
    "advanced/lyrebird_renew": "../../reference/lyrebird-renew/",
}

TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta http-equiv="refresh" content="0; url={target}">
<link rel="canonical" href="{target}"><title>Redirecting</title></head>
<body><p><a href="{target}">Continue to the new page</a>.</p></body></html>
"""


def main() -> None:
    if not SITE.is_dir():
        raise SystemExit("site/ does not exist; run 'zensical build --clean' first")
    for old_path, target in REDIRECTS.items():
        destination = SITE / old_path / "index.html"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(TEMPLATE.format(target=target))


if __name__ == "__main__":
    main()
