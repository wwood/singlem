#!/usr/bin/env python3
"""Draft (but do not publish) a new version of the SingleM metapackage Zenodo record.

Steps:
  1. Create a new-version draft from an existing published record.
  2. Remove the files carried over from the previous version.
  3. Upload the new zenodo_backpack archive (streamed, so large files are fine).
  4. Update version, publication date and prepend a changelog entry to the description.

The draft is left unpublished: review it on the Zenodo website and publish manually.

Example:
  python scripts/zenodo_draft_new_version.py \\
      --record-id 20150069 \\
      --file ~/s/S7.5.0.GTDB_r232.metapackage_20260805.smpkg.zb.tar.zst \\
      --version 7.5.0 \\
      --changelog 'Bundled a weebill database for joint SingleM+weebill profiling (singlem pipe/condense --joint).' \\
      --changelog 'The 7.x version number indicates the metapackage format. Old 6.x, 5.x and 4.x versions still work with the current software version.'
"""

import argparse
import datetime
import hashlib
import html
import logging
import os
import re
import sys

import requests
from tqdm import tqdm

ZENODO_API = "https://zenodo.org/api"
DEFAULT_TOKEN_FILE = os.path.expanduser("~/.zenodo_draft_release_api_token")


class ProgressFile:
    """File wrapper that shows a progress bar and computes md5 while streaming."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.size = os.path.getsize(path)
        self.md5 = hashlib.md5()
        self.pbar = tqdm(total=self.size, unit="B", unit_scale=True, unit_divisor=1024,
                         desc="Uploading", miniters=1, mininterval=1)

    def read(self, n=-1):
        chunk = self.f.read(n if n and n > 0 else 8 * 1024 * 1024)
        self.md5.update(chunk)
        self.pbar.update(len(chunk))
        return chunk

    def __len__(self):
        return self.size

    def close(self):
        self.pbar.close()
        self.f.close()


def check(resp, what):
    if not resp.ok:
        logging.error("%s failed: HTTP %s\n%s", what, resp.status_code, resp.text)
        sys.exit(1)
    return resp


def build_description(old_description, version, changelog_items):
    entry = "<p>version {}</p>\n<ul>\n{}\n</ul>\n".format(
        html.escape(version),
        "\n".join("<li>{}</li>".format(html.escape(c)) for c in changelog_items))
    marker = "<h2>Changelog</h2>\n"
    if marker in old_description:
        return old_description.replace(marker, marker + entry, 1)
    logging.warning("No '<h2>Changelog</h2>' heading found; appending changelog at end of description")
    return old_description + "\n" + marker + entry


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--record-id", default="20150069", help="Published record to create a new version of [default: %(default)s]")
    parser.add_argument("--file", required=True, help="zenodo_backpack archive to upload")
    parser.add_argument("--version", required=True, help="New version string, e.g. 7.5.0")
    parser.add_argument("--changelog", action="append", default=[], help="Changelog bullet point (repeatable)")
    parser.add_argument("--gtdb-version", help="If given, replace 'GTDB) version XX-RSNNN' in the description, e.g. '11-RS232'")
    parser.add_argument("--token-file", default=DEFAULT_TOKEN_FILE, help="File containing a Zenodo API token [default: %(default)s]")
    parser.add_argument("--keep-old-files", action="store_true", help="Do not delete files carried over from the previous version")
    parser.add_argument("--dry-run", action="store_true", help="Show the planned metadata changes without touching Zenodo")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    path = os.path.expanduser(args.file)
    if not os.path.isfile(path):
        sys.exit("File not found: {}".format(path))
    filename = os.path.basename(path)
    if not args.changelog:
        sys.exit("At least one --changelog entry is required")

    # Fetch the public metadata of the existing record (no auth needed)
    old = check(requests.get("{}/records/{}".format(ZENODO_API, args.record_id)), "Fetching record").json()
    old_meta = old["metadata"]
    if not old_meta["relations"]["version"][0].get("is_last", False):
        logging.warning("Record %s is not the latest version of its concept record", args.record_id)
    logging.info("Existing record %s is version %s", args.record_id, old_meta.get("version"))
    if old_meta.get("version") == args.version:
        sys.exit("New version is the same as the existing one ({})".format(args.version))

    description = build_description(old_meta["description"], args.version, args.changelog)
    if args.gtdb_version:
        description, n = re.subn(r"(GTDB\) version )\S+?(,)", r"\g<1>{}\2".format(args.gtdb_version), description, count=1)
        if n == 0:
            logging.warning("Could not find GTDB version text to replace in description")
    today = datetime.date.today().isoformat()

    if args.dry_run:
        logging.info("Dry run. Would upload %s (%.1f GiB) as new version %s, publication_date %s",
                     filename, os.path.getsize(path) / 1024**3, args.version, today)
        print(description)
        return

    with open(os.path.expanduser(args.token_file)) as f:
        token = f.read().strip()
    session = requests.Session()
    session.headers["Authorization"] = "Bearer {}".format(token)

    # 1. Create new version draft
    logging.info("Creating new version draft from record %s ..", args.record_id)
    resp = check(session.post("{}/deposit/depositions/{}/actions/newversion".format(ZENODO_API, args.record_id)),
                 "Creating new version")
    draft_url = resp.json()["links"]["latest_draft"]
    draft = check(session.get(draft_url), "Fetching draft").json()
    draft_id = draft["id"]
    bucket_url = draft["links"]["bucket"]
    logging.info("Draft deposition id %s: %s", draft_id, draft["links"].get("html"))

    # 2. Delete files carried over from the previous version
    if not args.keep_old_files:
        files = check(session.get("{}/deposit/depositions/{}/files".format(ZENODO_API, draft_id)), "Listing draft files").json()
        for fl in files:
            if fl["filename"] == filename:
                continue
            logging.info("Deleting carried-over file %s from draft ..", fl["filename"])
            check(session.delete("{}/deposit/depositions/{}/files/{}".format(ZENODO_API, draft_id, fl["id"])),
                  "Deleting {}".format(fl["filename"]))

    # 3. Upload the new file
    logging.info("Uploading %s (%.1f GiB) ..", filename, os.path.getsize(path) / 1024**3)
    pf = ProgressFile(path)
    try:
        resp = check(session.put("{}/{}".format(bucket_url, filename), data=pf,
                                 headers={"Content-Type": "application/octet-stream"}),
                     "Uploading file")
    finally:
        pf.close()
    remote_checksum = resp.json().get("checksum", "")
    local_checksum = "md5:" + pf.md5.hexdigest()
    if remote_checksum != local_checksum:
        logging.error("Checksum mismatch: local %s, Zenodo %s", local_checksum, remote_checksum)
        sys.exit(1)
    logging.info("Upload complete, checksum %s verified", local_checksum)

    # 4. Update metadata
    meta = draft["metadata"]
    meta["version"] = args.version
    meta["publication_date"] = today
    meta["description"] = description
    meta.pop("doi", None)  # let Zenodo assign a new DOI on publish
    check(session.put("{}/deposit/depositions/{}".format(ZENODO_API, draft_id), json={"metadata": meta}),
          "Updating metadata")

    logging.info("Draft ready (NOT published). Review and publish at: https://zenodo.org/uploads/%s", draft_id)


if __name__ == "__main__":
    main()
