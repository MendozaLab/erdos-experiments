#!/usr/bin/env python3
"""
scripts/zenodo/upload.py
Zenodo automation for MendozaLab preprints and experiment data.

Usage:
  # Upload a new restricted preprint, get DOI back
  python upload.py preprint \
      --file paper.pdf \
      --title "Lemniscate Duality and the Erdos #114 Channel Hypothesis" \
      --authors "Mendoza, Kenneth A." \
      --description "Preprint describing certified interval-arithmetic proofs..." \
      --keywords "Erdos problems,lemniscate,information theory,channel capacity" \
      --access restricted

  # Upload experiment data as a new version, get DOI back
  python upload.py data \
      --dir ../../results/erdos-114 \
      --title "EHP Channel Experiments n=3-11 certified results" \
      --related-doi 10.5281/zenodo.XXXXXXX

  # Publish (flip restricted -> open) an existing record
  python upload.py publish --record-id 1234567

  # List your records
  python upload.py list

Environment:
  ZENODO_TOKEN   Zenodo personal access token (required)
  ZENODO_SANDBOX Set to '1' to use sandbox.zenodo.org instead of zenodo.org

Dependencies:
  pip install requests
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

SANDBOX = os.environ.get("ZENODO_SANDBOX", "0") == "1"
BASE_URL = "https://sandbox.zenodo.org/api" if SANDBOX else "https://zenodo.org/api"


def get_token():
    token = os.environ.get("ZENODO_TOKEN", "")
    if not token:
        print("ERROR: ZENODO_TOKEN environment variable not set.", file=sys.stderr)
        print("  Get a token at: https://zenodo.org/account/settings/applications/", file=sys.stderr)
        sys.exit(1)
    return token


def headers(token):
    return {"Authorization": f"Bearer {token}"}


# ---------------------------------------------------------------------------
# Core API helpers
# ---------------------------------------------------------------------------

def create_deposition(token, metadata):
    """Create a new empty deposition and return its id + bucket URL."""
    r = requests.post(
        f"{BASE_URL}/deposit/depositions",
        json={"metadata": metadata},
        headers={**headers(token), "Content-Type": "application/json"},
    )
    r.raise_for_status()
    data = r.json()
    return data["id"], data["links"]["bucket"], data["links"]["html"]


def upload_file(token, bucket_url, filepath):
    """Upload a single file to the deposition bucket."""
    path = Path(filepath)
    with open(path, "rb") as f:
        r = requests.put(
            f"{bucket_url}/{path.name}",
            data=f,
            headers=headers(token),
        )
    r.raise_for_status()
    return r.json()


def upload_directory(token, bucket_url, dirpath):
    """Upload all files in a directory (non-recursive)."""
    uploaded = []
    for p in sorted(Path(dirpath).iterdir()):
        if p.is_file():
            print(f"  Uploading {p.name} ({p.stat().st_size} bytes)...")
            result = upload_file(token, bucket_url, p)
            uploaded.append({"file": p.name, "checksum": result.get("checksum", "")})
    return uploaded


def publish_deposition(token, deposition_id):
    """Publish (or re-publish) a deposition."""
    r = requests.post(
        f"{BASE_URL}/deposit/depositions/{deposition_id}/actions/publish",
        headers=headers(token),
    )
    r.raise_for_status()
    data = r.json()
    return data["doi"], data["links"]["html"]


def get_deposition(token, deposition_id):
    r = requests.get(
        f"{BASE_URL}/deposit/depositions/{deposition_id}",
        headers=headers(token),
    )
    r.raise_for_status()
    return r.json()


def list_depositions(token):
    r = requests.get(
        f"{BASE_URL}/deposit/depositions",
        headers=headers(token),
        params={"sort": "mostrecent", "size": 25},
    )
    r.raise_for_status()
    return r.json()


# ---------------------------------------------------------------------------
# Metadata builders
# ---------------------------------------------------------------------------

def build_preprint_metadata(args):
    authors = []
    for name in args.authors.split(";"):
        name = name.strip()
        if "," in name:
            family, given = [p.strip() for p in name.split(",", 1)]
        else:
            family, given = name, ""
        authors.append({"name": f"{family}, {given}", "affiliation": "MendozaLab"})

    keywords = [k.strip() for k in args.keywords.split(",")]

    meta = {
        "title": args.title,
        "upload_type": "publication",
        "publication_type": "preprint",
        "description": args.description,
        "creators": authors,
        "keywords": keywords,
        "access_right": args.access,
        "license": "cc-by-4.0",
    }

    if args.access == "restricted":
        meta["access_conditions"] = (
            "Contact author for access prior to public release. "
            "Will be made open after peer review or arXiv submission."
        )

    if hasattr(args, "related_doi") and args.related_doi:
        meta["related_identifiers"] = [
            {
                "identifier": args.related_doi,
                "relation": "isSupplementTo",
                "scheme": "doi",
            }
        ]

    return meta


def build_data_metadata(args):
    meta = {
        "title": args.title,
        "upload_type": "dataset",
        "description": (
            f"Certified interval-arithmetic experiment results from "
            f"MendozaLab/erdos-experiments. {args.title}"
        ),
        "creators": [{"name": "Mendoza, Kenneth A.", "affiliation": "MendozaLab"}],
        "keywords": [
            "Erdos problems",
            "lemniscate",
            "interval arithmetic",
            "certified proof",
            "channel capacity",
        ],
        "access_right": "open",
        "license": "cc-by-4.0",
    }

    if hasattr(args, "related_doi") and args.related_doi:
        meta["related_identifiers"] = [
            {
                "identifier": args.related_doi,
                "relation": "isSupplementTo",
                "scheme": "doi",
            }
        ]

    return meta


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_preprint(args):
    token = get_token()
    print(f"Creating {'SANDBOX' if SANDBOX else 'LIVE'} Zenodo preprint record...")
    print(f"  Title: {args.title}")
    print(f"  Access: {args.access}")

    meta = build_preprint_metadata(args)
    dep_id, bucket_url, html_url = create_deposition(token, meta)
    print(f"  Deposition ID: {dep_id}")
    print(f"  Draft URL: {html_url}")

    print(f"  Uploading {args.file}...")
    upload_file(token, bucket_url, args.file)

    # DOI is reserved but record stays as draft (not published yet)
    dep = get_deposition(token, dep_id)
    prereserved_doi = dep.get("metadata", {}).get("prereserve_doi", {}).get("doi", "pending")
    print(f"\n  Reserved DOI: {prereserved_doi}")
    print(f"  Status: DRAFT (restricted, not yet published)")
    print(f"  To publish: python upload.py publish --record-id {dep_id}")

    result = {
        "deposition_id": dep_id,
        "prereserved_doi": prereserved_doi,
        "html_url": html_url,
        "access": args.access,
        "status": "draft",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    out = Path("zenodo_record.json")
    out.write_text(json.dumps(result, indent=2))
    print(f"  Record saved to {out}")
    return result


def cmd_data(args):
    token = get_token()
    print(f"Creating {'SANDBOX' if SANDBOX else 'LIVE'} Zenodo dataset record...")
    print(f"  Title: {args.title}")
    print(f"  Directory: {args.dir}")

    meta = build_data_metadata(args)
    dep_id, bucket_url, html_url = create_deposition(token, meta)
    print(f"  Deposition ID: {dep_id}")

    uploaded = upload_directory(token, bucket_url, args.dir)
    print(f"  Uploaded {len(uploaded)} files.")

    dep = get_deposition(token, dep_id)
    prereserved_doi = dep.get("metadata", {}).get("prereserve_doi", {}).get("doi", "pending")
    print(f"\n  Reserved DOI: {prereserved_doi}")
    print(f"  To publish: python upload.py publish --record-id {dep_id}")

    result = {
        "deposition_id": dep_id,
        "prereserved_doi": prereserved_doi,
        "html_url": html_url,
        "files_uploaded": uploaded,
        "status": "draft",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    out = Path("zenodo_data_record.json")
    out.write_text(json.dumps(result, indent=2))
    print(f"  Record saved to {out}")
    return result


def cmd_publish(args):
    token = get_token()
    print(f"Publishing deposition {args.record_id}...")
    doi, html_url = publish_deposition(token, args.record_id)
    print(f"  Published DOI: {doi}")
    print(f"  URL: {html_url}")
    return {"doi": doi, "url": html_url}


def cmd_list(args):
    token = get_token()
    records = list_depositions(token)
    print(f"{'ID':<12} {'Status':<12} {'Title'}")
    print("-" * 70)
    for r in records:
        print(f"{r['id']:<12} {r['state']:<12} {r['title'][:44]}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Zenodo upload automation for MendozaLab preprints and data"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # preprint
    p_pre = sub.add_parser("preprint", help="Upload a paper PDF as a restricted preprint")
    p_pre.add_argument("--file", required=True, help="Path to PDF")
    p_pre.add_argument("--title", required=True)
    p_pre.add_argument("--authors", required=True,
                       help="Semicolon-separated 'Family, Given' names")
    p_pre.add_argument("--description", required=True)
    p_pre.add_argument("--keywords", default="Erdos problems,mathematics")
    p_pre.add_argument("--access", choices=["open", "restricted", "embargoed"],
                       default="restricted")
    p_pre.add_argument("--related-doi", default="", dest="related_doi")

    # data
    p_data = sub.add_parser("data", help="Upload a directory of result files as a dataset")
    p_data.add_argument("--dir", required=True, help="Directory to upload")
    p_data.add_argument("--title", required=True)
    p_data.add_argument("--related-doi", default="", dest="related_doi")

    # publish
    p_pub = sub.add_parser("publish", help="Publish a draft deposition (makes it live)")
    p_pub.add_argument("--record-id", required=True, dest="record_id")

    # list
    sub.add_parser("list", help="List recent depositions")

    args = parser.parse_args()

    dispatch = {
        "preprint": cmd_preprint,
        "data": cmd_data,
        "publish": cmd_publish,
        "list": cmd_list,
    }
    dispatch[args.command](args)


if __name__ == "__main__":
    main()
