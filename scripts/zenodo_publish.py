#!/usr/bin/env python3
"""
Zenodo new-version automation for the EHP preprint.
Record concept DOI: 10.5281/zenodo.19184468

Usage:
  # Dry run (creates draft, does NOT publish):
  python zenodo_publish.py --token $ZENODO_TOKEN --files results/erdos-114/ papers/ehp-conjecture-114/EHP_Erdos114_Preprint.pdf

  # Publish immediately:
  python zenodo_publish.py --token $ZENODO_TOKEN --files ... --publish

  # Use ZENODO_TOKEN env var instead of --token:
  ZENODO_TOKEN=xxx python zenodo_publish.py --files ...

API: Zenodo InvenioRDM REST v1 (https://developers.zenodo.org/)
"""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Install requests: pip install requests")

RECORD_ID = "19184468"
BASE = "https://zenodo.org/api"
TIMEOUT = 60

# Files to include in the new version (relative to repo root or absolute).
# Only files matching these patterns are uploaded.
INCLUDE_EXTENSIONS = {".pdf", ".json", ".sha256", ".md", ".rs", ".toml", ".tex"}
EXCLUDE_NAMES = {"Cargo.lock"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def api(method: str, path: str, token: str, **kwargs) -> requests.Response:
    headers = {"Authorization": f"Bearer {token}", "Accept": "application/json"}
    headers.update(kwargs.pop("headers", {}))
    url = f"{BASE}{path}"
    resp = getattr(requests, method)(url, headers=headers, timeout=TIMEOUT, **kwargs)
    if not resp.ok:
        print(f"  ERROR {resp.status_code}: {resp.text[:400]}", file=sys.stderr)
        resp.raise_for_status()
    return resp


def collect_files(paths: list[str]) -> list[Path]:
    """Expand directories and filter to uploadable files."""
    result = []
    for p in paths:
        fp = Path(p)
        if fp.is_dir():
            for child in sorted(fp.rglob("*")):
                if child.is_file() and child.suffix in INCLUDE_EXTENSIONS and child.name not in EXCLUDE_NAMES:
                    result.append(child)
        elif fp.is_file():
            result.append(fp)
        else:
            print(f"  WARN: {p} not found, skipping", file=sys.stderr)
    return result


def create_new_version(record_id: str, token: str) -> dict:
    print(f"Creating new version of record {record_id} …")
    r = api("post", f"/records/{record_id}/versions", token)
    draft = r.json()
    draft_id = draft["id"]
    print(f"  Draft ID: {draft_id}")
    return draft


def delete_existing_files(draft_id: str, token: str) -> None:
    print("Removing existing draft files …")
    r = api("get", f"/records/{draft_id}/draft/files", token)
    existing = r.json().get("entries", [])
    for entry in existing:
        key = entry["key"]
        api("delete", f"/records/{draft_id}/draft/files/{key}", token)
        print(f"  Deleted: {key}")


def upload_file(draft_id: str, token: str, local_path: Path, remote_name: str) -> None:
    size = local_path.stat().st_size
    print(f"  Uploading {remote_name} ({size:,} bytes) …", end=" ", flush=True)

    # 1. Reserve the slot
    api("post", f"/records/{draft_id}/draft/files", token,
        headers={"Content-Type": "application/json"},
        data=json.dumps([{"key": remote_name}]))

    # 2. Upload content (octet-stream)
    with open(local_path, "rb") as fh:
        api("put", f"/records/{draft_id}/draft/files/{remote_name}/content", token,
            headers={"Content-Type": "application/octet-stream"},
            data=fh)

    # 3. Commit
    r = api("post", f"/records/{draft_id}/draft/files/{remote_name}/commit", token)
    checksum = r.json().get("checksum", "?")
    print(f"ok  checksum={checksum}")


def update_metadata(draft_id: str, token: str, degree_max: int) -> None:
    print("Updating metadata …")
    # Fetch current metadata first
    r = api("get", f"/records/{draft_id}/draft", token)
    meta = r.json()["metadata"]

    meta["title"] = (
        f"Computational Verification of the Erdős-Herzog-Piranian Conjecture "
        f"for Degrees 3 ≤ n ≤ {degree_max}"
    )
    meta["description"] = (
        f"IEEE 1788 interval arithmetic certified branch-and-bound verification "
        f"of the EHP conjecture for monic polynomials of degrees 3 through {degree_max}. "
        f"The margin by which z^n−1 dominates grows from ~6% (n=3) to over 61% (n=11). "
        f"Includes closed-form formula L(z^n−1) = 2^{{1/n}}·√π·Γ(1/(2n))/Γ(1/(2n)+1/2). "
        f"Code + results: https://github.com/MendozaLab/erdos-experiments"
    )
    meta["version"] = f"n3-n{degree_max}"

    api("put", f"/records/{draft_id}/draft", token,
        headers={"Content-Type": "application/json"},
        data=json.dumps({"metadata": meta}))
    print(f"  Title updated to n={degree_max}")


def publish(draft_id: str, token: str) -> str:
    print("Publishing …")
    r = api("post", f"/records/{draft_id}/draft/actions/publish", token)
    rec = r.json()
    doi = rec.get("doi", rec.get("conceptdoi", "?"))
    links = rec.get("links", {})
    url = links.get("self_html", links.get("html", "?"))
    print(f"  Published! DOI={doi}  URL={url}")
    return url


def main() -> None:
    parser = argparse.ArgumentParser(description="Zenodo new-version publisher")
    parser.add_argument("--token", default=os.environ.get("ZENODO_TOKEN"),
                        help="Zenodo personal access token (or set ZENODO_TOKEN env var)")
    parser.add_argument("--record", default=RECORD_ID,
                        help=f"Zenodo record ID (default: {RECORD_ID})")
    parser.add_argument("--files", nargs="+", required=True,
                        help="Files or directories to upload")
    parser.add_argument("--degree-max", type=int, default=11,
                        help="Maximum degree proven (for metadata title)")
    parser.add_argument("--publish", action="store_true",
                        help="Publish immediately (default: leave as draft)")
    args = parser.parse_args()

    if not args.token:
        sys.exit("ERROR: provide --token or set ZENODO_TOKEN env var")

    files = collect_files(args.files)
    if not files:
        sys.exit("ERROR: no uploadable files found")

    print(f"\nFiles to upload ({len(files)}):")
    for f in files:
        print(f"  {f}  [{sha256_file(f)[:16]}…]")

    # Step 1: create draft
    draft = create_new_version(args.record, args.token)
    draft_id = draft["id"]

    # Step 2: clear old files
    delete_existing_files(draft_id, args.token)

    # Step 3: upload new files
    print(f"\nUploading {len(files)} files …")
    for local_path in files:
        upload_file(draft_id, args.token, local_path, local_path.name)

    # Step 4: update metadata
    update_metadata(draft_id, args.token, args.degree_max)

    # Step 5: publish or leave as draft
    if args.publish:
        url = publish(draft_id, args.token)
        print(f"\nDone. New version live at: {url}")
    else:
        edit_url = f"https://zenodo.org/uploads/{draft_id}"
        print(f"\nDraft ready for review: {edit_url}")
        print("Re-run with --publish to make it live.")


if __name__ == "__main__":
    main()
