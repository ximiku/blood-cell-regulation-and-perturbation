"""Download small authoritative annotations, retaining source and checksums."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "references"
SOURCES = {
    "HumanTFs_v1.01.csv": "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv",
    "GSE92742_Broad_LINCS_gene_info.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz",
    "GSE92742_Broad_LINCS_pert_info.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz",
    "hgnc_complete_set.txt": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
}


def fetch(item):
    name, url = item
    path = DEST / name
    record = {"file": str(path.relative_to(ROOT)), "url": url}
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "HW2-academic-analysis/1.0"})
        with urllib.request.urlopen(request, timeout=45) as response:
            body = response.read()
            record["final_url"] = response.url
        if len(body) < 100 or body[:100].lstrip().lower().startswith(b"<!doctype html"):
            raise ValueError("Expected an annotation file; response content failed validation")
        path.write_bytes(body)
        record.update(status="downloaded", bytes=len(body), sha256=hashlib.sha256(body).hexdigest(),
                      retrieved_utc=datetime.now(timezone.utc).isoformat())
    except Exception as error:
        record.update(status="error", error=f"{type(error).__name__}: {error}")
    print(json.dumps(record, ensure_ascii=False), flush=True)
    return record


if __name__ == "__main__":
    DEST.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(fetch, SOURCES.items()))
    (DEST / "source_manifest.json").write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n")
    if any(row["status"] != "downloaded" for row in records):
        raise SystemExit(1)
