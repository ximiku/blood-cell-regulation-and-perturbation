"""Snapshot stage-4 remote science files and verify their local byte/SHA copies."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json

ROOT = Path(__file__).resolve().parents[1]
REMOTE = ROOT / 'logs/stage4_remote_manifest.json'
LOCAL = ROOT / 'logs/stage4_artifact_manifest.json'
VALIDATION = ROOT / 'logs/stage4_remote_archive_validation.json'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(2**20), b''):
            digest.update(block)
    return digest.hexdigest()


def remote_files():
    patterns = [
        'results/stage4/**/*', 'results/stage4_exploratory/**/*',
        'tables/04_*', 'logs/04_*', 'logs/stage4_*',
        'data/splits/stage4_*', 'data/processed/stage4_*',
        'scripts/stage4_*.py', 'scripts/run_stage4_gpu.py',
        'scripts/check_stage4_solvers.py', 'scripts/verify_stage4*.py',
        'scripts/normalize_stage4_array_strings.py', 'config/stage4_config.json',
        'data/processed/analysis_matrices.h5', 'data/splits/gene_cv_splits.json',
    ]
    return sorted({p for pattern in patterns for p in ROOT.glob(pattern)})


def snapshot(remote):
    output = REMOTE if remote else LOCAL
    paths = remote_files() if remote else sorted(ROOT.rglob('*'))
    rows = []
    for path in paths:
        relative = path.relative_to(ROOT)
        if (not path.is_file() or path == output or path.suffix == '.tmp'
                or any(part.startswith('.') or part == '__pycache__' for part in relative.parts)
                or (remote and path in (LOCAL, VALIDATION))):
            continue
        rows.append({'path': str(relative), 'bytes': path.stat().st_size, 'sha256': sha(path)})
    result = {'generated_utc': datetime.now(timezone.utc).isoformat(), 'root': str(ROOT),
              'scope': ('stage4 remote scientific inputs, outputs, source, and logs; '
                        'earlier stages retained in their own snapshots') if remote else
                       'all visible solution artifacts, excluding this snapshot, temporary files and Python caches',
              'files': rows}
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'manifest': str(output), 'files': len(rows), 'bytes': sum(r['bytes'] for r in rows)}))


def verify_remote():
    manifest = json.loads(REMOTE.read_text())
    failures = []
    for row in manifest['files']:
        path = ROOT / row['path']
        if not path.is_file():
            failures.append({'path': row['path'], 'reason': 'missing'})
        elif path.stat().st_size != row['bytes'] or sha(path) != row['sha256']:
            failures.append({'path': row['path'], 'reason': 'bytes or SHA-256 mismatch'})
    result = {'generated_utc': datetime.now(timezone.utc).isoformat(),
              'passed': not failures, 'files_checked': len(manifest['files']),
              'bytes_checked': sum(r['bytes'] for r in manifest['files']),
              'failures': failures}
    VALIDATION.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--remote', action='store_true')
    parser.add_argument('--verify-remote', action='store_true')
    args = parser.parse_args()
    verify_remote() if args.verify_remote else snapshot(args.remote)
