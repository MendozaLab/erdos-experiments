#!/usr/bin/env python3
"""Recompile this package using an installed Lean and existing pinned dependencies.

Standard library only. Does not run Lake, install software, fetch dependencies,
or modify dependency caches. The output directory must not already exist.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lean', required=True, type=Path,
                        help='Installed Lean 4.33.1 executable')
    parser.add_argument('--packages', required=True, type=Path,
                        help='Existing directory containing the pinned dependency checkouts/caches')
    parser.add_argument('--output', required=True, type=Path,
                        help='New directory for compiler objects, logs and receipts')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    lean = args.lean.resolve(strict=True)
    packages = args.packages.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        parser.error('Output directory already exists; choose a fresh path.')
    output.mkdir(parents=True)
    lib = output / 'lib'
    logs = output / 'logs'
    lib.mkdir()
    logs.mkdir()
    manifest = json.loads((root / 'lake-manifest.json').read_text())
    provenance = json.loads((root / 'provenance/source-manifest.json').read_text())
    env = {'PATH': '/usr/bin:/bin', 'GIT_CONFIG_GLOBAL': os.devnull}
    version = subprocess.run([str(lean), '--version'], env=env, text=True,
                             capture_output=True, check=True).stdout.strip()
    if not re.search(r'version 4\.33\.1(?:,|\))', version):
        raise RuntimeError('Unexpected compiler: ' + version)
    paths = [str(lib)]
    dependencies = []
    for package in manifest['packages']:
        checkout = packages / package['name']
        revision = subprocess.run(['git', '-C', str(checkout), 'rev-parse', 'HEAD'],
                                  env=env, capture_output=True, text=True, check=True).stdout.strip()
        origin = subprocess.run(['git', '-C', str(checkout), 'remote', 'get-url', 'origin'],
                                env=env, capture_output=True, text=True, check=True).stdout.strip()
        if revision != package['rev']:
            raise RuntimeError('Dependency revision mismatch: ' + package['name'])
        if origin.removesuffix('.git') != package['url'].removesuffix('.git'):
            raise RuntimeError('Dependency origin mismatch: ' + package['name'])
        built = checkout / '.lake/build/lib/lean'
        if not built.is_dir():
            raise RuntimeError('Missing compiled cache: ' + package['name'])
        paths.append(str(built))
        dependencies.append({'name': package['name'], 'revision': revision, 'origin': origin})
    env['LEAN_PATH'] = ':'.join(paths)

    sources = {'.'.join(p.relative_to(root).with_suffix('').parts): p
               for folder in ('OAI', 'Erdos604Proof')
               for p in (root / folder).rglob('*.lean')}
    for name in ('OAI', 'Erdos604Proof'):
        sources[name] = root / (name + '.lean')
    source_hashes = {str(p.relative_to(root)): digest(p) for p in sources.values()}
    for p in sources.values():
        if re.search(r'\b(sorry|admit|axiom|native_decide|unsafe)\b', p.read_text()):
            raise RuntimeError('Disallowed proof token in ' + str(p.relative_to(root)))
    imports = {}
    for name, path in sources.items():
        imported = re.findall(r'^import\s+(\S+)\s*$', path.read_text(), re.M)
        missing = [x for x in imported if x.startswith(('OAI.', 'Erdos604Proof.')) and x not in sources]
        if missing:
            raise RuntimeError('Missing local import: ' + repr(missing))
        imports[name] = [x for x in imported if x in sources]
    ordered = []
    visiting = set()
    visited = set()

    def visit(name):
        if name in visited:
            return
        if name in visiting:
            raise RuntimeError('Import cycle at ' + name)
        visiting.add(name)
        for dep in imports[name]:
            visit(dep)
        visiting.remove(name)
        visited.add(name)
        ordered.append(name)

    for name in sorted(sources):
        visit(name)
    receipt = {'started_at_utc': datetime.now(timezone.utc).isoformat(),
               'compiler': version, 'dependencies': dependencies,
               'source_sha256': source_hashes,
               'lake_manifest_sha256': digest(root / 'lake-manifest.json'),
               'replay_script_sha256': digest(Path(__file__)),
               'modules': [], 'scope': 'Fresh direct Lean compilation against existing pinned caches'}
    for index, name in enumerate(ordered, 1):
        relative = sources[name].relative_to(root)
        obj = lib / relative.with_suffix('.olean')
        obj.parent.mkdir(parents=True, exist_ok=True)
        command = [str(lean), '-DautoImplicit=false', '-DrelaxedAutoImplicit=false',
                   '-DwarningAsError=true', '-R', '.', '-o', str(obj), str(relative)]
        log = logs / (name + '.log')
        started = time.monotonic()
        with log.open('w') as stream:
            result = subprocess.run(command, cwd=root, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT)
        receipt['modules'].append({'module': name, 'command': command,
                                   'exit_code': result.returncode,
                                   'seconds': round(time.monotonic() - started, 3),
                                   'log': str(log.relative_to(output)), 'log_sha256': digest(log)})
        (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        print(f'{index}/{len(ordered)} {name}: exit {result.returncode}', flush=True)
        if result.returncode:
            print(log.read_text(), file=sys.stderr)
            return result.returncode
    for relative, expected in source_hashes.items():
        if digest(root / relative) != expected:
            raise RuntimeError('Source changed during compilation: ' + relative)
    audit_text = (logs / 'Erdos604Proof.Audit.log').read_text()
    reports = re.findall(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", audit_text, re.S)
    allowed = {'propext', 'Classical.choice', 'Quot.sound'}
    if {name for name, _ in reports} != set(provenance['audit_declarations']):
        raise RuntimeError('Missing or unexpected axiom report')
    audited = {}
    for name, body in reports:
        axioms = {x.strip() for x in body.split(',') if x.strip()}
        if not axioms <= allowed:
            raise RuntimeError('Nonstandard axiom in ' + name)
        audited[name] = sorted(axioms)
    receipt.update({'finished_at_utc': datetime.now(timezone.utc).isoformat(),
                    'status': 'pass', 'module_count': len(ordered),
                    'source_hashes_unchanged': True, 'axiom_reports': audited})
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(f'PASS: {len(ordered)} modules and {len(audited)} axiom reports.', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
