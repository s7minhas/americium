"""Compile the two manuscript entry points and reject unresolved build problems."""

import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paper = root / 'paper'
assert {p.name for p in paper.glob('*.tex')} == {'main.tex', 'si.tex'}
checks = root / 'results' / 'checks'
checks.mkdir(parents=True, exist_ok=True)
receipt = {'built_utc': datetime.now(timezone.utc).isoformat(), 'documents': []}
for name in ['main', 'si']:
    source = paper / f'{name}.tex'
    text = source.read_text()
    assert not re.search(r'\\(?:input|include)\{', text), 'Keep prose and tables in the two entry points.'
    assert '\u2014' not in text and '---' not in text, 'Em dash found.'
    assert not re.search(r'\(\s*\$?p\s*[=<>]', text), 'Parenthetical p-value found.'
    with tempfile.TemporaryDirectory(prefix='americium_tex_') as temp:
        build = Path(temp)
        environment = os.environ.copy()
        environment['BIBINPUTS'] = str(paper) + os.pathsep + environment.get('BIBINPUTS', '')
        command = ['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error',
                   '-outdir=' + str(build), source.name]
        result = subprocess.run(command, cwd=paper, env=environment,
                                text=True, encoding='utf-8', errors='replace',
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(result.stdout[-12000:])
        log = (build / f'{name}.log').read_text(errors='replace')
        problems = re.findall(r'^.*(?:LaTeX Warning|Package .* Warning|pdfTeX warning|Overfull|undefined).*$',
                              log, flags=re.MULTILINE)
        if problems:
            raise RuntimeError('\n'.join(problems))
        pdf = build / f'{name}.pdf'
        info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
        subprocess.run(['pdftotext', '-layout', str(pdf), str(build / f'{name}.txt')], check=True)
        assert len((build / f'{name}.txt').read_text().strip()) > 1000
        shutil.copy2(pdf, paper / pdf.name)
        receipt['documents'].append({
            'name': name, 'pages': int(re.search(r'Pages:\s+(\d+)', info)[1]),
            'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
            'warnings': problems})
receipt['bibliography_sha256'] = hashlib.sha256((paper / 'references.bib').read_bytes()).hexdigest()
(checks / 'build.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
