"""Check manuscript PDFs, references, word counts and source consistency."""
import hashlib
import json
import re
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/revision'
receipt=json.loads((ROOT/'results/checks/build.json').read_text())
bibliography=(ROOT/'paper/references.bib').read_text()
keys=re.findall(r'@\w+\s*\{([^,]+),',bibliography)
assert len(keys)==len(set(keys)), 'Duplicate bibliography key'
metrics={}
for document in receipt['documents']:
    name=document['name']
    path=ROOT/'paper'/f'{name}.tex'
    text=path.read_text()
    assert document['source_sha256']==hashlib.sha256(path.read_bytes()).hexdigest(), 'Stale PDF source'
    assert document['pdf_sha256']==hashlib.sha256((ROOT/'paper'/f'{name}.pdf').read_bytes()).hexdigest(), 'PDF hash mismatch'
    assert not document['warnings']
    labels=re.findall(r'\\label\{([^}]+)\}',text)
    references=re.findall(r'\\ref\{([^}]+)\}',text)
    assert len(labels)==len(set(labels)),f'Duplicate label in {name}'
    assert set(references)<=set(labels),f'Missing reference in {name}'
    citations=[]
    for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\]){0,2}\{([^}]+)\}',text):
        citations.extend(x.strip() for x in group.split(','))
    assert set(citations)<=set(keys),f'Missing citation in {name}'
    result=subprocess.check_output(['texcount','-sum','-1',str(path)],cwd=ROOT,text=True).strip()
    words=int(result.splitlines()[-1])
    metrics[name]=dict(pages=document['pages'],texcount=words,source_sha256=document['source_sha256'])
    if name=='main':
        abstract=text.split('\\begin{abstract}',1)[1].split('\\vspace',1)[0]
        abstract=re.sub(r'\\(?:singlespacing|noindent)','',abstract)
        metrics[name]['abstract_words']=len(abstract.split())
        assert metrics[name]['abstract_words']<=150
        assert words<=12000
    else:
        assert document['pages']<=25
assert receipt['bibliography_sha256']==hashlib.sha256((ROOT/'paper/references.bib').read_bytes()).hexdigest()
claims=json.loads((OUT/'claim_check_summary.json').read_text())
assert claims['retained_numeric_checks']==claims['retained_passes']
assert not claims['unmatched_new_bindings']
for path,expected in claims['manuscripts'].items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,'Stale number bindings'
metrics['passed']=True
(OUT/'artifact_checks.json').write_text(json.dumps(metrics,indent=2)+'\n')
