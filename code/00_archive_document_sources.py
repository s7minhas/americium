"""Archive and verify documentary sources cited by the manuscripts."""
import argparse
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/revision_sources_20261101'
SOURCES={
    'kegley_hook_metadata.json':'https://api.crossref.org/works/10.2307/2600701',
    'jerusalem_2017_un.html':'https://www.un.org/unispal/document/ga-overwhelmingly-adopts-resolution-asking-nations-not-to-locate-diplomatic-missions-in-jerusalem-ga11995/',
    'lake_metadata.json':'https://api.crossref.org/works/10.1017/S1537592717003085',
    'eu_metadata.json':'https://api.crossref.org/works/10.1080/13501763.2014.919833',
    'crs_tariff_timeline.pdf':'https://www.congress.gov/crs_external_products/R/PDF/R48549/R48549.2.pdf',
    'canada_countertariffs.html':'https://gazette.gc.ca/rp-pr/p2/2025/2025-04-23/html/sor-dors119-eng.html',
    'lake_2018.pdf':'https://apsanet.org/Portals/54/goverance/international_legitimacy_lost_rule_and_resistance_when_america_is_first.pdf',
    'canada_tariff_order.pdf':'https://www.govinfo.gov/content/pkg/FR-2025-02-07/pdf/2025-02406.pdf',
    'mexico_tariff_order.pdf':'https://www.govinfo.gov/content/pkg/FR-2025-02-07/pdf/2025-02407.pdf',
    'aluminum_proclamation.pdf':'https://www.govinfo.gov/content/pkg/FR-2025-02-18/pdf/2025-02832.pdf',
    'steel_proclamation.pdf':'https://www.govinfo.gov/content/pkg/FR-2025-02-18/pdf/2025-02833.pdf',
    'automobile_proclamation.pdf':'https://www.govinfo.gov/content/pkg/FR-2025-04-03/pdf/2025-05930.pdf',
    'reciprocal_tariff_order.pdf':'https://www.govinfo.gov/content/pkg/FR-2025-04-07/pdf/2025-06063.pdf',
}
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('mode',choices=['fetch','verify'],default='verify',nargs='?')
args=parser.parse_args()
OUT.mkdir(parents=True,exist_ok=True)
manifest=OUT/'retrieval.json'
records=json.loads(manifest.read_text()) if manifest.exists() else []
if args.mode=='fetch':
    for filename,url in SOURCES.items():
        path=OUT/filename
        if path.exists():continue
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Research source preservation'}),timeout=60) as response:
            data=response.read();final_url=response.url
        if filename.endswith('.pdf') and not data.startswith(b'%PDF'):raise ValueError(f'Expected a PDF: {url}')
        path.write_bytes(data)
        records.append(dict(path=str(path.relative_to(ROOT)),url=url,final_url=final_url,
                            retrieved_utc=datetime.now(timezone.utc).isoformat(),sha256=hashlib.sha256(data).hexdigest()))
        manifest.write_text(json.dumps(records,indent=2)+'\n')
checks=[dict(path=r['path'],passed=hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256']) for r in records if 'path' in r]
(ROOT/'results/revision/source_checks.json').write_text(json.dumps(dict(checks=checks,passed=all(x['passed'] for x in checks)),indent=2)+'\n')
assert all(x['passed'] for x in checks)
