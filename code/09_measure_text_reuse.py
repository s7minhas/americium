"""Text carried forward between successive adopted General Assembly resolutions.

fetch:   python code/09_measure_text_reuse.py fetch FIRST_SESSION LAST_SESSION WORK_DIR
         downloads adopted resolutions as PDF, extracts text, and records every attempt in WORK_DIR/manifest.csv
letters: python code/09_measure_text_reuse.py letters WORK_DIR
         fetches lettered resolutions through combined symbols such as "A/RES/80/84 A-B"
measure: python code/09_measure_text_reuse.py measure TEXT_ARCHIVE OUT_DIR
         computes the measure from the archived texts alone; nothing is downloaded
terms:   python code/09_measure_text_reuse.py terms TEXT_ARCHIVE OUT_DIR
         counts the contested terms in the same-title pairs adopted in 2017 and 2025 and the operative
         paragraphs reproduced verbatim in the same-title pairs of every year, from OUT_DIR/pairs.csv,
         which the measure mode writes first

A predecessor is the most recent resolution with the same title within four sessions. An operative
paragraph is carried forward when at least 80 percent of its four-word sequences, digits removed,
occur anywhere in the successor. Each pair of resolutions has equal weight. No model judgment enters.
"""
import csv
import hashlib
import io
import json
import re
import subprocess
import sys
import threading
import time
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

project = Path(__file__).resolve().parents[1]
MONTHS = ['january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december']
SHINGLE = 4
MIN_WORDS = 8
MAX_GAP = 4
THRESHOLDS = (0.6, 0.8, 0.9)
TERM_YEARS = (2017, 2025)
TERMS = {
    'gender': r'\bgenders?\b',
    'sex': r'\bsex(?:es)?\b',
    'diversity': r'\bdiversit(?:y|ies)\b',
    'biological diversity': r'\bbiological diversit(?:y|ies)\b',
    'climate change': r'\bclimate change\b',
    'Sustainable Development Goals': r'\bsustainable development goals?\b',
    'Paris Agreement': r'\bparis agreement\b',
}
TERM_WORDS = ('agreement', 'biological', 'change', 'climate', 'development', 'diversities', 'diversity', 'gender', 'genders',
              'goal', 'goals', 'paris', 'sex', 'sexes', 'sustainable')
MAX_HEADER_LINES = 6

committee_names = {'first': 1, 'second': 2, 'third': 3, 'fourth': 4, 'fifth': 5, 'sixth': 6}
furniture = re.compile(r'^\s*(A/RES/\S+.*|\d+/\d+\s*|\S*\d{2}-\d{5}\S*.*|\*\d+\*|Distr\.:.*|United Nations\s+A/RES.*|General Assembly\s+Distr.*)$')
operative_start = re.compile(r'^\s*(\d{1,3})\.\s+(["“(]?[A-Z][a-z]+)')
title_line = re.compile(r'^\s*(\d{2})/(\d{1,3})\s*([A-Z]\s*)?\.\s+(\S.*)$')


def clean_pages(raw):
    pages = []
    for page in raw.split('\f'):
        lines = page.split('\n')
        cut = next((i for i, line in enumerate(lines) if re.match(r'^\s*_{5,}\s*$', line)), len(lines))
        pages.append([line for line in lines[:cut] if not furniture.match(line)])
    return pages


def split_parts(raw):
    """Return (letter, text) for each lettered resolution in a combined document, else one unlettered part."""
    starts = []
    offset = 0
    for line in raw.split('\n'):
        hit = title_line.match(line)
        if hit and hit.group(3):
            starts.append((offset, hit.group(3).strip()))
        offset += len(line) + 1
    if len(starts) < 2:
        return [('', raw, raw[:4000])]
    bounds = [p for p, _ in starts] + [len(raw)]
    return [(letter, raw[bounds[i]:bounds[i + 1]], raw[:4000]) for i, (_, letter) in enumerate(starts)]


def parse(path):
    raw = path.read_text(errors='replace')
    return parse_text(raw, raw[:4000])


def parse_text(raw, head):
    committee = 0
    found = re.search(r'on the report of the\s+(?:Special Political and Decolonization Committee\s*\(?)?(\w+)\s+Committee', head)
    if found and found.group(1).lower() in committee_names:
        committee = committee_names[found.group(1).lower()]
    elif re.search(r'Special Political', head):
        committee = 4
    pages = clean_pages(raw)
    title, started, paragraphs, current, kind = '', False, [], [], 'preamble'
    cited = []
    for lines in pages:
        indents = [len(l) - len(l.lstrip()) for l in lines if l.strip()]
        base = Counter(indents).most_common(1)[0][0] if indents else 0
        for index, line in enumerate(lines):
            if not title:
                hit = title_line.match(line)
                if hit:
                    parts, j = [hit.group(4).strip()], index + 1
                    while j < len(lines) and lines[j].strip() and not lines[j].strip().startswith('The General Assembly'):
                        parts.append(lines[j].strip())
                        j += 1
                    title = ' '.join(parts)
                continue
            if not started:
                if re.match(r'^\s*The General Assembly\s*,?\s*$', line):
                    started = True
                continue
            if re.match(r'^\s*Annex(es)?\s*(I|1)?\s*$', line):
                if current:
                    paragraphs.append((kind, ' '.join(current)))
                parse_text.cited = cited
                return title, committee, paragraphs
            if not line.strip():
                continue
            if kind == 'preamble' and len(paragraphs) < 4:
                cited += re.findall(r'\b(\d{2})/(\d{1,3})\b', line)
            indent = len(line) - len(line.lstrip())
            numbered = operative_start.match(line)
            if numbered or indent >= base + 3:
                if current:
                    paragraphs.append((kind, ' '.join(current)))
                current = []
                if numbered:
                    kind = 'operative'
            current.append(line.strip())
    if current:
        paragraphs.append((kind, ' '.join(current)))
    parse_text.cited = cited
    return title, committee, paragraphs


def normalize(text):
    text = text.lower().replace('­', '')
    text = re.sub(r'(\w)\s*-\s+(\w)', r'\1\2', text)
    text = re.sub(r'\d+', ' ', text)
    text = re.sub(r'[^a-z ]', ' ', text)
    return text.split()


def shingles(words, size=None):
    size = size or SHINGLE
    return {' '.join(words[i:i + size]) for i in range(len(words) - size + 1)}


def title_key(title):
    return ' '.join(normalize(title))



work = Path('.')
fields = ['symbol', 'session', 'status', 'http', 'bytes', 'sha256', 'text_chars', 'url', 'fetched_utc']
lock = threading.Lock()
agent = {'User-Agent': 'Mozilla/5.0 (academic research on resolution text reuse; low request rate)'}
endpoint = 'https://documents.un.org/api/symbol/access'


def load_manifest():
    if not (work / 'manifest.csv').exists():
        return {}
    with (work / 'manifest.csv').open(newline='') as handle:
        return {row['symbol']: row for row in csv.DictReader(handle)}


def record(row):
    with lock:
        new = not (work / 'manifest.csv').exists()
        with (work / 'manifest.csv').open('a', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            if new:
                writer.writeheader()
            writer.writerow(row)


def fetch(symbol):
    session = symbol.split('/')[2]
    name = symbol.replace('/', '_').replace(' ', '')
    pdf = work / 'pdf' / session / (name + '.pdf')
    text = work / 'text' / session / (name + '.txt')
    pdf.parent.mkdir(parents=True, exist_ok=True)
    text.parent.mkdir(parents=True, exist_ok=True)
    row = dict(symbol=symbol, session=session, status='error', http='', bytes=0, sha256='', text_chars=0, url='',
               fetched_utc=datetime.now(timezone.utc).isoformat())
    import requests  # needed only to download; the measure mode runs without it
    for attempt in range(3):
        try:
            reply = requests.get(endpoint + '?s=' + quote(symbol, safe='/') + '&l=en&t=pdf', headers=agent, timeout=90)
            row['http'] = reply.status_code
            row['url'] = reply.url
            if reply.content[:4] == b'%PDF':
                pdf.write_bytes(reply.content)
                subprocess.run(['pdftotext', '-layout', str(pdf), str(text)], check=False)
                row.update(status='ok', bytes=len(reply.content), sha256=hashlib.sha256(reply.content).hexdigest(),
                           text_chars=len(text.read_text(errors='replace').strip()) if text.exists() else 0)
                break
            if reply.status_code == 200:
                row['status'] = 'missing'
                break
        except requests.RequestException as error:
            row['status'] = 'error:' + type(error).__name__
        time.sleep(5 * (attempt + 1))
    time.sleep(0.4)
    record(row)
    return row


def listed_symbols():
    path = project / 'results/adoptions/history/all_adoption_units.csv'
    out = {}
    with path.open(newline='') as handle:
        for row in csv.DictReader(handle):
            out.setdefault(int(row['session']), []).append(row['canonical_resolution'])
    return out




def fetch_sessions(first, last, directory):
    global work
    work = directory
    work.mkdir(parents=True, exist_ok=True)
    listed = listed_symbols()
    for session in range(last, first - 1, -1):
        done = load_manifest()
        if session in listed:
            todo = [s for s in dict.fromkeys(listed[session]) if s not in done]
            with ThreadPoolExecutor(4) as pool:
                list(pool.map(fetch, todo))
            continue
        number, misses = 1, 0
        while misses < 25 and number < 450:
            batch = [f'A/RES/{session}/{n}' for n in range(number, number + 8)]
            with ThreadPoolExecutor(4) as pool:
                results = list(pool.map(fetch, [s for s in batch if s not in done]))
            for symbol in batch:
                status = next((r['status'] for r in results if r['symbol'] == symbol), done.get(symbol, {}).get('status'))
                misses = 0 if status == 'ok' else misses + 1
            number += 8


def fetch_letters(directory):
    global work
    work = directory
    rows = list(load_manifest().values())
    done = {r['symbol'] for r in rows}
    last_ok = {}
    for r in rows:
        hit = re.match(r'A/RES/(\d+)/(\d+)', r['symbol'])
        if r['status'] == 'ok':
            last_ok[hit.group(1)] = max(last_ok.get(hit.group(1), 0), int(hit.group(2)))
    bases = sorted({'/'.join(r['symbol'].split('/')[:3]) + '/' + re.match(r'A/RES/\d+/(\d+)', r['symbol']).group(1) for r in rows
                    if r['status'] == 'missing' and int(re.match(r'A/RES/\d+/(\d+)', r['symbol']).group(1)) < last_ok.get(r['symbol'].split('/')[2], 0)})

    def attempt(base):
        for end in 'BCDEFGHIJKLMNOPQRSTUVWX':
            symbol = f'{base} A-{end}'
            if symbol in done or fetch(symbol)['status'] == 'ok':
                return

    with ThreadPoolExecutor(4) as pool:
        list(pool.map(attempt, bases))


def adoption_year(raw, session):
    """Calendar year of adoption, read from the plenary-meeting line that closes each resolution."""
    found = re.findall(r'plenary\s+meeting\s*\n\s*(\d{1,2})\s+([A-Z][a-z]+)\s+(\d{4})', raw)
    years = [int(y) for _, m, y in found if m.lower() in MONTHS and session + 1944 <= int(y) <= session + 1947]
    return years[0] if years else None


def numbering_consistent(raw):
    """True when the numbered operative paragraphs found in the text run 1, 2, ... without a gap."""
    body = raw.split('Annex')[0]
    numbers = {int(n) for n in re.findall(r'(?m)^\s*(\d{1,3})\.\s+["“(]?[A-Z][a-z]+', body)}
    return bool(numbers) and numbers == set(range(1, max(numbers) + 1))


def measure(archive, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    bundle = zipfile.ZipFile(archive)
    requested = list(csv.DictReader(io.TextIOWrapper(bundle.open('manifest.csv'), encoding='utf-8', newline='')))
    manifest = [r for r in requested if r['status'] == 'ok']
    names = set(bundle.namelist())
    docs = {}
    for row in manifest:
        session = int(row['session'])
        member = f"text/{session}/{row['symbol'].replace('/', '_').replace(' ', '')}.txt"
        if member not in names or int(row['text_chars']) < 400:
            docs[row['symbol']] = dict(symbol=row['symbol'], session=session, usable=False, title='', committee=0,
                                       paragraphs=[], cited=[], year=None, numbering=False)
            continue
        raw = bundle.read(member).decode('utf-8', errors='replace')
        for letter, part, head in split_parts(raw):
            title, committee, paragraphs = parse_text(part, head)
            kept = [(k, normalize(t)) for k, t in paragraphs]
            kept = [(k, w) for k, w in kept if len(w) >= MIN_WORDS]
            symbol = row['symbol'].split(' ')[0] + (' ' + letter if letter else '') if ' ' in row['symbol'] else row['symbol']
            docs[symbol] = dict(symbol=symbol, session=session, usable=bool(title and kept), title=title, committee=committee,
                                paragraphs=kept, cited=[f'A/RES/{a}/{b}' for a, b in parse_text.cited],
                                year=adoption_year(part, session), numbering=numbering_consistent(part))

    def write(name, table):
        if table:
            with (out / name).open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(table[0].keys()))
                writer.writeheader()
                writer.writerows(table)

    write('manifest_summary.csv', [dict(documents=len(manifest), combined_lettered_documents=sum(' A-' in r['symbol'] for r in manifest),
                                        symbols_without_document=sum(r['status'] == 'missing' for r in requested))])
    write('documents.csv', [dict(symbol=d['symbol'], session=d['session'], adoption_year=d['year'] or '', usable=d['usable'],
                                 committee=d['committee'], numbering_consistent=d['numbering'],
                                 operative_paragraphs=sum(1 for k, _ in d['paragraphs'] if k == 'operative'),
                                 preamble_paragraphs=sum(1 for k, _ in d['paragraphs'] if k == 'preamble'), title=d['title'])
                            for d in docs.values()])

    by_title = defaultdict(list)
    for d in docs.values():
        if d['usable']:
            by_title[title_key(d['title'])].append(d)
    pairs = []
    for key, group in by_title.items():
        if not key:
            continue
        per_session = defaultdict(list)
        for d in group:
            per_session[d['session']].append(d)
        # a title used by more than one resolution in a session is ambiguous and skipped
        sessions = sorted(s for s, ds in per_session.items() if len(ds) == 1)
        for earlier, later in zip(sessions, sessions[1:]):
            if later - earlier <= MAX_GAP:
                pairs.append((per_session[earlier][0], per_session[later][0]))
    same_title = {(id(a), id(b)) for a, b in pairs}

    # sensitivity: link through resolutions cited in the opening preambular paragraphs when the title changed
    linked = {later['symbol'] for _, later in pairs}
    chain_pairs = []
    for d in docs.values():
        if not d['usable'] or d['symbol'] in linked:
            continue
        words, best = set(normalize(d['title'])), None
        for symbol in d['cited']:
            c = docs.get(symbol)
            if not c or not c['usable'] or not 0 < d['session'] - c['session'] <= MAX_GAP:
                continue
            other = set(normalize(c['title']))
            if len(words & other) / max(1, len(words | other)) >= 0.5 and (best is None or c['session'] > best['session']):
                best = c
        if best:
            chain_pairs.append((best, d))

    rows = []
    for prior, later in pairs + chain_pairs:
        prior_op = [w for k, w in prior['paragraphs'] if k == 'operative']
        later_op = [w for k, w in later['paragraphs'] if k == 'operative']
        if not prior_op or not later_op:
            continue
        row = dict(prior=prior['symbol'], later=later['symbol'], prior_session=prior['session'], later_session=later['session'],
                   session_year=later['session'] + 1945, adoption_year=later['year'] or '',
                   gap=later['session'] - prior['session'], committee=later['committee'], title=later['title'],
                   prior_operative=len(prior_op), later_operative=len(later_op),
                   link='same_title' if (id(prior), id(later)) in same_title else 'cited_predecessor')
        for size in (4, 3, 5):
            later_all = set().union(*(shingles(w, size) for _, w in later['paragraphs']))
            prior_all = set().union(*(shingles(w, size) for _, w in prior['paragraphs']))
            carried = [len(shingles(w, size) & later_all) / max(1, len(shingles(w, size))) for w in prior_op]
            if size == 4:
                inherited = [len(shingles(w, size) & prior_all) / max(1, len(shingles(w, size))) for w in later_op]
                row['mean_containment'] = sum(carried) / len(carried)
                for t in THRESHOLDS:
                    row[f'carried_{int(t * 100)}'] = sum(c >= t for c in carried) / len(carried)
                row['new_80'] = sum(c < 0.8 for c in inherited) / len(inherited)
            else:
                row[f'carried_80_k{size}'] = sum(c >= 0.8 for c in carried) / len(carried)
        rows.append(row)
    write('pairs.csv', rows)

    columns = ['carried_60', 'carried_80', 'carried_90', 'carried_80_k3', 'carried_80_k5', 'mean_containment', 'new_80']

    def summarize(selected, keys):
        groups = defaultdict(list)
        for r in selected:
            if all(r[k] != '' for k in keys):
                groups[tuple(r[k] for k in keys)].append(r)
        table = []
        for group_key, members in sorted(groups.items()):
            entry = dict(zip(keys, group_key))
            entry['pairs'] = len(members)
            for column in columns:
                entry[column] = round(100 * sum(m[column] for m in members) / len(members), 2)
            table.append(entry)
        return table

    main_rows = [r for r in rows if r['link'] == 'same_title']
    write('series_session.csv', summarize(main_rows, ['session_year']))
    write('series_adoption_year.csv', summarize(main_rows, ['adoption_year']))
    write('series_annual_gap.csv', summarize([r for r in main_rows if r['gap'] == 1], ['session_year']))
    write('series_by_committee.csv', summarize(main_rows, ['committee', 'session_year']))
    write('series_with_cited_predecessors.csv', summarize(rows, ['session_year']))
    coverage = defaultdict(Counter)
    for d in docs.values():
        c = coverage[d['session']]
        c['texts'] += 1
        c['usable'] += int(d['usable'])
        c['dated'] += int(d['year'] is not None)
        c['numbering_consistent'] += int(d['numbering'])
    paired = Counter(r['later_session'] for r in main_rows)
    write('coverage.csv', [dict(session=s, session_year=s + 1945, texts=c['texts'], usable=c['usable'], adoption_year_found=c['dated'],
                                numbering_consistent=c['numbering_consistent'], paired_as_successor=paired[s])
                           for s, c in sorted(coverage.items())])


def strip_running_headers(raw):
    """Remove from every page after the first the running header that repeats the symbol and the title.

    The header is the block of lines before the first blank line of the page. It is removed when it has at most
    MAX_HEADER_LINES lines and either holds the document symbol or consists of lines whose letters and digits recur
    on the first page, which the title does. Pages without such a block are left as they are.
    """
    pages = raw.split('\f')
    first = ' '.join(re.sub(r'[^a-z0-9]+', ' ', pages[0].lower()).split())
    kept = [pages[0]]
    for page in pages[1:]:
        lines = page.split('\n')
        block = lines[:next((i for i, line in enumerate(lines) if not line.strip()), len(lines))]
        words = [' '.join(re.sub(r'[^a-z0-9]+', ' ', line.lower()).split()) for line in block]
        if block and len(block) <= MAX_HEADER_LINES and (any('A/RES/' in line for line in block) or all(w in first for w in words)):
            lines = lines[len(block):]
        kept.append('\n'.join(lines))
    return '\f'.join(kept)


def term_text(text):
    """Lower-case letters and spaces, with a term word that pdftotext split by a hyphen at a line end joined again."""
    text = re.sub(r'\d+', ' ', text.lower().replace('­', ''))
    for word in TERM_WORDS:
        text = re.sub(r'\b' + r'(?:-\s+)?'.join(word) + r'\b', word, text)
    return re.sub(r'[^a-z]+', ' ', text)


def count_terms(paragraphs):
    """Occurrences of each term in the paragraphs joined into one text, and in the operative paragraphs joined alone."""
    # joining the paragraphs counts an occurrence that a page break split between two paragraphs
    whole = term_text(' '.join(text for _, text in paragraphs))
    operative = term_text(' '.join(text for kind, text in paragraphs if kind == 'operative'))
    return {term: [len(re.findall(pattern, whole)), len(re.findall(pattern, operative))] for term, pattern in TERMS.items()}


def term_outcome(before, after):
    if before == 0:
        return 'added' if after else 'absent'
    if after == 0:
        return 'dropped'
    return 'reduced' if after < before else 'unchanged' if after == before else 'increased'


def load_pair_texts(bundle, symbols):
    """Parsed paragraphs of the listed symbols, read from the archive as the measure mode reads them, running headers removed."""
    requested = csv.DictReader(io.TextIOWrapper(bundle.open('manifest.csv'), encoding='utf-8', newline=''))
    names = set(bundle.namelist())
    bases = {s.split(' ')[0] for s in symbols}
    docs = {}
    for row in requested:
        if row['status'] != 'ok' or row['symbol'].split(' ')[0] not in bases:
            continue
        session = int(row['session'])
        member = f"text/{session}/{row['symbol'].replace('/', '_').replace(' ', '')}.txt"
        if member not in names or int(row['text_chars']) < 400:
            continue
        raw = strip_running_headers(bundle.read(member).decode('utf-8', errors='replace'))
        for letter, part, head in split_parts(raw):
            title, committee, paragraphs = parse_text(part, head)
            symbol = row['symbol'].split(' ')[0] + (' ' + letter if letter else '') if ' ' in row['symbol'] else row['symbol']
            docs[symbol] = dict(year=adoption_year(part, session), paragraphs=paragraphs)
    return docs


def sha256_of(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_rows(path, table):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0].keys()))
        writer.writeheader()
        writer.writerows(table)


def terms(archive, out_dir):
    out = Path(out_dir)
    with (out / 'pairs.csv').open(newline='') as handle:
        pairs = [r for r in csv.DictReader(handle) if r['link'] == 'same_title' and r['adoption_year']]
    verbatim_years = sorted({int(r['adoption_year']) for r in pairs})
    bundle = zipfile.ZipFile(archive)
    docs = load_pair_texts(bundle, {r['prior'] for r in pairs} | {r['later'] for r in pairs})
    count_rows, verbatim_rows, unread, mismatched = [], [], [], []
    for r in pairs:
        if r['prior'] not in docs or r['later'] not in docs:
            unread.append((r['prior'], r['later']))
            continue
        prior, later = docs[r['prior']], docs[r['later']]
        year = int(r['adoption_year'])
        base = dict(prior=r['prior'], later=r['later'], prior_session=r['prior_session'], later_session=r['later_session'],
                    prior_adoption_year=prior['year'] or '', adoption_year=year, gap=r['gap'], committee=r['committee'],
                    title=r['title'], annual=prior['year'] == year - 1)
        if year in TERM_YEARS:
            before, after = count_terms(prior['paragraphs']), count_terms(later['paragraphs'])
            for term, pattern in TERMS.items():
                count_rows.append(dict(base, term=term, pattern=pattern, prior_count=before[term][0], later_count=after[term][0],
                                       prior_operative_count=before[term][1], later_operative_count=after[term][1],
                                       change=after[term][0] - before[term][0], outcome=term_outcome(before[term][0], after[term][0])))
        # operative paragraphs of at least MIN_WORDS words in each text, and all paragraphs of the successor, normalized as in the measure mode
        prior_op = [w for w in (normalize(t) for k, t in prior['paragraphs'] if k == 'operative') if len(w) >= MIN_WORDS]
        later_op = [w for w in (normalize(t) for k, t in later['paragraphs'] if k == 'operative') if len(w) >= MIN_WORDS]
        later_any = {' '.join(normalize(t)) for _, t in later['paragraphs']}
        later_op_keys = {' '.join(w) for w in later_op}
        if len(prior_op) != int(r['prior_operative']):
            mismatched.append((r['prior'], r['later'], year))
        if not prior_op:
            continue
        in_any = sum(' '.join(w) in later_any for w in prior_op)
        in_op = sum(' '.join(w) in later_op_keys for w in prior_op)
        verbatim_rows.append(dict(base, prior_operative=len(prior_op), later_operative=len(later_op),
                                  verbatim_paragraphs=in_any, verbatim_any=in_any / len(prior_op), verbatim_operative=in_op / len(prior_op),
                                  carried_80=float(r['carried_80'])))
    write_rows(out / 'term_counts.csv', count_rows)
    write_rows(out / 'verbatim_paragraphs.csv', verbatim_rows)

    def subsets(rows, years):
        for year in years:
            in_year = [x for x in rows if x['adoption_year'] == year]
            groups = [('all pairs', in_year), ('annual pairs', [x for x in in_year if x['annual']])]
            groups += [(f'committee {c}', [x for x in in_year if int(x['committee']) == c]) for c in range(7)]
            for subset, selected in groups:
                if selected:
                    yield year, subset, selected

    def share(part, whole):
        return 100 * part / whole if whole else ''

    term_summary = []
    for year, subset, rows in subsets(count_rows, TERM_YEARS):
        for term in TERMS:
            hits = [x for x in rows if x['term'] == term]
            tally = Counter(x['outcome'] for x in hits)
            with_term = sum(x['prior_count'] > 0 for x in hits)
            term_summary.append(dict(adoption_year=year, subset=subset, term=term, pairs=len(hits),
                                     predecessor_with_term=with_term, successor_with_term=sum(x['later_count'] > 0 for x in hits),
                                     dropped=tally['dropped'], reduced=tally['reduced'], unchanged=tally['unchanged'],
                                     increased=tally['increased'], added=tally['added'], absent=tally['absent'],
                                     dropped_share=share(tally['dropped'], with_term), reduced_share=share(tally['reduced'], with_term),
                                     dropped_or_reduced_share=share(tally['dropped'] + tally['reduced'], with_term),
                                     predecessor_mentions=sum(x['prior_count'] for x in hits),
                                     successor_mentions=sum(x['later_count'] for x in hits)))
    write_rows(out / 'term_summary.csv', term_summary)
    verbatim_summary = []
    for year, subset, rows in subsets(verbatim_rows, verbatim_years):
        paragraphs = sum(x['prior_operative'] for x in rows)
        verbatim_summary.append(dict(adoption_year=year, subset=subset, pairs=len(rows),
                                     verbatim_any=100 * sum(x['verbatim_any'] for x in rows) / len(rows),
                                     verbatim_operative=100 * sum(x['verbatim_operative'] for x in rows) / len(rows),
                                     carried_80=100 * sum(x['carried_80'] for x in rows) / len(rows),
                                     pairs_fully_verbatim=sum(x['verbatim_any'] == 1 for x in rows),
                                     pairs_none_verbatim=sum(x['verbatim_any'] == 0 for x in rows),
                                     predecessor_paragraphs=paragraphs, verbatim_paragraphs=sum(x['verbatim_paragraphs'] for x in rows),
                                     pooled_verbatim_any=100 * sum(x['verbatim_paragraphs'] for x in rows) / paragraphs))
    write_rows(out / 'verbatim_summary.csv', verbatim_summary)
    provenance = dict(
        script='code/09_measure_text_reuse.py', mode='terms', python=sys.version.split()[0],
        inputs=dict(archive=dict(path=str(archive), sha256=sha256_of(archive)),
                    pairs=dict(path=str(out / 'pairs.csv'), sha256=sha256_of(out / 'pairs.csv'))),
        parameters=dict(
            term_years=list(TERM_YEARS), verbatim_years=verbatim_years, link='same_title', terms=TERMS,
            rejoined_words=list(TERM_WORDS), min_words=MIN_WORDS, max_header_lines=MAX_HEADER_LINES,
            running_headers='on every page after the first, the block of lines before the first blank line is removed when it '
                            'has at most max_header_lines lines and holds the document symbol or repeats text of the first page; '
                            'the title that the header repeats is therefore not counted, and a paragraph that a page break '
                            'divided is joined',
            term_counting='the paragraphs after The General Assembly and before any annex, preamble and operative, are joined '
                          'into one text, lower-cased, a term word split by a hyphen at a line end is joined, digits and '
                          'punctuation become spaces, and whole-word matches of each pattern are counted; the operative count '
                          'repeats this on the operative paragraphs alone',
            verbatim='share of the predecessor operative paragraphs of at least min_words words whose normalized words '
                     '(lower case, digits removed, letters only) equal those of a paragraph of the successor; verbatim_any '
                     'matches any successor paragraph, verbatim_operative only operative ones; prior_operative and '
                     'later_operative count the operative paragraphs of at least min_words words in each text after header removal',
            annual='prior_adoption_year equals adoption_year minus one'),
        outputs=dict(term_counts=dict(path=str(out / 'term_counts.csv'), rows=len(count_rows)),
                     verbatim_paragraphs=dict(path=str(out / 'verbatim_paragraphs.csv'), rows=len(verbatim_rows)),
                     term_summary=dict(path=str(out / 'term_summary.csv'), rows=len(term_summary)),
                     verbatim_summary=dict(path=str(out / 'verbatim_summary.csv'), rows=len(verbatim_summary))),
        checks=dict(pairs_selected=len(pairs), pairs_in_term_years=sum(int(r['adoption_year']) in TERM_YEARS for r in pairs),
                    documents_read=len(docs), pairs_without_text=unread,
                    pairs_with_operative_count_differing_from_measure=dict(
                        count=len(mismatched),
                        note='the measure mode parses running headers as text, so its operative count can include a header '
                             'line; the pairs in the term years are listed',
                        term_years=[p[:2] for p in mismatched if p[2] in TERM_YEARS])))
    with (out / 'terms_provenance.json').open('w') as handle:
        json.dump(provenance, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'measure':
        measure(sys.argv[2], sys.argv[3])
    elif mode == 'fetch':
        fetch_sessions(int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4]))
    elif mode == 'letters':
        fetch_letters(Path(sys.argv[2]))
    elif mode == 'terms':
        terms(sys.argv[2], sys.argv[3])
