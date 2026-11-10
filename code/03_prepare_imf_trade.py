"""Prepare trade datasets from archived IMF responses; download only with --fetch."""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

BASE = 'https://api.imf.org/external/sdmx/3.0/'


def fetch(output, resource, parameters=None):
    import requests

    response = requests.get(BASE + resource, params=parameters,
                            headers={'Accept': 'application/json'}, timeout=180)
    response.raise_for_status()
    output.write_bytes(response.content)
    output.with_suffix('.manifest.json').write_text(json.dumps({
        'url': response.url, 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
        'headers': dict(response.headers),
        'sha256': hashlib.sha256(response.content).hexdigest()}, indent=2))


def member(dimension, index):
    if index is None:
        return None
    if 'values' not in dimension:
        return index
    value = dimension['values'][int(index)]
    return value.get('id', value.get('value'))


def parse(path):
    data = json.loads(path.read_text())['data']
    records = []
    for dataset in data['dataSets']:
        structure = data['structures'][dataset['structure']]
        dimensions = structure['dimensions']
        attributes = structure['attributes']
        for key, series in dataset['series'].items():
            record = {dim['id']: member(dim, code) for dim, code in
                      zip(dimensions['series'], key.split(':'))}
            record.update({attr['id']: member(attr, code) for attr, code in
                           zip(attributes['series'], series['attributes'])})
            for time, observation in series.get('observations', {}).items():
                row = record.copy()
                row.update({dim['id']: member(dim, code) for dim, code in
                            zip(dimensions['observation'], time.split(':'))})
                row['OBS_VALUE'] = float(observation[0]) if observation[0] is not None else np.nan
                row.update({attr['id']: member(attr, code) for attr, code in
                            zip(attributes['observation'], observation[1:])})
                records.append(row)
    frame = pd.DataFrame(records).rename(columns={
        'COUNTRY': 'reporter', 'COUNTERPART_COUNTRY': 'partner',
        'OBS_VALUE': 'exports_usd', 'TIME_PERIOD': 'period'})
    frame['year'] = frame.period.str[:4].astype(int)
    assert not frame.duplicated(['reporter', 'partner', 'period']).any()
    return frame


def verify_archive(work):
    manifest = json.loads((work / 'source_manifest.json').read_text())
    for name, entry in manifest['files'].items():
        actual = hashlib.sha256((work / name).read_bytes()).hexdigest()
        if actual != entry['sha256']:
            raise ValueError('Archived IMF source hash changed: ' + name)


def record_archive(work):
    paths = [p for p in work.glob('*.json') if p.name != 'source_manifest.json']
    (work / 'source_manifest.json').write_text(json.dumps({
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'purpose': 'Hashes recorded after an explicitly requested source refresh.',
        'files': {p.name: {'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                           'bytes': p.stat().st_size} for p in sorted(paths)}}, indent=2))


def audit(root):
    work = root / 'data' / 'imf'
    output = root / 'results' / 'trade' / 'imf'
    output.mkdir(parents=True, exist_ok=True)
    verify_archive(work)
    actors = set((root / 'data' / 'actors.txt').read_text().splitlines())
    annual = parse(work / 'imf_all_annual.json')
    # Descriptive trade shares also use recorded economies outside the UN panel.
    # The source's regional aggregates use other code formats and stay excluded.
    economies = annual.reporter.str.fullmatch(r'[A-Z]{3}') & annual.partner.str.fullmatch(r'[A-Z]{3}')
    all_annual = annual[economies & (annual.reporter != annual.partner) &
                        (annual.exports_usd > 0)].copy()
    all_annual.to_csv(output / 'imf_all_economies_annual.csv', index=False)
    annual = annual[(annual.reporter.isin(actors)) & (annual.partner.isin(actors)) &
                    (annual.reporter != annual.partner) & (annual.exports_usd > 0)].copy()
    annual.to_csv(output / 'imf_current_annual.csv', index=False)
    annual.groupby(['year', 'DERIVATION_TYPE', 'STATUS'], dropna=False).agg(
        records=('exports_usd', 'size'), nominal_exports=('exports_usd', 'sum')).reset_index().to_csv(
            output / 'imf_current_attributes.csv', index=False)
    cached = pd.read_csv(root / 'results' / 'trade' / 'baseline' / 'cached_terminal.csv')
    cached = cached[cached.reporter.isin(actors) & cached.partner.isin(actors)]
    match = cached.merge(annual, on=['reporter', 'partner', 'year'], how='outer',
                         suffixes=('_cached', '_current'), indicator=True)
    match['relative_revision'] = match.exports_usd_current / match.exports_usd_cached - 1
    match['abs_relative_revision'] = match.relative_revision.abs()
    match['identical_value'] = match.exports_usd_current == match.exports_usd_cached
    match.to_csv(output / 'imf_vintage_comparison.csv', index=False)
    match.groupby(['year', '_merge'], observed=True).agg(
        records=('reporter', 'size'),
        identical_values=('identical_value', 'sum'),
        median_abs_relative_revision=('abs_relative_revision', 'median'),
        nominal_total_cached=('exports_usd_cached', 'sum'),
        nominal_total_current=('exports_usd_current', 'sum')).reset_index().to_csv(
            output / 'imf_vintage_summary.csv', index=False)
    monthly_path = work / 'imf_all_monthly_2025.json'
    if monthly_path.exists():
        monthly = parse(monthly_path)
        all_monthly = monthly[monthly.reporter.str.fullmatch(r'[A-Z]{3}') &
                              monthly.partner.str.fullmatch(r'[A-Z]{3}') &
                              (monthly.reporter != monthly.partner) & monthly.exports_usd.notna()].copy()
        assert (all_monthly.exports_usd >= 0).all()
        all_sums = all_monthly.groupby(['reporter', 'partner', 'year']).agg(
            months=('period', 'nunique'), monthly_sum=('exports_usd', 'sum'),
            earliest_month=('period', 'min'), latest_month=('period', 'max')).reset_index()
        all_coverage = all_annual[all_annual.year == 2025].merge(
            all_sums, on=['reporter', 'partner', 'year'], how='left')
        all_coverage['monthly_annual_ratio'] = all_coverage.monthly_sum / all_coverage.exports_usd
        all_coverage['annual_matches_months'] = np.isclose(
            all_coverage.monthly_sum, all_coverage.exports_usd, rtol=1e-6, atol=12)
        all_coverage.to_csv(output / 'imf_all_economies_monthly_coverage.csv', index=False)
        monthly = monthly[monthly.reporter.isin(actors) & monthly.partner.isin(actors) &
                          (monthly.reporter != monthly.partner) & monthly.exports_usd.notna()]
        assert (monthly.exports_usd >= 0).all()
        sums = monthly.groupby(['reporter', 'partner', 'year']).agg(
            months=('period', 'nunique'), monthly_sum=('exports_usd', 'sum'),
            earliest_month=('period', 'min'), latest_month=('period', 'max')).reset_index()
        coverage = annual[annual.year == 2025].merge(sums, on=['reporter', 'partner', 'year'], how='left')
        coverage['monthly_annual_ratio'] = coverage.monthly_sum / coverage.exports_usd
        coverage['annual_matches_months'] = np.isclose(coverage.monthly_sum,
                                                       coverage.exports_usd, rtol=1e-6, atol=12)
        coverage.to_csv(output / 'imf_monthly_coverage.csv', index=False)
        coverage.groupby('months', dropna=False).agg(
            annual_directions=('reporter', 'size'),
            nominal_exports=('exports_usd', 'sum'),
            exact_annual_matches=('annual_matches_months', 'sum')).reset_index().to_csv(
                output / 'imf_monthly_summary.csv', index=False)
        monthly.groupby(['DERIVATION_TYPE', 'STATUS'], dropna=False).agg(
            records=('exports_usd', 'size')).reset_index().to_csv(
                output / 'imf_monthly_attributes.csv', index=False)
    source_files = [work / 'imf_all_annual.json', monthly_path,
                    work / 'imf_derivation_codes.json']
    (output / 'imf_audit_manifest.json').write_text(json.dumps({
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'inputs': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in source_files if p.exists()},
        'scope': 'Fresh source vintage. Current flags cannot recover metadata discarded by the earlier download.',
        'expanded_economy_scope': 'Additional descriptive files retain positive nonself observations for three-letter source economy codes, including recorded countries and territories outside the political model panel. Regional aggregates are excluded. Incoming flows are partner-reported exports, and an absent reporter, including Taiwan in this archive, is not imputed as zero.',
        'attribute_meanings': 'STATUS is preserved literally because the fetched scheme defines a String and supplies no authoritative mapping for e/r. A blank flag does not establish direct or complete reporting. DERIVATION_TYPE M is Mixed-type data. Twelve populated monthly observations do not imply unestimated values.',
        'mode': 'Archived offline audit; acquisition occurs only when --fetch is explicitly invoked.'
    }, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--fetch', action='store_true')
    args = parser.parse_args()
    args.root = args.root.resolve()
    work = args.root / 'data' / 'imf'
    if args.fetch:
        work.mkdir(parents=True, exist_ok=True)
        for frequency, filename, window in [
            ('A', 'imf_all_annual.json', 'ge:2024-A1+le:2025-A1'),
            ('M', 'imf_all_monthly_2025.json', 'ge:2025-M01+le:2025-M12')]:
            fetch(work / filename, 'data/dataflow/IMF.STA/IMTS/+/*.XG_FOB_USD.*.' + frequency,
                  {'dimensionAtObservation': 'TIME_PERIOD', 'attributes': 'dsd', 'measures': 'all',
                   'c[TIME_PERIOD]': window})
        fetch(work / 'imf_derivation_codes.json', 'structure/codelist/IMF/CL_DERIVATION_TYPE/+')
        record_archive(work)
    audit(args.root)
