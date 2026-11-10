"""Assemble manuscript tables and summary checks from retained analysis outputs."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/revision'
OUT.mkdir(parents=True, exist_ok=True)
SOURCES = {}

def rows(relative):
    path = ROOT / relative
    SOURCES[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return list(csv.DictReader(path.open()))

def table(caption, label, columns, headings, values, note=''):
    text = '\\begin{table}[ht]\n\\centering\\small\n'
    text += '\\caption{' + caption + '}\n\\label{' + label + '}\n'
    text += '\\begin{tabularx}{\\textwidth}{' + columns + '}\\toprule\n'
    text += ' & '.join(headings) + r'\\\midrule' + '\n'
    text += '\n'.join(' & '.join(map(str, row)) + r'\\' for row in values)
    text += '\n\\bottomrule\n\\end{tabularx}\n'
    if note:
        text += '\\begin{minipage}{\\textwidth}\\footnotesize ' + note + '\\end{minipage}\n'
    return text + '\\end{table}\n'

variants = [
    ('cold_full', 'Cold starts'), ('recurring_both_years', 'Recurring titles'),
    ('equal_weights', 'Equal resolution weights'), ('us_recurring_only_change', 'Only American recurring ties change'),
    ('rho_random_walk', 'Persistence fixed at 0.999'), ('rho_free_bound', 'Persistence bound 0.9999'),
    ('regular_voters_only', 'Regular voters only'), ('truncate_2020', 'Record ending in 2020'),
    ('partial_following_75', 'Independent following, 25 percent'), ('partial_following_50', 'Independent following, 50 percent'),
    ('allies_divide', 'Allied agreement reduced'), ('partial_following_75_joint', 'Joint following, 25 percent'),
    ('partial_following_50_joint', 'Joint following, 50 percent'),
    ('partial_following_75_joint_controlled', 'Joint following on controlled network'),
]
comparison = []
for variant, label in variants:
    folder = f'results/models/un_variants/{variant}'
    if not all((ROOT/folder/name).exists() for name in ['variant_summary.csv','terminal_probabilities.csv','provenance.json']):
        continue
    summary = {r['quantity']: r for r in rows(folder + '/variant_summary.csv')}
    terminal = rows(folder + '/terminal_probabilities.csv')
    us = next(r for r in terminal if r['actor'] == 'USA')
    allies = [r for r in terminal if r['actor'] != 'USA']
    item = dict(variant=variant, label=label, year=us['year'], us_probability=us['probability'],
                max_allied_probability=max(float(r['probability']) for r in allies),
                allies_above_half=sum(float(r['probability']) > .5 for r in allies))
    for name in ['sigma_uv', 'rho_uv', 'pi_snap', 'threshold_norm_1', 'max_rhat', 'min_ess_bulk', 'min_ess_tail']:
        item[name] = float(summary[name]['mean'])
    item['diagnostics_pass'] = item['max_rhat'] <= 1.01 and item['min_ess_bulk'] >= 400 and item['min_ess_tail'] >= 400
    comparison.append(item)
with (OUT / 'variant_comparison.csv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(comparison[0]))
    writer.writeheader()
    writer.writerows(comparison)

fragments = {}
reported = [r for r in comparison if r['diagnostics_pass'] and r['variant'] in {
    'cold_full', 'recurring_both_years', 'equal_weights', 'regular_voters_only', 'truncate_2020',
    'partial_following_75_joint', 'partial_following_50_joint', 'partial_following_75_joint_controlled', 'allies_divide'}]
values = [[r['label'], f"{float(r['us_probability']):.5f}", f"{r['max_allied_probability']:.5f}",
           r['allies_above_half'], 'Pass'] for r in reported]
fragments['variants'] = table('Terminal American and allied break probabilities under alternative inputs and priors.',
    'tab:variant_summaries', '>{\\raggedright\\arraybackslash}Xrrrr',
    ['Fit', 'US', 'Highest ally', 'Allies above 0.5', 'Diagnostics'], values,
    'Each new fit uses two cold-start chains with 3,000 warmup and 12,000 sampling sweeps, retaining every fifth draw. '
    'The terminal year is 2025 except for the record ending in 2020. Diagnostics require maximum split $\\widehat R\\leq1.01$ '
    'and minimum bulk and tail effective sample sizes of 400 over all checked quantities. Complete run records, including exploratory fits, are retained in the replication materials.')
fragments['parameters'] = table('Posterior parameter means and the implied threshold in the new fits.',
    'tab:variant_parameters', '>{\\raggedright\\arraybackslash}Xrrrrrr',
    ['Fit', '$\\sigma$', '$\\rho$', '$\\pi$', 'Threshold', 'Max $\\widehat R$', 'Min ESS'],
    [[r['label'], f"{r['sigma_uv']:.4f}", f"{r['rho_uv']:.5f}", f"{r['pi_snap']:.4f}",
      f"{r['threshold_norm_1']:.3f}", f"{r['max_rhat']:.3f}", f"{min(r['min_ess_bulk'],r['min_ess_tail']):.0f}"] for r in reported],
    'Thresholds are evaluated at a position with squared length one. Complete parameter intervals, diagnostic families, '
    'per-chain probabilities and runtimes are retained in each fit directory.')

selection = rows('results/observed/ally_selection/selection.csv')
# the output records the operational selection rule and each member's criterion values
selected = [r for r in selection if r['selected_ally'] == 'TRUE']
(OUT/'selection_columns.json').write_text(json.dumps(list(selection[0]),indent=2)+'\n')
allies = rows('results/observed/followers/ally_table.csv')
fragments['allies'] = table('The fourteen allies on recurring American reversals in 2025.', 'tab:allies_individual',
    '>{\\raggedright\\arraybackslash}Xrrrr', ['Ally', 'Paired votes', 'Kept', 'Break probability', 'Move / $\\sigma$'],
    [[r['member_name'].title().replace('Korea, Republic Of', 'South Korea').replace('United Kingdom Of Great Britain And Northern Ireland','United Kingdom'), r['switch_pairs'], r['votes_kept'], f"{float(r['break_probability_2025']):.4f}", f"{float(r['displacement_sigma_2025']):.1f}"] for r in allies],
    'Break probabilities and displacement are from the released weighted full-agenda fit. The recurring vote columns are observed counts.')

reading = [r for r in rows('results/proposal_comparison/reading_agreement/batch_agreement.csv') if r['scope'].startswith('common_families') or r['scope']=='all_substantive_pairs']
fragments['readings'] = table('Agreement between the archived first and second readings.', 'tab:reading_agreement',
    '>{\\raggedright\\arraybackslash}Xrrrr', ['Comparison', 'Outcomes paired', 'Agreement (percent)', 'Kappa', 'Continuation agreement'],
    [[r['scope'].replace('common_families_', '').replace('all_substantive_pairs', 'All substantive pairs').replace('_', ' to '),
      r['outcomes_compared'], f"{float(r['percent_agreement']):.1f}", f"{float(r['kappa']):.2f}", f"{float(r['percent_agreement_continued']):.1f}"] for r in reading],
    'Continuation agreement contrasts unchanged continuation with all other categories. Ambiguous outcomes are excluded from agreement, '
    'with counts retained in the replication tables. These are language-model readings with reused context, not independent human coding.')
(OUT / 'fragments.json').write_text(json.dumps(fragments, indent=2)+'\n')
(OUT / 'provenance.json').write_text(json.dumps(dict(script='code/19_assemble_tables.py', inputs=SOURCES), indent=2)+'\n')

terms = [r for r in rows('results/text_reuse/term_summary.csv') if r['adoption_year']=='2025' and r['subset']=='all pairs' and r['term']!='biological diversity']
term_names = {'gender':'Gender', 'sex':'Sex', 'diversity':'Diversity', 'climate_change':'Climate change',
              'climate change':'Climate change', 'Sustainable Development Goals':'Sustainable Development Goals', 'Paris Agreement':'Paris Agreement', 'sdgs':'Sustainable Development Goals', 'paris_agreement':'Paris Agreement'}
fragments['terms'] = table('Mentions of six terms in the 2025 successors and their predecessors.', 'tab:contested_terms',
    '>{\\raggedright\\arraybackslash}Xrrrr', ['Term', 'Predecessors naming it', 'Removed', 'Reduced', 'Added'],
    [[term_names.get(r['term'],r['term'].replace('_',' ')),r['predecessor_with_term'],r['dropped'],r['reduced'],r['added']] for r in terms],
    'The population is all 257 pairs with a successor adopted in 2025. Reduced excludes complete removal. '
    'Counts cover full text after running headers are stripped and words across page breaks are joined; '
    'they do not measure the strength of requests. Diversity includes biological diversity; its separate counts are retained in the replication tables.')
fragments['selection'] = table('States selected as the fourteen allies.', 'tab:ally_selection',
    '>{\\raggedright\\arraybackslash}Xrr', ['State', 'Treaty in force', '2024 GDP, billion dollars'],
    [[r['display_name'],r['treaty_in_force_year'],f"{float(r['gdp_2024_imf']):.1f}"] for r in selected],
    'All selected states have a mutual defense treaty with the United States continuously in force from 2001 through 2025 '
    'and a high-income classification in July 2026. GDP uses the October 2026 World Economic Outlook. '
    'The replication tables retain every member, exclusions and alternative criteria.')
trade_diagnostics = rows('results/models/trade/diagnostics.csv')
variance = next(float(r['mean']) for r in trade_diagnostics if r['specification']=='current_baseline' and r['parameter']=='ve')
relations = rows('results/models/trade/relations/us_ally_relationship.csv')
trade_checks = [dict(quantity='residual_sd',value=variance**.5),
    dict(quantity='max_abs_baseline_model_change',value=max(abs(float(r['model_mean'])) for r in relations if r['specification']=='current_baseline'))]
with (OUT/'trade_checks.csv').open('w') as handle:
    writer=csv.DictWriter(handle,fieldnames=['quantity','value']);writer.writeheader();writer.writerows(trade_checks)
(OUT / 'fragments.json').write_text(json.dumps(fragments, indent=2)+'\n')
(OUT / 'provenance.json').write_text(json.dumps(dict(script='code/19_assemble_tables.py', inputs=SOURCES), indent=2)+'\n')

selection_by_actor={r['iso3']:r for r in selected}
fragments['allies_combined'] = table('Selection and observed votes of the fourteen allies in 2025.', 'tab:allies_individual',
    '>{\\raggedright\\arraybackslash}Xrrrrr',
    ['Ally', 'Treaty year', '2024 GDP', 'Votes kept', 'Break probability', 'Move / $\\sigma$'],
    [[selection_by_actor[r['iso3']]['display_name'], selection_by_actor[r['iso3']]['treaty_in_force_year'],
      f"{float(selection_by_actor[r['iso3']]['gdp_2024_imf']):.1f}", r['votes_kept']+' / '+r['switch_pairs'],
      f"{float(r['break_probability_2025']):.4f}",f"{float(r['displacement_sigma_2025']):.1f}"] for r in allies],
    'All states have a treaty continuously in force during the selection window and high income in July 2026. '
    'GDP is in billion dollars. Votes kept concern recurring titles where Washington changed its vote. '
    'Break probabilities and displacement are from the released weighted full-agenda fit; the replication tables retain every selection criterion.')
(OUT/'fragments.json').write_text(json.dumps(fragments,indent=2)+'\n')

(OUT/'provenance.json').write_text(json.dumps(dict(
    script='code/19_assemble_tables.py',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    inputs=SOURCES),indent=2)+'\n')
