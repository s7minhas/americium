"""Check reported numbers against analysis outputs and manuscript passages."""
import csv
import hashlib
import json
import math
import operator
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/revision'
OUT.mkdir(parents=True, exist_ok=True)
CACHE = {}

def records(path):
    if path not in CACHE:
        with (ROOT/path).open() as handle:
            CACHE[path] = list(csv.DictReader(handle))
    return CACHE[path]

def scalar(value):
    if value in ['TRUE','True','true']:
        return True
    if value in ['FALSE','False','false']:
        return False
    try:
        return float(value)
    except (TypeError, ValueError):
        return value

def equal(a,b):
    a,b=scalar(a),scalar(b)
    return a==b

def matches(row,selector):
    for column,rule in selector.items():
        value=scalar(row.get(column,''))
        if isinstance(rule,dict):
            for op,target in rule.items():
                if op=='in':
                    if not any(equal(value,x) for x in target):return False
                elif op=='not_in':
                    if any(equal(value,x) for x in target):return False
                elif op in ['lt','le','gt','ge','eq','ne']:
                    if not {'lt':operator.lt,'le':operator.le,'gt':operator.gt,'ge':operator.ge,'eq':operator.eq,'ne':operator.ne}[op](value,scalar(target)):return False
                else:raise ValueError(f'Unknown selector operation {op}')
        elif not equal(value,rule):return False
    return True

def evaluate(source,selector,column,operation,scale,digits,expected):
    selected=[r for r in records(source) if matches(r,selector)]
    if not selected and operation!='count':raise ValueError('Empty selection')
    if operation=='count':actual=len(selected)
    elif operation=='distinct_count':actual=len({r[column] for r in selected})
    elif operation=='sum_ratio':
        numerator,denominator=column.split('|')
        actual=sum(float(r[numerator]) for r in selected)/sum(float(r[denominator]) for r in selected)
    else:
        values=[scalar(r[column]) for r in selected]
        if operation.startswith('all_'):
            compare={'all_eq':operator.eq,'all_lt':operator.lt,'all_le':operator.le,'all_gt':operator.gt,'all_ge':operator.ge}[operation]
            return all(compare(v if isinstance(v,str) else v*scale,scalar(expected)) for v in values),str(values[:4])
        if operation=='one':
            if len(values)!=1:raise ValueError(f'Expected one row, found {len(values)}')
            actual=values[0]
        elif operation=='sum':actual=sum(values)
        elif operation=='mean':actual=sum(values)/len(values)
        elif operation=='min':actual=min(values)
        elif operation=='max':actual=max(values)
        else:raise ValueError(f'Unknown operation {operation}')
    if isinstance(actual,str):
        return actual==expected,actual
    actual=float(actual)*scale
    rounded=round(actual,digits) if digits is not None else actual
    target=round(float(expected),digits) if digits is not None else float(expected)
    return math.isclose(rounded,target,abs_tol=1e-9,rel_tol=1e-9),actual

checks=[]
for claim in records('results/checks/claims.csv'):
    if claim['kind']!='numeric':continue
    try:
        passed,actual=evaluate(claim['source'],json.loads(claim['selector'] or '{}'),claim['column'],claim['operation'],float(claim['scale'] or 1),int(claim['digits']) if claim['digits'] else None,claim['expected'])
        error=''
    except Exception as exc:
        passed=False;actual='';error=str(exc)
    checks.append(dict(id=claim['id'],claim=claim['claim'],source=claim['source'],passed=passed,actual=actual,expected=claim['expected'],error=error))
with (OUT/'retained_claim_checks.csv').open('w') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(checks[0]));writer.writeheader();writer.writerows(checks)

# each binding records a displayed value, its exact source cell and the revised passage
new=[]
def bind(description,source,selector,column,document,anchor,digits=1,scale=1,operation='one'):
    selected=[r for r in records(source) if matches(r,selector)]
    if operation=='one':
        assert len(selected)==1,(description,len(selected))
        actual=float(selected[0][column])*scale
    elif operation=='min':actual=min(float(r[column]) for r in selected)*scale
    elif operation=='max':actual=max(float(r[column]) for r in selected)*scale
    else:raise ValueError(operation)
    value=f'{actual:.{digits}f}'
    text=(ROOT/document).read_text()
    passages=[(i,line) for i,line in enumerate(text.splitlines(),1) if anchor in line]
    assert passages,(description,'missing passage',anchor)
    pattern=r'(?<![\d.])'+re.escape(value)+r'(?!\d|\.\d)'
    # Main-text prose may spell out small counts, e.g. 'Eight members'.
    words='zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eighteen nineteen twenty'.split()
    if digits==0 and actual==int(actual) and 0<=actual<len(words):
        pattern='(?:'+pattern+r'|\b'+words[int(actual)]+r'\b)'
    matching=[(i,line) for i,line in passages if re.search(pattern,line.replace(',',''),re.IGNORECASE)]
    line,passage=(matching or passages)[0]
    number_present=bool(matching)
    new.append(dict(claim=description,value=value,document=document,line=line,passage=passage,source=source,selector=json.dumps(selector,sort_keys=True),column=column,operation=operation,scale=scale,digits=digits,script=producer(source),number_present=number_present))

def producer(source):
    rules=[('results/revision/','19_assemble_tables.py'),('un_variants/threshold','13_model_threshold.R'),('un_variants/decomposition','13_model_threshold.R'),('un_variants/ideal_points','15_ideal_points.R'),('un_variants/','11_model_variants.R'),('followers/','16_followers.R'),('ally_selection/','04_ally_selection.R'),('reading_agreement/','08_reading_agreement.R'),('text_reuse/','09_measure_text_reuse.py'),('trade/','03_economic.R')]
    return next(('code/'+script for pattern,script in rules if pattern in source),'code/01_observed_votes.R')

main='paper/main.tex';si='paper/si.tex'
bind('Recurring-title American probability','results/models/un_variants/recurring_both_years/terminal_probabilities.csv',{'actor':'USA'},'probability',main,'A fit of the last two years',5)
bind('Equal-weight allied probabilities above half','results/revision/variant_comparison.csv',{'variant':'equal_weights'},'allies_above_half',si,'Equal resolution weights',0)
for parameter,anchor,digits in [('sigma_uv','The estimated step size of the drift',3),('rho_uv','The estimated step size of the drift',5),('pi_snap','The estimated step size of the drift',3)]:
    bind('Released '+parameter,'results/models/un/diagnostics.csv',{'specification':'full','parameter':parameter,'family':'continuous'},'mean',main,anchor,digits)
for name,group,col,anchor,digits,scale in [
 ('Other-member unchanged','all','unchanged','Of the members outside the fourteen',0,1),
 ('Other-member paired','all','switch_pairs','Of the members outside the fourteen',0,1),
 ('Other-member unchanged share','all','unchanged_share_of_paired','Of the members outside the fourteen',1,100),
 ('Other-member moves toward','all','toward','Of the members outside the fourteen',0,1),
 ('Other-member moves away','all','away','Of the members outside the fourteen',0,1)]:
    bind(name,'results/observed/followers/response_by_group.csv',{'grouping':group},col,main,anchor,digits,scale)
for col,anchor in [('yes','The Second Committee provides'),('abstain','The Second Committee provides'),('no','The Second Committee provides'),('all_fourteen_yes','The Second Committee provides'),('resolutions','The Second Committee provides')]:
    bind('Second Committee '+col,'results/observed/followers/second_committee_totals.csv',{'cohort_year':2025},col,main,anchor,0)
for year in [2017,2024,2025]:
    bind('Verbatim continuity '+str(year),'results/text_reuse/verbatim_summary.csv',{'adoption_year':year,'subset':'all pairs'},'verbatim_any',si,'The exact-match check found',1)
for col in ['within_sector_pp','sector_composition_pp','change_pp']:
    # the paper reports absolute magnitudes of the share losses
    bind('Product decomposition '+col,'results/trade/census/detail/component_summary.csv',{'group':'fourteen','end_year':2025,'classification':'five_digit','horizon':1},col,si,'At the five-digit product level',2,-1)
# Observed trade supports the Discussion; these are descriptive flows, not model results.
trade_selector=dict(vintage='refreshed',support='paired_available_1989',omitted_partner='none',comparison='annual',horizon=1,baseline=2024,target=2025,direction='exports',group='allied_peers_without_us',weighting='value_weighted')
for column,digits,scale in [('share_before',0,100),('share_after',0,100),('allied_flow_growth_percent',2,1),('countries_share_declines',0,1)]:
    bind('Discussion observed trade '+column,'results/trade/all_economies/political_trade_terminal.csv',trade_selector,column,main,'The observed trade figures',digits,scale)
for column in ['share_before','share_after']:
    bind('SI observed trade '+column,'results/trade/all_economies/political_trade_terminal.csv',trade_selector,column,si,'The fourteen allies’ share',1,100)
for source,anchor in [('results/observed/followers/ally_table.csv','BEL &'),('results/observed/ally_selection/selection.csv','Germany &')]:
    if 'ally_table' in source:
        for row in records(source):
            for col,dig in [('switch_pairs',0),('votes_kept',0),('break_probability_2025',4),('displacement_sigma_2025',1)]:
                bind(row['iso3']+' '+col,source,{'iso3':row['iso3']},col,si,next(r['display_name'] for r in records('results/observed/ally_selection/selection.csv') if r['iso3']==row['iso3'])+' &',dig)
    else:
        for row in records(source):
            if row['selected_ally']=='TRUE':
                for col,dig in [('treaty_in_force_year',0),('gdp_2024_imf',1)]:
                    bind(row['iso3']+' '+col,source,{'iso3':row['iso3']},col,si,row['display_name']+' &',dig)

# direct source-cell bindings for the reported specification table
for variant in records('results/revision/variant_comparison.csv'):
    if variant['diagnostics_pass']!='True':continue
    text=(ROOT/si).read_text()
    if variant['label']+' &' not in text:continue
    for col,dig in [('us_probability',5),('max_allied_probability',5),('allies_above_half',0)]:
        bind(variant['label']+' '+col,'results/revision/variant_comparison.csv',{'variant':variant['variant']},col,si,variant['label']+' &',dig)

bind('Joint American first break','results/revision/joint_probability.csv',{},'probability',main,'The joint posterior probability',3)
for col in ['change_2025','prior_extreme']:
    bind('Yes-no-only '+col,'results/observed/allies/historical_ranks.csv',{'specification':'yes_no_only','horizon':1,'group':'us_allies','history_sample':'all_supported'},col,main,'Excluding abstentions gives',1,-100)

for row in records('results/text_reuse/term_summary.csv'):
    if row['adoption_year']!='2025' or row['subset']!='all pairs' or row['term']=='biological diversity':continue
    label=row['term'][0].upper()+row['term'][1:]
    for col in ['predecessor_with_term','dropped','reduced','added']:
        bind('Term '+row['term']+' '+col,'results/text_reuse/term_summary.csv',{'adoption_year':2025,'subset':'all pairs','term':row['term']},col,si,label+' &',0)
for col,dig in [('threshold_units',4),('length_square',4),('threshold_sigma_units',2)]:
    bind('Threshold '+col,'results/models/un_variants/threshold/threshold.csv',{'kappa':2,'length_reference':'usa_2025_mean'},col,si,'The displacement threshold',dig)
for col in ['adjusted_mean','adjusted_lower','adjusted_upper']:
    bind('American move '+col,'results/models/un_variants/threshold/displacements.csv',{'actor':'USA'},col,si,'The displacement threshold',4)
for operation in ['min','max']:
    bind('Allied move '+operation,'results/models/un_variants/threshold/displacements.csv',{'actor':{'ne':'USA'}},'adjusted_mean',si,'The displacement threshold',4,operation=operation)

with (OUT/'new_numbers.csv').open('w') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(new[0]));writer.writeheader();writer.writerows(new)
summary=dict(manuscripts={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in [main,si]},script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),retained_numeric_checks=len(checks),retained_passes=sum(r['passed'] for r in checks),retained_failures=[r for r in checks if not r['passed']],new_bindings=len(new),unmatched_new_bindings=[r['claim'] for r in new if not r['number_present']])
(OUT/'claim_check_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

assert summary["retained_passes"] == summary["retained_numeric_checks"], "Retained numeric claim mismatch"
assert not summary["unmatched_new_bindings"], "An added number is missing from its bound manuscript passage"
