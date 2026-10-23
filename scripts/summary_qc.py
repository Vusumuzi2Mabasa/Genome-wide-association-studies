"""Strict synthetic autosomal summary-statistics QC and allele alignment."""
import csv, math, re
from collections import defaultdict
from pathlib import Path
FIELDS='variant chr pos effect_allele other_allele beta se p eaf n info'.split()
NUMERIC='chr pos beta se p eaf n info'.split()

def read(path, required=()):
    with open(path, newline='') as f:
        reader=csv.DictReader((l for l in f if not l.startswith('##')),delimiter='\t')
        names=reader.fieldnames
        if not names or len(names)!=len(set(names)) or not set(required)<=set(names):
            raise ValueError(f'Invalid header: {path}')
        rows=list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError(f'Ragged table: {path}')
    return rows

def write(path, rows, fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def number(value):
    if not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?',str(value).strip()):return math.nan
    return float(value)

def clean_rows(source):
    rows=[]; reasons=[]; loci=defaultdict(list); ids=defaultdict(list);seen=set()
    for idx,original in enumerate(source):
        r={k:original[k] for k in FIELDS};r['variant']=r['variant'].strip()
        for k in ['effect_allele','other_allele']:r[k]=r[k].strip().upper()
        for k in NUMERIC:r[k]=number(r[k])
        why=[]
        if not r['variant']:why.append('invalid_variant_id')
        for k in NUMERIC:
            if not math.isfinite(r[k]):why.append('invalid_'+k)
        checks={'chr':r['chr'] in range(1,23),'pos':r['pos']>0 and r['pos'].is_integer(),
                'n':r['n']>0 and r['n'].is_integer(),'se':r['se']>0,'p':0<r['p']<=1,
                'eaf':0<r['eaf']<1,'info':.8<=r['info']<=1}
        for k,ok in checks.items():
            if not ok and 'invalid_'+k not in why:why.append('invalid_'+k)
        ea,oa=r['effect_allele'],r['other_allele']
        if ea not in ['A','C','G','T'] or oa not in ['A','C','G','T'] or ea==oa:why.append('invalid_alleles')
        if 0<r['eaf']<1 and min(r['eaf'],1-r['eaf'])<.01:why.append('low_maf')
        signature=tuple(str(r[k]) for k in FIELDS)
        if signature in seen:why.append('exact_duplicate')
        seen.add(signature)
        if not why or why==['exact_duplicate']:
            pair=tuple(sorted((ea,oa))); key=(r['chr'],r['pos'],pair)
            sig=(r['beta'] if ea==pair[0] else -r['beta'],r['se'],r['p'],r['eaf'] if ea==pair[0] else 1-r['eaf'],r['n'],r['info'])
            loci[key].append((idx,sig));ids[r['variant']].append((idx,key,sig))
        r['source_row']=idx+2;r['palindromic']=''.join(sorted((ea,oa))) in ['AT','CG']
        rows.append(r);reasons.append(why)
    def equal(a,b):return all(math.isclose(x,y,rel_tol=1e-8,abs_tol=1e-12) for x,y in zip(a,b))
    for group in loci.values():
        if any(not equal(group[0][1],sig) for _,sig in group):
            for i,_ in group:reasons[i].append('conflicting_locus')
        elif len(group)>1:
            for i,_ in group[1:]:
                if 'exact_duplicate' not in reasons[i]:reasons[i].append('redundant_locus')
    for group in ids.values():
        if any(key!=group[0][1] or not equal(sig,group[0][2]) for _,key,sig in group):
            for i,_,_ in group:reasons[i].append('conflicting_variant_id')
    for r,why in zip(rows,reasons):r.update(exclusion_reason=';'.join(why),keep=not why)
    return rows

def harmonise(primary,partner):
    reference={r['variant']:r for r in primary if r['keep']}
    audit=clean_rows(partner); out=[]; complement=str.maketrans('ACGT','TGCA')
    for original,p in zip(partner,audit):
        issues=p['exclusion_reason'].split(';') if p['exclusion_reason'] else []
        if original.get('build')!='SIM1':issues.append('build_mismatch')
        r=reference.get(p['variant']);orient='not_matched';beta=eaf='NA'
        if r is None:issues.append('not_in_clean_primary')
        elif (p['chr'],p['pos'])!=(r['chr'],r['pos']):issues.append('locus_mismatch')
        elif not issues:
            ea,oa=p['effect_allele'],p['other_allele'];target=(r['effect_allele'],r['other_allele'])
            if p['palindromic']:issues.append('palindromic_unresolved');orient='palindromic_unresolved'
            else:
                for alleles,label,swap in [((ea,oa),'same',False),((oa,ea),'swap',True),((ea.translate(complement),oa.translate(complement)),'strand',False),((oa.translate(complement),ea.translate(complement)),'strand_swap',True)]:
                    if alleles==target:
                        orient=label;beta=-p['beta'] if swap else p['beta'];eaf=1-p['eaf'] if swap else p['eaf'];break
                else:issues.append('allele_mismatch')
        if math.isfinite(p['beta']) and math.isfinite(p['se']) and p['se']>0 and 0<p['p']<=1:
            expected=math.erfc(abs(p['beta']/p['se'])/math.sqrt(2))
            if expected==0 or abs(math.log10(expected)-math.log10(p['p']))>.05:issues.append('wald_p_inconsistent')
        if issues:beta=eaf='NA'
        out.append(dict(variant=p['variant'],orientation=orient,aligned_beta=beta,aligned_eaf=eaf,issues=';'.join(issues),ready_for_combination=not issues))
    return out

def run(data,out):
    audit=clean_rows(read(data/'gwas_summary.tsv',FIELDS));clean=[r for r in audit if r['keep']]
    write(out/'audit.tsv',audit);write(out/'clean.tsv',clean,list(audit[0]))
    hits=sorted((r for r in clean if r['p']<5e-8),key=lambda r:r['p']);write(out/'hits.tsv',hits,list(audit[0]))
    aligned=harmonise(audit,read(data/'gwas_partner_for_harmonisation.tsv',FIELDS+['build']))
    write(out/'harmonisation.tsv',aligned)
    ann=read(data/'annotation.tsv',['variant','gene_label'])
    if len({r['variant'] for r in ann})!=len(ann):raise ValueError('Duplicate annotation key')
    labels={r['variant']:r['gene_label'] for r in ann}
    write(out/'annotated.tsv',[dict(r,gene_label=labels.get(r['variant'],'NA')) for r in clean],list(audit[0])+['gene_label'])
    return dict(input_rows=len(audit),retained=len(clean),excluded=len(audit)-len(clean),hits=len(hits),partner_rows=len(aligned),harmonised=sum(r['ready_for_combination'] for r in aligned),missing_annotation=sum(r['variant'] not in labels for r in clean))
