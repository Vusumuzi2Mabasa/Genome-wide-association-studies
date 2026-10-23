#!/usr/bin/env python3
"""Portable Linux/WSL workflow; no network access, installation or job submission."""
import argparse,csv,gzip,hashlib,itertools,json,math,os,shlex,shutil,statistics,subprocess,sys,time,zipfile
from pathlib import Path
from summary_qc import read,write,run as summary_run
ROOT=Path(__file__).resolve().parents[1]

def require(ok,message):
    if not ok:raise ValueError(message)

def fastq(path):
    with gzip.open(path,'rt') as f:
        while True:
            head=f.readline()
            if not head:break
            seq=f.readline().rstrip();plus=f.readline().rstrip();qual=f.readline().rstrip()
            require(head.startswith('@') and plus.startswith('+'),'Malformed FASTQ header')
            require(len(seq)>0 and len(seq)==len(qual),'Truncated FASTQ or length mismatch')
            require(set(seq.upper())<=set('ACGTN') and all(33<=ord(q)<=74 for q in qual),'Invalid DNA/Phred+33 characters')
            name=head[1:].split()[0]
            if name.endswith(('/1','/2')):name=name[:-2]
            yield name,seq,qual

def paired(r1,r2):
    count=0;seen=set()
    for a,b in itertools.zip_longest(fastq(r1),fastq(r2)):
        require(a is not None and b is not None and a[0]==b[0],'Mismatched paired reads')
        require(a[0] not in seen,'Duplicate read identifier');seen.add(a[0]);count+=1
    return count

def prepare_synthetic_pairs(r1,r2,out,sample):
    """Repair only the documented SIM1 fixture convention; never infer real pairs."""
    require(sample in ['sample_good','sample_adapter','sample_lowqual'],'Unknown synthetic fixture')
    out.mkdir(parents=True,exist_ok=True);paths=[out/f'{sample}_R{x}.fastq.gz' for x in [1,2]]
    with gzip.open(paths[0],'wt') as f1,gzip.open(paths[1],'wt') as f2:
        count=0
        for count,(a,b) in enumerate(itertools.zip_longest(fastq(r1),fastq(r2)),1):
            require(a is not None and b is not None,'Unequal synthetic mate counts')
            for mate,r,f in [(1,a,f1),(2,b,f2)]:
                expected=f'{sample}:{count}:SIM1:1:{1000+count}:{2000+mate}'
                require(r[0]==expected,'Unexpected or reordered synthetic identifier; refusing repair')
                f.write(f'@{sample}:{count}/{mate}\n{r[1]}\n+\n{r[2]}\n')
    require(paired(*paths)==count,'Synthetic identifier normalisation failed')
    return paths,count

FIXTURE_ADAPTER='AGATCGGAAGAGCACACGTCTGAACTCCAG'

def read_metrics(paths):
    bases=q30=adapters=low=0
    for path in paths:
        for _,seq,qual in fastq(path):
            bases+=len(seq);q30+=sum(ord(q)-33>=30 for q in qual)
            low+=sum(ord(q)-33<20 for q in qual)
            adapters+=FIXTURE_ADAPTER in seq
    require(bases>0,'Empty read input')
    return dict(bases=bases,q30_fraction=q30/bases,low_quality_bases=low,adapter_reads=adapters)

def reads_branch(data,out,run,tools,threads):
    raw=out/'raw_fastqc';post=out/'processed_fastqc';trim=out/'trimmed'
    for p in [raw,post,trim]:p.mkdir(parents=True)
    samples=sorted(p.name.replace('_R1.fastq.gz','') for p in (data/'fastq').glob('*_R1.fastq.gz'))
    require(samples and len(list((data/'fastq').glob('*.gz')))==2*len(samples),'Unpaired FASTQ filenames')
    prepared={s:prepare_synthetic_pairs(data/'fastq'/f'{s}_R1.fastq.gz',data/'fastq'/f'{s}_R2.fastq.gz',out/'paired_input',s) for s in samples}
    counts={s:prepared[s][1] for s in samples}
    run('fastqc_raw',[tools['fastqc'],'--threads',str(threads),'--outdir',str(raw),*map(str,sorted((data/'fastq').glob('*.gz')))])
    modules=[]
    def collect(folder,stage):
        for archive in folder.glob('*.zip'):
            with zipfile.ZipFile(archive) as z:
                name=next(n for n in z.namelist() if n.endswith('/summary.txt'))
                for line in z.read(name).decode().splitlines():
                    status,module,file=line.split('\t');modules.append(dict(stage=stage,file=file,module=module,status=status))
    collect(raw,'raw');records=[];processed=[]
    for s in samples:
        pair=prepared[s][0]
        before=read_metrics(pair)
        flagged=before['low_quality_bases']>0 or before['adapter_reads']>0
        if flagged:
            new=[trim/f'{s}_R{x}.fastq.gz' for x in [1,2]]
            run('fastp_'+s,[tools['fastp'],'-i',str(pair[0]),'-I',str(pair[1]),'-o',str(new[0]),'-O',str(new[1]),'--adapter_sequence',FIXTURE_ADAPTER,'--adapter_sequence_r2',FIXTURE_ADAPTER,'--cut_tail','--cut_tail_mean_quality','20','--length_required','30','--dont_eval_duplication','--thread',str(threads),'--json',str(out/(s+'_fastp.json')),'--html',str(out/(s+'_fastp.html'))])
            pair=new
        after=paired(*pair);require(0<after<=counts[s],'No usable pairs or increased read count')
        metrics=read_metrics(pair)
        require(metrics['adapter_reads']==0 and metrics['low_quality_bases']==0,'Residual fixture adapters or low-quality bases')
        records.append(dict(sample=s,input_pairs=counts[s],output_pairs=after,synthetic_ids_normalised=True,processed=flagged,**{'before_'+k:v for k,v in before.items()},**{'after_'+k:v for k,v in metrics.items()}));processed+=pair
    run('fastqc_processed',[tools['fastqc'],'--threads',str(threads),'--outdir',str(post),*map(str,processed)])
    collect(post,'processed');write(out/'fastqc_modules.tsv',modules);write(out/'read_retention.tsv',records)
    for stage,folder in [('raw',raw),('processed',post)]:
        run('multiqc_'+stage,[tools['multiqc'],str(folder),'--outdir',str(out/('multiqc_'+stage)),'--filename','multiqc_report.html','--force'])
    return records

def genotype_branch(data,out,run,tool,threads):
    out.mkdir(parents=True); variants=[];samples=None
    with open(data/'genotype/cohort_practice.vcf') as f:
        for line in f:
            if line.startswith('##'):continue
            row=line.rstrip().split('\t')
            if row[0]=='#CHROM':samples=row[9:];continue
            require(samples is not None and len(row)==len(samples)+9,'VCF column mismatch')
            gtindex=row[8].split(':').index('GT');gt=[r.split(':')[gtindex].replace('|','/') for r in row[9:]]
            require(set(gt)<=set(['0/0','0/1','1/0','1/1','./.']),'Unsupported GT')
            variants.append((row[:5],gt))
    require(len(set(samples))==len(samples) and len({v[0][2] for v in variants})==len(variants),'Duplicate sample/variant ID')
    cov=read(data/'genotype/covariates.tsv',['FID','IID','AGE','SEX','POP']);ph=read(data/'genotype/phenotype.tsv',['FID','IID','PHENO'])
    for table in [cov,ph]:
        require(len(table)==len(samples) and len({r['IID'] for r in table})==len(table) and {r['IID'] for r in table}==set(samples),'Metadata mismatch')
        require(all(r['FID']==r['IID'] for r in table),'Expected matching FID and IID')
    cv={r['IID']:r for r in cov};phen={r['IID']:r for r in ph}
    require(all(math.isfinite(float(r['AGE'])) and r['SEX'] in ['1','2'] for r in cov),'Invalid covariates')
    require(all(r['PHENO'] in ['1','2'] for r in ph),'Invalid binary phenotype')
    signatures={};audit=[]
    for i,s in enumerate(samples):
        g=tuple('/'.join(sorted(v[1][i].split('/'))) for v in variants)
        observed=[x for x in g if x!='./.'];het=sum(x=='0/1' for x in observed)/len(observed) if observed else math.nan
        duplicate=signatures.get(g,'');signatures.setdefault(g,s)
        audit.append(dict(IID=s,missing=1-len(observed)/len(g),heterozygosity=het,duplicate_of=duplicate,heterozygosity_outlier=False))
    for group in sorted({r['POP'] for r in cov}):
        eligible=[r for r in audit if cv[r['IID']]['POP']==group and r['missing']<=.1 and not r['duplicate_of']]
        center=statistics.median(r['heterozygosity'] for r in eligible)
        mad=1.4826*statistics.median(abs(r['heterozygosity']-center) for r in eligible)
        for r in eligible:r['heterozygosity_outlier']=abs(r['heterozygosity']-center)>6*max(mad,.005)
    write(out/'sample_audit.tsv',audit)
    remove=[r['IID'] for r in audit if r['duplicate_of'] or r['heterozygosity_outlier']]
    (out/'remove.tsv').write_text(''.join(f'{s}\t{s}\n' for s in remove))
    def plink(name,*args):run('plink_'+name,[tool,'--threads',str(threads),'--memory','1024',*map(str,args),'--out',str(out/name)])
    plink('raw','--vcf',data/'genotype/cohort_practice.vcf','--double-id','--make-pgen')
    plink('missing','--pfile',out/'raw','--missing')
    args=['--pfile',out/'raw','--mind','.10','--make-pgen']
    if remove:args+=['--remove',out/'remove.tsv']
    plink('samples',*args)
    kept=[r['IID'] for r in read(out/'samples.psam')]
    require(set(kept)=={r['IID'] for r in audit if r['IID'] not in remove and r['missing']<=.1+1e-12},'Sample QC mismatch')
    plink('base','--pfile',out/'samples','--geno','.10','--maf','.01','--mac','5','--make-pgen')
    fails=set();hwecounts=[]
    for i,group in enumerate(sorted({r['POP'] for r in cov})):
        controls=[s for s in kept if cv[s]['POP']==group and phen[s]['PHENO']=='1']
        require(len(controls)>=5,'Insufficient controls for within-population HWE report')
        keep=out/f'controls_{i}.tsv';keep.write_text(''.join(f'{s}\t{s}\n' for s in controls))
        plink(f'hwe_{i}','--pfile',out/'base','--keep',keep,'--hardy')
        for r in read(out/f'hwe_{i}.hardy'):
            p=float(r['P']);require(math.isfinite(p),'Nonfinite HWE p')
            if p<1e-6:fails.add(r['ID'])
        hwecounts.append(dict(population=group,controls=len(controls)))
    write(out/'hwe_groups.tsv',hwecounts)
    (out/'hwe_exclude.txt').write_text(''.join(s+'\n' for s in sorted(fails)))
    args=['--pfile',out/'base','--make-pgen']
    if fails:args+=['--exclude',out/'hwe_exclude.txt']
    plink('qc',*args);plink('qc_stats','--pfile',out/'qc','--freq','--missing')
    plink('ld','--pfile',out/'qc','--indep-pairwise','50','5','.2')
    plink('pcs','--pfile',out/'qc','--extract',out/'ld.prune.in','--pca','3')
    pcs=read(out/'pcs.eigenvec');require({r['IID'] for r in pcs}==set(kept),'PC sample mismatch')
    merged=[]
    for r in pcs:
        require(all(math.isfinite(float(r[k])) for k in ['PC1','PC2','PC3']),'Nonfinite PCs')
        merged.append(dict(FID=r['IID'],IID=r['IID'],AGE=cv[r['IID']]['AGE'],SEX=cv[r['IID']]['SEX'],**{k:r[k] for k in ['PC1','PC2','PC3']}))
    write(out/'covariates.tsv',merged)
    plink('assoc','--pfile',out/'qc','--pheno',data/'genotype/phenotype.tsv','--pheno-name','PHENO','--covar',out/'covariates.tsv','--covar-name','AGE,SEX,PC1,PC2,PC3','--covar-variance-standardize','AGE,PC1,PC2,PC3','--glm','hide-covar','firth-fallback')
    result=list(out.glob('assoc.*.glm.*'));require(len(result)==1,'Unexpected association outputs')
    assoc=read(result[0]);add=[r for r in assoc if r['TEST']=='ADD'];pvar=read(out/'qc.pvar')
    require(len(add)==len(pvar),'Missing additive tests')
    valid=[r for r in add if r.get('ERRCODE','.')=='.' and math.isfinite(float(r['P'])) and 0<float(r['P'])<=1]
    write(out/'association_valid.tsv',valid,list(add[0]));write(out/'association_audit.tsv',add)
    return dict(raw_samples=len(samples),raw_variants=len(variants),retained_samples=len(kept),retained_variants=len(pvar),pca_markers=len((out/'ld.prune.in').read_text().splitlines()),cases=sum(phen[s]['PHENO']=='2' for s in kept),controls=sum(phen[s]['PHENO']=='1' for s in kept),hwe_exclusions=len(fails),valid_tests=len(valid),failed_tests=len(add)-len(valid),bonferroni_hits=sum(float(r['P'])<.05/len(pvar) for r in valid),genomewide_hits=sum(float(r['P'])<5e-8 for r in valid),firth_tests=sum(r.get('FIRTH?')=='Y' for r in valid))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['all','summary','reads','genotype'],default='all');ap.add_argument('--out',default='runs/local');ap.add_argument('--threads',type=int,default=2);a=ap.parse_args()
    require(a.threads>0,'threads must be positive');os.chdir(ROOT)
    out=Path(a.out).resolve();require(not out.exists(),'Output already exists; choose a new --out directory')
    names=[] if a.stage=='summary' else (['fastqc','multiqc','fastp'] if a.stage=='reads' else ['plink2'] if a.stage=='genotype' else ['fastqc','multiqc','fastp','plink2'])
    tools={name:shutil.which(os.environ.get(name.upper(),name)) for name in names}
    require(all(tools.values()),'Missing tools: '+','.join(k for k,v in tools.items() if not v))
    out.mkdir(parents=True);logs=out/'logs';logs.mkdir()
    data=ROOT/'data';before={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in data.rglob('*') if p.is_file()}
    commands=[]
    def run(name,cmd):
        commands.append(shlex.join(cmd));(out/'commands.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\n'+'\n'.join(commands)+'\n')
        with open(logs/(name+'.log'),'w') as f:subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT)
    versions={'python':sys.version,'executables':tools}
    for name,exe in tools.items():
        v=subprocess.run([exe,'--version'],text=True,capture_output=True,check=True);versions[name]=(v.stdout+v.stderr).strip()
    (out/'versions.json').write_text(json.dumps(versions,indent=2)+'\n')
    result={'stage':a.stage,'hpc_execution':'not inferred; see scheduler logs','synthetic_only':True}
    if a.stage in ['all','summary']:result['summary']=summary_run(data,out/'summary')
    if a.stage in ['all','reads']:result['reads']=reads_branch(data,out/'reads',run,tools,a.threads)
    if a.stage in ['all','genotype']:result['genotype']=genotype_branch(data,out/'genotype',run,tools['plink2'],a.threads)
    after={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in data.rglob('*') if p.is_file()}
    require(before==after,'Input files changed during run')
    result['inputs_unchanged']=True;result['status']='PASS';result['input_sha256']=before
    (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
