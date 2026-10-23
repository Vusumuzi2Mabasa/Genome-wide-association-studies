import sys,unittest,math,tempfile,gzip
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from summary_qc import clean_rows,harmonise,read,FIELDS
from run_workflow import paired,prepare_synthetic_pairs,read_metrics,FIXTURE_ADAPTER

def row(**kw):
    r=dict(zip(FIELDS,['rs1','1','100','A','G','0.1','0.1','0.3173105079','0.3','100','0.95']));r.update(kw);return r

class QC(unittest.TestCase):
    def test_bad_numeric(self):
        for field,value in [('p','0'),('p','NA'),('p','inf'),('se','-1'),('n','2.5'),('chr','23'),('pos','-2'),('beta','NaN'),('eaf','1.1'),('info','.7'),('beta','1_0')]:
            with self.subTest(field=field,value=value):self.assertFalse(clean_rows([row(**{field:value})])[0]['keep'])
    def test_duplicate(self):
        result=clean_rows([row(),row()]);self.assertTrue(result[0]['keep']);self.assertFalse(result[1]['keep'])
    def test_locus_conflict(self):
        self.assertTrue(all(not x['keep'] for x in clean_rows([row(),row(variant='rs2',beta='.3')])))
    def test_id_conflict(self):
        self.assertTrue(all(not x['keep'] for x in clean_rows([row(),row(pos='200')])))
    def test_equivalent_swapped_locus(self):
        r=clean_rows([row(),row(variant='alias',effect_allele='G',other_allele='A',beta='-.1',eaf='.7')]);self.assertEqual(sum(x['keep'] for x in r),1)
    def test_swap(self):
        partner=row(effect_allele='G',other_allele='A',beta='-.1',eaf='.7',build='SIM1')
        x=harmonise(clean_rows([row()]),[partner])[0];self.assertTrue(x['ready_for_combination']);self.assertAlmostEqual(x['aligned_beta'],.1);self.assertAlmostEqual(x['aligned_eaf'],.3)
    def test_malformed_partner(self):
        x=harmonise(clean_rows([row()]),[row(se='0',build='SIM1')])[0];self.assertFalse(x['ready_for_combination']);self.assertEqual(x['aligned_beta'],'NA')
    def test_palindrome(self):
        r=row(effect_allele='A',other_allele='T');x=harmonise(clean_rows([r]),[dict(r,build='SIM1')])[0];self.assertFalse(x['ready_for_combination'])
    def test_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.tsv';p.write_text('a\tb\n1\t2\t3\n')
            with self.assertRaises(ValueError):read(p)
    def test_pairing(self):
        with tempfile.TemporaryDirectory() as d:
            paths=[Path(d)/f'{x}.gz' for x in [1,2]]
            for p in paths:
                with gzip.open(p,'wt') as f:f.write('@x\nACGT\n+\nIIII\n')
            self.assertEqual(paired(*paths),1)
            with gzip.open(paths[1],'wt') as f:f.write('@y\nACGT\n+\nIIII\n')
            with self.assertRaises(ValueError):paired(*paths)
    def test_fixture_repair_and_reordering(self):
        with tempfile.TemporaryDirectory() as d:
            paths=[Path(d)/f'{x}.gz' for x in [1,2]]
            for mate,p in enumerate(paths,1):
                with gzip.open(p,'wt') as f:f.write(f'@sample_good:1:SIM1:1:1001:{2000+mate}\nACGT\n+\nIIII\n')
            repaired,n=prepare_synthetic_pairs(*paths,Path(d)/'out','sample_good')
            self.assertEqual(n,1);self.assertEqual(paired(*repaired),1)
            self.assertEqual(read_metrics(repaired)['q30_fraction'],1)
            with gzip.open(paths[1],'wt') as f:f.write('@sample_good:2:SIM1:1:1002:2002\nACGT\n+\nIIII\n')
            with self.assertRaises(ValueError):prepare_synthetic_pairs(*paths,Path(d)/'bad','sample_good')
    def test_quality_length(self):
        with tempfile.TemporaryDirectory() as d:
            paths=[Path(d)/f'{x}.gz' for x in [1,2]]
            for p in paths:
                with gzip.open(p,'wt') as f:f.write('@x\nACGT\n+\nIII\n')
            with self.assertRaises(ValueError):paired(*paths)

if __name__=='__main__':unittest.main()
