"""Regenerate Report/overview/sample_table.tex from figures/samples/*/*/*_params.dat.
Run from the repo root. Verdicts are the judgements recorded in PROGRESS.md, step S3."""
import glob, os
verdict = {
 ('bulge','both_007'):'good',('bulge','both_008'):'good',('bulge','both_009'):'good (showcase)',
 ('bulge','gap_filler_010'):'good',('bulge','gap_filler_011'):'weak',('bulge','gap_filler_012'):'good',
 ('bulge','roman_only_003'):'acceptable',('bulge','roman_only_004'):'acceptable',
 ('bulge','any_005'):'acceptable',('bulge','any_006'):'acceptable',
 ('bulge','astrometric_001'):'weak',('bulge','astrometric_002'):'weak',
 ('bulge','astrometric_b001'):'good',('bulge','astrometric_b002'):'good',('bulge','astrometric_b003'):'good (showcase)',
 ('bulge','rubin_only_c001'):'good',('bulge','rubin_only_c002'):'good',
 ('ns','astrometric_001'):'good',('ns','astrometric_002'):'weak',('ns','both_008'):'good',('ns','both_009'):'good',('ns','both_010'):'good',
 ('ns','gap_filler_011'):'good',('ns','gap_filler_012'):'weak',('ns','ns_typical_004'):'weak',('ns','ns_typical_006'):'good',
 ('ns','ns_typical_007'):'weak',('ns','roman_only_003'):'acceptable',('ns','roman_only_005'):'acceptable',('ns','rubin_only_b001'):'good',
 ('bh','both_006'):'good',('bh','both_007'):'good',('bh','both_011'):'good',('bh','astrometric_001'):'good',('bh','astrometric_005'):'good',
 ('bh','bh_long_009'):'acceptable',('bh','bh_long_010'):'acceptable',('bh','roman_only_002'):'acceptable',('bh','roman_only_003'):'acceptable',
 ('bh','bh_short_004'):'weak',('bh','bh_short_008'):'weak',('bh','bh_short_b001'):'weak',('bh','bh_short_b002'):'good',
}
zone={'0':'season','1':'gap','2':'off'}
def pct(x):
    x=float(x)
    if x<0: return 'n/m'
    p=100*x
    if p>=1000: return '$>$1000'
    if p>=10: return f'{p:.0f}'
    return f'{p:.1g}' if p<1 else f'{p:.1f}'
rows=[]
for pop in ['bulge','ns','bh']:
    for f in sorted(glob.glob(f'figures/samples/{pop}/*/*_params.dat')):
        k={}
        for line in open(f):
            if line.startswith('#'): continue
            t=line.split()
            if len(t)>=2:
                for i in range(0,len(t)-1,2): k[t[i]]=t[i+1]
                if t[0] in ('magb','blend'): k[t[0]]=t[1:]
        ev=os.path.basename(f).replace('_params.dat','')
        tE=float(k['tE'])
        det=('L' if k['detL']=='1' else '')+('R' if k['detR']=='1' else '')
        def rel(key): 
            v=float(k[key]); return pct(v/tE) if v>=0 else 'n/m'
        rows.append(f"{pop} & \\texttt{{{ev.replace('_','\\_')}}} & {tE:.1f} & {float(k['u0']):.3f} & {float(k['Ml']):.3g} & {float(k['Dl']):.2f} & {det or '--'} & {zone[k['t0zone']]} & {float(k['max_shift']):.2f} & {rel('sigtE_J')} / {rel('sigtE_L')} / {rel('sigtE_R')} & {pct(k['relMl_J'])} / {pct(k['relMl_R'])} & {verdict[(pop,ev)]} \\\\")
open('Report/overview/sample_table.tex','w').write('\n'.join(rows)+'\n\\bottomrule\n')
print(len(rows))
