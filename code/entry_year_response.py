"""Response rate in the year each population first became eligible, by approval era,
and reconciliation of the 2006 estimate with Marquart et al (2018)."""
import importlib.util,io,contextlib,os
here=os.path.dirname(os.path.abspath(__file__))
src=os.path.join(here,'full_series.py') if os.path.exists(os.path.join(here,'full_series.py')) else '/home/claude/full.py'
spec=importlib.util.spec_from_file_location('full',src); F=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(F)
pt=F.ppt(); T=F.D['TOT']
def cohort(a,b,drop_unrated=False):
    e=r=0
    for name,s,en,c,base,prev,sched,note in F.POP:
        if c!='T' or not (a<=s<=b): continue
        pv=pt[prev[1]] if prev[0]=='ph' else prev[0]/prev[1]; el=float(base(pt)*pv)
        o,_=F.best_orr(sched,s)
        if o is None and drop_unrated: continue
        e+=el; r+=el*(o or 0)
    return r/e
for a,b in ((2006,2014),(2015,2020),(2021,2026)):
    print(f"{a}-{b}: response at entry {cohort(a,b):.1%} (excluding populations with no labeled rate at entry: {cohort(a,b,True):.1%})")
E06=F.Ep['T'][2006][0]
egfr=float(F.D['LUNG']*pt['ns'])*5301/34613
her2=F.D['BREAST']*(0.20-1841/15267)
print(f"2006 eligible: {E06/T:.1%}; + EGFR-mutant NSCLC from 2006: {(E06+egfr)/T:.1%}; + HER2 breast at 20%: {(E06+egfr+her2)/T:.1%}")
