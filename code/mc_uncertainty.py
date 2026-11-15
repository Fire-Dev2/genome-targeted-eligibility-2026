import json, numpy as np
rng=np.random.default_rng(20261114)
M=json.load(open('/home/claude/met.json')); M.append([1,1,1,1])  # idx 34: VHL disease (prevalence 1 by definition)
D={'NHL':19970,'LUNG':124990,'LIVER':30980,'GB':4590,'SOFT':5400,'BREAST':42670,'BRAIN':18350,'CRC':55230,'AML':11500,'ALL':1600,'PANC':52740,'OVARY':12450,'TOT':626140,'KID':15160,'MEL':8510,'THY':2320,'HEME':55830}
BRF=D['TOT']-D['HEME']-D['MEL']-D['LUNG']-D['CRC']; RETP=D['TOT']-D['HEME']-D['LUNG']-D['THY']
# populations: name, year, base fn(params), genie idx, CNS?, orr, orr_n, benefit tags
P=[
("ALK+ ALCL",2021,lambda p:D['NHL']*p['alcl'],25,False,.90,20,"Regular","Single-arm ORR","—"),
("EGFR exon 20 NSCLC",2021,lambda p:D['LUNG']*p['ns'],2,False,.67,153,"Accelerated","Single-arm ORR","Verified 3/2024 (PFS)"),
("KRAS G12C NSCLC",2021,lambda p:D['LUNG']*p['ns'],1,False,.43,112,"Accelerated","Single-arm ORR","Unconfirmed"),
("IDH1 cholangiocarcinoma",2021,lambda p:D['LIVER']*p['icc']+D['GB']*p['ecc'],12,False,.024,126,"Regular","Randomized PFS (OS NS)","—"),
("BRAF V600E other solid tumors",2022,lambda p:BRF,10,False,.41,131,"Accelerated","Single-arm ORR","Unconfirmed"),
("ALK+ IMT",2022,lambda p:D['SOFT']*p['imt'],24,False,.86,14,"Regular","Single-arm ORR","—"),
("HER2-mutant NSCLC",2022,lambda p:D['LUNG']*p['ns'],5,False,.76,72,"Accelerated","Single-arm ORR","Unconfirmed"),
("FGFR1 MLN",2022,lambda p:0.0,20,False,.78,18,"Regular","Single-arm CR","—"),
("RET fusion other solid tumors",2022,lambda p:RETP,11,False,.44,41,"Accelerated","Single-arm ORR","Verified 7/2026 (ORR)"),
("ESR1-mutant breast",2023,lambda p:D['BREAST'],15,False,.19,97,"Regular","Randomized PFS","—"),
("IDH1 MDS",2023,lambda p:p['mds'],13,False,.389,18,"Regular","Single-arm CR","—"),
("AKT1/PTEN breast (no PIK3CA)",2023,lambda p:D['BREAST'],16,False,.26,132,"Regular","Randomized PFS","—"),
("BRAF-fusion low-grade glioma",2024,lambda p:D['BRAIN']*p['lgg'],27,True,.51,76,"Accelerated","Single-arm ORR","Unconfirmed"),
("KRAS G12C colorectal",2024,lambda p:D['CRC'],0,False,.26,53,"Accelerated (withdrawn 9/2026); sotorasib+pani regular 2025","Single-arm ORR; then randomized PFS","Adagrasib withdrawn; still covered"),
("KMT2A AML",2024,lambda p:D['AML'],18,False,.212,104,"Regular","Single-arm CR+CRh","—"),
("KMT2A ALL",2024,lambda p:D['ALL'],19,False,.212,104,"Regular","Single-arm CR+CRh","—"),
("NRG1 NSCLC",2024,lambda p:D['LUNG']*p['ns'],6,False,.33,64,"Accelerated","Single-arm ORR","Unconfirmed"),
("NRG1 pancreatic",2024,lambda p:D['PANC'],7,False,.40,30,"Accelerated","Single-arm ORR","Unconfirmed"),
("KRAS low-grade serous ovarian",2025,lambda p:D['OVARY']*p['lgs'],23,False,.44,57,"Accelerated","Single-arm ORR","Unconfirmed"),
("NPM1 AML",2025,lambda p:D['AML'],17,False,.231,65,"Regular","Single-arm CR+CRh","—"),
("NRG1 cholangiocarcinoma",2026,lambda p:D['LIVER']*p['icc']+D['GB']*p['ecc'],8,False,.368,19,"Regular","Single-arm ORR","—"),
("VHL disease tumors",2021,lambda p:D['KID']*p['vhl'],34,False,.49,61,"Regular","Single-arm ORR","—"),
("HER2-amplified colorectal",2023,lambda p:D['CRC'],28,False,.38,84,"Accelerated","Single-arm ORR","Unconfirmed"),
("HER2-amplified biliary tract",2024,lambda p:D['LIVER']*p['icc']+D['GB'],29,False,.52,62,"Accelerated","Single-arm ORR","Unconfirmed"),
]
point=dict(ns=.85,icc=.30,ecc=.28,alcl=.02,imt=.005,lgg=.01,lgs=.05,mds=10000,vhl=.003)
def draw(n):
    u=lambda c,lo,hi: rng.uniform(c*lo,c*hi,n)
    return dict(ns=rng.uniform(.80,.90,n),icc=rng.uniform(.25,.35,n),ecc=rng.uniform(.20,.36,n),alcl=rng.uniform(.01,.03,n),imt=rng.uniform(.0025,.0075,n),lgg=rng.uniform(.005,.02,n),lgs=rng.uniform(.03,.08,n),mds=rng.uniform(7000,13000,n),vhl=rng.uniform(.0015,.0045,n))
def run(metastatic, n=20000):
    pars=draw(n); E=np.zeros(n); R=np.zeros(n); byyr={}
    pe=pr=0; rows=[]
    for name,yr,base,gi,cns,orr,on,ap,basis,status in P:
        a,d,ma,md=M[gi]; (k,t)=(ma,md) if (metastatic and not cns and md>=30) else (a,d)
        prev_pt=k/t if t else 0; prev=rng.beta(k+.5,t-k+.5,n)
        b=base(pars) if callable(base) else base
        b=np.broadcast_to(np.asarray(b,dtype=float),(n,))
        e=b*prev; E+=e; byyr.setdefault(yr,[np.zeros(n),np.zeros(n)]); byyr[yr][0]+=e
        e_pt=base(point)*prev_pt; pe+=e_pt
        if orr is not None:
            o=rng.beta(orr*on+.5,(1-orr)*on+.5,n); R+=e*o; byyr[yr][1]+=e*o; pr+=e_pt*orr
        rows.append((name,yr,e_pt,None if orr is None else e_pt*orr,prev_pt,k,t,ap,basis,status))
    q=lambda x:(np.percentile(x,2.5),np.percentile(x,97.5))
    out=dict(E=pe,Ec=q(E),R=pr,Rc=q(R),Ep=pe/D['TOT'],Epc=tuple(v/D['TOT'] for v in q(E)),Rp=pr/D['TOT'],Rpc=tuple(v/D['TOT'] for v in q(R)),rows=rows)
    cumE=np.zeros(n); cumR=np.zeros(n); yrs=[]
    for y in sorted(byyr): cumE+=byyr[y][0]; cumR+=byyr[y][1]; yrs.append((y,np.median(cumE)/D['TOT'],*[v/D['TOT'] for v in q(cumE)],np.median(cumR)/D['TOT'],*[v/D['TOT'] for v in q(cumR)]))
    out['years']=yrs; return out
A=run(False); B=run(True)
json.dump({'all':{k:v for k,v in A.items()},'met':{k:v for k,v in B.items()}},open('/home/claude/mc_out.json','w'),default=float)
for lab,o in (("All samples",A),("Metastatic samples",B)):
    print(lab, f"eligible {o['E']:.0f} ({o['Ec'][0]:.0f}-{o['Ec'][1]:.0f}) = {o['Ep']*100:.2f}% ({o['Epc'][0]*100:.2f}-{o['Epc'][1]*100:.2f}); respond {o['R']:.0f} ({o['Rc'][0]:.0f}-{o['Rc'][1]:.0f}) = {o['Rp']*100:.2f}% ({o['Rpc'][0]*100:.2f}-{o['Rpc'][1]*100:.2f})")
    for y in o['years']: print('  ',y[0],' '.join(f'{v*100:.2f}' for v in y[1:]))
# benefit summary
rows=A['rows']; tot=sum(r[2] for r in rows)
import collections
for key,idx in (("approval",7),("basis",8)):
    c=collections.defaultdict(float)
    for r in rows:
        k=r[idx].split(';')[0].split(' (')[0] if key=='approval' else ('Randomized' if r[8].startswith('Randomized') else 'Single-arm')
        c[k]+=r[2]
    print(key,{k:f'{v:.0f} ({v/tot*100:.0f}%)' for k,v in c.items()}, 'n pops', collections.Counter((r[7].split(' (')[0] if key=='approval' else ('Randomized' if r[8].startswith('Randomized') else 'Single-arm')) for r in rows))
