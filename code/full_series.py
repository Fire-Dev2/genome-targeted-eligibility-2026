import json, numpy as np
rng=np.random.default_rng(20261114); N=20000
D=dict(NHL=19970,LUNG=124990,LIVER=30980,GB=4590,SOFT=5400,BREAST=42670,BRAIN=18350,CRC=55230,AML=11500,ALL=1600,PANC=52740,OVARY=12450,TOT=626140,MEL=8510,THY=2320,HEME=55830,CML=1170,ESO=16290,STOM=10740,BLAD=17870,PROS=36320,KID=15160)
SOLID=D['TOT']-D['HEME']; BRF=SOLID-D['MEL']-D['LUNG']-D['CRC']; RETP=SOLID-D['LUNG']-D['THY']
# literature-sourced shares / counts: (point, lo, hi) -- sources in SOURCES dict
PH=dict(ns=(.85,.80,.90),icc=(.30,.25,.35),ecc=(.28,.20,.36),alcl=(.02,.01,.03),imt=(.005,.0025,.0075),
        lgg=(.01,.005,.02),lgg2=(.05,.03,.10),lgs=(.05,.03,.08),mds=(10000,7000,13000),gist=(1000,500,1500),
        phall=(.25,.20,.30),atc=(.40,.30,.50),mtc=(.13,.08,.18),fl=(.11,.08,.14),msi=(.035,.027,.045),
        tmb=(.11,.08,.14),vhl=(.003,.0015,.0045),flt3=(.30,.25,.35))
SOURCES=dict(
 ns="NSCLC = 85% of lung cancer (Haslam 2021 method)",
 icc="CDC WONDER 2023: 9,370 ICC (C22.1) deaths ~30% of liver+IBD deaths (Clin Mol Hepatol 2024, PMC11540356)",
 ecc="Same source: ECC 1,004 of 3,533 ECC+gallbladder+ampulla deaths (28%)",
 alcl="PTCL = 7% of NHL incidence-based mortality (Howlader, CEBP 2016); ALCL ~1/3 of PTCL -> ~2% (assumption)",
 imt="Assumption (IMT extremely rare); wide range",
 lgg="Assumption: pediatric-type low-grade glioma deaths ~1% of brain cancer deaths",
 lgg2="Assumption: grade 2 IDH-mutant diffuse glioma ~5% (3-10%) of brain cancer deaths; derived from CBTRUS incidence",
 lgs="LGSOC ~5% of ovarian carcinomas (Wohlmuth, Cancers 2022; WHO 2020)",
 mds="~10,000 incident MDS/yr in US (Ma, Cancer 2007; ACS), poor survival -> deaths assumed ~10,000 (7,000-13,000)",
 gist="Assumption from 4,000-6,000 incident GIST/yr (ACS) with ~20% disease-specific mortality",
 phall="Ph+ = ~25% of adult ALL (LLS; Yilmaz, Clin Adv Hematol Oncol 2018)",
 atc="ATC ~40% of thyroid cancer deaths (Yamazaki, Eur Thyroid J 2025)",
 mtc="MTC ~13% of thyroid cancer deaths (Randle, Surgery 2017 PMC5164945)",
 fl="FL = 11% of NHL incidence-based mortality, 2011 (Howlader, CEBP 2016)",
 msi="MSI-H 3.8% across 39 cancer types (Bonneville, JCO PO 2017); pooled 2.7% (scoping review 2022)",
 tmb="TMB>=10 pooled 14.0% pan-cancer (PMC9705554) minus MSI overlap -> ~11%; KEYNOTE-158 13%",
 vhl="Assumption", flt3="FLT3 mutations ~30% of AML (Daver, Blood Cancer J 2021)")
def pdraw(): return {k:rng.uniform(lo,hi,N) for k,(p,lo,hi) in PH.items()}
def ppt(): return {k:p for k,(p,lo,hi) in PH.items()}
CCA=lambda p:D['LIVER']*p['icc']+D['GB']*p['ecc']; NS=lambda p:D['LUNG']*p['ns']
# ORR schedule: list of (drug_start, drug_end, orr, n, source) -- current FDA label values;
# population response in year y = highest ORR among drugs on the market for it in y (Haslam 2021 rule).
def O(*a): return list(a)
# (name, start, end, class, base(p), prev:(k,n) or ('ph',key), orr_schedule, note)
POP=[
("BCR-ABL CML",2006,None,"T",lambda p:D['CML'],(1,1),O((2006,None,.747,553,"imatinib CCyR IRIS"),(2010,None,.80,282,"nilotinib CCyR by 12 mo ENESTnd")),""),
("Ph+ ALL",2006,None,"T",lambda p:D['ALL'],('ph','phall'),O((2006,None,.31,43,"imatinib MaHR R/R Ph+ ALL"),(2012,None,.41,32,"ponatinib MaHR PACE")),"published prevalence"),
("KIT-mutant GIST",2006,None,"T",lambda p:p['gist'],(1550,2155),O((2006,None,.514,818,"imatinib 400 mg CR+PR")),""),
("HER2-amplified breast",2006,None,"T",lambda p:D['BREAST'],(1841,15267),O((2006,None,.45,235,"trastuzumab+chemo H0648g"),(2012,None,.802,343,"pertuzumab CLEOPATRA"),(2019,None,.603,184,"T-DXd DESTINY-Breast01"),(2025,None,.87,374,"T-DXd+pertuzumab DESTINY-Breast09")),""),
("HER2-amplified esophagogastric",2010,None,"T",lambda p:D['ESO']+D['STOM'],(652,5513),O((2010,None,.47,294,"trastuzumab ToGA (Bang 2010)"),(2021,None,.73,298,"pembrolizumab+trastuzumab KEYNOTE-811")),""),
("BRAF V600 melanoma",2011,None,"T",lambda p:D['MEL'],(2668,9134),O((2011,None,.484,219,"vemurafenib BRIM-3"),(2014,None,.66,211,"dabrafenib+trametinib COMBI-d")),""),
("ALK NSCLC",2011,None,"T",NS,(600,25209),O((2011,None,.74,172,"crizotinib PROFILE 1014"),(2017,None,.79,152,"alectinib ALEX")),""),
("EGFR common NSCLC",2013,None,"T",NS,(5301,34613),O((2013,None,.65,86,"erlotinib EURTAC"),(2018,None,.77,279,"osimertinib FLAURA"),(2024,None,.78,429,"amivantamab+lazertinib MARIPOSA")),""),
("ROS1 NSCLC",2016,None,"T",NS,(251,25209),O((2016,None,.66,50,"crizotinib IRR"),(2019,None,.74,92,"entrectinib"),(2023,None,.79,71,"repotrectinib TKI-naive"),(2025,None,.90,103,"taletrectinib TRUST-I")),""),
("BRAF V600E NSCLC",2017,None,"T",NS,(525,34617),O((2017,None,.61,57,"dabrafenib+trametinib"),(2023,None,.75,59,"encorafenib+binimetinib PHAROS naive")),""),
("IDH2 AML",2017,None,"T",lambda p:D['AML'],(458,4424),O((2017,None,.23,199,"enasidenib CR/CRh")),""),
("FLT3 AML",2017,None,"T",lambda p:D['AML'],('ph','flt3'),O((2018,None,.226,247,"gilteritinib CR/CRh ADMIRAL"),(2023,None,.55,268,"quizartinib+chemo CR QuANTUM-First (equal to placebo arm)")),"midostaurin (2017) OS-based, no response rate"),
("IDH1 AML",2018,None,"T",lambda p:D['AML'],(330,4670),O((2018,None,.328,174,"ivosidenib R/R CR/CRh"),(2019,None,.429,28,"ivosidenib ND unfit CR/CRh")),""),
("BRAF V600E anaplastic thyroid",2018,None,"T",lambda p:D['THY']*p['atc'],(152,379),O((2018,None,.53,36,"dabrafenib+trametinib")),""),
("NTRK fusion solid tumors",2018,None,"T",lambda p:SOLID,(484,172014),O((2018,None,.60,339,"larotrectinib current label")),""),
("EGFR uncommon NSCLC",2018,None,"T",NS,(722,34613),O((2018,None,.656,32,"afatinib 21/32")),""),
("FGFR3 urothelial",2019,None,"T",lambda p:D['BLAD'],(1230,6411),O((2019,None,.322,87,"erdafitinib BLC2001"),(2024,None,.353,136,"erdafitinib THOR")),""),
("PIK3CA breast",2019,None,"T",lambda p:D['BREAST'],(5612,19620),O((2019,None,.357,169,"alpelisib SOLAR-1"),(2024,None,.58,161,"inavolisib INAVO120")),""),
("FGFR2 cholangiocarcinoma",2020,None,"T",CCA,(194,2580),O((2020,None,.36,107,"pemigatinib"),(2022,None,.42,103,"futibatinib")),""),
("MET exon 14 NSCLC",2020,None,"T",NS,(692,34538),O((2020,None,.68,60,"capmatinib treatment-naive")),""),
("RET fusion NSCLC",2020,None,"T",NS,(269,25209),O((2020,None,.84,69,"selpercatinib treatment-naive")),""),
("RET fusion thyroid",2020,None,"T",lambda p:D['THY']*(1-p['atc']-p['mtc']),(147,2700),O((2020,None,.79,19,"selpercatinib")),""),
("RET-mutant medullary thyroid",2020,None,"T",lambda p:D['THY']*p['mtc'],(268,378),O((2020,None,.69,55,"selpercatinib")),""),
("BRAF V600E colorectal",2020,None,"T",lambda p:D['CRC'],(1934,22771),O((2020,None,.20,220,"encorafenib+cetuximab BEACON"),(2024,None,.61,110,"encorafenib+cetuximab+mFOLFOX6 BREAKWATER")),""),
("PDGFRA exon 18 GIST",2020,None,"T",lambda p:p['gist'],(156,1950),O((2020,None,.84,43,"avapritinib")),""),
("EZH2 follicular lymphoma",2020,2026,"T",lambda p:D['NHL']*p['fl'],(209,1058),O((2020,2026,.69,45,"tazemetostat")),"tazemetostat withdrawn 6/2026"),
("KRAS/NRAS wild-type colorectal",2009,None,"I",lambda p:D['CRC'],(12187,22775),O((2009,None,.10,231,"panitumumab monotherapy")),""),
("BRCA ovarian",2014,None,"I",lambda p:D['OVARY'],(846,9119),O((2014,None,.34,137,"olaparib")),""),
("MSI-H/dMMR solid tumors",2017,None,"I",lambda p:SOLID,('ph','msi'),O((2017,None,.40,149,"pembrolizumab")),""),
("BRCA breast",2018,None,"I",lambda p:D['BREAST'],(748,18258),O((2018,None,.52,205,"olaparib OlympiAD")),""),
("BRCA pancreatic",2019,None,"I",lambda p:D['PANC'],(282,11071),O((2019,None,.23,78,"olaparib POLO")),""),
("HRR prostate",2020,None,"I",lambda p:D['PROS'],(1193,8506),O((2020,None,.33,162,"olaparib/rucaparib")),""),
("TMB-H solid tumors (non-MSI)",2020,None,"I",lambda p:SOLID,('ph','tmb'),O((2020,None,.29,102,"pembrolizumab KEYNOTE-158")),""),
("ALK+ ALCL",2021,None,"T",lambda p:D['NHL']*p['alcl'],(21,101),O((2021,None,.90,20,"crizotinib")),""),
("EGFR exon 20 NSCLC",2021,None,"T",NS,(585,34613),O((2021,None,.40,81,"amivantamab CHRYSALIS"),(2024,None,.67,153,"amivantamab+chemo PAPILLON")),""),
("KRAS G12C NSCLC",2021,None,"T",NS,(4032,34619),O((2021,None,.36,124,"sotorasib"),(2022,None,.43,112,"adagrasib")),""),
("IDH1 cholangiocarcinoma",2021,None,"T",CCA,(443,3264),O((2021,None,.024,126,"ivosidenib ClarIDHy (PFS-based)")),"PFS-based"),
("VHL disease tumors",2021,None,"T",lambda p:D['KID']*p['vhl'],(1,1),O((2021,None,.49,61,"belzutifan")),""),
("BRAF V600E other solid tumors",2022,None,"T",lambda p:BRF,(2826,148901),O((2022,None,.41,131,"dabrafenib+trametinib")),""),
("ALK+ IMT",2022,None,"T",lambda p:D['SOFT']*p['imt'],(24,55),O((2022,None,.86,14,"crizotinib")),""),
("HER2-mutant NSCLC",2022,None,"T",NS,(746,34523),O((2022,None,.58,52,"T-DXd DESTINY-Lung02"),(2025,None,.75,71,"zongertinib Beamion LUNG-1 pretreated"),(2026,None,.76,72,"zongertinib Beamion LUNG-1 untreated")),""),
("RET fusion other solid tumors",2022,None,"T",lambda p:RETP,(112,144970),O((2022,None,.44,41,"selpercatinib")),""),
("ESR1 breast",2023,None,"T",lambda p:D['BREAST'],(1338,16262),O((2025,None,.14,112,"imlunestrant EMBER-3 (PFS-based)"),(2026,None,.19,97,"vepdegestrant VERITAC-2 (PFS-based)")),"PFS-based"),
("IDH1 MDS",2023,None,"T",lambda p:p['mds'],(63,2371),O((2023,None,.389,18,"ivosidenib")),""),
("AKT1/PTEN breast (no PIK3CA)",2023,None,"T",lambda p:D['BREAST'],(1650,19593),O((2023,None,.26,132,"capivasertib CAPItello-291 (PFS-based)")),"PFS-based"),
("HER2-amplified colorectal",2023,None,"T",lambda p:D['CRC'],(343,14571),O((2023,None,.38,84,"tucatinib+trastuzumab")),"decision A"),
("BRAF-fusion low-grade glioma",2024,None,"T",lambda p:D['BRAIN']*p['lgg'],(569,1885),O((2024,None,.51,76,"tovorafenib")),""),
("KRAS G12C colorectal",2024,None,"T",lambda p:D['CRC'],(699,22775),O((2024,2026,.34,94,"adagrasib+cetuximab (withdrawn 9/2026)"),(2025,None,.26,53,"sotorasib+panitumumab (rate used from 2026, after adagrasib withdrawal)")),""),
("KMT2A AML",2024,None,"T",lambda p:D['AML'],(67,2085),O((2024,None,.212,104,"revumenib CR/CRh")),""),
("KMT2A ALL",2024,None,"T",lambda p:D['ALL'],(24,955),O((2024,None,.212,104,"revumenib CR/CRh")),""),
("NRG1 NSCLC",2024,None,"T",NS,(40,10389),O((2024,None,.33,64,"zenocutuzumab")),""),
("NRG1 pancreatic",2024,None,"T",lambda p:D['PANC'],(10,3766),O((2024,None,.40,30,"zenocutuzumab")),""),
("HER2-amplified biliary tract",2024,None,"T",lambda p:CCA(p)+D['GB']*(1-p['ecc']),(113,2753),O((2024,None,.52,62,"zanidatamab")),"decision A"),
("KRAS low-grade serous ovarian",2025,None,"T",lambda p:D['OVARY']*p['lgs'],(137,514),O((2025,None,.44,57,"avutometinib+defactinib")),""),
("H3 K27M diffuse midline glioma",2025,None,"I",lambda p:D['BRAIN'],(470,14525),O((2025,None,.22,50,"dordaviprone")),""),
("NPM1 AML",2025,None,"T",lambda p:D['AML'],(584,4669),O((2025,None,.231,65,"revumenib CR/CRh")),""),
("NRG1 cholangiocarcinoma",2026,None,"T",CCA,(1,962),O((2026,None,.368,19,"zenocutuzumab")),""),
("FGFR1 MLN",2022,None,"T",lambda p:0.0,(1,1),O((2022,None,.78,18,"pemigatinib CR")),"deaths not estimable; contributes 0 patients"),
]
yrs=list(range(2006,2027))
def simulate(draw=True):
    p=pdraw() if draw else ppt()
    n=N if draw else 1
    E={c:{y:np.zeros(n) for y in yrs} for c in ("T","I")}; R={c:{y:np.zeros(n) for y in yrs} for c in ("T","I")}
    for name,s,e,c,base,prev,sched,note in POP:
        b=np.broadcast_to(np.asarray(base(p),dtype=float),(n,))
        if prev[0]=='ph': pv=np.broadcast_to(np.asarray(p[prev[1]],float),(n,))
        else:
            k,t=prev; pv=rng.beta(k+.5,t-k+.5,n) if (draw and t>1) else np.full(n,k/t)
        el=b*pv
        od=[(ds,de,(rng.beta(o*on+.5,(1-o)*on+.5,n) if draw else np.full(n,o))) for ds,de,o,on,_ in sched]
        for y in yrs:
            if y>=s and (e is None or y<e):
                E[c][y]+=el
                act=[x for ds,de,x in od if y>=ds and (de is None or y<de)]
                if act: R[c][y]+=el*np.max(np.vstack(act),axis=0)
    return E,R
def best_orr(sched,y):
    act=[(o,src) for ds,de,o,on,src in sched if y>=ds and (de is None or y<de)]
    return max(act) if act else (None,None)
Ep,Rp=simulate(False); Ed,Rd=simulate(True)
q=lambda x:(np.percentile(x,2.5),np.percentile(x,97.5))
out=[]
for y in yrs:
    row=dict(year=y)
    for lab,c in (("T",("T",)),("TI",("T","I"))):
        ep=sum(Ep[k][y][0] for k in c); rp=sum(Rp[k][y][0] for k in c)
        ed=sum(Ed[k][y] for k in c); rd=sum(Rd[k][y] for k in c)
        row[lab]=dict(E=ep/D['TOT'],Elo=q(ed)[0]/D['TOT'],Ehi=q(ed)[1]/D['TOT'],R=rp/D['TOT'],Rlo=q(rd)[0]/D['TOT'],Rhi=q(rd)[1]/D['TOT'],En=ep,Rn=rp)
    out.append(row)
json.dump(out,open('/home/claude/full_out.json','w'))
# population table (point, 2026 and 2020)
pt=ppt(); tab=[]
for name,s,e,c,base,prev,sched,note in POP:
    pv=pt[prev[1]] if prev[0]=='ph' else prev[0]/prev[1]
    el=float(base(pt)*pv); yy=2026 if (e is None or e>2026) else e-1
    o,srcn=best_orr(sched,yy)
    tab.append(dict(name=name,start=s,end=e,cls=c,eligible=round(el),orr=o,orr_src=srcn,responders=None if o is None else round(el*o),prev=pv,note=note,
                    sched=[list(x) for x in sched]))
json.dump(tab,open('/home/claude/full_pops.json','w'))
json.dump({k:[*PH[k],SOURCES.get(k,'')] for k in PH},open('/home/claude/full_params.json','w'))
print("Year | Targeted elig% (CI) | resp% | Targeted+informed elig% | resp%")
for r in out:
    T,TI=r['T'],r['TI']; print(r['year'],f"{T['E']*100:.2f} ({T['Elo']*100:.2f}-{T['Ehi']*100:.2f}) | {T['R']*100:.2f} ({T['Rlo']*100:.2f}-{T['Rhi']*100:.2f}) | {TI['E']*100:.2f} ({TI['Elo']*100:.2f}-{TI['Ehi']*100:.2f}) | {TI['R']*100:.2f}")
