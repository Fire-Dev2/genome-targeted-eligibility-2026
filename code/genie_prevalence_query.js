// GENIE biomarker prevalence for the genome-targeted eligibility update.
// Run in the browser console while signed in at https://genie.cbioportal.org
// Dataset: AACR Project GENIE v20.0-public (study id: genie_public). Queried 2026-11-14.
// Output: patient-level counts (altered / profiled) per biomarker and cancer cohort.
// GENIE terms: do not redistribute the underlying data; share this code and the summary counts only.

const STUDY = 'genie_public';
const MU = `${STUDY}_mutations`, SVP = `${STUDY}_structural_variants`, CN = `${STUDY}_cna`;
const post = async (url, body) => (await fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body)})).json();
const get = async url => (await fetch(url)).json();

// 1. Sample clinical data
const clin = async a => get(`/api/studies/${STUDY}/clinical-data?clinicalDataType=SAMPLE&attributeId=${a}&projection=SUMMARY&pageSize=10000000`);
const S = {};
for (const r of await clin('CANCER_TYPE')) S[r.sampleId] = {p: r.patientId, ct: r.value};
for (const r of await clin('ONCOTREE_CODE')) if (S[r.sampleId]) S[r.sampleId].ot = r.value;

// 2. Gene panels and which samples were profiled in each molecular profile
const GP = {};
for (const p of await get('/api/gene-panels?projection=DETAILED&pageSize=100000')) GP[p.genePanelId] = new Set(p.genes.map(g => g.hugoGeneSymbol));
const ids = Object.keys(S), PROF = {};
for (const prof of [MU, SVP, CN]) {
  const m = new Map();
  for (let i = 0; i < ids.length; i += 40000)
    for (const x of await post(`/api/molecular-profiles/${prof}/gene-panel-data/fetch`, {sampleIds: ids.slice(i, i + 40000)}))
      if (x.profiled) m.set(x.sampleId, x.genePanelId);
  PROF[prof] = m;
}
// A sample is profiled for a gene if its panel includes the gene (whole-genome assays without a panel list count as profiled)
const profiled = (s, genes, p) => { const pid = PROF[p].get(s); if (!pid) return false; const set = GP[pid]; return set ? genes.some(g => set.has(g)) : true; };

// 3. Alterations
const GENES = {KRAS:3845,EGFR:1956,ERBB2:2064,BRAF:673,IDH1:3417,IDH2:3418,ESR1:2099,AKT1:207,PTEN:5728,NPM1:4869,'H3-3A':3020,H3C2:8358,
  BRCA2:675,BRCA1:672,ATM:472,CHEK2:11200,PALB2:79728,CDK12:51755,ATR:545,FANCA:2175,MLH1:4292,MRE11:4361,NBN:4683,RAD51C:5889,PIK3CA:5290};
const MUT = [];
for (const [g, id] of Object.entries(GENES))
  for (const m of await post(`/api/molecular-profiles/${MU}/mutations/fetch?projection=SUMMARY`, {sampleListId: `${STUDY}_all`, entrezGeneIds: [id]}))
    MUT.push({s: m.sampleId, g, pc: m.proteinChange, t: m.mutationType});
const SV = (await post('/api/structural-variant/fetch', {molecularProfileIds: [SVP], entrezGeneIds: [3084,5979,238,2260,4297,673,2263,9817]}))
  .map(v => ({s: v.sampleId, g1: v.site1HugoSymbol, g2: v.site2HugoSymbol}));
const CNA = [];
for (const [g, id] of Object.entries({ERBB2:2064,BRCA2:675,PTEN:5728,BRCA1:672,ATM:472,CHEK2:11200,PALB2:79728,CDK12:51755,ATR:545,FANCA:2175,MLH1:4292,MRE11:4361,NBN:4683,RAD51C:5889}))
  for (const x of await post(`/api/molecular-profiles/${CN}/discrete-copy-number/fetch?discreteCopyNumberEventType=HOMDEL_AND_AMP&projection=SUMMARY`, {sampleListId: `${STUDY}_all`, entrezGeneIds: [id]}))
    CNA.push({s: x.sampleId, g, v: x.alteration});

// 4. Helpers
const by = (arr) => arr.reduce((o, x) => ((o[x.s] = o[x.s] || []).push(x), o), {});
const mutBy = by(MUT), svBy = by(SV), cnaBy = by(CNA);
const codon = pc => { const m = (pc || '').match(/^[A-Z*]?(\d+)/); return m ? +m[1] : null; };
const TRUNC = new Set(['Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','Translation_Start_Site','Nonstop_Mutation']);
const muts = (s, g, f) => (mutBy[s] || []).some(m => m.g === g && f(m));
const fusion = (s, g) => (svBy[s] || []).some(v => (v.g1 === g || v.g2 === g) && v.g1 && v.g2 && v.g1 !== v.g2);
const cna = (s, g, v) => (cnaBy[s] || []).some(c => c.g === g && c.v === v);
const HEME = /Leukemia|Lymphoma|B-Cell|T and NK|Myelo|Hodgkin|Histiocytosis|Mastocytosis|Blood Cancer|Lymphatic|Posttransplant|Blastic Plasmacytoid/;
const solid = x => !HEME.test(x.ct) && x.ct !== 'UNKNOWN';
const ct = (...t) => x => t.includes(x.ct), ot = (...c) => x => c.includes(x.ot);
const CHOL = ot('IHCH','CHOL','EHCH','PHCH'), BTC = ot('IHCH','CHOL','EHCH','PHCH','GBC','GBAD');
const AML = x => x.ct === 'Leukemia' && (/^AML/.test(x.ot) || ['AM','AMKL','AMOL','AMML','AWM','TAML','APLPMLRARA','MS'].includes(x.ot));
const HRR = ['ATM','ATR','BRCA1','BRCA2','CDK12','CHEK2','FANCA','MLH1','MRE11','NBN','PALB2','RAD51C'];

// 5. Queries: [biomarker, cohort, cohort filter, panel gene(s) for profiling, profile, alteration test]
const Q = [
  ['KRAS G12C','Colorectal',ct('Colorectal Cancer'),['KRAS'],MU,s=>muts(s,'KRAS',m=>m.pc==='G12C')],
  ['KRAS G12C','NSCLC',ct('Non-Small Cell Lung Cancer'),['KRAS'],MU,s=>muts(s,'KRAS',m=>m.pc==='G12C')],
  ['EGFR exon 20 insertion','NSCLC',ct('Non-Small Cell Lung Cancer'),['EGFR'],MU,s=>muts(s,'EGFR',m=>/ins|dup/.test(m.pc||'')&&codon(m.pc)>=762&&codon(m.pc)<=775)],
  ['EGFR ex19del or L858R','NSCLC',ct('Non-Small Cell Lung Cancer'),['EGFR'],MU,s=>muts(s,'EGFR',m=>m.pc==='L858R'||(m.t==='In_Frame_Del'&&codon(m.pc)>=745&&codon(m.pc)<=753))],
  ['ERBB2 TKD mutation','NSCLC',ct('Non-Small Cell Lung Cancer'),['ERBB2'],MU,s=>muts(s,'ERBB2',m=>['Missense_Mutation','In_Frame_Ins','In_Frame_Del'].includes(m.t)&&codon(m.pc)>=720&&codon(m.pc)<=987)],
  ['ERBB2 activating (TKD + S310F/Y)','NSCLC',ct('Non-Small Cell Lung Cancer'),['ERBB2'],MU,s=>muts(s,'ERBB2',m=>(['Missense_Mutation','In_Frame_Ins','In_Frame_Del'].includes(m.t)&&codon(m.pc)>=720&&codon(m.pc)<=987)||/^S310[FY]$/.test(m.pc))],
  ['NRG1 fusion','NSCLC',ct('Non-Small Cell Lung Cancer'),['NRG1'],SVP,s=>fusion(s,'NRG1')],
  ['NRG1 fusion','Pancreatic',ct('Pancreatic Cancer'),['NRG1'],SVP,s=>fusion(s,'NRG1')],
  ['NRG1 fusion','Cholangiocarcinoma',CHOL,['NRG1'],SVP,s=>fusion(s,'NRG1')],
  ['BRAF V600E','All solid tumors',solid,['BRAF'],MU,s=>muts(s,'BRAF',m=>m.pc==='V600E')],
  ['BRAF V600E','Solid tumors not previously covered',x=>solid(x)&&!['Melanoma','Non-Small Cell Lung Cancer','Colorectal Cancer'].includes(x.ct)&&x.ot!=='THAP',['BRAF'],MU,s=>muts(s,'BRAF',m=>m.pc==='V600E')],
  ['RET fusion','Solid tumors excl. NSCLC and thyroid',x=>solid(x)&&!['Non-Small Cell Lung Cancer','Thyroid Cancer'].includes(x.ct),['RET'],SVP,s=>fusion(s,'RET')],
  ['IDH1 R132','Cholangiocarcinoma',CHOL,['IDH1'],MU,s=>muts(s,'IDH1',m=>codon(m.pc)===132)],
  ['IDH1 R132','MDS',ct('Myelodysplastic Syndromes'),['IDH1'],MU,s=>muts(s,'IDH1',m=>codon(m.pc)===132)],
  ['IDH1 R132 or IDH2 R172','Lower-grade glioma',ot('ASTR','DASTR','ODG','OAST'),['IDH1'],MU,s=>muts(s,'IDH1',m=>codon(m.pc)===132)||muts(s,'IDH2',m=>codon(m.pc)===172)],
  ['ESR1 ligand-binding domain','Breast',ct('Breast Cancer'),['ESR1'],MU,s=>muts(s,'ESR1',m=>m.t==='Missense_Mutation'&&[380,463,536,537,538].includes(codon(m.pc)))],
  ['AKT1 E17K or PTEN loss, no PIK3CA hotspot','Breast',ct('Breast Cancer'),['AKT1'],MU,s=>(muts(s,'AKT1',m=>m.pc==='E17K')||muts(s,'PTEN',m=>TRUNC.has(m.t))||cna(s,'PTEN',-2))&&!muts(s,'PIK3CA',m=>[420,542,545,546,1047].includes(codon(m.pc)))],
  ['NPM1 exon 12 insertion','AML',AML,['NPM1'],MU,s=>muts(s,'NPM1',m=>/Frame_Shift|In_Frame_Ins/.test(m.t)&&codon(m.pc)>=280)],
  ['KMT2A fusion','AML',AML,['KMT2A'],SVP,s=>fusion(s,'KMT2A')],
  ['KMT2A fusion','B-lymphoblastic leukemia',ct('B-Lymphoblastic Leukemia/Lymphoma'),['KMT2A'],SVP,s=>fusion(s,'KMT2A')],
  ['FGFR1 fusion','MLNER',ot('MLNER'),['FGFR1'],SVP,s=>fusion(s,'FGFR1')],
  ['H3 K27M','DIPG',ot('DIPG'),['H3F3A','HIST1H3B'],MU,s=>muts(s,'H3-3A',m=>/^K2[78]M$/.test(m.pc))||muts(s,'H3C2',m=>/^K2[78]M$/.test(m.pc))],
  ['H3 K27M','All glioma',ct('Glioma'),['H3F3A','HIST1H3B'],MU,s=>muts(s,'H3-3A',m=>/^K2[78]M$/.test(m.pc))||muts(s,'H3C2',m=>/^K2[78]M$/.test(m.pc))],
  ['KRAS codon 12/13/61','Low-grade serous ovarian',ot('LGSOC'),['KRAS'],MU,s=>muts(s,'KRAS',m=>[12,13,61].includes(codon(m.pc)))],
  ['ALK fusion','Inflammatory myofibroblastic tumor',ot('IMT'),['ALK'],SVP,s=>fusion(s,'ALK')],
  ['ALK fusion','ALCL',ot('ALCL','ALCLALKP','ALCLALKN'),['ALK'],SVP,s=>fusion(s,'ALK')],
  ['BRAF fusion or V600','Low-grade glioma histologies',ot('PAST','LGGNOS','GNG','DNT','PXA'),['BRAF'],MU,s=>fusion(s,'BRAF')||muts(s,'BRAF',m=>codon(m.pc)===600)],
  ['BRAF fusion only','Low-grade glioma histologies',ot('PAST','LGGNOS','GNG','DNT','PXA'),['BRAF'],SVP,s=>fusion(s,'BRAF')],
  ['ERBB2 amplification','Colorectal',ct('Colorectal Cancer'),['ERBB2'],CN,s=>cna(s,'ERBB2',2)],
  ['ERBB2 amplification','Biliary tract',BTC,['ERBB2'],CN,s=>cna(s,'ERBB2',2)],
  ['ERBB2 amplification','Esophagogastric',ct('Esophagogastric Cancer'),['ERBB2'],CN,s=>cna(s,'ERBB2',2)],
  ['BRCA2 truncating or deep deletion','Prostate',ct('Prostate Cancer'),['BRCA2'],MU,s=>muts(s,'BRCA2',m=>TRUNC.has(m.t))||cna(s,'BRCA2',-2)],
  ['HRR (12 genes) truncating or deep deletion','Prostate',ct('Prostate Cancer'),['BRCA2'],MU,s=>HRR.some(g=>muts(s,g,m=>TRUNC.has(m.t))||cna(s,g,-2))],
  ['PTEN truncating or deep deletion','Prostate',ct('Prostate Cancer'),['PTEN'],MU,s=>muts(s,'PTEN',m=>TRUNC.has(m.t))||cna(s,'PTEN',-2)],
];

// 6. Patient-level counts
const results = Q.map(([biomarker, cohort, f, genes, prof, alt]) => {
  const denom = new Set(), num = new Set();
  for (const [s, x] of Object.entries(S)) {
    if (!f(x) || !profiled(s, genes, prof)) continue;
    denom.add(x.p); if (alt(s)) num.add(x.p);
  }
  return {biomarker, cohort, altered: num.size, profiled: denom.size, pct: denom.size ? +(100 * num.size / denom.size).toFixed(2) : null};
});
console.table(results);
