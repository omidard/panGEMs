/* panGEMs Analytics — compare models, compare groups, clustermap, species tree.
   Globals used: jQuery ($), Chart, d3. Data: gems_metadata.json + assets/presence_* */
(function () {
'use strict';
const ECO = '#2c6fbb', LACTO = '#c26a10';
let META = [], ORDER = [], N = 0;
const MAT = {};                 // 'reactions'|'metabolites' -> {buf,rowbytes,V,vocab,idx}
const rowByFile = {};
let charts = {};

// ---- bit helpers ----
const POP = new Uint8Array(256);
const BYTE_BITS = [];
for (let b = 0; b < 256; b++) { let c = 0, arr = []; for (let k = 0; k < 8; k++) if (b & (1 << k)) { c++; arr.push(k); } POP[b] = c; BYTE_BITS[b] = arr; }
function countRow(row, m) { const rb = m.rowbytes, base = row * rb, buf = m.buf; let c = 0; for (let i = 0; i < rb; i++) c += POP[buf[base + i]]; return c; }
function interCount(a, b, m) { const rb = m.rowbytes, ba = a * rb, bb = b * rb, buf = m.buf; let c = 0; for (let i = 0; i < rb; i++) c += POP[buf[ba + i] & buf[bb + i]]; return c; }
function has(row, col, m) { return (m.buf[row * m.rowbytes + (col >> 3)] >> (col & 7)) & 1; }
function prevalence(rows, m) { // Int32 counts per column
  const rb = m.rowbytes, buf = m.buf, cnt = new Int32Array(m.V);
  for (let r = 0; r < rows.length; r++) { const base = rows[r] * rb; for (let i = 0; i < rb; i++) { let v = buf[base + i]; if (!v) continue; const bits = BYTE_BITS[v], off = i << 3; for (let k = 0; k < bits.length; k++) cnt[off + bits[k]]++; } }
  return cnt;
}
function colsOfRow(row, m) { const rb = m.rowbytes, base = row * rb, out = []; for (let i = 0; i < rb; i++) { let v = m.buf[base + i]; if (!v) continue; const bits = BYTE_BITS[v], off = i << 3; for (let k = 0; k < bits.length; k++) out.push(off + bits[k]); } return out; }
function mean(a){ return a.reduce((s,x)=>s+x,0)/a.length; }
function sd(a){ const mu=mean(a); return Math.sqrt(a.reduce((s,x)=>s+(x-mu)*(x-mu),0)/a.length); }
function esc(s){ return String(s==null?'':s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c])); }

// ---- load ----
async function load() {
  const msg = document.getElementById('loading-msg');
  msg.textContent = 'Loading model metadata…';
  META = await fetch('gems_metadata.json').then(r => r.json());
  const man = await fetch('assets/presence_manifest.json').then(r => r.json());
  ORDER = man.order; N = man.n_models;
  META.forEach((r, i) => { rowByFile[r.gem_file] = i; });
  msg.textContent = 'Loading reaction matrix…';
  const [rv, mv, rb, mb] = await Promise.all([
    fetch('assets/reactions_vocab.json').then(r => r.json()),
    fetch('assets/metabolites_vocab.json').then(r => r.json()),
    fetch('assets/presence_reactions.bin').then(r => r.arrayBuffer()),
    fetch('assets/presence_metabolites.bin').then(r => r.arrayBuffer()),
  ]);
  MAT.reactions = { buf: new Uint8Array(rb), rowbytes: man.reaction_rowbytes, V: man.n_reactions, vocab: rv, idx: null };
  MAT.metabolites = { buf: new Uint8Array(mb), rowbytes: man.metabolite_rowbytes, V: man.n_metabolites, vocab: mv, idx: null };
  document.getElementById('loading').style.display = 'none';
  document.getElementById('app').style.display = 'block';
  initUI();
}

function organisms() { return [...new Set(META.map(r => r.organism).filter(Boolean))].sort(); }

// ---- overlap note ----
function overlapNote() {
  const eco = META.map((r,i)=>[r,i]).filter(x=>x[0].dataset==='EcopanGEM').map(x=>x[1]);
  const lac = META.map((r,i)=>[r,i]).filter(x=>x[0].dataset==='LactoPanGEM').map(x=>x[1]);
  const pe = prevalence(eco, MAT.reactions), pl = prevalence(lac, MAT.reactions);
  let e=0,l=0,sh=0; for (let c=0;c<MAT.reactions.V;c++){ const a=pe[c]>0,b=pl[c]>0; if(a)e++; if(b)l++; if(a&&b)sh++; }
  document.getElementById('overlap-note').innerHTML =
    `<strong>Tip:</strong> comparisons are richest <em>within</em> a collection (E.&nbsp;coli&nbsp;↔&nbsp;E.&nbsp;coli, or Lactobacillaceae species&nbsp;↔&nbsp;species). Across collections the reactomes overlap only partially — <strong>${sh.toLocaleString()}</strong> reactions are shared of ${e.toLocaleString()} (E.&nbsp;coli) and ${l.toLocaleString()} (Lactobacillaceae).`;
}

// ---- tabs ----
function initUI() {
  overlapNote();
  document.querySelectorAll('#tabs .nav-link').forEach(btn => btn.addEventListener('click', () => {
    document.querySelectorAll('#tabs .nav-link').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    document.querySelectorAll('.tabpane').forEach(p => p.style.display = 'none');
    document.getElementById('tab-' + btn.dataset.tab).style.display = 'block';
  }));
  initModels(); initGroups(); initCluster(); initTree();
}

/* ================= TOOL 1: COMPARE MODELS ================= */
let mSel = [], mMetric = 'reactions';
function initModels() {
  const dl = document.getElementById('m-list');
  META.forEach(r => { const o = document.createElement('option'); o.value = r.gem_file; o.label = `${r.organism || ''} ${r.strain || ''}`.trim(); dl.appendChild(o); });
  document.getElementById('m-add').addEventListener('click', addModel);
  document.getElementById('m-search').addEventListener('keydown', e => { if (e.key === 'Enter') addModel(); });
  document.getElementById('m-clear').addEventListener('click', () => { mSel = []; renderModels(); });
  document.getElementById('m-examples').addEventListener('click', () => {
    // three E. coli + two related lacto to show within/across contrast
    const pick = [];
    const eco = META.filter(r => r.dataset === 'EcopanGEM'); const lac = META.filter(r => r.organism === 'Lactobacillus gasseri');
    [eco[0], eco[1], lac[0], lac[1]].forEach(r => r && pick.push(rowByFile[r.gem_file]));
    mSel = [...new Set(pick)].slice(0, 4); renderModels();
  });
  document.querySelectorAll('#m-metric button').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('#m-metric button').forEach(x => x.classList.remove('active')); b.classList.add('active');
    mMetric = b.dataset.metric; document.querySelectorAll('.m-metric-label').forEach(e => e.textContent = mMetric); renderModels();
  }));
  document.getElementById('m-csv').addEventListener('click', () => downloadDiffCSV());
}
function addModel() {
  const q = document.getElementById('m-search').value.trim(); if (!q) return;
  let row = rowByFile[q];
  if (row == null) { const ql = q.toLowerCase(); const hit = META.find(r => (r.gem_file + ' ' + r.organism + ' ' + r.strain).toLowerCase().includes(ql)); if (hit) row = rowByFile[hit.gem_file]; }
  if (row == null) { alert('No model matches “' + q + '”'); return; }
  if (!mSel.includes(row)) mSel.push(row);
  if (mSel.length > 12) mSel = mSel.slice(-12);
  document.getElementById('m-search').value = ''; renderModels();
}
function renderModels() {
  const chips = document.getElementById('m-chips');
  chips.innerHTML = mSel.map(r => { const m = META[r]; const c = m.dataset === 'EcopanGEM' ? ECO : LACTO;
    return `<span class="chip" style="background:${c}18;color:${c}"><em>${esc(m.organism)}</em> ${esc(m.strain||m.gem_file)} <span class="x" data-r="${r}">×</span></span>`; }).join('');
  chips.querySelectorAll('.x').forEach(x => x.addEventListener('click', () => { mSel = mSel.filter(r => r != x.dataset.r); renderModels(); }));
  const res = document.getElementById('m-results');
  if (mSel.length < 2) { res.style.display = 'none'; return; }
  res.style.display = 'block';
  const m = MAT[mMetric];
  const sets = mSel.map(r => new Set(colsOfRow(r, m)));
  const union = new Set(); sets.forEach(s => s.forEach(x => union.add(x)));
  const core = [...union].filter(c => sets.every(s => s.has(c)));
  // stats
  document.getElementById('m-stats').innerHTML =
    stat(mSel.length, 'Models') + stat(union.size, 'Union (pan)') + stat(core.length, 'Core (in all)') +
    stat((union.size - core.length), 'Variable');
  // pairwise Jaccard heatmap (HTML table)
  let jh = '<table class="mini" style="width:auto"><tr><th></th>' + mSel.map(r => `<th title="${esc(META[r].organism)} ${esc(META[r].strain||'')}">${shortLabel(r)}</th>`).join('') + '</tr>';
  for (let i = 0; i < mSel.length; i++) { jh += `<tr><th style="text-align:right">${shortLabel(mSel[i])}</th>`;
    for (let j = 0; j < mSel.length; j++) { const inter = interCount(mSel[i], mSel[j], m); const uni = sets[i].size + sets[j].size - inter; const J = uni ? inter / uni : 1; const bg = `rgba(44,111,187,${(0.12 + 0.85 * J).toFixed(2)})`; jh += `<td style="background:${bg};text-align:center;color:${J>0.6?'#fff':'#334'}" title="${(J*100).toFixed(1)}% shared">${(J*100).toFixed(0)}</td>`; } jh += '</tr>'; }
  jh += '</table>'; document.getElementById('m-jaccard').innerHTML = jh;
  // set breakdown chart: per model unique vs shared
  const uniqueOf = mSel.map((r, i) => [...sets[i]].filter(c => sets.every((s, j) => j === i || !s.has(c))).length);
  if (charts.mset) charts.mset.destroy();
  charts.mset = new Chart(document.getElementById('m-setchart'), {
    type: 'bar',
    data: { labels: mSel.map(r => shortLabel(r)), datasets: [
      { label: 'Shared with ≥1 other', data: mSel.map((r,i)=>sets[i].size - uniqueOf[i]), backgroundColor: '#9ec3ea' },
      { label: 'Unique to this model', data: uniqueOf, backgroundColor: LACTO } ] },
    options: { responsive:true, maintainAspectRatio:false, plugins:{legend:{labels:{boxWidth:10,font:{size:10}}}},
      scales:{ x:{stacked:true,grid:{display:false},ticks:{font:{size:9}}}, y:{stacked:true,title:{display:true,text:mMetric,font:{size:10}}} } }
  });
  // differential table
  renderDiffTable(sets, union, core);
}
function shortLabel(r){ const m=META[r]; const g=(m.organism||'').split(' '); const ab=g.length>1?g[0][0]+'. '+g[1]:m.organism; return (m.strain? m.strain : (m.gem_file.replace('.json.json','').replace('.json',''))); }
function stat(v,l){ return `<div class="statbox"><div class="v">${typeof v==='number'?v.toLocaleString():v}</div><div class="l">${l}</div></div>`; }
let mDiffRows = [];
function renderDiffTable(sets, union, core){
  const m = MAT[mMetric]; const coreSet = new Set(core);
  const diff = [...union].filter(c => !coreSet.has(c));
  diff.sort((a,b)=> (m.vocab[a].id<m.vocab[b].id?-1:1));
  mDiffRows = diff;
  const cap = 1200;
  let h = '<thead><tr><th>ID</th><th>Name</th>' + mSel.map(r=>`<th title="${esc(META[r].organism)}">${esc(shortLabel(r))}</th>`).join('') + '</tr></thead><tbody>';
  diff.slice(0,cap).forEach(c => { h += `<tr><td>${esc(m.vocab[c].id)}</td><td>${esc((m.vocab[c].name||'').slice(0,42))}</td>` +
    sets.map(s => s.has(c)?'<td class="present">✓</td>':'<td class="absent">·</td>').join('') + '</tr>'; });
  h += '</tbody>';
  document.getElementById('m-diff').innerHTML = h;
  if (diff.length > cap) document.getElementById('m-diff').insertAdjacentHTML('afterend','');
}
function downloadDiffCSV(){
  if (!mSel.length) return; const m = MAT[mMetric];
  let csv = 'id,name,' + mSel.map(r=>'"'+ (META[r].gem_file) +'"').join(',') + '\n';
  const sets = mSel.map(r=>new Set(colsOfRow(r,m)));
  mDiffRows.forEach(c => { csv += `"${m.vocab[c].id}","${(m.vocab[c].name||'').replace(/"/g,'""')}",` + sets.map(s=>s.has(c)?1:0).join(',') + '\n'; });
  dl(csv, `panGEMs_compare_${mMetric}_${mSel.length}models.csv`);
}
function dl(text, name){ const b=new Blob([text],{type:'text/csv'}); const a=document.createElement('a'); a.href=URL.createObjectURL(b); a.download=name; a.click(); URL.revokeObjectURL(a.href); }

/* ================= TOOL 2: COMPARE GROUPS ================= */
let gMetric='reactions', groupCache={A:[],B:[]};
function groupBuilder(id, letter, defaults){
  const orgs = organisms();
  const countries = [...new Set(META.map(r=>r.country).filter(Boolean))].sort();
  const hosts = [...new Set(META.map(r=>r.host_name).filter(Boolean))].sort();
  const el = document.getElementById(id);
  el.innerHTML = `<div style="font-weight:700;color:${letter==='A'?ECO:LACTO};margin-bottom:0.4rem">Group ${letter}</div>
    <div class="row g-2">
      <div class="col-6"><label>Dataset</label><select class="form-select form-select-sm" data-f="dataset"><option value="">All</option><option>EcopanGEM</option><option>LactoPanGEM</option></select></div>
      <div class="col-6"><label>Country</label><select class="form-select form-select-sm" data-f="country"><option value="">All</option>${countries.map(c=>`<option>${esc(c)}</option>`).join('')}</select></div>
      <div class="col-12"><label>Species (pick one or more; empty = all)</label><select multiple class="form-select form-select-sm" data-f="organism">${orgs.map(o=>`<option>${esc(o)}</option>`).join('')}</select></div>
      <div class="col-12"><label>Host</label><select class="form-select form-select-sm" data-f="host_name"><option value="">All</option>${hosts.map(h=>`<option>${esc(h)}</option>`).join('')}</select></div>
    </div>
    <div class="mt-1" style="font-size:0.8rem;color:#556">Matches: <strong class="gcount">0</strong> models</div>`;
  const resolve = () => {
    const ds = el.querySelector('[data-f=dataset]').value;
    const ct = el.querySelector('[data-f=country]').value;
    const hn = el.querySelector('[data-f=host_name]').value;
    const orgSel = [...el.querySelector('[data-f=organism]').selectedOptions].map(o=>o.value);
    const rows = [];
    for (let i=0;i<META.length;i++){ const r=META[i];
      if (ds && r.dataset!==ds) continue;
      if (ct && r.country!==ct) continue;
      if (hn && r.host_name!==hn) continue;
      if (orgSel.length && !orgSel.includes(r.organism)) continue;
      rows.push(i); }
    el.querySelector('.gcount').textContent = rows.length.toLocaleString();
    return rows;
  };
  el.addEventListener('change', ()=>{ groupCache[letter]=resolve(); });
  if (defaults){ if(defaults.dataset) el.querySelector('[data-f=dataset]').value=defaults.dataset;
    if(defaults.organism){ [...el.querySelector('[data-f=organism]').options].forEach(o=>{ if(defaults.organism.includes(o.value)) o.selected=true; }); } }
  groupCache[letter]=resolve();
  return resolve;
}
function initGroups(){
  const rA = groupBuilder('grpA','A',{dataset:'EcopanGEM'});
  const rB = groupBuilder('grpB','B',{organism:['Lactiplantibacillus plantarum']});
  document.querySelectorAll('#g-metric button').forEach(b=>b.addEventListener('click',()=>{ document.querySelectorAll('#g-metric button').forEach(x=>x.classList.remove('active')); b.classList.add('active'); gMetric=b.dataset.metric; document.querySelectorAll('.g-metric-label').forEach(e=>e.textContent=gMetric[0].toUpperCase()+gMetric.slice(1)); }));
  document.getElementById('g-run').addEventListener('click', ()=>runGroups(rA(), rB()));
  document.getElementById('g-csv').addEventListener('click', ()=>gDownloadCSV());
}
let gDiffData=null;
function runGroups(A, B){
  if (!A.length || !B.length){ alert('Both groups need at least one model.'); return; }
  const m = MAT[gMetric];
  const cntA = prevalence(A,m), cntB = prevalence(B,m);
  const perA = A.map(r=>countRow(r,m)), perB = B.map(r=>countRow(r,m));
  let panA=0,panB=0,coreA=0,coreB=0,shared=0,Aonly=0,Bonly=0,bothCore=0;
  for (let c=0;c<m.V;c++){ const a=cntA[c]>0,b=cntB[c]>0; if(a)panA++; if(b)panB++; if(cntA[c]===A.length)coreA++; if(cntB[c]===B.length)coreB++; if(a&&b)shared++; else if(a)Aonly++; else if(b)Bonly++; if(cntA[c]===A.length&&cntB[c]===B.length)bothCore++; }
  document.getElementById('g-results').style.display='block';
  document.getElementById('g-stats').innerHTML =
    stat(A.length,'Group A models')+stat(B.length,'Group B models')+
    stat(Math.round(mean(perA)),'A mean '+gMetric)+stat(Math.round(mean(perB)),'B mean '+gMetric)+
    stat(coreA,'A core')+stat(coreB,'B core')+stat(shared,'Shared (pan)')+stat(bothCore,'Shared core');
  // distribution chart (overlaid histograms)
  distChart(perA, perB);
  overlapChart(Aonly, shared, Bonly);
  // differential prevalence table
  const rows=[]; for (let c=0;c<m.V;c++){ const pa=cntA[c]/A.length, pb=cntB[c]/B.length; if (pa===0&&pb===0) continue; rows.push([c, pa, pb, pa-pb]); }
  rows.sort((x,y)=>Math.abs(y[3])-Math.abs(x[3]));
  gDiffData={rows,m,A,B};
  const cap=600; let h='<thead><tr><th>ID</th><th>Name</th><th>% in A</th><th>% in B</th><th>Δ</th><th>Enriched</th></tr></thead><tbody>';
  rows.slice(0,cap).forEach(([c,pa,pb,d])=>{ const en=d>0?`<span class="ds-eco">A</span>`:`<span class="ds-lacto">B</span>`; h+=`<tr><td>${esc(m.vocab[c].id)}</td><td>${esc((m.vocab[c].name||'').slice(0,40))}</td><td>${(pa*100).toFixed(0)}%</td><td>${(pb*100).toFixed(0)}%</td><td style="font-weight:600;color:${d>0?ECO:LACTO}">${d>0?'+':''}${(d*100).toFixed(0)}</td><td>${Math.abs(d)>0.001?en:'—'}</td></tr>`; });
  h+='</tbody>'; document.getElementById('g-diff').innerHTML=h;
}
function distChart(perA, perB){
  const all=perA.concat(perB); const lo=Math.min(...all), hi=Math.max(...all); const w=Math.max(20,Math.round((hi-lo)/22)||20);
  const start=Math.floor(lo/w)*w, nb=Math.max(1,Math.ceil((hi-start)/w)+1);
  const labels=[], A=new Array(nb).fill(0), B=new Array(nb).fill(0);
  for(let i=0;i<nb;i++) labels.push((start+i*w));
  perA.forEach(v=>A[Math.min(nb-1,Math.floor((v-start)/w))]++); perB.forEach(v=>B[Math.min(nb-1,Math.floor((v-start)/w))]++);
  if(charts.gdist) charts.gdist.destroy();
  charts.gdist=new Chart(document.getElementById('g-dist'),{type:'bar',
    data:{labels,datasets:[{label:'Group A',data:A,backgroundColor:ECO+'bb'},{label:'Group B',data:B,backgroundColor:LACTO+'bb'}]},
    options:{responsive:true,maintainAspectRatio:false,plugins:{legend:{labels:{boxWidth:10,font:{size:10}}}},scales:{x:{grid:{display:false},title:{display:true,text:gMetric+' per model',font:{size:10}},ticks:{font:{size:9}}},y:{title:{display:true,text:'models',font:{size:10}}}}}});
}
function overlapChart(Aonly, shared, Bonly){
  if(charts.gov) charts.gov.destroy();
  charts.gov=new Chart(document.getElementById('g-overlap'),{type:'bar',
    data:{labels:['Repertoire (pan)'],datasets:[
      {label:'A only',data:[Aonly],backgroundColor:ECO},
      {label:'Shared',data:[shared],backgroundColor:'#7ea9d8'},
      {label:'B only',data:[Bonly],backgroundColor:LACTO}]},
    options:{indexAxis:'y',responsive:true,maintainAspectRatio:false,plugins:{legend:{labels:{boxWidth:10,font:{size:10}}},tooltip:{callbacks:{label:c=>`${c.dataset.label}: ${c.parsed.x.toLocaleString()} ${gMetric}`}}},scales:{x:{stacked:true,title:{display:true,text:gMetric,font:{size:10}}},y:{stacked:true}}}});
}
function gDownloadCSV(){ if(!gDiffData) return; const {rows,m,A,B}=gDiffData; let csv=`id,name,pct_in_A,pct_in_B,delta\n`; rows.forEach(([c,pa,pb,d])=>{ csv+=`"${m.vocab[c].id}","${(m.vocab[c].name||'').replace(/"/g,'""')}",${(pa*100).toFixed(1)},${(pb*100).toFixed(1)},${(d*100).toFixed(1)}\n`; }); dl(csv,`panGEMs_groups_${gMetric}_diff.csv`); }

/* ================= TOOL 3: CLUSTERMAP ================= */
let cMetric='reactions';
function initCluster(){
  const sp=document.getElementById('c-species'); organisms().forEach(o=>{ const opt=document.createElement('option'); opt.value=o; opt.textContent=`${o} (${META.filter(r=>r.organism===o).length})`; sp.appendChild(opt); });
  document.getElementById('c-source').addEventListener('change',e=>{ document.getElementById('c-species-wrap').style.display = e.target.value==='species'?'block':'none'; });
  document.querySelectorAll('#c-metric button').forEach(b=>b.addEventListener('click',()=>{ document.querySelectorAll('#c-metric button').forEach(x=>x.classList.remove('active')); b.classList.add('active'); cMetric=b.dataset.metric; }));
  document.getElementById('c-run').addEventListener('click', buildCluster);
}
function clusterRows(){
  const src=document.getElementById('c-source').value;
  if (src==='species'){ const o=document.getElementById('c-species').value; return META.map((r,i)=>[r,i]).filter(x=>x[0].organism===o).map(x=>x[1]); }
  if (src==='groupA') return groupCache.A.slice();
  if (src==='groupB') return groupCache.B.slice();
  if (src==='selected') return mSel.slice();
  return [];
}
function sampleRows(rows, cap){ if (rows.length<=cap) return rows.slice(); const step=rows.length/cap, out=[]; for(let i=0;i<cap;i++) out.push(rows[Math.floor(i*step)]); return out; }
function buildCluster(){
  let rows=clusterRows();
  if (rows.length<3){ alert('Pick a set with at least 3 models (choose a species, or run Compare groups first).'); return; }
  const capM=120, capF=140;
  const sampled = sampleRows(rows, capM);
  const m=MAT[cMetric];
  const cnt=prevalence(sampled,m);
  // variable features
  const varc=[]; for(let c=0;c<m.V;c++){ const k=cnt[c]; if(k>0&&k<sampled.length){ const p=k/sampled.length; varc.push([c, p*(1-p)]); } }
  varc.sort((a,b)=>b[1]-a[1]);
  const cols=varc.slice(0,capF).map(x=>x[0]);
  if (!cols.length){ alert('These models are (nearly) identical in content — no variable '+cMetric+' to cluster.'); return; }
  // binary submatrix
  const R=sampled.length, C=cols.length; const bin=new Uint8Array(R*C);
  for(let i=0;i<R;i++){ const row=sampled[i]; for(let j=0;j<C;j++) bin[i*C+j]=has(row,cols[j],m); }
  // distances (Hamming) + cluster rows and cols
  const rowOrder=hclust(bin,R,C,true);
  const colOrder=hclust(bin,R,C,false);
  drawHeat(sampled,cols,bin,rowOrder,colOrder,m,rows.length);
}
function hclust(bin,R,C,byRow){
  const n=byRow?R:C;
  const dist=[]; for(let i=0;i<n;i++) dist.push(new Float32Array(n));
  const get=(i,k)=> byRow? bin[i*C+k] : bin[k*C+i];
  const len=byRow?C:R;
  for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){ let d=0; for(let k=0;k<len;k++) if(get(i,k)!==get(j,k)) d++; dist[i][j]=d; dist[j][i]=d; }
  // average linkage agglomerative
  const nodes={}; let active=[]; for(let i=0;i<n;i++){ nodes[i]={m:[i],l:null,r:null}; active.push(i); }
  const D={}; active.forEach(a=>{ D[a]={}; }); for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){ D[i][j]=dist[i][j]; D[j][i]=dist[i][j]; }
  let next=n;
  while(active.length>1){ let bi=0,bj=1,best=Infinity; for(let x=0;x<active.length;x++) for(let y=x+1;y<active.length;y++){ const a=active[x],b=active[y],d=D[a][b]; if(d<best){best=d;bi=a;bj=b;} }
    const na=nodes[bi].m.length, nb=nodes[bj].m.length, id=next++; nodes[id]={m:nodes[bi].m.concat(nodes[bj].m),l:bi,r:bj}; D[id]={};
    active.forEach(k=>{ if(k===bi||k===bj) return; const d=(na*D[bi][k]+nb*D[bj][k])/(na+nb); D[id][k]=d; D[k][id]=d; });
    active=active.filter(k=>k!==bi&&k!==bj); active.push(id); }
  const order=[]; (function trav(id){ const nd=nodes[id]; if(nd.l==null){order.push(id);return;} trav(nd.l); trav(nd.r); })(active[0]);
  return order;
}
function drawHeat(rows,cols,bin,rowOrder,colOrder,m,totalN){
  const cv=document.getElementById('cluster-canvas'), ctx=cv.getContext('2d');
  const R=rows.length,C=cols.length; const cell=Math.max(6,Math.min(11,Math.floor(900/C)));
  const labelW=190, stripW=10, topH=8, dpr=window.devicePixelRatio||1;
  const W=labelW+stripW+C*cell, H=topH+R*cell;
  cv.width=W*dpr; cv.height=H*dpr; cv.style.width=W+'px'; cv.style.height=H+'px'; ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H); ctx.font='9px Segoe UI'; ctx.textBaseline='middle';
  const spCol={}; const palette=d3.schemeTableau10;
  const uniqOrg=[...new Set(rows.map(r=>META[r].organism))]; uniqOrg.forEach((o,i)=>spCol[o]=palette[i%palette.length]);
  rowOrder.forEach((ri,y)=>{ const row=rows[ri], mrow=META[row];
    // species strip
    ctx.fillStyle=spCol[mrow.organism]; ctx.fillRect(labelW,topH+y*cell,stripW-1,cell-0.5);
    // label
    ctx.fillStyle='#334'; const lab=(mrow.strain||mrow.gem_file.replace('.json.json','').replace('.json','')).slice(0,30); ctx.textAlign='right'; ctx.fillText(lab, labelW-4, topH+y*cell+cell/2);
    colOrder.forEach((ci,x)=>{ if(bin[ri*C+ci]){ ctx.fillStyle=mrow.dataset==='EcopanGEM'?ECO:LACTO; ctx.fillRect(labelW+stripW+x*cell, topH+y*cell, cell-0.6, cell-0.6); } });
  });
  document.getElementById('c-results').style.display='block';
  const shown = rows.length<totalN? ` (showing ${rows.length} of ${totalN}, evenly sampled)`:'';
  document.getElementById('c-legend').innerHTML=`<strong>${rows.length}</strong> models${shown} × <strong>${cols.length}</strong> most-variable ${cMetric} · filled cell = present · rows clustered by shared content · left strip = species.`;
  // hover
  let tip=document.querySelector('.cluster-tip'); if(!tip){ tip=document.createElement('div'); tip.className='cluster-tip'; document.body.appendChild(tip); }
  cv.onmousemove=(ev)=>{ const rect=cv.getBoundingClientRect(); const x=ev.clientX-rect.left-labelW-stripW, y=ev.clientY-rect.top-topH; const cx=Math.floor(x/cell), cy=Math.floor(y/cell); if(cx<0||cy<0||cx>=C||cy>=R){tip.style.opacity=0;return;} const row=rows[rowOrder[cy]], col=cols[colOrder[cx]]; const pres=bin[rowOrder[cy]*C+colOrder[cx]]; tip.innerHTML=`<strong>${esc(META[row].strain||META[row].gem_file)}</strong><br><em>${esc(META[row].organism)}</em><br>${esc(m.vocab[col].id)} — ${pres?'present':'absent'}`; tip.style.opacity=1; tip.style.left=(ev.pageX+12)+'px'; tip.style.top=(ev.pageY-10)+'px'; };
  cv.onmouseleave=()=>{ tip.style.opacity=0; };
}

/* ================= TOOL 4: SPECIES TREE ================= */
let tMode='taxonomy';
const FAMILY = org => org==='Escherichia coli' ? 'Enterobacteriaceae' : 'Lactobacillaceae';
function initTree(){
  document.querySelectorAll('#t-mode button').forEach(b=>b.addEventListener('click',()=>{ document.querySelectorAll('#t-mode button').forEach(x=>x.classList.remove('active')); b.classList.add('active'); tMode=b.dataset.mode; drawTree(); }));
  drawTree();
}
function speciesMembers(o){ return META.map((r,i)=>[r,i]).filter(x=>x[0].organism===o).map(x=>x[1]); }
function taxonomyRoot(){
  const orgs=organisms(); const fam={};
  orgs.forEach(o=>{ const f=FAMILY(o), g=o.split(' ')[0]; (fam[f]=fam[f]||{})[g]=fam[f][g]||[]; fam[f][g].push(o); });
  return { name:'Bacteria', children: Object.entries(fam).map(([f,gs])=>({ name:f, children:Object.entries(gs).map(([g,ss])=>({ name:g, children:ss.map(s=>({name:s,leaf:true,count:speciesMembers(s).length,organism:s})) })) })) };
}
function contentRoot(){
  const orgs=organisms(); const m=MAT.reactions;
  const vecs=orgs.map(o=>{ const rows=speciesMembers(o); const cnt=prevalence(rows,m); const v=new Float32Array(m.V); for(let c=0;c<m.V;c++) v[c]=cnt[c]/rows.length; return v; });
  const n=orgs.length; const D={}; for(let i=0;i<n;i++) D[i]={};
  for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){ let dot=0,na=0,nb=0; const a=vecs[i],b=vecs[j]; for(let c=0;c<m.V;c++){ dot+=a[c]*b[c]; na+=a[c]*a[c]; nb+=b[c]*b[c]; } const cos=dot/(Math.sqrt(na*nb)||1); const d=1-cos; D[i][j]=d; D[j][i]=d; }
  // average-linkage tree -> nested
  const nodes={}; let active=[]; for(let i=0;i<n;i++){ nodes[i]={leaf:true,name:orgs[i],organism:orgs[i],count:speciesMembers(orgs[i]).length,m:[i]}; active.push(i); }
  const DD={}; active.forEach(a=>DD[a]={}); for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){ DD[i][j]=D[i][j]; DD[j][i]=D[i][j]; }
  let next=n;
  while(active.length>1){ let bi=0,bj=1,best=Infinity; for(let x=0;x<active.length;x++) for(let y=x+1;y<active.length;y++){ const a=active[x],b=active[y],d=DD[a][b]; if(d<best){best=d;bi=a;bj=b;} }
    const na=nodes[bi].m.length,nb=nodes[bj].m.length,id=next++; nodes[id]={name:'',children:[nodes[bi],nodes[bj]],m:nodes[bi].m.concat(nodes[bj].m)}; DD[id]={};
    active.forEach(k=>{ if(k===bi||k===bj) return; const d=(na*DD[bi][k]+nb*DD[bj][k])/(na+nb); DD[id][k]=d; DD[k][id]=d; });
    active=active.filter(k=>k!==bi&&k!==bj); active.push(id); }
  return nodes[active[0]];
}
function drawTree(){
  const root = tMode==='taxonomy'? taxonomyRoot() : contentRoot();
  const holder=document.getElementById('tree-holder'); holder.innerHTML='';
  const hier=d3.hierarchy(root); const leaves=hier.leaves().length;
  const W=Math.max(760, holder.clientWidth||880), rowH=26, H=leaves*rowH+40;
  const svg=d3.select(holder).append('svg').attr('width',W).attr('height',H);
  const layout=d3.cluster().size([H-30, W-480]); layout(hier);
  const g=svg.append('g').attr('transform','translate(150,15)');
  g.selectAll('path.link').data(hier.links()).join('path').attr('class','link').attr('fill','none').attr('stroke','#cdd6e2').attr('stroke-width',1)
    .attr('d',d=>`M${d.source.y},${d.source.x} C${(d.source.y+d.target.y)/2},${d.source.x} ${(d.source.y+d.target.y)/2},${d.target.x} ${d.target.y},${d.target.x}`);
  let tip=document.querySelector('.tree-tip'); if(!tip){ tip=document.createElement('div'); tip.className='tree-tip'; document.body.appendChild(tip); }
  const node=g.selectAll('g.node').data(hier.descendants()).join('g').attr('transform',d=>`translate(${d.y},${d.x})`);
  node.append('circle').attr('r',d=>d.data.leaf?4:3).attr('fill',d=>{ if(!d.data.leaf) return '#9fb0c4'; return FAMILY(d.data.organism)==='Enterobacteriaceae'?ECO:LACTO; })
    .attr('class',d=>d.data.leaf?'node-species':'').style('cursor',d=>d.data.leaf?'pointer':'default');
  node.append('text').attr('dy','0.32em').attr('x',d=>d.data.leaf?8:-6).attr('text-anchor',d=>d.data.leaf?'start':'end')
    .style('font-size',d=>d.data.leaf?'11px':'10px').style('font-style',d=>d.data.leaf?'italic':'normal')
    .style('fill',d=>d.data.leaf?'#223':'#7a8698').style('cursor',d=>d.data.leaf?'pointer':'default')
    .text(d=>d.data.leaf?`${d.data.name} (${d.data.count})`:d.data.name)
    .attr('class',d=>d.data.leaf?'node-species':'');
  node.filter(d=>d.data.leaf).style('cursor','pointer')
    .on('mousemove',(ev,d)=>{ tip.innerHTML=`<strong><em>${esc(d.data.name)}</em></strong><br>${d.data.count} models · click to open in clustermap`; tip.style.opacity=1; tip.style.left=(ev.pageX+12)+'px'; tip.style.top=(ev.pageY-10)+'px'; })
    .on('mouseleave',()=>tip.style.opacity=0)
    .on('click',(ev,d)=>openSpecies(d.data.organism));
  d3.select(holder).append('div').attr('style','font-size:0.78rem;color:#8792a3;margin-top:0.4rem')
    .html(tMode==='taxonomy'? 'NCBI-style taxonomy: family → genus → species. Circle colour = collection.' : 'Dendrogram from mean reaction-presence per species (1 − cosine, average linkage): species that are metabolically similar sit together.');
}
function openSpecies(o){
  // jump to clustermap for this species
  document.querySelectorAll('#tabs .nav-link').forEach(b=>b.classList.remove('active'));
  document.querySelector('#tabs .nav-link[data-tab=cluster]').classList.add('active');
  document.querySelectorAll('.tabpane').forEach(p=>p.style.display='none');
  document.getElementById('tab-cluster').style.display='block';
  document.getElementById('c-source').value='species'; document.getElementById('c-species-wrap').style.display='block';
  document.getElementById('c-species').value=o;
  buildCluster();
}

load().catch(e=>{ document.getElementById('loading-msg').textContent='Failed to load analytics data: '+e.message; console.error(e); });
})();
