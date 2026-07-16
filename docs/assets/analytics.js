/* panGEMs Analytics — comparative genomics of 4,659 genome-scale metabolic models.
   Client-side. Globals: jQuery ($), Chart, d3, Plotly.
   Data: gems_metadata.json + assets/presence_{manifest,reactions,metabolites}.{json,bin} */
(function () {
'use strict';

/* =============================== DESIGN TOKENS =============================== */
const ECO = '#2563EB', LACTO = '#D97706';
const INK1 = '#152238', INK2 = '#4a5567', INK3 = '#8b95a5', LINE = '#e6eaf0';
const PLOT_FONT = { family: "'Fira Sans', system-ui, -apple-system, sans-serif", color: INK2, size: 11 };
// validated 8-hue categorical palette (dataviz default; WARN-band CVD -> used only with legend + hover)
const CAT8 = ['#2a78d6', '#1baf7a', '#eda100', '#008300', '#4a3aa7', '#e34948', '#e87ba4', '#eb6834'];
const OTHER_COL = '#9aa4b2';
const BLUES = d3.interpolateBlues;

let META = [], ORDER = [], N = 0;
const MAT = {};                 // 'reactions'|'metabolites' -> {buf,rowbytes,V,vocab}
const rowByFile = {};
let charts = {};
let GENUS_COL = {};             // genus -> color (top-8 + Other)
const plotDivs = new Set();

/* =============================== BIT HELPERS =============================== */
const POP = new Uint8Array(256);
const BYTE_BITS = [];
for (let b = 0; b < 256; b++) { let c = 0, arr = []; for (let k = 0; k < 8; k++) if (b & (1 << k)) { c++; arr.push(k); } POP[b] = c; BYTE_BITS[b] = arr; }
function countRow(row, m) { const rb = m.rowbytes, base = row * rb, buf = m.buf; let c = 0; for (let i = 0; i < rb; i++) c += POP[buf[base + i]]; return c; }
function interCount(a, b, m) { const rb = m.rowbytes, ba = a * rb, bb = b * rb, buf = m.buf; let c = 0; for (let i = 0; i < rb; i++) c += POP[buf[ba + i] & buf[bb + i]]; return c; }
function has(row, col, m) { return (m.buf[row * m.rowbytes + (col >> 3)] >> (col & 7)) & 1; }
function prevalence(rows, m) {
  const rb = m.rowbytes, buf = m.buf, cnt = new Int32Array(m.V);
  for (let r = 0; r < rows.length; r++) { const base = rows[r] * rb; for (let i = 0; i < rb; i++) { let v = buf[base + i]; if (!v) continue; const bits = BYTE_BITS[v], off = i << 3; for (let k = 0; k < bits.length; k++) cnt[off + bits[k]]++; } }
  return cnt;
}
function colsOfRow(row, m) { const rb = m.rowbytes, base = row * rb, out = []; for (let i = 0; i < rb; i++) { let v = m.buf[base + i]; if (!v) continue; const bits = BYTE_BITS[v], off = i << 3; for (let k = 0; k < bits.length; k++) out.push(off + bits[k]); } return out; }
function meanVecOf(rows, m) { const cnt = prevalence(rows, m), v = new Float32Array(m.V), inv = 1 / rows.length; for (let c = 0; c < m.V; c++) v[c] = cnt[c] * inv; return v; }
function mean(a) { return a.reduce((s, x) => s + x, 0) / a.length; }
function esc(s) { return String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c])); }
function fmt(n) { return typeof n === 'number' ? n.toLocaleString() : n; }
function genusOf(org) { return (org || '').split(' ')[0]; }
function abbr(org) { const g = (org || '').split(' '); return g.length > 1 ? g[0][0] + '. ' + g.slice(1).join(' ') : org; }
function shortLabel(r) { const m = META[r]; return (m.strain ? m.strain : m.gem_file.replace(/\.json(\.json)?$/, '')).slice(0, 30); }

/* =============================== TOOLTIP =============================== */
let TIP;
function tip() { if (!TIP) { TIP = document.createElement('div'); TIP.className = 'vtip'; document.body.appendChild(TIP); } return TIP; }
function showTip(html, x, y) { const t = tip(); t.innerHTML = html; t.style.opacity = 1; t.style.left = (x + 14) + 'px'; t.style.top = (y + 14) + 'px'; }
function hideTip() { if (TIP) TIP.style.opacity = 0; }

/* =============================== LOAD =============================== */
async function load() {
  const msg = document.getElementById('loading-msg'), bar = document.getElementById('load-i');
  const setp = p => { if (bar) bar.style.width = (p * 100) + '%'; };
  msg.textContent = 'Loading model metadata…'; setp(0.08);
  META = await fetch('gems_metadata.json').then(r => r.json());
  const man = await fetch('assets/presence_manifest.json').then(r => r.json());
  ORDER = man.order; N = man.n_models;
  META.forEach((r, i) => { rowByFile[r.gem_file] = i; }); setp(0.2);
  msg.textContent = 'Loading reaction & metabolite matrices…';
  const [rv, mv, rb, mb] = await Promise.all([
    fetch('assets/reactions_vocab.json').then(r => r.json()),
    fetch('assets/metabolites_vocab.json').then(r => r.json()),
    fetch('assets/presence_reactions.bin').then(r => r.arrayBuffer()),
    fetch('assets/presence_metabolites.bin').then(r => r.arrayBuffer()),
  ]); setp(0.85);
  MAT.reactions = { buf: new Uint8Array(rb), rowbytes: man.reaction_rowbytes, V: man.n_reactions, vocab: rv };
  MAT.metabolites = { buf: new Uint8Array(mb), rowbytes: man.metabolite_rowbytes, V: man.n_metabolites, vocab: mv };
  buildGenusPalette(); setp(1);
  document.getElementById('loading').style.display = 'none';
  document.getElementById('app').style.display = 'flex';
  initUI();
}
function organisms() { return [...new Set(META.map(r => r.organism).filter(Boolean))].sort(); }
function speciesMembers(o) { const out = []; for (let i = 0; i < META.length; i++) if (META[i].organism === o) out.push(i); return out; }
function rowsOfDataset(ds) { const out = []; for (let i = 0; i < META.length; i++) if (META[i].dataset === ds) out.push(i); return out; }

function buildGenusPalette() {
  const cnt = {}; META.forEach(r => { const g = genusOf(r.organism); cnt[g] = (cnt[g] || 0) + 1; });
  const top = Object.entries(cnt).sort((a, b) => b[1] - a[1]).slice(0, 8).map(x => x[0]);
  GENUS_COL = {}; top.forEach((g, i) => GENUS_COL[g] = CAT8[i]); GENUS_COL.__other__ = OTHER_COL;
}
function genusColor(org) { const g = genusOf(org); return GENUS_COL[g] || GENUS_COL.__other__; }
function genusLabel(org) { const g = genusOf(org); return GENUS_COL[g] ? g : 'Other'; }

/* =============================== SHARED MATH =============================== */
// PCoA (classical MDS) via power iteration with deflation. D = n×n 2D array of distances.
function pcoa(D, n, k = 2) {
  const A = new Float64Array(n * n);
  const rowm = new Float64Array(n); let gm = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { const v = -0.5 * D[i][j] * D[i][j]; A[i * n + j] = v; rowm[i] += v; }
  for (let i = 0; i < n; i++) { rowm[i] /= n; gm += rowm[i]; } gm /= n;
  const B = A; // reuse
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) B[i * n + j] = B[i * n + j] - rowm[i] - rowm[j] + gm;
  const matvec = (v) => { const w = new Float64Array(n); for (let i = 0; i < n; i++) { let s = 0, base = i * n; for (let j = 0; j < n; j++) s += B[base + j] * v[j]; w[i] = s; } return w; };
  const dot = (a, b) => { let s = 0; for (let i = 0; i < n; i++) s += a[i] * b[i]; return s; };
  const vecs = [], eigs = [];
  for (let e = 0; e < k; e++) {
    let v = new Float64Array(n); for (let i = 0; i < n; i++) v[i] = Math.sin((i + 1) * (e + 1) * 2.399) + 1e-3;
    let nv = Math.sqrt(dot(v, v)); for (let i = 0; i < n; i++) v[i] /= nv;
    for (let it = 0; it < 200; it++) {
      let w = matvec(v);
      for (const pv of vecs) { const d = dot(w, pv); for (let i = 0; i < n; i++) w[i] -= d * pv[i]; }
      const nw = Math.sqrt(dot(w, w)) || 1; for (let i = 0; i < n; i++) w[i] /= nw;
      let diff = 0; for (let i = 0; i < n; i++) diff += Math.abs(Math.abs(w[i]) - Math.abs(v[i]));
      v = w; if (diff < 1e-9) break;
    }
    const lam = dot(v, matvec(v)); eigs.push(lam); vecs.push(v);
  }
  const total = (() => { let t = 0; for (let i = 0; i < n; i++) t += B[i * n + i]; return t; })();
  const x = new Float64Array(n), y = new Float64Array(n);
  const s1 = Math.sqrt(Math.max(eigs[0], 0)), s2 = Math.sqrt(Math.max(eigs[1], 0));
  for (let i = 0; i < n; i++) { x[i] = vecs[0][i] * s1; y[i] = vecs[1][i] * s2; }
  return { x, y, varexp: [eigs[0] / total, eigs[1] / total] };
}

// hierarchical clustering (UPGMA) from a distance matrix -> {order, root} with node heights
function hclustTree(D, n) {
  if (n === 1) return { order: [0], root: { isLeaf: true, idx: 0, height: 0 } };
  const nodes = {}; let active = [];
  for (let i = 0; i < n; i++) { nodes[i] = { isLeaf: true, idx: i, height: 0, size: 1 }; active.push(i); }
  const dd = {}; active.forEach(a => dd[a] = {});
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { dd[i][j] = D[i][j]; dd[j][i] = D[i][j]; }
  let next = n;
  while (active.length > 1) {
    let bi = active[0], bj = active[1], best = Infinity;
    for (let x = 0; x < active.length; x++) for (let y = x + 1; y < active.length; y++) { const a = active[x], b = active[y], d = dd[a][b]; if (d < best) { best = d; bi = a; bj = b; } }
    const na = nodes[bi].size, nb = nodes[bj].size, id = next++;
    nodes[id] = { isLeaf: false, left: nodes[bi], right: nodes[bj], height: best, size: na + nb };
    dd[id] = {};
    active.forEach(kk => { if (kk === bi || kk === bj) return; const d = (na * dd[bi][kk] + nb * dd[bj][kk]) / (na + nb); dd[id][kk] = d; dd[kk][id] = d; });
    active = active.filter(kk => kk !== bi && kk !== bj); active.push(id);
  }
  const root = nodes[active[0]], order = [];
  (function trav(nd) { if (nd.isLeaf) { order.push(nd.idx); return; } trav(nd.left); trav(nd.right); })(root);
  return { order, root };
}

// binary Hamming distance matrix for a sub-matrix (rows over given columns)
function hammingRowsD(rowIdx, cols, m) {
  const R = rowIdx.length, C = cols.length, D = [];
  const bits = []; // pack each row's presence over cols into a Uint8Array
  for (let i = 0; i < R; i++) { const a = new Uint8Array(C); for (let j = 0; j < C; j++) a[j] = has(rowIdx[i], cols[j], m); bits.push(a); }
  for (let i = 0; i < R; i++) { D.push(new Float64Array(R)); }
  for (let i = 0; i < R; i++) for (let j = i + 1; j < R; j++) { let d = 0; const a = bits[i], b = bits[j]; for (let k = 0; k < C; k++) if (a[k] !== b[k]) d++; D[i][j] = d; D[j][i] = d; }
  return { D, bits };
}
// hamming distance matrix over columns of a packed bit set (transpose)
function hammingColsD(bits, R, C) {
  const D = []; for (let j = 0; j < C; j++) D.push(new Float64Array(C));
  for (let a = 0; a < C; a++) for (let b = a + 1; b < C; b++) { let d = 0; for (let i = 0; i < R; i++) if (bits[i][a] !== bits[i][b]) d++; D[a][b] = d; D[b][a] = d; }
  return D;
}

// distance between two mean-presence vectors
function vecDist(a, b, V, kind) {
  if (kind === 'cosine') { let dot = 0, na = 0, nb = 0; for (let c = 0; c < V; c++) { dot += a[c] * b[c]; na += a[c] * a[c]; nb += b[c] * b[c]; } return 1 - dot / (Math.sqrt(na * nb) || 1); }
  // weighted Jaccard (Ruzicka): 1 - sum(min)/sum(max)
  let mn = 0, mx = 0; for (let c = 0; c < V; c++) { const x = a[c], y = b[c]; mn += Math.min(x, y); mx += Math.max(x, y); } return 1 - (mx ? mn / mx : 1);
}

// erfc + chi-square(1df) survival as -log10(p)
function erfc(x) {
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const ans = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? ans : 2 - ans;
}
function chi2NegLog10P(chi2) {
  if (chi2 <= 0) return 0;
  const x = Math.sqrt(chi2 / 2), p = erfc(x);
  if (p > 1e-300) return -Math.log10(p);
  const lnErfc = -x * x - 0.5 * Math.log(Math.PI) - Math.log(x); // asymptotic
  return -lnErfc / Math.LN10;
}
function yatesChi2(a, b, c, d) { // 2x2: a=A+,b=A-,c=B+,d=B-
  const n = a + b + c + d, r1 = a + b, r2 = c + d, c1 = a + c, c2 = b + d;
  const den = r1 * r2 * c1 * c2; if (!den) return 0;
  const num = Math.max(0, Math.abs(a * d - b * c) - n / 2);
  return n * num * num / den;
}

/* =============================== INIT / NAV =============================== */
const shown = {};
const onShow = {};
function initUI() {
  computeGlobals();
  document.getElementById('navToggle').addEventListener('click', () => document.getElementById('sidebar').classList.toggle('open'));
  document.querySelectorAll('.nav-item').forEach(btn => btn.addEventListener('click', () => switchView(btn.dataset.view)));
  window.addEventListener('resize', () => { plotDivs.forEach(id => { const el = document.getElementById(id); if (el && el.offsetParent) try { Plotly.Plots.resize(el); } catch (e) {} }); });
  initOverview(); initLandscape(); initOpenness(); initDistance(); initTree();
  initProps(); initGeo();
  initModels(); initGroups(); initCluster();
  const start = (location.hash || '').replace(/^#/, '');
  switchView(document.getElementById('view-' + start) ? start : 'overview');
  window.addEventListener('hashchange', () => { const n = (location.hash || '').replace(/^#/, ''); if (document.getElementById('view-' + n)) switchView(n); });
}
function switchView(name) {
  document.querySelectorAll('.nav-item').forEach(b => b.classList.toggle('active', b.dataset.view === name));
  document.querySelectorAll('.view').forEach(v => v.classList.toggle('active', v.id === 'view-' + name));
  try { history.replaceState(null, '', '#' + name); } catch (e) {}
  document.getElementById('sidebar').classList.remove('open');
  if (!shown[name]) { shown[name] = true; if (onShow[name]) onShow[name](); }
  // resize plotly in the now-visible view
  requestAnimationFrame(() => plotDivs.forEach(id => { const el = document.getElementById(id); if (el && el.offsetParent) try { Plotly.Plots.resize(el); } catch (e) {} }));
}
/* Shared crisp axis styling — merged into every axis the caller defines, so one edit
   lifts every plot: hairline grid, outside ticks, a readable 11.5px tick font and a
   12px semibold title, with a subtle hover spike line. Caller's title/type/range win. */
const AX = {
  gridcolor: '#EEF2F8', griddash: 'solid', zeroline: true, zerolinecolor: '#D4DCE8', zerolinewidth: 1.4,
  linecolor: '#CBD5E1', linewidth: 1, showline: true, mirror: false,
  ticks: 'outside', tickcolor: '#CBD5E1', ticklen: 4, tickwidth: 1,
  tickfont: { size: 11.5, family: PLOT_FONT.family, color: INK2 },
  titlefont: { size: 12.5, family: PLOT_FONT.family, color: INK1 },
  showspikes: true, spikecolor: '#B7C1D1', spikethickness: 1, spikedash: 'dot', spikemode: 'across', spikesnap: 'cursor',
  automargin: true,
};
function styleAxes(layout) {
  Object.keys(layout).forEach(k => {
    if (/^[xy]axis\d*$/.test(k)) layout[k] = Object.assign({}, AX, layout[k]);
  });
  return layout;
}
function newPlot(id, data, layout, extra) {
  plotDivs.add(id);
  const cfg = { responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'].filter(() => !(extra && extra.keepSelect)), toImageButtonOptions: { format: 'png', scale: 3, filename: 'panGEMs_' + id } };
  if (extra && extra.keepSelect) cfg.modeBarButtonsToRemove = [];
  const base = {
    font: PLOT_FONT, margin: { l: 52, r: 18, t: 14, b: 44 },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    hovermode: 'closest', hoverdistance: 24,
    hoverlabel: { bgcolor: INK1, bordercolor: INK1, font: { color: '#fff', family: PLOT_FONT.family, size: 12 }, align: 'left' },
    colorway: CAT8,
    legend: { bgcolor: 'rgba(255,255,255,.65)', bordercolor: '#E3E8EF', borderwidth: 1, font: { size: 11 } },
  };
  return Plotly.newPlot(id, data, styleAxes(Object.assign(base, layout)), cfg);
}

/* globals: shared reactions/metabolites, per-collection totals */
let G = {};
function computeGlobals() {
  const eco = rowsOfDataset('EcopanGEM'), lac = rowsOfDataset('LactoPanGEM');
  const out = {};
  ['reactions', 'metabolites'].forEach(k => {
    const m = MAT[k], pe = prevalence(eco, m), pl = prevalence(lac, m);
    let e = 0, l = 0, sh = 0; for (let c = 0; c < m.V; c++) { const a = pe[c] > 0, b = pl[c] > 0; if (a) e++; if (b) l++; if (a && b) sh++; }
    out[k] = { eco: e, lac: l, shared: sh, jac: sh / (e + l - sh) };
  });
  G = { eco, lac, r: out.reactions, m: out.metabolites };
}
function stat(v, l, cls) { return `<div class="statbox ${cls || ''}"><div class="v tabular">${fmt(v)}</div><div class="l">${esc(l)}</div></div>`; }

/* =============================== OVERVIEW =============================== */
function initOverview() {
  onShow.overview = () => {}; // rendered eagerly below
  // metrics
  document.getElementById('ov-metrics').innerHTML =
    `<div class="metric hero"><div class="v tabular">4,659</div><div class="l">genome-scale metabolic models · reconstructed, harmonised &amp; comparable</div></div>` +
    metric('2', 'collections', ECO) + metric('27', 'species', LACTO) +
    metric('3,840', 'reactions', ECO) + metric(fmt(G.r.shared), 'shared reactions', '#5a9') +
    metric('2,976', 'metabolites', LACTO);
  // hero embedding (species) + spectrum + species bar + rxn dist
  renderEmbedSpecies('ov-embed', 'jaccard', 'dataset', 'ov-embed-cap');
  document.querySelectorAll('#ov-embed-mode button').forEach(b => b.addEventListener('click', () => {
    seg(b, '#ov-embed-mode'); if (b.dataset.mode === 'species') renderEmbedSpecies('ov-embed', 'jaccard', 'dataset', 'ov-embed-cap');
    else renderEmbedModels('ov-embed', 600, 'jaccard', 'dataset', 'ov-embed-cap');
  }));
  renderSpectrum('ov-spectrum', G.eco, MAT.reactions, 'E. coli', 'ov-spec-cap');
  document.querySelectorAll('#ov-spec-scope button').forEach(b => b.addEventListener('click', () => {
    seg(b, '#ov-spec-scope'); const eco = b.dataset.scope === 'eco';
    renderSpectrum('ov-spectrum', eco ? G.eco : G.lac, MAT.reactions, eco ? 'E. coli' : 'Lactobacillaceae', 'ov-spec-cap');
  }));
  renderSpeciesBar();
  renderRxnDist();
}
function metric(v, l, col) { return `<div class="metric"><div class="accent-bar" style="background:${col}"></div><div class="v tabular">${v}</div><div class="l">${esc(l)}</div></div>`; }
function seg(btn, sel) { document.querySelectorAll(sel + ' button').forEach(b => b.classList.remove('active')); btn.classList.add('active'); }

function renderSpeciesBar() {
  const orgs = organisms().map(o => [o, speciesMembers(o).length]).sort((a, b) => b[1] - a[1]).slice(0, 14).reverse();
  const ds = orgs.map(o => META[speciesMembers(o[0])[0]].dataset);
  const yEco = orgs.map((o, i) => ds[i] === 'EcopanGEM' ? o[1] : null);
  const yLac = orgs.map((o, i) => ds[i] === 'LactoPanGEM' ? o[1] : null);
  const labels = orgs.map(o => abbr(o[0]));
  newPlot('ov-species-bar', [
    { x: yEco, y: labels, type: 'bar', orientation: 'h', name: 'EcopanGEM', marker: { color: ECO }, hovertemplate: '%{y}<br>%{x} models<extra>EcopanGEM</extra>' },
    { x: yLac, y: labels, type: 'bar', orientation: 'h', name: 'LactoPanGEM', marker: { color: LACTO }, hovertemplate: '%{y}<br>%{x} models<extra>LactoPanGEM</extra>' }
  ], {
    barmode: 'stack', height: 340, margin: { l: 128, r: 14, t: 8, b: 34 },
    xaxis: { title: { text: 'models', font: PLOT_FONT }, gridcolor: LINE, zeroline: false },
    yaxis: { tickfont: { size: 10, family: PLOT_FONT.family, color: INK2 }, automargin: true },
    legend: { orientation: 'h', y: 1.08, x: 0, font: { size: 10 } }
  });
}
function renderRxnDist() {
  const m = MAT.reactions;
  const eA = G.eco.map(r => countRow(r, m)), lB = G.lac.map(r => countRow(r, m));
  const all = eA.concat(lB), lo = Math.min(...all), hi = Math.max(...all), w = Math.max(20, Math.round((hi - lo) / 30));
  const start = Math.floor(lo / w) * w, nb = Math.max(1, Math.ceil((hi - start) / w) + 1);
  const labels = [], A = new Array(nb).fill(0), B = new Array(nb).fill(0);
  for (let i = 0; i < nb; i++) labels.push(start + i * w);
  eA.forEach(v => A[Math.min(nb - 1, Math.floor((v - start) / w))]++);
  lB.forEach(v => B[Math.min(nb - 1, Math.floor((v - start) / w))]++);
  if (charts.ovrxn) charts.ovrxn.destroy();
  charts.ovrxn = new Chart(document.getElementById('ov-rxn-canvas'), {
    type: 'bar', data: { labels, datasets: [
      { label: 'EcopanGEM', data: A, backgroundColor: ECO + 'cc', borderRadius: 2 },
      { label: 'LactoPanGEM', data: B, backgroundColor: LACTO + 'cc', borderRadius: 2 }] },
    options: chartOpts('reactions per model', 'models')
  });
}
function chartOpts(xt, yt) {
  return { responsive: true, maintainAspectRatio: false, animation: { duration: 350 },
    plugins: { legend: { labels: { boxWidth: 11, font: { size: 11, family: PLOT_FONT.family }, color: INK2, usePointStyle: true, pointStyle: 'rectRounded' } },
      tooltip: { backgroundColor: INK1, padding: 9, cornerRadius: 8, titleFont: { size: 11 }, bodyFont: { size: 12 } } },
    scales: { x: { grid: { display: false }, title: { display: true, text: xt, font: { size: 11 }, color: INK3 }, ticks: { font: { size: 9 }, color: INK3, maxRotation: 0 } },
      y: { grid: { color: LINE }, border: { display: false }, title: { display: true, text: yt, font: { size: 11 }, color: INK3 }, ticks: { font: { size: 10 }, color: INK3 } } } };
}

/* ---- reaction frequency spectrum ---- */
function renderSpectrum(id, rows, m, label, capId) {
  const cnt = prevalence(rows, m), Ng = rows.length;
  const nbins = 24, bins = new Array(nbins).fill(0), edges = [];
  let core = 0, unique = 0, present = 0;
  for (let c = 0; c < m.V; c++) { const k = cnt[c]; if (k === 0) continue; present++; const frac = k / Ng; const bi = Math.min(nbins - 1, Math.floor(frac * nbins)); bins[bi]++; if (k === Ng) core++; if (k === 1) unique++; }
  const centers = []; for (let i = 0; i < nbins; i++) centers.push(((i + 0.5) / nbins * 100).toFixed(0));
  const colors = bins.map((_, i) => BLUES(0.25 + 0.72 * (i / (nbins - 1))));
  newPlot(id, [{
    x: centers, y: bins, type: 'bar', marker: { color: colors, line: { width: 0 } },
    hovertemplate: '%{x}% of genomes<br>%{y} features<extra></extra>'
  }], {
    height: id === 'ov-spectrum' ? 300 : 330, margin: { l: 52, r: 14, t: 22, b: 42 }, bargap: 0.06,
    xaxis: { title: { text: '% of ' + label + ' genomes carrying the feature', font: PLOT_FONT }, gridcolor: LINE, zeroline: false, tickmode: 'array', tickvals: ['0', '25', '50', '75', '100'] },
    yaxis: { title: { text: 'number of features (log)', font: PLOT_FONT }, type: 'log', gridcolor: LINE, zeroline: false },
    shapes: [
      { type: 'rect', xref: 'paper', yref: 'paper', x0: 0, x1: 0.16, y0: 0, y1: 1, fillcolor: 'rgba(217,119,6,.07)', line: { width: 0 }, layer: 'below' },
      { type: 'rect', xref: 'paper', yref: 'paper', x0: 0.84, x1: 1, y0: 0, y1: 1, fillcolor: 'rgba(37,99,235,.08)', line: { width: 0 }, layer: 'below' }
    ],
    annotations: [
      { x: 0.15, y: 0.97, xref: 'paper', yref: 'paper', text: '◀ <b>cloud</b>', showarrow: false, font: { size: 10.5, color: '#B45309', family: PLOT_FONT.family }, xanchor: 'right', yanchor: 'top' },
      { x: 0.99, y: 0.97, xref: 'paper', yref: 'paper', text: '<b>core</b> ▶', showarrow: false, font: { size: 10.5, color: ECO, family: PLOT_FONT.family }, xanchor: 'right', yanchor: 'top' }
    ]
  });
  if (capId) document.getElementById(capId).innerHTML = `<b>${label}</b> (${fmt(Ng)} genomes): <b>${fmt(core)}</b> core features (in every genome), <b>${fmt(unique)}</b> unique to a single genome, <b>${fmt(present)}</b> in the pan-repertoire. The bimodal U shape is the panreactome signature.`;
}

/* =============================== METABOLIC LANDSCAPE =============================== */
// species-level embedding
let speciesPCoACache = {};
function speciesEmbedding(metric, dist) {
  const key = metric + dist; if (speciesPCoACache[key]) return speciesPCoACache[key];
  const orgs = organisms(), m = MAT[metric], vecs = orgs.map(o => meanVecOf(speciesMembers(o), m));
  const n = orgs.length, D = []; for (let i = 0; i < n; i++) D.push(new Float64Array(n));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = vecDist(vecs[i], vecs[j], m.V, dist === 'cosine' ? 'cosine' : 'jaccard'); D[i][j] = d; D[j][i] = d; }
  const p = pcoa(D, n); const res = { orgs, p, counts: orgs.map(o => speciesMembers(o).length) };
  speciesPCoACache[key] = res; return res;
}
function renderEmbedSpecies(id, dist, colorBy, capId) {
  const { orgs, p, counts } = speciesEmbedding('reactions', dist);
  const maxC = Math.max(...counts);
  const groups = {}; // legend group -> indices
  orgs.forEach((o, i) => { const key = colorBy === 'dataset' ? META[speciesMembers(o)[0]].dataset : genusLabel(o); (groups[key] = groups[key] || []).push(i); });
  const traces = Object.entries(groups).map(([key, idx]) => ({
    x: idx.map(i => p.x[i]), y: idx.map(i => p.y[i]), text: idx.map(i => `<b>${esc(abbr(orgs[i]))}</b><br>${fmt(counts[i])} models · ${esc(genusLabel(orgs[i]))}`),
    customdata: idx.map(i => orgs[i]),
    mode: 'markers', type: 'scatter', name: key,
    marker: { size: idx.map(i => 11 + 30 * Math.sqrt(counts[i] / maxC)), sizemode: 'diameter', color: colorBy === 'dataset' ? (key === 'EcopanGEM' ? ECO : LACTO) : (GENUS_COL[key] || OTHER_COL), opacity: 0.82, line: { color: '#fff', width: 1.5 } },
    hovertemplate: '%{text}<extra></extra>'
  }));
  drawEmbed(id, traces, p, capId, orgs, true);
}
function renderEmbedModels(id, nSample, dist, colorBy, capId) {
  const cap = document.getElementById(capId); if (cap) cap.innerHTML = '<span style="color:var(--ink-3)">Computing PCoA on sampled models…</span>';
  setTimeout(() => {
    const m = MAT.reactions;
    const sample = stratifiedSample(nSample);
    const n = sample.length, D = []; for (let i = 0; i < n; i++) D.push(new Float64Array(n));
    for (let i = 0; i < n; i++) { const ci = countRow(sample[i], m); for (let j = i + 1; j < n; j++) { const inter = interCount(sample[i], sample[j], m); const uni = ci + countRow(sample[j], m) - inter; const d = uni ? 1 - inter / uni : 0; D[i][j] = d; D[j][i] = d; } }
    const p = pcoa(D, n);
    const groups = {}; sample.forEach((r, i) => { const key = colorBy === 'dataset' ? META[r].dataset : genusLabel(META[r].organism); (groups[key] = groups[key] || []).push(i); });
    const traces = Object.entries(groups).map(([key, idx]) => ({
      x: idx.map(i => p.x[i]), y: idx.map(i => p.y[i]),
      customdata: idx.map(i => sample[i]),
      text: idx.map(i => `<b>${esc(META[sample[i]].strain || META[sample[i]].gem_file)}</b><br><i>${esc(META[sample[i]].organism)}</i>`),
      mode: 'markers', type: 'scatter', name: key,
      marker: { size: 6, color: colorBy === 'dataset' ? (key === 'EcopanGEM' ? ECO : LACTO) : (GENUS_COL[key] || OTHER_COL), opacity: 0.6, line: { color: '#fff', width: 0.4 } },
      hovertemplate: '%{text}<extra></extra>'
    }));
    drawEmbed(id, traces, p, capId, null, false, sample.length);
  }, 30);
}
function drawEmbed(id, traces, p, capId, orgs, isSpecies, nModels) {
  // on-plot cloud labels so the two collections read at a glance
  const acc = { eco: { sx: 0, sy: 0, n: 0 }, lac: { sx: 0, sy: 0, n: 0 } };
  traces.forEach(t => { for (let i = 0; i < t.x.length; i++) { const cd = t.customdata[i]; const eco = isSpecies ? FAMILY(cd) === 'Enterobacteriaceae' : META[cd].dataset === 'EcopanGEM'; const g = eco ? acc.eco : acc.lac; g.sx += t.x[i]; g.sy += t.y[i]; g.n++; } });
  const anns = [];
  if (acc.lac.n) anns.push({ x: acc.lac.sx / acc.lac.n, y: acc.lac.sy / acc.lac.n, text: '<b>Lactobacillaceae</b>', showarrow: false, font: { size: 12.5, color: LACTO, family: PLOT_FONT.family }, bgcolor: 'rgba(255,255,255,.72)', bordercolor: 'rgba(217,119,6,.35)', borderpad: 3, borderwidth: 1, opacity: 0.95 });
  if (acc.eco.n) anns.push({ x: acc.eco.sx / acc.eco.n, y: acc.eco.sy / acc.eco.n, ax: -4, ay: -40, text: '<b><i>E. coli</i></b> — metabolic outgroup', showarrow: true, arrowhead: 2, arrowsize: 1, arrowwidth: 1.3, arrowcolor: ECO, font: { size: 11, color: ECO, family: PLOT_FONT.family }, bgcolor: 'rgba(255,255,255,.85)', bordercolor: 'rgba(37,99,235,.4)', borderpad: 3, borderwidth: 1 });
  newPlot(id, traces, {
    height: id === 'ov-embed' ? 430 : 560, margin: { l: 48, r: 14, t: 10, b: 44 },
    xaxis: { title: { text: `PCo1 (${(p.varexp[0] * 100).toFixed(1)}%)`, font: PLOT_FONT }, gridcolor: LINE, zeroline: true, zerolinecolor: '#eef1f5' },
    yaxis: { title: { text: `PCo2 (${(p.varexp[1] * 100).toFixed(1)}%)`, font: PLOT_FONT }, gridcolor: LINE, zeroline: true, zerolinecolor: '#eef1f5' },
    legend: { font: { size: 11 }, itemsizing: 'constant', bgcolor: 'rgba(255,255,255,.6)' },
    showlegend: true, annotations: anns
  }, { keepSelect: true });
  const el = document.getElementById(id);
  el.removeAllListeners && el.removeAllListeners('plotly_click');
  el.on('plotly_click', ev => { const pt = ev.points[0]; if (!pt) return; const cd = pt.data.customdata[pt.pointNumber]; if (isSpecies && typeof cd === 'string') openSpecies(cd); });
  if (capId && isSpecies) document.getElementById(capId).innerHTML = `Look for: the compact Lactobacillaceae neighbourhood vs the distant <b><i>E. coli</i></b> island. Point size ∝ √(model count). Click a species to open its clustermap. First two axes explain <b>${((p.varexp[0] + p.varexp[1]) * 100).toFixed(0)}%</b> of the metabolic variance.`;
  else if (capId) document.getElementById(capId).innerHTML = `Look for: same-species models landing together (metabolic cohesion) and the <b><i>E. coli</i></b> / Lactobacillaceae separation. <b>${fmt(nModels)}</b> models, stratified by species; first two axes explain <b>${((p.varexp[0] + p.varexp[1]) * 100).toFixed(0)}%</b> of variance.`;
}
function stratifiedSample(target) {
  const orgs = organisms(), total = META.length, out = [];
  orgs.forEach(o => { const mem = speciesMembers(o); const k = Math.max(1, Math.round(target * mem.length / total)); const step = mem.length / k; for (let i = 0; i < k; i++) out.push(mem[Math.floor(i * step)]); });
  return out;
}
function initLandscape() {
  onShow.landscape = () => runLandscape();
  document.querySelectorAll('#ls-unit button').forEach(b => b.addEventListener('click', () => { seg(b, '#ls-unit'); document.getElementById('ls-n-wrap').style.display = b.dataset.unit === 'models' ? 'flex' : 'none'; }));
  ['#ls-color', '#ls-dist'].forEach(sel => document.querySelectorAll(sel + ' button').forEach(b => b.addEventListener('click', () => seg(b, sel))));
  document.getElementById('ls-run').addEventListener('click', runLandscape);
}
function runLandscape() {
  const unit = document.querySelector('#ls-unit button.active').dataset.unit;
  const colorBy = document.querySelector('#ls-color button.active').dataset.c;
  const dist = document.querySelector('#ls-dist button.active').dataset.d;
  renderLandscapeLegend(colorBy);
  if (unit === 'species') renderEmbedSpecies('ls-plot', dist, colorBy, 'ls-cap');
  else renderEmbedModels('ls-plot', +document.getElementById('ls-n').value, dist, colorBy, 'ls-cap');
}
function renderLandscapeLegend(colorBy) {
  const el = document.getElementById('ls-legend');
  if (colorBy === 'dataset') el.innerHTML = `<span class="it"><span class="dot" style="background:${ECO}"></span>EcopanGEM (E. coli)</span><span class="it"><span class="dot" style="background:${LACTO}"></span>LactoPanGEM (Lactobacillaceae)</span>`;
  else el.innerHTML = Object.keys(GENUS_COL).filter(g => g !== '__other__').map(g => `<span class="it"><span class="dot" style="background:${GENUS_COL[g]}"></span>${esc(g)}</span>`).join('') + `<span class="it"><span class="dot" style="background:${OTHER_COL}"></span>Other</span>`;
}

/* =============================== PANGENOME OPENNESS =============================== */
function initOpenness() {
  const sel = document.getElementById('op-scope');
  sel.innerHTML = `<option value="ds:EcopanGEM">E. coli — all (2,313)</option><option value="ds:LactoPanGEM">Lactobacillaceae — all (2,346)</option>` +
    organisms().map(o => `<option value="sp:${esc(o)}">${esc(o)} (${speciesMembers(o).length})</option>`).join('');
  document.querySelectorAll('#op-metric button').forEach(b => b.addEventListener('click', () => seg(b, '#op-metric')));
  document.getElementById('op-run').addEventListener('click', runOpenness);
  onShow.openness = () => runOpenness();
}
function opScopeRows() {
  const v = document.getElementById('op-scope').value;
  if (v.startsWith('ds:')) return { rows: rowsOfDataset(v.slice(3)), label: v.slice(3) === 'EcopanGEM' ? 'E. coli' : 'Lactobacillaceae' };
  const o = v.slice(3); return { rows: speciesMembers(o), label: o };
}
function runOpenness() {
  const { rows, label } = opScopeRows();
  const metric = document.querySelector('#op-metric button.active').dataset.m;
  const m = MAT[metric];
  // spectrum
  renderSpectrum('op-spectrum', rows, m, label, 'op-spec-cap');
  // stats
  const cnt = prevalence(rows, m), Ng = rows.length; let core = 0, unique = 0, pan = 0, soft = 0;
  for (let c = 0; c < m.V; c++) { const k = cnt[c]; if (k === 0) continue; pan++; if (k === Ng) core++; else if (k === 1) unique++; if (k >= 0.95 * Ng) soft++; }
  // rarefaction
  const rar = rarefaction(rows, m, 25, 250);
  document.getElementById('op-stats').innerHTML =
    stat(Ng, 'genomes') + stat(pan, 'pan (' + metric + ')') + stat(core, 'strict core') + stat(soft, 'soft core (≥95%)') +
    stat(unique, 'unique (n=1)') + `<div class="statbox"><div class="v tabular">${rar.gamma.toFixed(3)}</div><div class="l">Heaps' γ (openness)</div></div>`;
  renderRarefaction('op-rarefaction', rar, metric, 'op-rare-cap');
}
function rarefaction(rows, m, perms, cap) {
  const rb = m.rowbytes, buf = m.buf;
  const pool = rows.slice(); const Ncap = Math.min(pool.length, cap);
  const panSum = new Float64Array(Ncap), coreSum = new Float64Array(Ncap);
  for (let p = 0; p < perms; p++) {
    // Fisher-Yates partial shuffle for first Ncap
    for (let i = 0; i < Ncap; i++) { const j = i + Math.floor(Math.random() * (pool.length - i)); const t = pool[i]; pool[i] = pool[j]; pool[j] = t; }
    const uni = new Uint8Array(rb), core = new Uint8Array(rb); core.fill(0xff);
    for (let i = 0; i < Ncap; i++) { const base = pool[i] * rb; let uc = 0, cc = 0; for (let b = 0; b < rb; b++) { uni[b] |= buf[base + b]; core[b] &= buf[base + b]; uc += POP[uni[b]]; cc += POP[core[b]]; } panSum[i] += uc; coreSum[i] += cc; }
  }
  const Ns = [], pan = [], corev = [];
  for (let i = 0; i < Ncap; i++) { Ns.push(i + 1); pan.push(panSum[i] / perms); corev.push(coreSum[i] / perms); }
  // Heaps' law fit pan = k * N^gamma on log-log (skip first 3)
  let sx = 0, sy = 0, sxx = 0, sxy = 0, nn = 0;
  for (let i = 3; i < Ncap; i++) { const lx = Math.log(Ns[i]), ly = Math.log(pan[i]); sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly; nn++; }
  const gamma = nn > 1 ? (nn * sxy - sx * sy) / (nn * sxx - sx * sx) : 0;
  const lk = nn > 1 ? (sy - gamma * sx) / nn : 0; const kappa = Math.exp(lk);
  return { Ns, pan, core: corev, gamma, kappa, perms, Ncap, total: rows.length };
}
function renderRarefaction(id, rar, metric, capId) {
  const fitY = rar.Ns.map(n => rar.kappa * Math.pow(n, rar.gamma));
  newPlot(id, [
    { x: rar.Ns, y: rar.pan, mode: 'lines', name: 'Pan (union)', line: { color: ECO, width: 2.5 }, hovertemplate: 'N=%{x}<br>pan %{y:.0f}<extra></extra>' },
    { x: rar.Ns, y: fitY, mode: 'lines', name: `Heaps' fit (γ=${rar.gamma.toFixed(3)})`, line: { color: ECO, width: 1, dash: 'dot' }, hoverinfo: 'skip' },
    { x: rar.Ns, y: rar.core, mode: 'lines', name: 'Core (intersection)', line: { color: LACTO, width: 2.5 }, hovertemplate: 'N=%{x}<br>core %{y:.0f}<extra></extra>' }
  ], {
    height: 330, margin: { l: 52, r: 14, t: 8, b: 42 },
    xaxis: { title: { text: 'genomes sampled (N)', font: PLOT_FONT }, gridcolor: LINE, zeroline: false },
    yaxis: { title: { text: metric + ' (log)', font: PLOT_FONT }, type: 'log', gridcolor: LINE, zeroline: false },
    legend: { orientation: 'h', y: 1.12, x: 0, font: { size: 10 } }
  });
  const open = rar.gamma > 0.05 ? 'an <b>open</b>' : 'a <b>closed</b>';
  document.getElementById(capId).innerHTML = `Averaged over <b>${rar.perms}</b> random permutations of up to <b>${rar.Ncap}</b> genomes${rar.total > rar.Ncap ? ' (of ' + fmt(rar.total) + ')' : ''}. Heaps' law γ&nbsp;=&nbsp;<b>${rar.gamma.toFixed(3)}</b> → ${open} metabolic panreactome (pan keeps growing as N^γ; core plateaus).`;
}

/* =============================== SPECIES DISTANCE HEATMAP =============================== */
function initDistance() {
  ['#dm-metric', '#dm-dist'].forEach(sel => document.querySelectorAll(sel + ' button').forEach(b => b.addEventListener('click', () => seg(b, sel))));
  document.getElementById('dm-run').addEventListener('click', runDistance);
  onShow.distance = () => runDistance();
}
function runDistance() {
  const metric = document.querySelector('#dm-metric button.active').dataset.m;
  const dist = document.querySelector('#dm-dist button.active').dataset.d;
  const orgs = organisms(), m = MAT[metric];
  const vecs = orgs.map(o => meanVecOf(speciesMembers(o), m));
  const n = orgs.length, D = []; for (let i = 0; i < n; i++) D.push(new Float64Array(n));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = vecDist(vecs[i], vecs[j], m.V, dist === 'cosine' ? 'cosine' : 'jaccard'); D[i][j] = d; D[j][i] = d; }
  const { order } = hclustTree(D, n);
  const labels = order.map(i => abbr(orgs[i]));
  const z = order.map(i => order.map(j => D[i][j]));
  const isEco = orgs.map(o => META[speciesMembers(o)[0]].dataset === 'EcopanGEM');
  // Dynamic range: one lone E. coli row/col (~0.7) otherwise swamps the whole
  // Lactobacillaceae block (~0.1), flattening it. Cap the colour scale at the 97th
  // percentile of the within-family (non-outgroup) distances so that block gets the
  // full ramp; E. coli then reads as a saturated off-scale outgroup, as it should.
  const sameFam = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (isEco[i] === isEco[j]) sameFam.push(D[i][j]);
  sameFam.sort((a, b) => a - b);
  const cap = sameFam.length ? sameFam[Math.min(sameFam.length - 1, Math.floor(0.97 * sameFam.length))] : Math.max(...z.flat());
  const capR = cap < 0.1 ? cap.toFixed(3) : cap.toFixed(2);
  // find E. coli display position for the outgroup annotation
  const ecoPos = order.findIndex(i => isEco[i]);
  const tickText = order.map(i => { const ds = META[speciesMembers(orgs[i])[0]].dataset; const col = ds === 'EcopanGEM' ? ECO : LACTO; return `<span style="color:${col}">●</span> ${esc(abbr(orgs[i]))}`; });
  const unitTxt = dist === 'cosine' ? '1 − cosine' : 'weighted Jaccard';
  newPlot('dm-plot', [{
    z, x: labels, y: labels, type: 'heatmap',
    colorscale: [[0, '#FCFBFF'], [0.15, '#EDE6FB'], [0.4, '#C9B3F2'], [0.7, '#9B72E6'], [1, '#5B21B6']],
    zmin: 0, zmax: cap, zsmooth: false,
    colorbar: { title: { text: unitTxt + ' &nbsp;(capped)', side: 'right', font: { size: 11, family: PLOT_FONT.family, color: INK2 } }, thickness: 13, len: 0.6, outlinewidth: 0, tickfont: { size: 9.5, color: INK2 }, ticks: 'outside', ticklen: 3, tickcolor: '#CBD5E1' },
    xgap: 1.4, ygap: 1.4,
    hovertemplate: '<b>%{y}</b> ↔ <b>%{x}</b><br>' + unitTxt + ' = %{z:.3f}<extra></extra>'
  }], {
    height: 660, margin: { l: 158, r: 24, t: 30, b: 158 },
    xaxis: { tickangle: -55, tickfont: { size: 10, family: PLOT_FONT.family, color: INK2 }, ticktext: tickText, tickvals: labels, automargin: true, showgrid: false, showspikes: false, zeroline: false, showline: false, ticks: '' },
    yaxis: { autorange: 'reversed', tickfont: { size: 10, family: PLOT_FONT.family, color: INK2 }, ticktext: tickText, tickvals: labels, automargin: true, showgrid: false, showspikes: false, zeroline: false, showline: false, ticks: '' },
    annotations: (ecoPos >= 0 ? [{
      x: labels[ecoPos], y: labels[Math.min(ecoPos + 2, n - 1)], xref: 'x', yref: 'y',
      text: '<b>E. coli</b><br>outgroup<br>(off-scale)', showarrow: true, arrowhead: 2, arrowsize: 0.9, arrowwidth: 1.4, arrowcolor: '#5B21B6', ax: 46, ay: 26,
      font: { size: 9.5, color: '#5B21B6', family: PLOT_FONT.family }, align: 'center', bgcolor: 'rgba(255,255,255,.9)', bordercolor: '#C9B3F2', borderpad: 3, borderwidth: 1
    }] : []).concat([{
      x: 0, y: 1.045, xref: 'paper', yref: 'paper', text: 'colour scale capped at the 97th-percentile within-family distance (' + capR + ') so the Lactobacillaceae block is legible', showarrow: false, font: { size: 10, color: INK3, family: PLOT_FONT.family }, align: 'left', xanchor: 'left'
    }])
  });
  // computed interpretation: within-Lacto cohesion vs the E. coli outgroup gap
  let ll = 0, lln = 0, el = 0, eln = 0, best = Infinity, bi = 0, bj = 1;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const d = D[i][j];
    if (!isEco[i] && !isEco[j]) { ll += d; lln++; if (d < best) { best = d; bi = i; bj = j; } }
    else if (isEco[i] !== isEco[j]) { el += d; eln++; }
  }
  const mLL = lln ? ll / lln : 0, mEL = eln ? el / eln : 0, unit = dist === 'cosine' ? '1−cosine' : 'weighted-Jaccard';
  document.getElementById('dm-interp').innerHTML =
    ico() + 'Within the <span class="lm">Lactobacillaceae</span>, species sit a mean <b>' + mLL.toFixed(3) + '</b> ' + unit + ' apart — a metabolically coherent block. '
    + '<span class="em">E. coli</span> sits <b>' + mEL.toFixed(3) + '</b> from them on average, <b>' + (mLL > 0 ? (mEL / mLL).toFixed(1) : '∞') + '×</b> farther than the lactobacilli are from each other: the bright cross that marks it as the metabolic outgroup. '
    + 'The most metabolically similar pair here is <b><i>' + esc(abbr(orgs[bi])) + '</i></b> and <b><i>' + esc(abbr(orgs[bj])) + '</i></b> (' + best.toFixed(3) + '). '
    + 'Blocks of low distance are clades that could share media and engineering strategies.';
}

/* =============================== SPECIES TREE (polished) =============================== */
let tMode = 'taxonomy';
const FAMILY = org => org === 'Escherichia coli' ? 'Enterobacteriaceae' : 'Lactobacillaceae';
function initTree() {
  document.querySelectorAll('#t-mode button').forEach(b => b.addEventListener('click', () => { seg(b, '#t-mode'); tMode = b.dataset.mode; drawTree(); }));
  onShow.tree = () => drawTree();
}
function taxonomyRoot() {
  const orgs = organisms(), fam = {};
  orgs.forEach(o => { const f = FAMILY(o), g = genusOf(o); (fam[f] = fam[f] || {})[g] = fam[f][g] || []; fam[f][g].push(o); });
  return { name: 'Bacteria', children: Object.entries(fam).map(([f, gs]) => ({ name: f, children: Object.entries(gs).map(([g, ss]) => ({ name: g, children: ss.map(s => ({ name: s, leaf: true, count: speciesMembers(s).length, organism: s })) })) })) };
}
function contentRoot() {
  const orgs = organisms(), m = MAT.reactions, vecs = orgs.map(o => meanVecOf(speciesMembers(o), m));
  const n = orgs.length, D = []; for (let i = 0; i < n; i++) D.push(new Float64Array(n));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = vecDist(vecs[i], vecs[j], m.V, 'cosine'); D[i][j] = d; D[j][i] = d; }
  // convert hclustTree to d3 hierarchy
  const { root } = hclustTree(D, n);
  const conv = nd => nd.isLeaf ? { name: orgs[nd.idx], organism: orgs[nd.idx], leaf: true, count: speciesMembers(orgs[nd.idx]).length } : { name: '', children: [conv(nd.left), conv(nd.right)], height: nd.height };
  return conv(root);
}
function drawTree() {
  const root = tMode === 'taxonomy' ? taxonomyRoot() : contentRoot();
  const holder = document.getElementById('tree-holder'); holder.innerHTML = '';
  const hier = d3.hierarchy(root); const leaves = hier.leaves().length;
  const W = Math.max(760, holder.clientWidth || 880), rowH = 26, H = leaves * rowH + 40;
  const svg = d3.select(holder).append('svg').attr('width', W).attr('height', H).style('font-family', PLOT_FONT.family);
  const layout = d3.cluster().size([H - 30, W - 470]); layout(hier);
  const g = svg.append('g').attr('transform', 'translate(150,15)');
  g.selectAll('path.link').data(hier.links()).join('path').attr('fill', 'none').attr('stroke', '#d3dbe6').attr('stroke-width', 1.2)
    .attr('d', d => `M${d.source.y},${d.source.x} C${(d.source.y + d.target.y) / 2},${d.source.x} ${(d.source.y + d.target.y) / 2},${d.target.x} ${d.target.y},${d.target.x}`);
  const node = g.selectAll('g.node').data(hier.descendants()).join('g').attr('transform', d => `translate(${d.y},${d.x})`);
  node.append('circle').attr('r', d => d.data.leaf ? 4.5 : 3).attr('fill', d => !d.data.leaf ? '#aab6c6' : (FAMILY(d.data.organism) === 'Enterobacteriaceae' ? ECO : LACTO)).attr('stroke', '#fff').attr('stroke-width', 1);
  node.append('text').attr('dy', '0.32em').attr('x', d => d.data.leaf ? 9 : -6).attr('text-anchor', d => d.data.leaf ? 'start' : 'end')
    .style('font-size', d => d.data.leaf ? '11.5px' : '10px').style('font-style', d => d.data.leaf ? 'italic' : 'normal')
    .style('fill', d => d.data.leaf ? INK1 : INK3).style('font-weight', d => d.data.leaf ? '500' : '600')
    .text(d => d.data.leaf ? `${d.data.name} (${d.data.count})` : d.data.name);
  node.filter(d => d.data.leaf).style('cursor', 'pointer')
    .on('mousemove', (ev, d) => showTip(`<b><i>${esc(d.data.name)}</i></b><br><span class="k">${d.data.count} models · click to open in clustermap</span>`, ev.clientX, ev.clientY))
    .on('mouseleave', hideTip).on('click', (ev, d) => openSpecies(d.data.organism));
  d3.select(holder).append('div').attr('class', 'interp').style('margin-top', '12px')
    .html(ico() + (tMode === 'taxonomy'
      ? 'This is the <b>NCBI taxonomy</b> — family → genus → species, a fixed reference, not inferred from the models. Circle colour marks the collection. Switch to <b>Metabolic content</b> to see whether reaction repertoires alone rebuild this same shape.'
      : 'This dendrogram is built <b>only from metabolism</b> — mean reaction presence per species, 1−cosine distance, average linkage. It independently recovers the <span class="lm">Lactobacillaceae</span> genera and isolates <span class="em">E. coli</span> as a lone outgroup, so metabolic content carries the same signal as the taxonomy above: closely related species really do build similar networks. Click any leaf to open its clustermap.'));
}

/* =============================== COMPARE MODELS (+ UpSet) =============================== */
let mSel = [], mMetric = 'reactions';
function initModels() {
  const dl = document.getElementById('m-list');
  META.forEach(r => { const o = document.createElement('option'); o.value = r.gem_file; o.label = `${r.organism || ''} ${r.strain || ''}`.trim(); dl.appendChild(o); });
  document.getElementById('m-add').addEventListener('click', addModel);
  document.getElementById('m-search').addEventListener('keydown', e => { if (e.key === 'Enter') addModel(); });
  document.getElementById('m-clear').addEventListener('click', () => { mSel = []; renderModels(); });
  document.getElementById('m-examples').addEventListener('click', () => {
    const pick = []; const eco = rowsOfDataset('EcopanGEM'), lac = speciesMembers('Lactobacillus gasseri');
    [eco[0], eco[1], lac[0], lac[1]].forEach(r => r != null && pick.push(r)); mSel = [...new Set(pick)].slice(0, 4); renderModels();
  });
  document.querySelectorAll('#m-metric button').forEach(b => b.addEventListener('click', () => { seg(b, '#m-metric'); mMetric = b.dataset.metric; document.querySelectorAll('.m-metric-label').forEach(e => e.textContent = mMetric); renderModels(); }));
  document.getElementById('m-csv').addEventListener('click', downloadDiffCSV);
  onShow.models = () => {};
}
function addModel() {
  const q = document.getElementById('m-search').value.trim(); if (!q) return;
  let row = rowByFile[q];
  if (row == null) { const ql = q.toLowerCase(); const hit = META.find(r => (r.gem_file + ' ' + r.organism + ' ' + r.strain).toLowerCase().includes(ql)); if (hit) row = rowByFile[hit.gem_file]; }
  if (row == null) { alert('No model matches "' + q + '"'); return; }
  if (!mSel.includes(row)) mSel.push(row);
  if (mSel.length > 12) mSel = mSel.slice(-12);
  document.getElementById('m-search').value = ''; renderModels();
}
function renderModels() {
  const chips = document.getElementById('m-chips');
  chips.innerHTML = mSel.map(r => { const m = META[r], c = m.dataset === 'EcopanGEM' ? ECO : LACTO; return `<span class="chip" style="background:${c}18;color:${c}"><i>${esc(abbr(m.organism))}</i> ${esc(m.strain || m.gem_file)} <span class="x" data-r="${r}">×</span></span>`; }).join('');
  chips.querySelectorAll('.x').forEach(x => x.addEventListener('click', () => { mSel = mSel.filter(r => r != x.dataset.r); renderModels(); }));
  const res = document.getElementById('m-results');
  if (mSel.length < 2) { res.style.display = 'none'; return; }
  res.style.display = 'block';
  const m = MAT[mMetric], sets = mSel.map(r => new Set(colsOfRow(r, m)));
  const union = new Set(); sets.forEach(s => s.forEach(x => union.add(x)));
  const core = [...union].filter(c => sets.every(s => s.has(c)));
  document.getElementById('m-stats').innerHTML = stat(mSel.length, 'Models') + stat(union.size, 'Union (pan)') + stat(core.length, 'Core (in all)') + stat(union.size - core.length, 'Variable');
  // pairwise Jaccard
  let jh = '<table class="mini" style="width:auto"><tr><th></th>' + mSel.map(r => `<th title="${esc(META[r].organism)}">${esc(shortLabel(r))}</th>`).join('') + '</tr>';
  for (let i = 0; i < mSel.length; i++) { jh += `<tr><th style="text-align:right">${esc(shortLabel(mSel[i]))}</th>`;
    for (let j = 0; j < mSel.length; j++) { const inter = interCount(mSel[i], mSel[j], m), uni = sets[i].size + sets[j].size - inter, J = uni ? inter / uni : 1; const bg = BLUES(0.1 + 0.85 * J); jh += `<td style="background:${bg};text-align:center;color:${J > 0.55 ? '#fff' : INK1};font-variant-numeric:tabular-nums" title="${(J * 100).toFixed(1)}% shared">${(J * 100).toFixed(0)}</td>`; } jh += '</tr>'; }
  jh += '</table>'; document.getElementById('m-jaccard').innerHTML = jh;
  // shared vs unique bar
  const uniqueOf = mSel.map((r, i) => [...sets[i]].filter(c => sets.every((s, j) => j === i || !s.has(c))).length);
  if (charts.mset) charts.mset.destroy();
  charts.mset = new Chart(document.getElementById('m-setchart'), {
    type: 'bar', data: { labels: mSel.map(r => shortLabel(r)), datasets: [
      { label: 'Shared with ≥1 other', data: mSel.map((r, i) => sets[i].size - uniqueOf[i]), backgroundColor: '#9ec3ea', borderRadius: 2 },
      { label: 'Unique to this model', data: uniqueOf, backgroundColor: LACTO, borderRadius: 2 }] },
    options: Object.assign(chartOpts('', mMetric), { scales: Object.assign(chartOpts('', mMetric).scales, { x: { stacked: true, grid: { display: false }, ticks: { font: { size: 9 }, color: INK3 } }, y: { stacked: true, grid: { color: LINE }, border: { display: false }, title: { display: true, text: mMetric, font: { size: 11 }, color: INK3 } } }) })
  });
  drawUpSet(sets, union);
  renderDiffTable(sets, union, core);
}
function drawUpSet(sets, union) {
  const K = Math.min(sets.length, 6);
  const card = document.getElementById('m-upset-card');
  const note = sets.length > 6 ? `<div class="sub" style="color:var(--lacto)">Showing the first 6 of ${sets.length} models (UpSet caps at 6).</div>` : '';
  // tally by membership bitmask over first K sets
  const tally = new Map();
  union.forEach(c => { let mask = 0; for (let i = 0; i < K; i++) if (sets[i].has(c)) mask |= (1 << i); if (mask) tally.set(mask, (tally.get(mask) || 0) + 1); });
  let inters = [...tally.entries()].sort((a, b) => b[1] - a[1]).slice(0, 16);
  const setSize = []; for (let i = 0; i < K; i++) setSize.push(sets[i].size);
  const maxInter = Math.max(...inters.map(x => x[1])), maxSet = Math.max(...setSize);
  // geometry
  const padL = 150, cellW = 34, barMaxH = 130, gap = 8, matrixRowH = 22, dotR = 6;
  const W = padL + inters.length * cellW + 20, matrixH = K * matrixRowH + 8;
  const H = barMaxH + 26 + matrixH + 10;
  let svg = `<svg width="${W}" height="${H}" font-family="${PLOT_FONT.family}">`;
  // intersection bars
  inters.forEach(([mask, sz], ci) => {
    const x = padL + ci * cellW + cellW / 2, h = maxInter ? (sz / maxInter) * barMaxH : 0, y = barMaxH - h;
    const single = (mask & (mask - 1)) === 0;
    svg += `<rect x="${x - 10}" y="${y}" width="20" height="${h}" rx="3" fill="${single ? LACTO : ECO}"><title>${sz} ${mMetric} present in exactly this set</title></rect>`;
    svg += `<text x="${x}" y="${y - 4}" text-anchor="middle" font-size="10" fill="${INK2}" font-variant-numeric="tabular-nums">${sz}</text>`;
  });
  svg += `<text x="6" y="14" font-size="10" fill="${INK3}">Intersection size</text>`;
  svg += `<line x1="${padL}" y1="${barMaxH + 1}" x2="${W - 16}" y2="${barMaxH + 1}" stroke="${LINE}"/>`;
  // matrix
  const my0 = barMaxH + 26;
  for (let i = 0; i < K; i++) {
    const cy = my0 + i * matrixRowH + matrixRowH / 2;
    if (i % 2 === 0) svg += `<rect x="${padL - 4}" y="${cy - matrixRowH / 2}" width="${inters.length * cellW + 8}" height="${matrixRowH}" fill="#f7f9fc"/>`;
    const lab = shortLabel(mSel[i]); const col = META[mSel[i]].dataset === 'EcopanGEM' ? ECO : LACTO;
    svg += `<circle cx="8" cy="${cy}" r="4" fill="${col}"/>`;
    svg += `<text x="18" y="${cy}" dy="0.32em" font-size="10.5" fill="${INK1}">${esc(lab.slice(0, 20))}</text>`;
    // set-size mini bar
    const sw = maxSet ? (setSize[i] / maxSet) * 40 : 0;
    svg += `<rect x="${padL - 50}" y="${cy - 4}" width="${sw}" height="8" rx="2" fill="#c9d6e8"><title>${setSize[i]} total</title></rect>`;
  }
  inters.forEach(([mask, sz], ci) => {
    const x = padL + ci * cellW + cellW / 2;
    // connect line for filled members
    const filled = []; for (let i = 0; i < K; i++) if (mask & (1 << i)) filled.push(i);
    if (filled.length > 1) svg += `<line x1="${x}" y1="${my0 + filled[0] * matrixRowH + matrixRowH / 2}" x2="${x}" y2="${my0 + filled[filled.length - 1] * matrixRowH + matrixRowH / 2}" stroke="${INK1}" stroke-width="2"/>`;
    for (let i = 0; i < K; i++) { const cy = my0 + i * matrixRowH + matrixRowH / 2; const on = mask & (1 << i); svg += `<circle cx="${x}" cy="${cy}" r="${dotR}" fill="${on ? INK1 : '#dfe5ee'}"/>`; }
  });
  svg += '</svg>';
  document.getElementById('m-upset').innerHTML = note + svg;
}
let mDiffRows = [];
function renderDiffTable(sets, union, core) {
  const m = MAT[mMetric], coreSet = new Set(core);
  const diff = [...union].filter(c => !coreSet.has(c)).sort((a, b) => (m.vocab[a].id < m.vocab[b].id ? -1 : 1));
  mDiffRows = diff; const cap = 1200;
  let h = '<thead><tr><th>ID</th><th>Name</th>' + mSel.map(r => `<th title="${esc(META[r].organism)}">${esc(shortLabel(r))}</th>`).join('') + '</tr></thead><tbody>';
  diff.slice(0, cap).forEach(c => { h += `<tr><td class="tabular">${esc(m.vocab[c].id)}</td><td>${esc((m.vocab[c].name || '').slice(0, 44))}</td>` + sets.map(s => s.has(c) ? '<td class="present">✓</td>' : '<td class="absent">·</td>').join('') + '</tr>'; });
  h += '</tbody>'; document.getElementById('m-diff').innerHTML = h;
}
function downloadDiffCSV() {
  if (!mSel.length) return; const m = MAT[mMetric], sets = mSel.map(r => new Set(colsOfRow(r, m)));
  let csv = 'id,name,' + mSel.map(r => '"' + META[r].gem_file + '"').join(',') + '\n';
  mDiffRows.forEach(c => { csv += `"${m.vocab[c].id}","${(m.vocab[c].name || '').replace(/"/g, '""')}",` + sets.map(s => s.has(c) ? 1 : 0).join(',') + '\n'; });
  dl(csv, `panGEMs_compare_${mMetric}_${mSel.length}models.csv`);
}
function dl(text, name) { const b = new Blob([text], { type: 'text/csv' }); const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = name; a.click(); URL.revokeObjectURL(a.href); }

/* =============================== COMPARE GROUPS (+ volcano) =============================== */
let gMetric = 'reactions', groupCache = { A: [], B: [] };
function groupBuilder(id, letter, defaults) {
  const orgs = organisms(), countries = [...new Set(META.map(r => r.country).filter(Boolean))].sort(), hosts = [...new Set(META.map(r => r.host_name).filter(Boolean))].sort();
  const el = document.getElementById(id); const col = letter === 'A' ? ECO : LACTO;
  el.className = 'card2'; el.style.borderTop = `3px solid ${col}`; el.style.margin = 0;
  el.innerHTML = `<div style="font-weight:750;color:${col};margin-bottom:8px;font-size:14px">Group ${letter}</div>
    <div class="ctl-row" style="gap:8px">
      <div class="ctl" style="flex:1"><label>Dataset</label><select data-f="dataset"><option value="">All</option><option>EcopanGEM</option><option>LactoPanGEM</option></select></div>
      <div class="ctl" style="flex:1"><label>Country</label><select data-f="country"><option value="">All</option>${countries.map(c => `<option>${esc(c)}</option>`).join('')}</select></div>
    </div>
    <div class="ctl" style="margin-top:8px"><label>Species (pick one or more; empty = all)</label><select multiple data-f="organism" style="width:100%">${orgs.map(o => `<option>${esc(o)}</option>`).join('')}</select></div>
    <div class="ctl" style="margin-top:8px"><label>Host</label><select data-f="host_name" style="width:100%"><option value="">All</option>${hosts.map(h => `<option>${esc(h)}</option>`).join('')}</select></div>
    <div style="margin-top:8px;font-size:12.5px;color:var(--ink-2)">Matches: <b class="gcount tabular">0</b> models</div>`;
  const resolve = () => {
    const ds = el.querySelector('[data-f=dataset]').value, ct = el.querySelector('[data-f=country]').value, hn = el.querySelector('[data-f=host_name]').value;
    const orgSel = [...el.querySelector('[data-f=organism]').selectedOptions].map(o => o.value);
    const rows = [];
    for (let i = 0; i < META.length; i++) { const r = META[i]; if (ds && r.dataset !== ds) continue; if (ct && r.country !== ct) continue; if (hn && r.host_name !== hn) continue; if (orgSel.length && !orgSel.includes(r.organism)) continue; rows.push(i); }
    el.querySelector('.gcount').textContent = rows.length.toLocaleString(); return rows;
  };
  el.addEventListener('change', () => { groupCache[letter] = resolve(); });
  if (defaults) { if (defaults.dataset) el.querySelector('[data-f=dataset]').value = defaults.dataset; if (defaults.organism) [...el.querySelector('[data-f=organism]').options].forEach(o => { if (defaults.organism.includes(o.value)) o.selected = true; }); }
  groupCache[letter] = resolve(); return resolve;
}
function initGroups() {
  const rA = groupBuilder('grpA', 'A', { dataset: 'EcopanGEM' });
  const rB = groupBuilder('grpB', 'B', { organism: ['Lactiplantibacillus plantarum'] });
  document.querySelectorAll('#g-metric button').forEach(b => b.addEventListener('click', () => { seg(b, '#g-metric'); gMetric = b.dataset.metric; document.querySelectorAll('.g-metric-label').forEach(e => e.textContent = gMetric[0].toUpperCase() + gMetric.slice(1)); }));
  document.getElementById('g-run').addEventListener('click', () => runGroups(rA(), rB()));
  document.getElementById('g-csv').addEventListener('click', gDownloadCSV);
  onShow.groups = () => {};
}
let gDiffData = null;
function runGroups(A, B) {
  if (!A.length || !B.length) { alert('Both groups need at least one model.'); return; }
  const m = MAT[gMetric], cntA = prevalence(A, m), cntB = prevalence(B, m);
  const perA = A.map(r => countRow(r, m)), perB = B.map(r => countRow(r, m));
  let coreA = 0, coreB = 0, shared = 0, Aonly = 0, Bonly = 0, bothCore = 0;
  for (let c = 0; c < m.V; c++) { const a = cntA[c] > 0, b = cntB[c] > 0; if (cntA[c] === A.length) coreA++; if (cntB[c] === B.length) coreB++; if (a && b) shared++; else if (a) Aonly++; else if (b) Bonly++; if (cntA[c] === A.length && cntB[c] === B.length) bothCore++; }
  document.getElementById('g-results').style.display = 'block';
  document.getElementById('g-stats').innerHTML =
    stat(A.length, 'A models', 'eco') + stat(B.length, 'B models', 'lacto') +
    stat(Math.round(mean(perA)), 'A mean ' + gMetric) + stat(Math.round(mean(perB)), 'B mean ' + gMetric) +
    stat(coreA, 'A core') + stat(coreB, 'B core') + stat(shared, 'Shared (pan)') + stat(bothCore, 'Shared core');
  distChart(perA, perB); overlapChart(Aonly, shared, Bonly);
  // differential + volcano
  const rows = [], volc = [];
  const nTest = (() => { let t = 0; for (let c = 0; c < m.V; c++) if (cntA[c] > 0 || cntB[c] > 0) t++; return t; })();
  for (let c = 0; c < m.V; c++) {
    const a = cntA[c], cc = cntB[c]; if (a === 0 && cc === 0) continue;
    const pa = a / A.length, pb = cc / B.length, d = pa - pb;
    rows.push([c, pa, pb, d]);
    const chi = yatesChi2(a, A.length - a, cc, B.length - cc);
    volc.push([c, d, chi2NegLog10P(chi), pa, pb]);
  }
  rows.sort((x, y) => Math.abs(y[3]) - Math.abs(x[3]));
  gDiffData = { rows, m, A, B };
  renderVolcano(volc, nTest, A.length, B.length);
  const cap = 600; let h = '<thead><tr><th>ID</th><th>Name</th><th>% in A</th><th>% in B</th><th>Δ</th><th>Enriched</th></tr></thead><tbody>';
  rows.slice(0, cap).forEach(([c, pa, pb, d]) => { const en = d > 0 ? `<span style="color:${ECO};font-weight:700">A</span>` : `<span style="color:${LACTO};font-weight:700">B</span>`; h += `<tr><td class="tabular">${esc(m.vocab[c].id)}</td><td>${esc((m.vocab[c].name || '').slice(0, 40))}</td><td class="tabular">${(pa * 100).toFixed(0)}%</td><td class="tabular">${(pb * 100).toFixed(0)}%</td><td class="tabular" style="font-weight:700;color:${d > 0 ? ECO : LACTO}">${d > 0 ? '+' : ''}${(d * 100).toFixed(0)}</td><td>${Math.abs(d) > 0.001 ? en : '—'}</td></tr>`; });
  h += '</tbody>'; document.getElementById('g-diff').innerHTML = h;
}
function renderVolcano(volc, nTest, nA, nB) {
  const m = MAT[gMetric];
  const bonf = -Math.log10(0.05 / Math.max(1, nTest));
  const sigLine = -Math.log10(0.05);
  const groups = { A: [], B: [], ns: [] };
  volc.forEach(v => { const sig = v[2] >= sigLine && Math.abs(v[1]) >= 0.05; if (!sig) groups.ns.push(v); else if (v[1] > 0) groups.A.push(v); else groups.B.push(v); });
  const mk = (arr, color, name) => ({
    x: arr.map(v => v[1] * 100), y: arr.map(v => v[2]),
    customdata: arr.map(v => [esc(m.vocab[v[0]].id), esc((m.vocab[v[0]].name || '').slice(0, 40)), (v[3] * 100).toFixed(0), (v[4] * 100).toFixed(0)]),
    mode: 'markers', type: 'scattergl', name,
    marker: { size: name === 'n.s.' ? 4 : 6, color, opacity: name === 'n.s.' ? 0.35 : 0.75, line: { width: 0 } },
    hovertemplate: '<b>%{customdata[0]}</b> %{customdata[1]}<br>Δ %{x:.0f} pts · A %{customdata[2]}% / B %{customdata[3]}%<br>−log₁₀p %{y:.1f}<extra></extra>'
  });
  const yMax = Math.max(bonf + 1, ...volc.map(v => v[2]).filter(v => isFinite(v)));
  newPlot('g-volcano', [mk(groups.ns, '#c2cad6', 'n.s.'), mk(groups.B, LACTO, 'Enriched in B'), mk(groups.A, ECO, 'Enriched in A')], {
    height: 440, margin: { l: 56, r: 16, t: 12, b: 46 },
    xaxis: { title: { text: 'Δ prevalence (A − B), percentage points', font: PLOT_FONT }, gridcolor: LINE, zeroline: true, zerolinecolor: '#dfe5ee', range: [-105, 105] },
    yaxis: { title: { text: '−log₁₀ p (χ², Yates)', font: PLOT_FONT }, gridcolor: LINE, zeroline: false, rangemode: 'tozero' },
    legend: { orientation: 'h', y: 1.1, x: 0, font: { size: 10 } },
    shapes: [
      { type: 'line', x0: -105, x1: 105, y0: sigLine, y1: sigLine, line: { color: '#9aa4b2', width: 1, dash: 'dash' } },
      { type: 'line', x0: -105, x1: 105, y0: bonf, y1: bonf, line: { color: '#d0563b', width: 1, dash: 'dot' } }
    ],
    annotations: [
      { x: -103, y: sigLine, xanchor: 'left', yanchor: 'bottom', text: 'p=0.05', showarrow: false, font: { size: 9, color: INK3 } },
      { x: 103, y: bonf, xanchor: 'right', yanchor: 'bottom', text: 'Bonferroni', showarrow: false, font: { size: 9, color: '#d0563b' } }
    ]
  });
  document.getElementById('g-volcano-cap').innerHTML = `<b>${fmt(groups.A.length)}</b> features significantly enriched in A, <b>${fmt(groups.B.length)}</b> in B (p&lt;0.05 &amp; |Δ|≥5pts). Dashed = p=0.05; dotted red = Bonferroni over ${fmt(nTest)} tests. χ² with Yates' correction on a 2×2 present/absent table.`;
}
function distChart(perA, perB) {
  const all = perA.concat(perB), lo = Math.min(...all), hi = Math.max(...all), w = Math.max(20, Math.round((hi - lo) / 22) || 20);
  const start = Math.floor(lo / w) * w, nb = Math.max(1, Math.ceil((hi - start) / w) + 1), labels = [], A = new Array(nb).fill(0), B = new Array(nb).fill(0);
  for (let i = 0; i < nb; i++) labels.push(start + i * w);
  perA.forEach(v => A[Math.min(nb - 1, Math.floor((v - start) / w))]++); perB.forEach(v => B[Math.min(nb - 1, Math.floor((v - start) / w))]++);
  if (charts.gdist) charts.gdist.destroy();
  charts.gdist = new Chart(document.getElementById('g-dist'), { type: 'bar', data: { labels, datasets: [{ label: 'Group A', data: A, backgroundColor: ECO + 'cc', borderRadius: 2 }, { label: 'Group B', data: B, backgroundColor: LACTO + 'cc', borderRadius: 2 }] }, options: chartOpts(gMetric + ' per model', 'models') });
}
function overlapChart(Aonly, shared, Bonly) {
  if (charts.gov) charts.gov.destroy();
  charts.gov = new Chart(document.getElementById('g-overlap'), {
    type: 'bar', data: { labels: ['Repertoire (pan)'], datasets: [
      { label: 'A only', data: [Aonly], backgroundColor: ECO }, { label: 'Shared', data: [shared], backgroundColor: '#8fb4dd' }, { label: 'B only', data: [Bonly], backgroundColor: LACTO }] },
    options: { indexAxis: 'y', responsive: true, maintainAspectRatio: false, plugins: { legend: { labels: { boxWidth: 11, font: { size: 11 }, color: INK2, usePointStyle: true, pointStyle: 'rectRounded' } }, tooltip: { backgroundColor: INK1, padding: 9, cornerRadius: 8, callbacks: { label: c => `${c.dataset.label}: ${c.parsed.x.toLocaleString()} ${gMetric}` } } }, scales: { x: { stacked: true, grid: { color: LINE }, border: { display: false }, title: { display: true, text: gMetric, font: { size: 11 }, color: INK3 } }, y: { stacked: true, grid: { display: false } } } }
  });
}
function gDownloadCSV() { if (!gDiffData) return; const { rows, m } = gDiffData; let csv = 'id,name,pct_in_A,pct_in_B,delta\n'; rows.forEach(([c, pa, pb, d]) => { csv += `"${m.vocab[c].id}","${(m.vocab[c].name || '').replace(/"/g, '""')}",${(pa * 100).toFixed(1)},${(pb * 100).toFixed(1)},${(d * 100).toFixed(1)}\n`; }); dl(csv, `panGEMs_groups_${gMetric}_diff.csv`); }

/* =============================== CLUSTERMAP (dendrograms + strip) =============================== */
let cMetric = 'reactions', cColTree = true, lastHeat = null;
function initCluster() {
  const sp = document.getElementById('c-species'); organisms().forEach(o => { const opt = document.createElement('option'); opt.value = o; opt.textContent = `${o} (${speciesMembers(o).length})`; sp.appendChild(opt); });
  document.getElementById('c-source').addEventListener('change', e => { document.getElementById('c-species-wrap').style.display = e.target.value === 'species' ? 'flex' : 'none'; });
  document.querySelectorAll('#c-metric button').forEach(b => b.addEventListener('click', () => { seg(b, '#c-metric'); cMetric = b.dataset.metric; }));
  document.querySelectorAll('#c-coltree button').forEach(b => b.addEventListener('click', () => { seg(b, '#c-coltree'); cColTree = b.dataset.v === '1'; }));
  document.getElementById('c-run').addEventListener('click', buildCluster);
  document.getElementById('c-png').addEventListener('click', () => { if (!lastHeat) return; const a = document.createElement('a'); a.href = document.getElementById('cluster-canvas').toDataURL('image/png'); a.download = 'panGEMs_clustermap.png'; a.click(); });
  onShow.cluster = () => {};
}
function clusterRows() {
  const src = document.getElementById('c-source').value;
  if (src === 'species') { const o = document.getElementById('c-species').value; return speciesMembers(o); }
  if (src === 'groupA') return groupCache.A.slice();
  if (src === 'groupB') return groupCache.B.slice();
  if (src === 'selected') return mSel.slice();
  return [];
}
function sampleRows(rows, cap) { if (rows.length <= cap) return rows.slice(); const step = rows.length / cap, out = []; for (let i = 0; i < cap; i++) out.push(rows[Math.floor(i * step)]); return out; }
function buildCluster() {
  let rows = clusterRows();
  if (rows.length < 3) { alert('Pick a set with at least 3 models (choose a species, or run Compare groups first).'); return; }
  const capM = 120, capF = 140, sampled = sampleRows(rows, capM), m = MAT[cMetric], cnt = prevalence(sampled, m);
  const varc = []; for (let c = 0; c < m.V; c++) { const k = cnt[c]; if (k > 0 && k < sampled.length) { const p = k / sampled.length; varc.push([c, p * (1 - p)]); } }
  varc.sort((a, b) => b[1] - a[1]);
  const cols = varc.slice(0, capF).map(x => x[0]);
  if (!cols.length) { alert('These models are (nearly) identical in content — no variable ' + cMetric + ' to cluster.'); return; }
  const { D: rowD, bits } = hammingRowsD(sampled, cols, m);
  const rowT = hclustTree(rowD, sampled.length);
  let colT = { order: cols.map((_, i) => i), root: null };
  if (cColTree && cols.length <= 160) { const colD = hammingColsD(bits, sampled.length, cols.length); colT = hclustTree(colD, cols.length); }
  drawHeat(sampled, cols, bits, rowT, colT, m, rows.length);
}
function drawHeat(rows, cols, bits, rowT, colT, m, totalN) {
  const cv = document.getElementById('cluster-canvas'), ctx = cv.getContext('2d');
  const R = rows.length, C = cols.length;
  const cellW = Math.max(5, Math.min(13, Math.floor(820 / C))), cellH = Math.max(8, Math.min(16, Math.floor(760 / R)));
  const rowDendW = 66, stripW = 13, labelW = 178, colDendH = cColTree && colT.root ? 48 : 6, topPad = 6;
  const gridW = C * cellW, gridH = R * cellH;
  const W = rowDendW + stripW + gridW + labelW, H = topPad + colDendH + gridH + 6;
  const dpr = window.devicePixelRatio || 1;
  cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + 'px'; cv.style.height = H + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  ctx.font = '9px ' + PLOT_FONT.family; ctx.textBaseline = 'middle';
  const gridX = rowDendW + stripW, gridY = topPad + colDendH;
  const rowOrder = rowT.order, colOrder = colT.order;
  const rowPos = {}; rowOrder.forEach((ri, y) => rowPos[ri] = y); // ri = index into rows
  const colPos = {}; colOrder.forEach((ci, x) => colPos[ci] = x);
  // ---- row dendrogram (left, horizontal) ----
  if (rowT.root && R > 2) drawDendro(ctx, rowT.root, rowPos, { horiz: true, x0: 2, x1: rowDendW - 3, cross0: gridY, cellCross: cellH, len: rowDendW - 5 });
  // ---- col dendrogram (top, vertical) ----
  if (cColTree && colT.root && C > 2) drawDendro(ctx, colT.root, colPos, { horiz: false, x0: topPad, x1: topPad + colDendH - 3, cross0: gridX, cellCross: cellW, len: colDendH - 5 });
  // ---- annotation strip + rows ----
  rowOrder.forEach((ri, y) => {
    const row = rows[ri], mrow = META[row], yy = gridY + y * cellH;
    ctx.fillStyle = genusColor(mrow.organism); ctx.fillRect(rowDendW, yy + 0.5, stripW - 1.5, cellH - 1);
    // cells
    for (let x = 0; x < C; x++) { const ci = colOrder[x]; if (bits[ri][ci]) { const v = 0.4 + 0.5 * (x / C); ctx.fillStyle = ECO; ctx.fillRect(gridX + x * cellW, yy, cellW - (cellW > 6 ? 1 : 0.4), cellH - (cellH > 9 ? 1 : 0.4)); } }
    // label (right)
    ctx.fillStyle = INK1; ctx.textAlign = 'left'; ctx.fillText((mrow.strain || mrow.gem_file.replace(/\.json(\.json)?$/, '')).slice(0, 26), gridX + gridW + 6, yy + cellH / 2);
  });
  document.getElementById('c-results').style.display = 'block';
  // legend: genera present + present/absent key
  const orgsHere = [...new Set(rows.map(r => META[r].organism))];
  const genPresent = [...new Set(orgsHere.map(o => genusLabel(o)))];
  document.getElementById('c-legend').innerHTML =
    `<span class="it"><span class="sw" style="background:${ECO}"></span>present</span><span class="it"><span class="sw" style="background:#fff;border:1px solid ${LINE}"></span>absent</span>` +
    '<span style="margin:0 4px;color:var(--ink-4)">|</span>' +
    genPresent.slice(0, 8).map(g => `<span class="it"><span class="sw" style="background:${GENUS_COL[g] || OTHER_COL}"></span>${esc(g)}</span>`).join('');
  const shownTxt = rows.length < totalN ? ` (showing ${rows.length} of ${fmt(totalN)}, evenly sampled)` : '';
  document.getElementById('c-caption').innerHTML = `<b>${rows.length}</b> models${shownTxt} × <b>${cols.length}</b> most-variable ${cMetric}. Rows${cColTree && colT.root ? ' and columns' : ''} ordered by average-linkage hierarchical clustering (Hamming distance). Left strip = genus. Filled = feature present. Hover for detail; download PNG above.`;
  lastHeat = { rows, cols, bits, rowT, colT, m, gridX, gridY, cellW, cellH, R, C };
  // hover
  cv.onmousemove = ev => {
    const rect = cv.getBoundingClientRect(); const x = Math.floor((ev.clientX - rect.left - gridX) / cellW), y = Math.floor((ev.clientY - rect.top - gridY) / cellH);
    if (x < 0 || y < 0 || x >= C || y >= R) { hideTip(); return; }
    const ri = rowOrder[y], ci = colOrder[x], row = rows[ri], col = cols[ci], pres = bits[ri][ci];
    showTip(`<b>${esc(META[row].strain || META[row].gem_file)}</b><br><i>${esc(META[row].organism)}</i><br><span class="k">${esc(m.vocab[col].id)}</span> — ${pres ? 'present' : 'absent'}`, ev.clientX, ev.clientY);
  };
  cv.onmouseleave = hideTip;
}
// draw a dendrogram from a UPGMA tree; leafPos maps leaf-idx -> ordinal position
function drawDendro(ctx, root, leafPos, o) {
  let maxH = 0; (function mx(nd) { if (nd.isLeaf) return; maxH = Math.max(maxH, nd.height); mx(nd.left); mx(nd.right); })(root);
  if (maxH <= 0) maxH = 1;
  const crossOf = pos => o.cross0 + (pos + 0.5) * o.cellCross;
  const depthOf = h => o.horiz ? (o.x1 - (h / maxH) * o.len) : (o.x1 - (h / maxH) * o.len); // near grid = small h
  ctx.strokeStyle = '#aab6c6'; ctx.lineWidth = 1;
  (function rec(nd) {
    if (nd.isLeaf) return { cross: crossOf(leafPos[nd.idx]), depth: o.x1 };
    const a = rec(nd.left), b = rec(nd.right), d = depthOf(nd.height);
    ctx.beginPath();
    if (o.horiz) { // depth is x, cross is y
      ctx.moveTo(a.depth, a.cross); ctx.lineTo(d, a.cross); ctx.lineTo(d, b.cross); ctx.lineTo(b.depth, b.cross);
    } else { // depth is y, cross is x
      ctx.moveTo(a.cross, a.depth); ctx.lineTo(a.cross, d); ctx.lineTo(b.cross, d); ctx.lineTo(b.cross, b.depth);
    }
    ctx.stroke();
    return { cross: (a.cross + b.cross) / 2, depth: d };
  })(root);
}

/* =============================== GENOME PROPERTIES =============================== */
function pearson(xs, ys) { const n = xs.length; let sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0; for (let i = 0; i < n; i++) { const x = xs[i], y = ys[i]; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; } const cov = sxy - sx * sy / n, vx = sxx - sx * sx / n, vy = syy - sy * sy / n; return cov / Math.sqrt(vx * vy); }
function linfit(xs, ys) { const n = xs.length; let sx = 0, sy = 0, sxx = 0, sxy = 0; for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; } const b = (n * sxy - sx * sy) / (n * sxx - sx * sx), a = (sy - b * sx) / n; return { a, b }; }
function median(a) { if (!a.length) return NaN; const s = a.slice().sort((x, y) => x - y), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
const PR_YLAB = { n_reactions: 'reaction', n_metabolites: 'metabolite', n_genes: 'gene' };

function initProps() {
  onShow.props = renderProps;
  document.querySelectorAll('#pr-color button').forEach(b => b.addEventListener('click', () => { seg(b, '#pr-color'); renderProps(); }));
  document.querySelectorAll('#pr-ymetric button').forEach(b => b.addEventListener('click', () => { seg(b, '#pr-ymetric'); renderProps(); }));
}
function propTraces(rows, colorBy, xf, yf, hover) {
  const groups = {};
  rows.forEach(r => { const key = colorBy === 'dataset' ? (r.dataset === 'EcopanGEM' ? 'E. coli' : 'Lactobacillaceae') : genusLabel(r.organism); (groups[key] = groups[key] || []).push(r); });
  return Object.entries(groups).sort((a, b) => b[1].length - a[1].length).map(([name, rs]) => ({
    type: 'scatter', mode: 'markers', name,
    x: rs.map(xf), y: rs.map(yf),
    marker: { size: 5, color: colorBy === 'dataset' ? (name === 'E. coli' ? ECO : LACTO) : genusColor(rs[0].organism), opacity: 0.5, line: { width: 0 } },
    text: rs.map(r => abbr(r.organism) + (r.strain ? ' ' + r.strain : '')),
    hovertemplate: '%{text}<br>' + hover + '<extra></extra>',
  }));
}
function renderProps() {
  const colorBy = document.querySelector('#pr-color button.active').dataset.c;
  const ym = document.querySelector('#pr-ymetric button.active').dataset.y;
  document.querySelectorAll('.pr-ylabel').forEach(e => e.textContent = PR_YLAB[ym]);
  const eco = rowsOfDataset('EcopanGEM').map(i => META[i]), lac = rowsOfDataset('LactoPanGEM').map(i => META[i]);
  const med = (rows, f) => median(rows.map(r => +r[f]).filter(v => v > 0));
  document.getElementById('pr-stats').innerHTML =
    stat((med(eco, 'genome_length') / 1e6).toFixed(2) + ' Mb', 'E. coli median genome', 'eco') +
    stat((med(lac, 'genome_length') / 1e6).toFixed(2) + ' Mb', 'Lacto median genome', 'lacto') +
    stat(med(eco, 'gc_content').toFixed(1) + '%', 'E. coli median GC', 'eco') +
    stat(med(lac, 'gc_content').toFixed(1) + '%', 'Lacto median GC', 'lacto') +
    stat(fmt(Math.round(med(eco, 'n_reactions'))), 'E. coli median reactions', 'eco') +
    stat(fmt(Math.round(med(lac, 'n_reactions'))), 'Lacto median reactions', 'lacto');

  // scatter 1: GC vs genome length
  const gcRows = META.filter(r => +r.gc_content > 0 && +r.genome_length > 0);
  newPlot('pr-scatter1', propTraces(gcRows, colorBy, r => +r.gc_content, r => +r.genome_length / 1e6,
    'GC %{x:.1f}%<br>genome %{y:.2f} Mb'),
    { xaxis: { title: 'GC content (%)' }, yaxis: { title: 'genome length (Mb)' }, height: 380, legend: { orientation: 'h', y: 1.12, font: { size: 10 } }, margin: { l: 52, r: 12, t: 8, b: 42 } });

  // scatter 2: genome length vs model size + trend
  const szRows = META.filter(r => +r.genome_length > 0 && +r[ym] > 0);
  const xs = szRows.map(r => +r.genome_length / 1e6), ys = szRows.map(r => +r[ym]);
  const r = pearson(xs, ys), f = linfit(xs, ys);
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const trend = { type: 'scatter', mode: 'lines', name: 'trend', x: [xmin, xmax], y: [f.a + f.b * xmin, f.a + f.b * xmax], line: { color: INK1, width: 2, dash: 'dash' }, hoverinfo: 'skip', showlegend: false };
  newPlot('pr-scatter2', propTraces(szRows, colorBy, r2 => +r2.genome_length / 1e6, r2 => +r2[ym], 'genome %{x:.2f} Mb<br>' + PR_YLAB[ym] + 's %{y}').concat([trend]),
    { xaxis: { title: 'genome length (Mb)' }, yaxis: { title: PR_YLAB[ym] + 's per model' }, height: 380, legend: { orientation: 'h', y: 1.12, font: { size: 10 } }, margin: { l: 56, r: 12, t: 8, b: 42 },
      annotations: [{ x: 0.02, y: 0.97, xref: 'paper', yref: 'paper', text: '<b>Pearson r = ' + r.toFixed(2) + '</b>', showarrow: false, font: { size: 12, color: INK1 }, bgcolor: 'rgba(255,255,255,.82)', bordercolor: LINE, borderpad: 4, align: 'left' }] });

  // box: distribution of model size by collection
  const mets = [['n_reactions', 'Reactions'], ['n_metabolites', 'Metabolites'], ['n_genes', 'Genes']];
  const boxTrace = (rows, name, color) => {
    const xcat = [], yv = [];
    mets.forEach(([f, lbl]) => rows.forEach(rr => { const v = +rr[f]; if (v > 0) { xcat.push(lbl); yv.push(v); } }));
    return { type: 'box', name, x: xcat, y: yv, marker: { color }, line: { width: 1.2 }, boxpoints: false, fillcolor: color + '33' };
  };
  newPlot('pr-box', [boxTrace(eco, 'E. coli', ECO), boxTrace(lac, 'Lactobacillaceae', LACTO)],
    { boxmode: 'group', yaxis: { title: 'count per model' }, height: 340, legend: { orientation: 'h', y: 1.1 }, margin: { l: 56, r: 12, t: 8, b: 34 } });

  const dGC = med(eco, 'gc_content') - med(lac, 'gc_content');
  document.getElementById('pr-interp').innerHTML =
    ico() + 'The two collections separate cleanly before any modelling. '
    + '<span class="em">E. coli</span> genomes are larger (median <b>' + (med(eco, 'genome_length') / 1e6).toFixed(2) + ' Mb</b>, GC <b>' + med(eco, 'gc_content').toFixed(1) + '%</b>) than the '
    + '<span class="lm">Lactobacillaceae</span> (<b>' + (med(lac, 'genome_length') / 1e6).toFixed(2) + ' Mb</b>, GC <b>' + med(lac, 'gc_content').toFixed(1) + '%</b>) — a ' + Math.abs(dGC).toFixed(0) + '-point GC gap that is the classic low-GC signature of the lactobacilli. '
    + 'Across all ' + fmt(szRows.length) + ' genomes, genome length and ' + PR_YLAB[ym] + ' count correlate at <b>r = ' + r.toFixed(2) + '</b>: '
    + (r > 0.6 ? 'a bigger genome does buy a richer model' : r > 0.3 ? 'genome size explains part of model richness, but reconstruction and curation add scatter' : 'model richness is only weakly set by genome size — curation dominates') + '. '
    + 'The box plot shows the spread each collection contributes to the pan-model.';
}

/* =============================== GEOGRAPHY & ECOLOGY =============================== */
let GEO_JSON = null;
function ico() { return '<span class="lead"><svg viewBox="0 0 24 24" fill="none" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M9 18h6M10 21h4M12 3a6 6 0 0 1 4 10.5c-.7.7-1 1.2-1 2.5H9c0-1.3-.3-1.8-1-2.5A6 6 0 0 1 12 3Z"/></svg>Interpretation</span> ';}
function initGeo() {
  onShow.geo = renderGeo;
  document.querySelectorAll('#ge-scope button').forEach(b => b.addEventListener('click', () => { seg(b, '#ge-scope'); renderGeo(); }));
  document.querySelectorAll('#ge-scale button').forEach(b => b.addEventListener('click', () => { seg(b, '#ge-scale'); renderGeo(); }));
}
function cleanStr(s) { return String(s == null ? '' : s).trim(); }
function topCounts(rows, field, k) {
  const c = {}; rows.forEach(r => { const v = cleanStr(r[field]).toLowerCase(); if (v && v !== 'nan' && v !== 'none' && v !== 'missing' && v !== 'not collected') c[v] = (c[v] || 0) + 1; });
  return Object.entries(c).sort((a, b) => b[1] - a[1]).slice(0, k);
}
async function renderGeo() {
  if (!GEO_JSON) { try { GEO_JSON = await fetch('assets/world.geojson').then(r => r.json()); } catch (e) { GEO_JSON = { type: 'FeatureCollection', features: [] }; } }
  const scope = document.querySelector('#ge-scope button.active').dataset.s;
  const logScale = document.querySelector('#ge-scale button.active').dataset.v === 'log';
  const rows = scope === 'all' ? META : META.filter(r => r.dataset === scope);
  // country counts
  const cc = {}; rows.forEach(r => { const iso = cleanStr(r.country_iso).toUpperCase(); if (iso && iso !== 'NAN') cc[iso] = (cc[iso] || 0) + 1; });
  const isos = Object.keys(cc), counts = isos.map(i => cc[i]);
  const withCountry = counts.reduce((a, b) => a + b, 0);
  const top = Object.entries(cc).sort((a, b) => b[1] - a[1]);
  const allIsos = (GEO_JSON.features || []).map(f => f.properties && f.properties.iso).filter(Boolean);

  document.getElementById('ge-stats').innerHTML =
    stat(fmt(withCountry), 'genomes with a country') +
    stat(Math.round(100 * withCountry / rows.length) + '%', 'of the selection') +
    stat(fmt(isos.length), 'countries') +
    (top[0] ? stat(top[0][0] + ' · ' + fmt(top[0][1]), 'top country') : '');

  // base (all countries grey) + data overlay
  const base = { type: 'choropleth', geojson: GEO_JSON, featureidkey: 'properties.iso', locations: allIsos, z: allIsos.map(() => 0), showscale: false, colorscale: [[0, '#E9EDF4'], [1, '#E9EDF4']], marker: { line: { color: '#FFFFFF', width: 0.4 } }, hoverinfo: 'skip' };
  const data = { type: 'choropleth', geojson: GEO_JSON, featureidkey: 'properties.iso', locations: isos, z: logScale ? counts.map(c => Math.log10(c)) : counts, text: isos.map(i => i + ': ' + cc[i] + ' genomes'), colorscale: [[0, '#EDE7FB'], [0.5, '#B79CF0'], [1, '#6D28D9']], marker: { line: { color: '#FFFFFF', width: 0.4 } }, colorbar: { title: { text: logScale ? 'log₁₀ genomes' : 'genomes', font: { size: 10 } }, thickness: 12, len: 0.75, x: 0.98 }, hovertemplate: '%{text}<extra></extra>' };
  newPlot('ge-map', [base, data], { geo: { projection: { type: 'natural earth' }, showframe: false, showcoastlines: false, showland: false, showocean: false, bgcolor: 'rgba(0,0,0,0)', lataxis: { range: [-56, 82] }, lonaxis: { range: [-170, 190] } }, margin: { l: 0, r: 0, t: 4, b: 0 }, height: 460 });
  const topStr = top.slice(0, 5).map(([i, n]) => '<b>' + i + '</b> ' + fmt(n)).join(' · ');
  document.getElementById('ge-map-cap').innerHTML = 'Top isolation countries: ' + topStr + '. ' + Math.round(100 * (rows.length - withCountry) / rows.length) + '% of the selection has no recorded country.';

  // isolation sources
  const src = topCounts(rows, 'isolation_source', 15).reverse();
  newPlot('ge-source', [{ type: 'bar', orientation: 'h', x: src.map(s => s[1]), y: src.map(s => s[0]), marker: { color: ECO }, hovertemplate: '%{y}: %{x} genomes<extra></extra>' }],
    { xaxis: { title: 'genomes' }, yaxis: { automargin: true, tickfont: { size: 10 } }, height: 360, margin: { l: 8, r: 12, t: 8, b: 36 } });
  // hosts
  const host = topCounts(rows, 'host_name', 12).reverse();
  newPlot('ge-host', [{ type: 'bar', orientation: 'h', x: host.map(s => s[1]), y: host.map(s => s[0]), marker: { color: LACTO }, hovertemplate: '%{y}: %{x} genomes<extra></extra>' }],
    { xaxis: { title: 'genomes' }, yaxis: { automargin: true, tickfont: { size: 10 } }, height: 360, margin: { l: 8, r: 12, t: 8, b: 36 } });

  const top3 = top.slice(0, 3), top3sum = top3.reduce((a, b) => a + b[1], 0);
  const srcTop = topCounts(rows, 'isolation_source', 6).map(s => s[0]);
  document.getElementById('ge-interp').innerHTML =
    ico() + 'Of the ' + fmt(rows.length) + ' models in this selection, <b>' + fmt(withCountry) + '</b> (' + Math.round(100 * withCountry / rows.length) + '%) carry an isolation country, spread over <b>' + fmt(isos.length) + '</b> nations. '
    + (top3.length ? 'Sampling is heavily concentrated: <b>' + top3.map(t => t[0]).join(', ') + '</b> alone account for <b>' + Math.round(100 * top3sum / withCountry) + '%</b> of the geolocated genomes. ' : '')
    + (srcTop.length ? 'The habitats — <b>' + srcTop.slice(0, 5).join(', ') + '</b> — map onto the gut, dairy/fermented-food and clinical niches these species occupy. ' : '')
    + 'Read this as the <b>sampling frame</b>, not the true distribution: it reflects where genomes were sequenced and deposited, so absence from a country means "not sampled", not "not present".';
}

/* =============================== CROSS-VIEW NAV =============================== */
function openSpecies(o) {
  switchView('cluster');
  document.getElementById('c-source').value = 'species'; document.getElementById('c-species-wrap').style.display = 'flex';
  document.getElementById('c-species').value = o; buildCluster();
}

load().catch(e => { document.getElementById('loading-msg').textContent = 'Failed to load analytics data: ' + e.message; console.error(e); });
})();
