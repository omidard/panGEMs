/* panGEMs landing dashboard — size distribution, species breakdown, isolation world map.
   Tailor-made for panGEMs. Uses globals: Chart (chart.js) + d3 (v7). Data: assets/dashboard.json + assets/world.geojson */
(function () {
  'use strict';
  const DS_COLORS = { EcopanGEM: '#2563EB', LactoPanGEM: '#D97706' };
  const INK = '#0F1B2D', MUTED = '#64748B', GRID = '#E3E8EF';
  const FONT = "'Fira Sans', system-ui, -apple-system, 'Segoe UI', sans-serif";

  Promise.all([
    fetch('assets/dashboard.json').then(r => r.json()),
    fetch('assets/world.geojson').then(r => r.json())
  ]).then(([db, world]) => {
    if (window.Chart) {
      Chart.defaults.font.family = FONT;
      Chart.defaults.color = MUTED;
      Chart.defaults.font.size = 12;
    }
    renderSize(db);
    renderSpecies(db);
    renderMap(db, world);
  }).catch(e => {
    console.error('dashboard load failed', e);
    const s = document.getElementById('dashboard-section');
    if (s) s.style.display = 'none';
  });

  // ---- 1. Model size distribution (metric toggle) ----
  let sizeChart = null;
  function renderSize(db) {
    const ctx = document.getElementById('chart-size');
    if (!ctx || !window.Chart) return;
    const METRIC_LABEL = { n_reactions: 'Reactions', n_genes: 'Genes', n_metabolites: 'Metabolites' };
    function build(metric) {
      const h = db.size[metric];
      const order = ['LactoPanGEM', 'EcopanGEM'].filter(d => h.series[d]);
      return {
        labels: h.labels,
        datasets: order.map(ds => ({
          label: ds,
          data: h.series[ds],
          backgroundColor: DS_COLORS[ds] + 'cc',
          borderColor: DS_COLORS[ds],
          borderWidth: 1,
          borderRadius: 3,
          categoryPercentage: 0.98,
          barPercentage: 0.96
        }))
      };
    }
    function draw(metric) {
      const data = build(metric);
      if (sizeChart) { sizeChart.data = data; sizeChart.options.scales.x.title.text = METRIC_LABEL[metric] + ' per model'; sizeChart.update(); return; }
      sizeChart = new Chart(ctx, {
        type: 'bar',
        data,
        options: {
          responsive: true, maintainAspectRatio: false,
          plugins: {
            legend: { position: 'top', labels: { usePointStyle: true, pointStyle: 'rectRounded', boxWidth: 10, boxHeight: 10, padding: 14, color: INK } },
            tooltip: { callbacks: { title: items => METRIC_LABEL[metric] + ': ' + items[0].label, label: c => `${c.dataset.label}: ${c.parsed.y.toLocaleString()} models` } }
          },
          scales: {
            x: { stacked: false, title: { display: true, text: 'Reactions per model', color: MUTED, font: { size: 11 } }, grid: { display: false }, ticks: { maxRotation: 60, minRotation: 60, autoSkip: true, maxTicksLimit: 12, font: { size: 10 } } },
            y: { title: { display: true, text: 'Models', color: MUTED, font: { size: 11 } }, grid: { color: GRID }, ticks: { precision: 0 } }
          }
        }
      });
    }
    draw('n_reactions');
    document.querySelectorAll('#size-toggle .metric-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        document.querySelectorAll('#size-toggle .metric-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        draw(btn.dataset.metric);
      });
    });
  }

  // ---- 2. Lactobacillaceae genomes per species ----
  function renderSpecies(db) {
    const ctx = document.getElementById('chart-species');
    if (!ctx || !window.Chart) return;
    const rows = db.species; // sorted desc
    new Chart(ctx, {
      type: 'bar',
      data: {
        labels: rows.map(r => r.name),
        datasets: [{
          label: 'Genomes',
          data: rows.map(r => r.count),
          backgroundColor: DS_COLORS.LactoPanGEM + 'dd',
          borderColor: DS_COLORS.LactoPanGEM,
          borderWidth: 1,
          borderRadius: 3,
          barPercentage: 0.86, categoryPercentage: 0.86
        }]
      },
      options: {
        indexAxis: 'y', responsive: true, maintainAspectRatio: false,
        plugins: {
          legend: { display: false },
          tooltip: { callbacks: { label: c => `${c.parsed.x.toLocaleString()} genomes` } }
        },
        scales: {
          x: { title: { display: true, text: 'Genome-scale models', color: MUTED, font: { size: 11 } }, grid: { color: GRID }, ticks: { precision: 0 } },
          y: { grid: { display: false }, ticks: { font: { size: 10, style: 'italic' }, color: INK, autoSkip: false } }
        }
      }
    });
  }

  // ---- 3. Isolation world map (choropleth) ----
  function renderMap(db, world) {
    const host = document.getElementById('map');
    if (!host || !window.d3) return;
    const geo = db.geo, max = db.geo_max || 1;
    const color = c => c > 0 ? d3.interpolateBlues(0.18 + 0.82 * Math.sqrt(c / max)) : '#EEF2F7';

    let tip = document.getElementById('map-tip');
    if (!tip) { tip = document.createElement('div'); tip.id = 'map-tip'; tip.className = 'map-tip'; document.body.appendChild(tip); }

    function paint() {
      const w = host.clientWidth || 900;
      const h = Math.max(320, Math.round(w * 0.5));
      host.innerHTML = '';
      const svg = d3.select(host).append('svg')
        .attr('width', w).attr('height', h)
        .attr('viewBox', `0 0 ${w} ${h}`).style('display', 'block');
      const fc = { type: 'FeatureCollection', features: world.features.filter(f => f.properties.iso !== 'ATA') };
      const proj = d3.geoNaturalEarth1().fitSize([w, h - 8], fc);
      const path = d3.geoPath(proj);
      svg.append('g').selectAll('path').data(fc.features).join('path')
        .attr('d', path)
        .attr('fill', d => color(geo[d.properties.iso] || 0))
        .attr('stroke', '#FFFFFF').attr('stroke-width', 0.5)
        .style('cursor', d => (geo[d.properties.iso] ? 'pointer' : 'default'))
        .on('mousemove', function (ev, d) {
          const n = geo[d.properties.iso] || 0;
          tip.innerHTML = `<strong>${d.properties.name}</strong><br>${n ? n.toLocaleString() + ' model' + (n === 1 ? '' : 's') : 'no isolates'}`;
          tip.style.opacity = 1;
          tip.style.left = (ev.pageX + 14) + 'px';
          tip.style.top = (ev.pageY - 10) + 'px';
          d3.select(this).attr('stroke', '#0F1B2D').attr('stroke-width', 1.2).raise();
        })
        .on('mouseleave', function () { tip.style.opacity = 0; d3.select(this).attr('stroke', '#fff').attr('stroke-width', 0.4); })
        .on('click', function (ev, d) {
          const n = geo[d.properties.iso] || 0;
          if (!n) return;
          const sel = document.getElementById('filter-country');
          if (sel && [...sel.options].some(o => o.value === d.properties.name)) {
            sel.value = d.properties.name;
            sel.dispatchEvent(new Event('change'));
            tip.style.opacity = 0;
            const tw = document.querySelector('.table-wrap');
            if (tw) tw.scrollIntoView({ behavior: 'smooth', block: 'start' });
          }
        });
    }
    paint();
    let t; window.addEventListener('resize', () => { clearTimeout(t); t = setTimeout(paint, 200); });

    // legend
    const leg = document.getElementById('map-legend');
    if (leg) {
      const stops = [];
      for (let i = 0; i <= 6; i++) stops.push(color(Math.round(max * (i / 6) * (i / 6))) );
      leg.querySelector('.grad').style.background = `linear-gradient(90deg, ${['#EEF2F7'].concat(stops.slice(1)).join(',')})`;
      leg.querySelector('.max').textContent = max.toLocaleString();
    }
  }
})();
