/* The scale page: database-wide conflict rates. */

const $ = (id) => document.getElementById(id);

function text(tag, cls, content) {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (content !== undefined && content !== null) n.textContent = content;
  return n;
}

function clear(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
}

const fmt = (n) => (n === null || n === undefined ? "—" : n.toLocaleString());

function renderTiles(data) {
  const box = $("tiles");
  clear(box);

  const r = data.rates_pct || {};
  const c = data.counts || {};

  const tiles = [
    {
      figure: `${r.any_label_disagreement_among_multi_submitter}%`,
      label: "of variants with more than one record are classified differently by them",
      sub: `${fmt(c.any_label_disagreement)} variants`,
      tone: "bad",
    },
    {
      figure: `${r.stale_5y}%`,
      label: "have not been re-evaluated in five years",
      sub: `${fmt(c.stale_5y)} variants`,
      tone: "warn",
    },
    {
      figure: `${r.every_record_lacks_criteria}%`,
      label: "have no record stating any assertion criteria",
      sub: `${fmt(c.every_record_lacks_criteria)} variants`,
      tone: "warn",
    },
    {
      figure: `${r.pathogenic_vs_benign}%`,
      label: "carry an outright pathogenic-versus-benign conflict",
      sub: `${fmt(c.pathogenic_vs_benign)} variants`,
      tone: "bad",
    },
  ];

  for (const t of tiles) {
    const card = text("div", `tile ${t.tone}`);
    card.appendChild(text("div", "tile-figure", t.figure));
    card.appendChild(text("div", "tile-label", t.label));
    card.appendChild(text("div", "tile-sub", t.sub));
    box.appendChild(card);
  }

  const note = text("p", "tiles-source",
    `Measured ${data.generated} across ${fmt(data.submissions)} submissions ` +
    `covering ${fmt(data.variants)} variants.`);
  box.parentElement.appendChild(note);
}

/* Horizontal bars, drawn here rather than in charts.js because this one needs
   a rank column and a reference line at the database-wide rate. */
function geneChart(genes, overallRate) {
  const SVG = "http://www.w3.org/2000/svg";
  const el = (name, attrs = {}, content = null) => {
    const n = document.createElementNS(SVG, name);
    for (const [k, v] of Object.entries(attrs)) {
      if (v !== null && v !== undefined) n.setAttribute(k, String(v));
    }
    if (content !== null) n.textContent = content;
    return n;
  };

  const data = genes.slice(0, 20);
  const rowH = 38;
  const padL = 128;
  // Two fixed right-hand columns: the rate, then the counts. Floating the
  // rate at each bar's end put it under the reference line.
  const padR = 230;
  const padT = 46;
  const padB = 30;
  const plotW = 560;
  const width = padL + plotW + padR;
  const height = padT + data.length * rowH + padB;
  const maxRate = Math.max(...data.map((g) => g.rate_pct), overallRate) * 1.04;

  const svg = el("svg", {
    viewBox: `0 0 ${width} ${height}`,
    role: "img",
    preserveAspectRatio: "xMinYMin meet",
  });

  const x = (v) => padL + (v / maxRate) * plotW;

  // The database-wide rate, so each gene can be read against it.
  const refX = x(overallRate);
  svg.appendChild(el("line", {
    x1: refX, y1: padT - 16, x2: refX, y2: padT + data.length * rowH,
    stroke: "var(--neutral)", "stroke-width": 2, "stroke-dasharray": "5 4",
  }));
  svg.appendChild(el("text", {
    x: refX, y: padT - 22, "text-anchor": "middle", class: "thresh-label",
    fill: "var(--neutral)",
  }, `all genes  ${overallRate}%`));

  data.forEach((g, i) => {
    const y = padT + i * rowH;
    const barY = y + 5;
    const barH = rowH - 14;
    const w = Math.max(x(g.rate_pct) - padL, 2);

    svg.appendChild(el("text", {
      x: padL - 14, y: barY + barH / 2 + 6, "text-anchor": "end", class: "bar-name",
    }, g.gene));
    svg.appendChild(el("rect", {
      x: padL, y: barY, width: w, height: barH, rx: 3,
      fill: "var(--path)", "fill-opacity": 0.82,
    }));
    svg.appendChild(el("text", {
      x: padL + plotW + 56, y: barY + barH / 2 + 6, "text-anchor": "end",
      class: "ax-value",
    }, `${g.rate_pct.toFixed(1)}%`));
    svg.appendChild(el("text", {
      x: width - 8, y: barY + barH / 2 + 6, "text-anchor": "end", class: "ax-label",
    }, `${g.discordant.toLocaleString()} of ${g.variants.toLocaleString()}`));
  });

  return svg;
}

function submittersChart(histogram) {
  const SVG = "http://www.w3.org/2000/svg";
  const el = (name, attrs = {}, content = null) => {
    const n = document.createElementNS(SVG, name);
    for (const [k, v] of Object.entries(attrs)) {
      if (v !== null && v !== undefined) n.setAttribute(k, String(v));
    }
    if (content !== null) n.textContent = content;
    return n;
  };

  const entries = Object.entries(histogram)
    .map(([k, v]) => [Number(k), Number(v)])
    .sort((a, b) => a[0] - b[0]);
  if (!entries.length) return null;

  const total = entries.reduce((s, [, v]) => s + v, 0);
  const max = Math.max(...entries.map(([, v]) => v));

  const barW = 62;
  const gap = 14;
  const padT = 34;
  const padB = 64;
  const plotH = 250;
  const width = entries.length * (barW + gap);
  const height = padT + plotH + padB;

  const svg = el("svg", {
    viewBox: `0 0 ${width} ${height}`,
    role: "img",
    preserveAspectRatio: "xMinYMin meet",
  });

  entries.forEach(([n, count], i) => {
    const h = Math.max((count / max) * plotH, 2);
    const x = i * (barW + gap);
    const y = padT + plotH - h;
    // One record means nothing to disagree with, so it is drawn neutral.
    const color = n === 1 ? "var(--neutral)" : "var(--path)";

    svg.appendChild(el("rect", {
      x, y, width: barW, height: h, rx: 3, fill: color, "fill-opacity": 0.82,
    }));
    svg.appendChild(el("text", {
      x: x + barW / 2, y: y - 9, "text-anchor": "middle", class: "ax-value",
    }, `${Math.round((count / total) * 100)}%`));
    svg.appendChild(el("text", {
      x: x + barW / 2, y: padT + plotH + 26, "text-anchor": "middle", class: "ax-label",
    }, n >= 10 ? "10+" : String(n)));
  });

  svg.appendChild(el("text", {
    x: 0, y: height - 14, class: "ax-title",
  }, "LABORATORIES CLASSIFYING ONE VARIANT"));

  return svg;
}

async function lookupGene(symbol) {
  const box = $("gene-result");
  clear(box);
  if (!symbol.trim()) return;

  try {
    const res = await fetch(`/api/baseline/gene?symbol=${encodeURIComponent(symbol.trim())}`);
    if (!res.ok) {
      const body = await res.json().catch(() => ({}));
      box.appendChild(text("div", "error", body.detail || "Lookup failed."));
      return;
    }
    const g = await res.json();
    const card = text("div", "gene-card");

    // Express the position as a percentile from the top. Rounding
    // "more discordant than 99.93% of genes" gave 100%, which cannot be true
    // of a gene that is in the set being compared against.
    const topPct = (g.rank / g.of_genes) * 100;
    const band =
      topPct < 1 ? `the most discordant ${topPct.toFixed(2)}%`
      : topPct < 10 ? `the most discordant ${topPct.toFixed(1)}%`
      : `the most discordant ${Math.round(topPct)}%`;

    card.appendChild(text("div", "gene-name", g.gene));
    card.appendChild(text("div", "gene-rate", `${g.rate_pct}%`));
    card.appendChild(text("div", "gene-detail",
      `${g.discordant.toLocaleString()} of ${g.variants.toLocaleString()} classified ` +
      `variants are classified differently by different laboratories.`));
    card.appendChild(text("div", "gene-rank",
      `Rank ${g.rank.toLocaleString()} of ${g.of_genes.toLocaleString()} genes, in ` +
      `${band} of them. The figure across all genes is ${g.overall_rate_pct}%.`));
    box.appendChild(card);
  } catch (err) {
    box.appendChild(text("div", "error", String(err)));
  }
}

const LABEL_CLASS = {
  pathogenic: "p",
  "likely pathogenic": "p",
  "pathogenic/likely pathogenic": "p",
  "uncertain significance": "u",
  "likely benign": "b",
  benign: "b",
  "benign/likely benign": "b",
};

async function loadOpposed() {
  const expertOnly = $("expert-only").checked;
  const box = $("opposed");
  clear(box);

  const res = await fetch(`/api/baseline/opposed?limit=25&min_stars=${expertOnly ? 3 : 0}`);
  if (!res.ok) {
    box.appendChild(text("div", "error", "Could not load the list."));
    return;
  }
  const data = await res.json();

  $("opposed-count").textContent =
    `${data.total.toLocaleString()} variants` +
    (expertOnly ? " with an expert-panel review" : " in total") +
    `, showing ${data.variants.length}`;

  for (const v of data.variants) {
    const row = document.createElement("a");
    row.className = "opposed-row";
    // The interpreter accepts a ClinVar variation id directly.
    row.href = `/?q=${v.id}`;
    row.title = `Open the full worksheet for ClinVar variation ${v.id}`;

    const left = text("div", "opposed-left");
    left.appendChild(text("span", "opposed-gene", v.gene || "—"));
    left.appendChild(text("span", "opposed-id", `ClinVar ${v.id}`));
    if (v.best_stars >= 3) left.appendChild(text("span", "opposed-badge", "expert panel"));
    if (v.last) left.appendChild(text("span", "opposed-date", `last evaluated ${v.last}`));
    row.appendChild(left);

    const labels = text("div", "opposed-labels");
    for (const [label, n] of v.labels) {
      const cls = LABEL_CLASS[label] || "n";
      labels.appendChild(text("span", `pill ${cls}`, `${n}× ${label}`));
    }
    row.appendChild(labels);
    box.appendChild(row);
  }
}

function initTheme() {
  const stored = (() => {
    try { return localStorage.getItem("dissent-theme"); } catch (_) { return null; }
  })();
  if (stored === "light" || stored === "dark") {
    document.documentElement.setAttribute("data-theme", stored);
  }
  $("theme-toggle").addEventListener("click", () => {
    const current = document.documentElement.getAttribute("data-theme");
    const prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
    const effective = current || (prefersDark ? "dark" : "light");
    const next = effective === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    try { localStorage.setItem("dissent-theme", next); } catch (_) { /* private mode */ }
  });
}

async function load() {
  const res = await fetch("/api/baseline");
  if (!res.ok) {
    $("tiles").appendChild(text("div", "error", "Baseline data is not installed."));
    return;
  }
  const data = await res.json();
  renderTiles(data);

  const overall = data.rates_pct.any_label_disagreement_among_multi_submitter;
  const chart = geneChart(data.most_discordant_genes, overall);
  if (chart) $("gene-chart").appendChild(chart);
  $("gene-chart-note").textContent =
    `Limited to genes with at least ${data.leaderboard_min_variants} classified ` +
    `variants, so the ordering is not driven by small-sample noise. A rate is quoted ` +
    `for any gene with at least ${data.gene_min_variants}, of which there are ` +
    `${data.genes_with_a_rate.toLocaleString()} — a larger and noisier population, ` +
    `and the one the lookup below ranks against.`;

  const hist = submittersChart(data.submitters_per_variant);
  if (hist) $("submitters-chart").appendChild(hist);

  await loadOpposed();
}

$("expert-only").addEventListener("change", loadOpposed);

$("gene-form").addEventListener("submit", (e) => {
  e.preventDefault();
  lookupGene($("gene").value);
});

initTheme();
load();

// Arriving from a variant report: look that gene up immediately.
const requested = new URLSearchParams(location.search).get("gene");
if (requested) {
  $("gene").value = requested;
  lookupGene(requested);
  document.addEventListener("DOMContentLoaded", () =>
    $("gene-form").scrollIntoView({ behavior: "smooth", block: "center" })
  );
}
