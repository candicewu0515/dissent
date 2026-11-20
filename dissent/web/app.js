/* Dissent — front end.
 *
 * Renders a report. Deliberately plain DOM: no framework, no build step, one
 * fewer thing to break while recording a demo, and the repository stays
 * runnable from a clone with no npm install.
 */

const $ = (id) => document.getElementById(id);

const THRESHOLDS = { ba1: 0.05, bs1: 0.01, pm2: 0.0001 };

const STATUS_LABEL = {
  met: "met",
  not_met: "argues against",
  not_applicable: "not applicable",
  insufficient: "data missing",
  manual: "needs a person",
};

const SIG_CLASS = {
  pathogenic: "p",
  likely_pathogenic: "p",
  uncertain_significance: "u",
  conflicting: "u",
  likely_benign: "b",
  benign: "b",
};

let currentReport = null;

/* ---------------- helpers ---------------- */

function text(tag, cls, content) {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (content !== undefined && content !== null) n.textContent = content;
  return n;
}

function clear(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
}

function callClass(call) {
  if (!call) return "v-none";
  const c = call.toLowerCase();
  if (c.startsWith("pathogenic") || c.startsWith("likely pathogenic")) return "v-pathogenic";
  if (c.startsWith("benign") || c.startsWith("likely benign")) return "v-benign";
  return "v-uncertain";
}

function prettySig(sig) {
  return (sig || "").replace(/_/g, " ");
}

function starString(n) {
  const filled = "★".repeat(n);
  const empty = "☆".repeat(4 - n);
  return filled + empty;
}

/* ---------------- sections ---------------- */

function renderIdentity(report) {
  const box = $("identity");
  clear(box);

  box.appendChild(text("span", "name", report.display_name || report.query.raw));

  const q = report.query;
  if (q.chrom && q.pos) {
    box.appendChild(text("span", "meta",
      `${q.chrom}:${q.pos}${q.ref || ""}>${q.alt || ""}  (${q.assembly})`));
  }
  if (report.gene_symbol) {
    const gene = text("span", "meta");
    gene.appendChild(document.createTextNode("gene "));
    const a = document.createElement("a");
    a.href = `/scale?gene=${encodeURIComponent(report.gene_symbol)}`;
    a.textContent = report.gene_symbol;
    a.title = `How often laboratories disagree about ${report.gene_symbol} variants`;
    gene.appendChild(a);
    box.appendChild(gene);
  }

  const links = document.createElement("span");
  links.className = "meta";
  const cv = report.bundle && report.bundle.clinvar;
  if (cv && cv.variant_id) {
    const a = document.createElement("a");
    a.href = `https://www.ncbi.nlm.nih.gov/clinvar/variation/${cv.variant_id}/`;
    a.target = "_blank";
    a.rel = "noopener noreferrer";
    a.textContent = "ClinVar";
    links.appendChild(a);
    links.appendChild(document.createTextNode("  "));
  }
  if (q.chrom && q.pos && q.ref && q.alt && q.assembly === "hg38") {
    const a = document.createElement("a");
    const vid = `${q.chrom.replace(/^chr/, "")}-${q.pos}-${q.ref}-${q.alt}`;
    a.href = `https://gnomad.broadinstitute.org/variant/${vid}?dataset=gnomad_r4`;
    a.target = "_blank";
    a.rel = "noopener noreferrer";
    a.textContent = "gnomAD";
    links.appendChild(a);
    links.appendChild(document.createTextNode("  "));
  }
  if (q.rsid) {
    const a = document.createElement("a");
    a.href = `https://www.ncbi.nlm.nih.gov/snp/${q.rsid}`;
    a.target = "_blank";
    a.rel = "noopener noreferrer";
    a.textContent = q.rsid;
    links.appendChild(a);
  }
  if (links.childNodes.length) box.appendChild(links);
}

/* ClinVar's majority label, with the counts that produced it. A majority is
 * reported as a majority, never as a consensus. */
function clinvarSummary(report) {
  const cv = report.bundle && report.bundle.clinvar;
  if (!cv || !cv.assertions || !cv.assertions.length) return null;

  const actionable = cv.assertions.filter(
    (a) => a.significance !== "not_provided" && a.significance !== "other"
  );
  if (!actionable.length) return null;

  const counts = {};
  for (const a of actionable) counts[a.significance] = (counts[a.significance] || 0) + 1;
  const ordered = Object.entries(counts).sort((x, y) => y[1] - x[1]);
  const [topSig, topN] = ordered[0];
  // -1 marks a review status this build does not recognise; it must not be
  // read as zero stars.
  const known = cv.assertions.filter((a) => a.stars >= 0);
  const maxStars = known.length ? Math.max(...known.map((a) => a.stars)) : -1;

  return {
    label: prettySig(topSig),
    detail: `${topN} of ${actionable.length} records` +
            (ordered.length > 1 ? `, ${ordered.length} distinct labels` : ", unanimous") +
            (maxStars >= 0 ? ` · best review status ${maxStars}/4 stars` : ""),
    unanimous: ordered.length === 1,
  };
}

function renderSummary(report) {
  const section = $("summary-section");
  const box = $("summary");
  if (!report.summary) { section.classList.add("hidden"); return; }
  section.classList.remove("hidden");
  box.textContent = report.summary;
}

function renderVerdicts(report) {
  const box = $("verdicts");
  clear(box);

  const mine = document.createElement("div");
  mine.className = `verdict ${callClass(report.classification.call)}`;
  mine.appendChild(text("div", "who", "Dissent draft worksheet"));
  mine.appendChild(text("div", "call", report.classification.call));
  const met = report.classification.met_codes || [];
  mine.appendChild(text("div", "sub",
    met.length ? `criteria met: ${met.join(", ")}` : "no criteria met from public data"));
  box.appendChild(mine);

  const cvs = clinvarSummary(report);
  const theirs = document.createElement("div");
  if (cvs) {
    theirs.className = `verdict ${callClass(cvs.label)}`;
    theirs.appendChild(text("div", "who", "ClinVar, majority of records"));
    theirs.appendChild(text("div", "call", cvs.label));
    theirs.appendChild(text("div", "sub", cvs.detail));
  } else {
    theirs.className = "verdict v-none";
    theirs.appendChild(text("div", "who", "ClinVar"));
    theirs.appendChild(text("div", "call", "No interpretation on record"));
    theirs.appendChild(text("div", "sub",
      "Nothing to compare against. The worksheet stands on its own evidence."));
  }
  box.appendChild(theirs);

  // State the disagreement between the two cards rather than leaving the
  // reader to notice it. It is the product, not a glitch.
  const mineDir = callClass(report.classification.call);
  const theirsDir = cvs ? callClass(cvs.label) : null;
  const banner = $("verdict-banner");
  if (theirsDir && mineDir !== theirsDir) {
    const opposed =
      (mineDir === "v-benign" && theirsDir === "v-pathogenic") ||
      (mineDir === "v-pathogenic" && theirsDir === "v-benign");
    banner.textContent = opposed
      ? "These two readings point in opposite directions. That is the finding."
      : "These two readings do not agree. The gap is what this report is about.";
    banner.className = `verdict-banner${opposed ? " opposed" : ""}`;
    banner.classList.remove("hidden");
  } else {
    banner.classList.add("hidden");
  }

  $("reasoning").textContent = report.classification.reasoning;
}

function renderConflicts(report) {
  const box = $("conflicts");
  clear(box);

  if (!report.conflicts.length) {
    const quiet = text("div", "quiet-panel");
    quiet.appendChild(text("strong", null, "No conflicts detected. "));
    quiet.appendChild(document.createTextNode(
      "Submitters agree, the records state their criteria, the frequency data is " +
      "consistent with the reported classification, and the predictors are not split. " +
      "That is a meaningful result, not an empty one."
    ));
    box.appendChild(quiet);
    return;
  }

  for (const c of report.conflicts) {
    const card = text("div", `conflict ${c.severity}`);
    const top = text("div", "top");
    top.appendChild(text("span", `sev ${c.severity}`, c.severity));
    top.appendChild(text("h3", null, c.headline));
    card.appendChild(top);
    card.appendChild(text("p", null, c.detail));

    // The measured base rate, set apart so it reads as a number rather than
    // more prose. It is the one line a reader without a genetics background
    // can act on.
    if (c.baseline) {
      card.appendChild(text("p", "baseline", c.baseline));
    }

    if (c.sources && c.sources.length) {
      const srcs = text("div", "srcs");
      for (const s of c.sources) {
        if (!s.url) continue;
        const a = document.createElement("a");
        a.href = s.url;
        a.target = "_blank";
        a.rel = "noopener noreferrer";
        a.textContent = s.label;
        srcs.appendChild(a);
      }
      if (srcs.childNodes.length) card.appendChild(srcs);
    }
    box.appendChild(card);
  }
}

function chartCard(title, caption, svg) {
  if (!svg) return null;
  const card = text("div", "chart-card");
  card.appendChild(text("h3", null, title));
  card.appendChild(text("p", "cap", caption));
  card.appendChild(svg);
  return card;
}

function renderCharts(report) {
  const box = $("charts");
  clear(box);
  const b = report.bundle;
  if (!b) return;

  const C = window.DissentCharts;

  const pop = chartCard(
    "Population frequency",
    "Dashed lines are the ACMG frequency thresholds. A bar crossing BA1 classifies the " +
    "variant benign on its own.",
    C.populationChart(b.frequencies, THRESHOLDS)
  );
  if (pop) box.appendChild(pop);

  const pred = chartCard(
    "Computational predictors",
    "Faded dots were derived from a numeric threshold rather than the tool's own " +
    "categorical call. Only REVEL feeds PP3, because only it has published strength " +
    "calibration.",
    C.predictorChart(b.predictors)
  );
  if (pred) box.appendChild(pred);

  const cv = b.clinvar;
  const timeline = cv ? C.timelineChart(cv.assertions) : null;
  const tl = chartCard(
    "When each assertion was last evaluated",
    "Old records sitting under a recent one are easy to miss in a list and obvious here.",
    timeline
  );
  if (tl) box.appendChild(tl);

  if (!box.childNodes.length) {
    box.appendChild(text("div", "quiet-panel",
      "No frequency, predictor or assertion data was returned for this variant."));
  }
}

function renderClinvarTable(report) {
  const section = $("clinvar-section");
  const box = $("clinvar-table");
  clear(box);

  const cv = report.bundle && report.bundle.clinvar;
  if (!cv || !cv.assertions || !cv.assertions.length) {
    section.classList.add("hidden");
    return;
  }
  section.classList.remove("hidden");

  const table = text("table", "data");
  const thead = document.createElement("thead");
  const hrow = document.createElement("tr");
  for (const h of ["Classification", "Review status", "Stars", "Submitters",
                   "Last evaluated", "Condition", "Origin"]) {
    hrow.appendChild(text("th", null, h));
  }
  thead.appendChild(hrow);
  table.appendChild(thead);

  const tbody = document.createElement("tbody");
  const sorted = [...cv.assertions].sort((a, b) => {
    if ((b.stars || 0) !== (a.stars || 0)) return (b.stars || 0) - (a.stars || 0);
    return (b.last_evaluated || "").localeCompare(a.last_evaluated || "");
  });

  for (const a of sorted) {
    const tr = document.createElement("tr");

    const sig = document.createElement("td");
    sig.appendChild(text("span", `pill ${SIG_CLASS[a.significance] || "n"}`,
      prettySig(a.significance)));
    tr.appendChild(sig);

    tr.appendChild(text("td", "mono", a.review_status || "—"));
    const st = text("td", "stars", starString(a.stars || 0));
    tr.appendChild(st);
    tr.appendChild(text("td", "mono", String(a.number_submitters || 0)));
    tr.appendChild(text("td", "mono", a.last_evaluated || "—"));
    tr.appendChild(text("td", null, a.condition || "—"));
    tr.appendChild(text("td", "mono", a.origin || "—"));
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);
  box.appendChild(table);
}

function critClass(c) {
  if (c.status === "met") return c.direction === "pathogenic" ? "met-path" : "met-ben";
  return c.status;
}

function renderCriteria(report) {
  const grid = $("crit-grid");
  clear(grid);
  clear($("crit-detail"));

  for (const c of report.criteria) {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = `crit ${critClass(c)}`;
    btn.appendChild(text("div", "code", c.code));
    const label = c.status === "met" && c.strength
      ? c.strength.replace(/_/g, " ")
      : STATUS_LABEL[c.status];
    btn.appendChild(text("div", "st", label));
    // The code is the field's vocabulary; the plain reading is what makes the
    // grid legible to a reader who does not have it.
    if (c.plain) btn.appendChild(text("div", "plain", c.plain));
    btn.title = c.plain || c.name;
    btn.addEventListener("click", () => showCriterion(c));
    grid.appendChild(btn);
  }

  const met = report.criteria.filter((c) => c.status === "met");
  showCriterion(met.length ? met[0] : report.criteria[0]);
}

function showCriterion(c) {
  const box = $("crit-detail");
  clear(box);
  const card = text("div", "crit-detail");
  card.appendChild(text("h4", null, c.plain || c.name));
  card.appendChild(text("div", "formal", `${c.code} · ${c.name}`));

  const bits = [`toward ${c.direction}`, STATUS_LABEL[c.status]];
  if (c.strength) bits.push(`${c.strength.replace(/_/g, " ")} strength`);
  card.appendChild(text("div", "codeline", bits.join("  ·  ")));

  card.appendChild(text("p", null, c.rationale));

  if (c.what_would_change) {
    const needs = text("p", "needs");
    needs.appendChild(text("strong", null, "What would settle this: "));
    needs.appendChild(document.createTextNode(c.what_would_change));
    card.appendChild(needs);
  }

  if (c.sources && c.sources.length) {
    const srcs = text("div", "srcs");
    for (const s of c.sources) {
      if (!s.url) continue;
      const a = document.createElement("a");
      a.href = s.url;
      a.target = "_blank";
      a.rel = "noopener noreferrer";
      a.textContent = s.label;
      srcs.appendChild(a);
    }
    if (srcs.childNodes.length) card.appendChild(srcs);
  }
  box.appendChild(card);
}

function renderSources(report) {
  const box = $("sources");
  clear(box);
  const b = report.bundle;
  if (!b) return;

  for (const s of b.sources) {
    const card = text("div", `src${s.ok ? "" : " bad"}`);
    card.appendChild(text("div", "nm", s.source));
    const parts = [s.ok ? "ok" : "failed"];
    if (s.cached) parts.push("cached");
    else if (s.elapsed_ms !== null && s.elapsed_ms !== undefined) parts.push(`${s.elapsed_ms} ms`);
    if (s.detail) parts.push(s.detail);
    card.appendChild(text("div", "dt", parts.join(" · ")));
    box.appendChild(card);
  }
}

function renderDisclaimer(report) {
  const box = $("disclaimer");
  clear(box);
  const manual = report.criteria.filter((c) => c.status === "manual").map((c) => c.code);
  const missing = report.criteria.filter((c) => c.status === "insufficient").map((c) => c.code);
  const na = report.criteria.filter((c) => c.status === "not_applicable").length;

  box.appendChild(text("strong", null, "This is a draft worksheet, not a clinical determination. "));
  box.appendChild(document.createTextNode(
    `${manual.length} of 28 criteria cannot be decided from public data and were not ` +
    `assessed (${manual.join(", ")}). `
  ));
  if (missing.length) {
    box.appendChild(document.createTextNode(
      `A further ${missing.length} lacked the data needed (${missing.join(", ")}). `
    ));
  }
  if (na) {
    box.appendChild(document.createTextNode(
      `${na} do not apply to a variant of this kind, which is not the same as arguing against it. `
    ));
  }
  box.appendChild(document.createTextNode(
    "Frequency thresholds are defaults, not gene-specific ones. A classification for " +
    "clinical use requires a qualified person applying gene-specific guidance to the " +
    "patient's own phenotype and family history."
  ));
}

/* ---------------- flow ---------------- */

function showStatus(msg) {
  const box = $("status");
  clear(box);
  if (msg) {
    box.appendChild(text("span", "spinner"));
    box.appendChild(document.createTextNode(msg));
    box.classList.remove("hidden");
  } else {
    box.classList.add("hidden");
  }
}

function showError(msg) {
  const box = $("error");
  if (!msg) { box.classList.add("hidden"); return; }
  box.textContent = msg;
  box.classList.remove("hidden");
}

function render(report) {
  currentReport = report;
  renderIdentity(report);
  renderSummary(report);
  renderVerdicts(report);
  renderConflicts(report);
  renderCharts(report);
  renderClinvarTable(report);
  renderCriteria(report);
  renderSources(report);
  renderDisclaimer(report);

  $("intro").classList.add("hidden");
  $("report").classList.remove("hidden");
}

async function interpret(query) {
  if (!query || !query.trim()) return;
  const assembly = $("assembly").value;

  $("q").value = query;
  $("go").disabled = true;
  showError(null);
  $("report").classList.add("hidden");
  showStatus(`Querying five sources for ${query}…`);

  try {
    const res = await fetch(
      `/api/interpret?q=${encodeURIComponent(query)}&assembly=${assembly}`
    );
    if (!res.ok) {
      let detail = `Request failed (${res.status})`;
      try {
        const body = await res.json();
        if (body.detail) detail = body.detail;
      } catch (_) { /* response was not JSON */ }
      throw new Error(detail);
    }
    render(await res.json());
    history.replaceState(null, "", `?q=${encodeURIComponent(query)}`);
  } catch (err) {
    showError(err.message || String(err));
    $("intro").classList.remove("hidden");
  } finally {
    showStatus(null);
    $("go").disabled = false;
  }
}

async function loadExamples() {
  try {
    const res = await fetch("/api/examples");
    const examples = await res.json();
    const chips = $("chips");
    clear(chips);
    for (const ex of examples) {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "chip";
      b.textContent = ex.label;
      b.title = ex.why;
      b.addEventListener("click", () => interpret(ex.query));
      chips.appendChild(b);
    }
  } catch (_) {
    /* Examples are a convenience; the search box works without them. */
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
    if (currentReport) renderCharts(currentReport);
  });
}

$("search-form").addEventListener("submit", (e) => {
  e.preventDefault();
  interpret($("q").value);
});

initTheme();
loadExamples();

const initial = new URLSearchParams(location.search).get("q");
if (initial) interpret(initial);
