/* The panel page: interpret a list of variants, then order it.

   Nodes are built rather than assembled from HTML strings. Every string on
   this page -- a gene symbol, a ClinVar condition name, an error message --
   arrives from a public API, and textContent cannot be talked into executing
   anything. */

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

/* How many variants go in one request. Each fans out to six public APIs, so a
   whole panel in one call would sit silent for a minute; a chunk lands every
   few seconds and the table fills in front of the reader. */
const CHUNK = 5;

const TIER_TONE = {
  contradicted: "bad",
  split: "warn",
  thin: "warn",
  concordant: "good",
  unresolved: "neutral",
};

let TIERS = [];
let ENP_TEXT = "";
let rows = [];
let running = false;

// --------------------------------------------------------------------------
// Rendering
// --------------------------------------------------------------------------

function tierName(key) {
  const t = TIERS.find((x) => x.key === key);
  return t ? t.name : key;
}

function stars(n) {
  if (n === null || n === undefined) return "";
  if (n < 0) return "unknown";
  return "★".repeat(n) + "☆".repeat(Math.max(0, 4 - n));
}

function renderTiles() {
  const box = $("tier-tiles");
  clear(box);
  const counts = {};
  rows.forEach((r) => {
    counts[r.tier] = (counts[r.tier] || 0) + 1;
  });

  TIERS.forEach((t) => {
    const n = counts[t.key] || 0;
    if (!n && t.key === "unresolved") return;
    const tile = text("div", `tile ${TIER_TONE[t.key] || "neutral"}`);
    tile.appendChild(text("div", "tile-figure", String(n)));
    tile.appendChild(text("div", "tile-label", t.name));
    box.appendChild(tile);
  });

  const enp = rows.filter((r) => r.evidence_not_public).length;
  $("enp-note").textContent = enp
    ? `${enp} of ${rows.length} rows are reported confidently while the computed call ` +
      `stays uncertain. ${ENP_TEXT}`
    : "";
}

function callClass(call) {
  if (!call) return "pill p-none";
  if (call.startsWith("Pathogenic") || call.startsWith("Likely pathogenic")) return "pill p-path";
  if (call.startsWith("Benign") || call.startsWith("Likely benign")) return "pill p-benign";
  return "pill p-uncertain";
}

function renderRow(r) {
  const card = text("div", `prow prow-${r.tier}`);

  const head = text("div", "prow-head");
  const name = text("div", "prow-name");
  if (r.resolved) {
    const a = text("a", "prow-link", r.display_name || r.query);
    a.href = `/?q=${encodeURIComponent(r.query)}`;
    name.appendChild(a);
  } else {
    name.appendChild(text("span", "prow-link", r.query));
  }
  if (r.display_name && r.display_name !== r.query) {
    name.appendChild(text("span", "prow-asked", `asked as ${r.query}`));
  }
  head.appendChild(name);

  const badge = text("span", `prow-tier tone-${TIER_TONE[r.tier] || "neutral"}`, tierName(r.tier));
  head.appendChild(badge);
  card.appendChild(head);

  const meta = text("div", "prow-meta");
  const add = (label, node) => {
    const cell = text("div", "prow-cell");
    cell.appendChild(text("span", "prow-key", label));
    cell.appendChild(node);
    meta.appendChild(cell);
  };

  add("computed", text("span", callClass(r.computed_call), r.computed_call || "—"));
  add("reported", text("span", callClass(r.reported_call), r.reported_call || "no record"));
  if (r.best_stars !== null && r.best_stars !== undefined) {
    add("review", text("span", "stars", stars(r.best_stars)));
  }
  if (r.records) add("records", text("span", "prow-val", String(r.records)));
  if (r.conflicts) {
    add(
      "conflicts",
      text("span", "prow-val", r.high_severity ? `${r.conflicts} (${r.high_severity} high)` : String(r.conflicts))
    );
  }
  if (r.met_codes && r.met_codes.length) {
    add("codes", text("span", "prow-codes", r.met_codes.join(" ")));
  }
  card.appendChild(meta);

  const why = text("ul", "prow-reasons");
  (r.reasons || []).forEach((reason) => why.appendChild(text("li", null, reason)));
  if (r.evidence_not_public) {
    why.appendChild(text("li", "prow-normal", "Submitter's evidence is not public — expected, not a flag"));
  }
  if (r.error) why.appendChild(text("li", "prow-err", r.error));
  card.appendChild(why);

  return card;
}

function renderRows() {
  const box = $("rows");
  clear(box);
  let lastTier = null;
  rows.forEach((r) => {
    if (r.tier !== lastTier) {
      lastTier = r.tier;
      const t = TIERS.find((x) => x.key === r.tier);
      const head = text("div", "prow-group");
      head.appendChild(text("h3", "prow-group-name", t ? t.name : r.tier));
      if (t) head.appendChild(text("p", "prow-group-why", t.meaning));
      box.appendChild(head);
    }
    box.appendChild(renderRow(r));
  });
}

function sortRows() {
  // Must match dissent/panel.py sort_rows, so the page, the CSV and the CLI
  // present one order.
  rows.sort((a, b) => {
    if (a.tier_rank !== b.tier_rank) return a.tier_rank - b.tier_rank;
    if (a.high_severity !== b.high_severity) return b.high_severity - a.high_severity;
    if (a.conflicts !== b.conflicts) return b.conflicts - a.conflicts;
    const sa = a.best_stars === null || a.best_stars === undefined ? 9 : a.best_stars;
    const sb = b.best_stars === null || b.best_stars === undefined ? 9 : b.best_stars;
    if (sa !== sb) return sa - sb;
    return (a.display_name || a.query).localeCompare(b.display_name || b.query);
  });
}

function renderTierDefs() {
  const box = $("tier-defs");
  clear(box);
  TIERS.forEach((t) => {
    const d = text("div", "tier-def");
    d.appendChild(text("span", `prow-tier tone-${TIER_TONE[t.key] || "neutral"}`, t.name));
    d.appendChild(text("p", "tier-def-text", t.meaning));
    box.appendChild(d);
  });
  const d = text("div", "tier-def");
  d.appendChild(text("span", "prow-tier tone-neutral", "Not a tier"));
  d.appendChild(text("p", "tier-def-text", ENP_TEXT));
  box.appendChild(d);
  $("rule").classList.remove("hidden");
}

// --------------------------------------------------------------------------
// Running
// --------------------------------------------------------------------------

function parseInput(raw) {
  const seen = new Set();
  const out = [];
  raw
    .replace(/,/g, "\n")
    .split("\n")
    .forEach((line) => {
      const q = line.trim().replace(/^["']|["']$/g, "").trim();
      if (!q || q.startsWith("#")) return;
      const key = q.toUpperCase();
      if (seen.has(key)) return;
      seen.add(key);
      out.push(q);
    });
  return out;
}

function progress(done, total) {
  $("progress-wrap").classList.remove("hidden");
  const pct = total ? Math.round((done / total) * 100) : 0;
  $("progress-fill").style.width = `${pct}%`;
  $("progress-text").textContent =
    done < total ? `${done} of ${total} interpreted` : `${total} interpreted`;
}

async function run() {
  if (running) return;
  const queries = parseInput($("queries").value);
  if (!queries.length) {
    $("hint").textContent = "Nothing to interpret.";
    return;
  }

  running = true;
  rows = [];
  $("run").disabled = true;
  $("hint").textContent = "";
  $("results").classList.remove("hidden");
  clear($("rows"));
  clear($("tier-tiles"));
  progress(0, queries.length);

  const assembly = $("assembly").value;
  let done = 0;

  for (let i = 0; i < queries.length; i += CHUNK) {
    const chunk = queries.slice(i, i + CHUNK);
    try {
      const res = await fetch("/api/panel", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ queries: chunk, assembly }),
      });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      rows = rows.concat(data.rows || []);
    } catch (err) {
      // A failed chunk costs that chunk, not the panel. The rows still appear,
      // marked unresolved, so the reader can see what was not covered rather
      // than reading a short table as a complete one.
      chunk.forEach((q) =>
        rows.push({
          query: q,
          resolved: false,
          tier: "unresolved",
          tier_rank: 5,
          conflicts: 0,
          high_severity: 0,
          records: 0,
          reasons: ["Request failed"],
          error: String(err.message || err),
        })
      );
    }
    done += chunk.length;
    sortRows();
    renderTiles();
    renderRows();
    progress(done, queries.length);
  }

  running = false;
  $("run").disabled = false;
}

async function download() {
  if (!rows.length) return;
  const res = await fetch("/api/panel/csv", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ rows }),
  });
  if (!res.ok) {
    $("hint").textContent = "Could not build the file.";
    return;
  }
  const blob = await res.blob();
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "dissent-panel.csv";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// --------------------------------------------------------------------------
// Setup
// --------------------------------------------------------------------------

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

async function loadMeta() {
  try {
    const res = await fetch("/api/panel/tiers");
    const data = await res.json();
    TIERS = data.tiers || [];
    ENP_TEXT = data.evidence_not_public || "";
    renderTierDefs();
    $("hint").textContent = `Up to ${data.max_per_request} variants per request.`;
  } catch (_) {
    $("hint").textContent = "";
  }

  try {
    const res = await fetch("/api/panel/examples");
    const panels = await res.json();
    const box = $("panel-chips");
    clear(box);
    panels.forEach((p) => {
      const b = text("button", "chip", `${p.name} (${p.queries.length})`);
      b.type = "button";
      b.addEventListener("click", () => {
        $("queries").value = p.queries.join("\n");
        $("panel-note").textContent = p.note;
      });
      box.appendChild(b);
    });
  } catch (_) { /* examples are a convenience, not a requirement */ }
}

initTheme();
loadMeta();
$("panel-form").addEventListener("submit", (e) => {
  e.preventDefault();
  run();
});
$("download").addEventListener("click", download);
