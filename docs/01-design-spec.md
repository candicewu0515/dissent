# Dissent — Design Spec

**A variant interpretation engine that shows you the argument, not just the answer.**

Date: 2026-11-19
Track: AI, Agentic Systems, Cybersecurity, and Dev Tools
Deadline: working deadline 2026-12-24 (official close 2026-12-26 00:00 ET)

---

## 1. The problem

Interpreting a single genetic variant means manually opening ClinVar, gnomAD,
dbNSFP/predictor sites, Open Targets and PubMed, then reconciling them by hand.
This engine replaces eight such lookups with one.

(An earlier draft put a figure on that — "20-60 minutes per variant". It was
not measured and not cited, and a project whose argument is that unsourced
numbers mislead should not open with one. The count of lookups is checkable;
the minutes were not.)

But the real problem is not the time. It is that the reconciliation step is where
errors happen, and existing tools hide it:

- ClinVar submitters frequently **disagree** with each other, and the aggregate
  label buries the dissent.
- Many ClinVar records carry `no assertion criteria provided` — no better than an
  opinion — displayed with the same visual weight as expert-panel reviews.
- Records go **stale**. Observed in testing: a BRAF V600E record last evaluated
  2009-05-29, seventeen years old, sitting under a 2025 record and presented
  with the same weight.
- Computational predictors **split**, and tools typically show only the 2-3 that
  agree with the headline.
- **Somatic variants get germline ACMG criteria applied to them**, which is a
  category error. Nothing warns you.

Existing tools (VarSome, Franklin, OpenCRAVAT, wANNOVAR) aggregate well and
output a label. None of them treat *disagreement* as the primary output.

## 2. What Dissent does

Input: a variant, as rsID, HGVS, or `chr-pos-ref-alt`.

Output, in three panels:

### Panel 1 — ACMG/AMP criteria worksheet
Each criterion Dissent can evaluate from public data is listed with:
status (`met` / `not met` / `insufficient data`), the strength applied, the
rationale in one sentence, and the source record it rests on.

Criteria **not** automatable are listed explicitly as `requires manual review`
rather than silently omitted. Honesty here is a feature: the user needs to know
what the machine did not look at.

### Panel 2 — Conflict panel *(the signature feature)*
Actively surfaced, ranked by severity:

| Conflict type | Detection |
|---|---|
| Submitter discordance | spread of `clinical_significance` across ClinVar RCVs |
| Weak evidence base | share of RCVs with `no assertion criteria provided` |
| Stale assertion | `last_evaluated` age vs. today |
| Frequency vs. pathogenicity | ClinVar P/LP but gnomAD AF above BA1/BS1 threshold |
| Predictor split | disagreement across ~20 dbNSFP predictors |
| Computed vs. reported | Dissent's ACMG call differs from ClinVar's label |
| Somatic/germline misuse | `origin: somatic` while germline criteria are being applied |

### Panel 3 — What would change this
For each `insufficient data` criterion and each conflict, state the specific
evidence that would resolve it (e.g. "a segregation study in ≥3 affected
relatives would allow PP1").

## 3. Data sources (all verified reachable 2026-11-19, all keyless)

| Source | Endpoint | Supplies |
|---|---|---|
| MyVariant.info | `myvariant.info/v1` | ClinVar RCVs (review_status, number_submitters, last_evaluated, origin, conditions), dbNSFP 20+ predictors, gnomAD AF, CADD, SpliceAI |
| gnomAD | `gnomad.broadinstitute.org/api` (GraphQL, r4) | allele/homozygote counts, per-population AF, gene constraint (pLI, oe_lof, mis_z) |
| Ensembl VEP | `rest.ensembl.org/vep/human` | consequence, transcript, protein position, canonical flag |
| Open Targets | `api.platform.opentargets.org/api/v4/graphql` | gene function, disease associations |
| NCBI LitVar2 | `ncbi.nlm.nih.gov/research/litvar2-api` | variant → PMID links |
| ~~NCBI E-utilities~~ | ~~`eutils.ncbi.nlm.nih.gov`~~ | *planned, never built — see below* |

**Not built: NCBI E-utilities.** The plan was to pull PubMed abstracts so the
literature panel could show what a paper says rather than only that it exists.
Nothing in the engine calls it. The row is struck through rather than deleted
because a design document that quietly matches whatever got built is no longer
a record of anything, and because the submission materials listed E-utilities
as a data source on the strength of this table — a disclosure naming a service
the project never contacts.

So the count is **eight databases through five endpoints**: ClinVar, dbNSFP,
dbSNP, snpEff and gnomAD's exome frequencies all arrive inside the MyVariant
document, and gnomAD, Ensembl VEP, Open Targets and LitVar2 are queried
directly.

No API key is required for any source. This is deliberate: a judge must be able
to clone the repo and run it without signing up for anything.

## 4. Architecture

```
dissent/
  sources/        one client per external API, each returns a typed model
  normalize/      collapse multi-transcript arrays into one canonical record
  acmg/           one pure function per criterion -> CriterionResult
  conflicts/      pure functions over the evidence bundle -> Conflict
  api/            FastAPI routes
  web/            single-page UI, no build step
  cache/          SQLite response cache, TTL'd
tests/
  fixtures/       recorded API responses; tests never touch the network
```

**Design rules**

- Source clients are fetched **in parallel** (httpx.AsyncClient + gather). Total
  latency is bounded by the slowest source, not the sum. Target < 8s cold.
- Every criterion and every conflict is a **pure function** over the normalized
  evidence bundle. No network, no state. This is what makes them testable and is
  why the test suite can run offline.
- Responses are **cached in SQLite** with a TTL. Two reasons: NCBI rate-limits at
  3 req/s unkeyed, and a demo must not depend on live network conditions.
- Every claim rendered in the UI carries a **provenance link** back to the source
  record. No unsourced assertions.
- No build step on the frontend. One less thing to break while recording the demo.

## 5. Scope

**In scope**
- Single-variant interpretation, the three panels, 11 computable ACMG criteria
  of 28, 12 conflict detectors, provenance links, SQLite cache, offline test
  suite.
- *(added 2026-11-20)* Panel triage: a list of variants interpreted and sorted
  into four tiers. See the amendment below.

**Out of scope (stated so it is not ambiguous)**
- VCF parsing. A VCF carries coordinates this engine already accepts, but also
  genotypes, filters and sample columns that imply a per-sample analysis this
  is not. The panel takes a list of variants, not a sequencing output.
- User accounts, saved sessions.
- Returning a final clinical classification as if authoritative. Dissent returns
  a *draft worksheet*. The UI says so.

### Amendment, 2026-11-20 — batch was out of scope, and that was wrong

The original line read: *"Batch/VCF upload. Attractive but it dilutes the demo;
the story is one variant interpreted deeply."*

That reasoning was about the demo, not about the user. One variant interpreted
deeply replaces eight lookups; a panel changes which
variants get read at all, which is the larger claim and the one worth making.
The demo does not suffer — a panel filling in tiers on camera shows the engine
doing something a reader cannot do by hand, which a single lookup does not.

What was genuinely learned by building it is recorded here because it changed
the design twice:

- The second tier was originally *"the computed call differs from the reported
  one"*. On a thirty-variant test panel it captured twenty-one rows. That count
  is correct, and it makes the tier useless: a laboratory classifies from
  segregation, functional assays and internal case series, none of which is in
  a public API, so the difference is the expected output of a public-data tool
  rather than a finding. It became a flag labelled "expected", and the tier was
  replaced with one that measures disagreement *inside* ClinVar.
- The split test counted every ClinVar record as a claim about pathogenicity,
  which reported MTHFR A222V as *"8 records disagree: benign vs uncertain"*
  when several of those records were about warfarin dosing. Only the five
  labels on the pathogenic/benign axis count; drug response, risk factor and
  low penetrance answer a different question and are named and excluded.

Design rule carried forward: **no score.** A single number invites trust in a
precision that is not present. Four named tiers, each row carrying the reasons
that put it in one, can be argued with — and the tier is computed from
structured fields only, never from the wording of a conflict message.

## 6. Risks

| Risk | Mitigation |
|---|---|
| An upstream API changes or goes down mid-demo | SQLite cache + fixtures; demo variants pre-warmed |
| Assembly confusion (hg19 vs hg38) | assembly is explicit in every internal model; never inferred |
| Multi-transcript arrays mis-collapsed | normalization layer is separately unit-tested against recorded fixtures |
| Looks like medical advice | prominent non-diagnostic notice; "draft worksheet" wording throughout |
| Judge cannot run it | zero API keys, single command to start, pinned deps; a Dockerfile verified to run under `--read-only` |
| Deployment has no writable disk | the cache disables itself rather than failing a request, and `/api/health` reports that it did |

## 7. Honesty commitments (rules require disclosure)

- Project started 2026-11-19, after the Nov 1 window opened.
- Built with AI assistance (Claude). To be listed in the Devpost description, as
  the rules require.
- All external libraries, APIs and data sources to be listed in the submission.
- Criteria that cannot be automated are shown as such, not hidden.
