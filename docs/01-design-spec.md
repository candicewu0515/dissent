# Dissent — Design Spec

**A variant interpretation engine that shows you the argument, not just the answer.**

Date: 2026-11-19
Track: AI, Agentic Systems, Cybersecurity, and Dev Tools
Deadline: working deadline 2026-12-24 (official close 2026-12-26 00:00 ET)

---

## 1. The problem

Interpreting a single genetic variant means manually opening ClinVar, gnomAD,
dbNSFP/predictor sites, Open Targets, and PubMed, then reconciling them by hand.
That costs 20-60 minutes per variant.

But the real problem is not the time. It is that the reconciliation step is where
errors happen, and existing tools hide it:

- ClinVar submitters frequently **disagree** with each other, and the aggregate
  label buries the dissent.
- Many ClinVar records carry `no assertion criteria provided` — no better than an
  opinion — displayed with the same visual weight as expert-panel reviews.
- Records go **stale**. Observed in testing: a BRAF V600E record last evaluated
  2014-09-04, twelve years old, still presented as current.
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
| NCBI E-utilities | `eutils.ncbi.nlm.nih.gov` | PubMed abstracts |

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
- Single-variant interpretation, the three panels, ~16 computable ACMG criteria,
  all 7 conflict detectors, provenance links, SQLite cache, offline test suite.

**Out of scope (stated so it is not ambiguous)**
- Batch/VCF upload. Attractive but it dilutes the demo; the story is one variant
  interpreted deeply.
- User accounts, saved sessions.
- Returning a final clinical classification as if authoritative. Dissent returns
  a *draft worksheet*. The UI says so.

## 6. Risks

| Risk | Mitigation |
|---|---|
| An upstream API changes or goes down mid-demo | SQLite cache + fixtures; demo variants pre-warmed |
| Assembly confusion (hg19 vs hg38) | assembly is explicit in every internal model; never inferred |
| Multi-transcript arrays mis-collapsed | normalization layer is separately unit-tested against recorded fixtures |
| Looks like medical advice | prominent non-diagnostic notice; "draft worksheet" wording throughout |
| Judge cannot run it | zero API keys, single command to start, pinned deps |

## 7. Honesty commitments (rules require disclosure)

- Project started 2026-11-19, after the Nov 1 window opened.
- Built with AI assistance (Claude). To be listed in the Devpost description, as
  the rules require.
- All external libraries, APIs and data sources to be listed in the submission.
- Criteria that cannot be automated are shown as such, not hidden.
