# Dissent

[![tests](https://github.com/OWNER/REPO/actions/workflows/tests.yml/badge.svg)](https://github.com/OWNER/REPO/actions/workflows/tests.yml)

**A variant interpretation engine that shows you the argument, not just the answer.**

Built for [GenHacks 2026](https://genhacks-open.devpost.com/).
Track: AI, Agentic Systems, Cybersecurity, and Dev Tools.

---

## The problem

Interpreting one genetic variant means opening ClinVar, gnomAD, a predictor
site, Open Targets and PubMed, then reconciling them by hand. That costs 20–60
minutes per variant.

The time is not the real problem. The reconciliation step is where errors
happen, and the tools that automate it hide exactly that step behind a single
label. Observed while building this, in real records:

- ClinVar records **contradict each other**, and the aggregate label buries it.
  Measured across the whole database: **34.5% of variants with more than one
  record are classified differently by them**, and the rate is highest in the
  genes where the answer matters most — BRCA2 33.4%, TP53 33.1%, BRCA1 29.6%.
- Many records carry `no assertion criteria provided` — an opinion with no
  stated method — rendered with the same visual weight as an expert review.
  144,092 variants have no record above that level.
- Records go **stale**. 696,673 variants have not been re-evaluated in five
  years. BRAF V600E carries one last evaluated in 2009, under a 2025 record
  that makes the whole entry look current.
- Predictors **split**. HFE C282Y: 14 call it damaging, 5 call it benign.
- **Germline criteria get applied to somatic variants**, which is a category
  error, and nothing warns you.

Those figures are measured, not estimated: see
[`analysis/clinvar_baseline.py`](analysis/clinvar_baseline.py), which streams
ClinVar's own submission release — 6,895,437 submissions across 4,470,690
variants.

VarSome, Franklin, OpenCRAVAT and wANNOVAR aggregate well and output a label.
None of them treats *disagreement* as the primary output.

## What Dissent does

```bash
python -m dissent "HFE C282Y"
```

Three things, for any variant:

**1. A draft ACMG/AMP worksheet.** All 28 criteria, each with a status, the
strength applied, a one-sentence rationale, and a link to the record it rests
on. 12 are computed from public data; the other 16 say what evidence would
settle them. Nothing is omitted — a criterion that cannot be assessed is shown
as such, so an absent row is never mistaken for absence of evidence.

**2. A conflict panel.** Twelve detectors, ranked by severity, each carrying a
measured base rate so a finding can be placed rather than merely stated:

| Detector | What it catches |
|---|---|
| Record discordance | ClinVar records disagreeing, and whether within one condition or across several |
| Weak evidence base | Share of records stating no assertion criteria |
| Stale assertion | Old records sitting underneath a recent one |
| Frequency vs claim | Pathogenic records above the BA1/BS1 thresholds |
| Predictor split | Disagreement across ~20 dbNSFP predictors |
| Computed vs reported | This worksheet landing opposite ClinVar, *or* uncertain where ClinVar is confident |
| Somatic/germline misuse | Germline ACMG applied to a somatic variant |
| Internal contradiction | Pathogenic and benign criteria both met |
| Ambiguous input | An identifier covering several alleles at one position |
| Annotator disagreement | snpEff and VEP differing on the consequence |
| Unrecognised review status | A ClinVar vocabulary this build does not know |
| Unreviewed literature | Substantial literature, weak curation |

Each report also opens with one sentence that needs no vocabulary. For
HFE C282Y: *"7 laboratories call this variant disease-causing. It is carried by
about 1 in 8 people of Non-Finnish European descent. Those two facts are hard to
reconcile, and nothing in either database mentions the other."*

**3. What would change this.** For every unresolved criterion and every
conflict, the specific evidence that would settle it.

## Try these

| Variant | Why it is interesting |
|---|---|
| `HFE C282Y` | Asserted pathogenic at ~7% allele frequency — above the stand-alone benign threshold. Hereditary haemochromatosis has very low penetrance, and the framework has no clean answer. |
| `MTHFR A222V` | 49% allele frequency in one population, yet 12 of 19 predictors call it damaging. BA1 and PP3 are met at once. |
| `BRAF V600E` | Somatic and germline assertions mixed in one record; 13 of 17 state no assertion criteria. |
| `CFTR G542X` | Nonsense variant in a recessive gene, where gnomAD constraint is uninformative by construction. |
| `NF1 R1947*` | PVS1 met, but the evidence falls into a gap in the published combining rules. |
| `BRCA1 C61G` | Well curated and unanimous. Included so the conflict panel can be seen staying quiet. |

## Running it

Requires Python 3.11+. **No API keys, no accounts, no configuration.**

```bash
git clone <this repo> && cd dissent
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
```

Web interface:

```bash
python -m dissent serve
```

Terminal:

```bash
python -m dissent "BRAF V600E"
```

JSON, for piping into something else:

```bash
python -m dissent "rs1800562" --json
```

Tests — 308 of them, none touching the network, so they run anywhere:

```bash
python -m pytest
```

CI runs them with the socket constructor replaced, so "offline" is enforced
rather than asserted:

```bash
python -m pytest --no-network
```

## How it works

```
dissent/
  sources/      one client per external API
  normalize/    input parsing; transcript selection from snpEff
  acmg/         one pure function per criterion, plus the combining rules
  conflicts/    one pure function per detector
  api/          FastAPI routes
  web/          single-page UI, no build step
  cache/        SQLite response cache
```

Six sources are queried concurrently in two phases, because the gene symbol is
not known until the first phase answers:

| Source | Supplies |
|---|---|
| MyVariant.info | ClinVar records, ~20 dbNSFP predictors, frequencies, snpEff annotations |
| gnomAD v4 (GraphQL) | Allele and homozygote counts, per-population frequency, gene constraint |
| Ensembl VEP | Consequence, as an independent cross-check |
| Open Targets | Gene function and disease associations |
| NCBI LitVar2 | Variant-to-publication links |
| dbSNP | Citations, carried inside the MyVariant document |

### Decisions worth stating

**Every criterion and detector is a pure function** over the evidence bundle.
No I/O, no state. That is what lets the 308 tests run offline against recorded
fixtures, and what lets a reviewer audit a rule by reading one function. The
fixtures come in two layers: normalised bundles for the scoring rules, and
verbatim API payloads for the code that produces them.

**PP3/BP4 use one calibrated predictor at its calibrated strength**
([Pejaver 2022](https://pubmed.ncbi.nlm.nih.gov/36413997/)), not a tally of
tools that happen to agree. Predictors trained on overlapping data are not
independent votes, which is why the conflict panel reports the split separately
rather than scoring it.

**PM2 is applied at Supporting**, per ClinGen SVI, not Moderate.

**BS2 is reported, never scored.** It asks whether a variant is seen in a
*healthy adult* in the genotype state the disorder requires, and public data
supplies neither inheritance mode, age of onset, nor phenotype. Scoring it from
a homozygote count inverts the answer on the clearest case: HFE c.845G>A has
3,248 gnomAD homozygotes, and homozygosity for that variant *is* what causes
hereditary haemochromatosis. It was also a second frequency criterion wearing a
genotype criterion's name — anything above 1% has homozygotes, so BS1 and BS2
fired together and two Strong criteria reached a terminal "Benign" from one row
of data.

**Frequency, cohort and allele number travel together.** Reading the frequency
from the highest-popmax record and the allele-number guard from the
largest-cohort record let an 8% estimate from 50 alleles satisfy stand-alone
BA1 while the rationale credited an 800,000-allele cohort.

**PP5 and BP6 are reported but never applied.** ClinGen retired them: using
another laboratory's conclusion as evidence double-counts the evidence they
already used.

**PVS1 is not gated on gnomAD constraint.** pLI and LOEUF measure selection
against *heterozygous* loss of function, so recessive disease genes score as
unconstrained because carriers are healthy. CFTR is pLI 3.8e-65 while loss of
function is certainly the mechanism of cystic fibrosis. Gating on constraint
would have discarded every recessive disease gene.

**Where the framework runs out, it says so** rather than inventing a rule.
Richards 2015 Table 5 has no row for one very strong criterion plus a single
supporting one — a gap created by the later PM2 downgrade. Many laboratories
do call that combination likely pathogenic; this tool reports the gap instead.

**Timeouts are a total budget per source, not per attempt.** A 20s timeout with
two retries is a 60s worst case; one slow gnomAD query consumed 63s and the
surrounding `asyncio.gather` discarded a sibling source's finished result.

**A source failure never fails the report.** The UI shows which sources
answered, which were cached, and which did not respond.

**A source failure is never reported as an empty result.** GraphQL endpoints
answer schema errors with HTTP 200 and an `errors` array. Without checking it, a
renamed gnomAD field reads as "variant not in database" — which would flip PM2
to met for every common variant while the sources panel still said `ok`. The
same applies to ClinVar review-status vocabulary: an unrecognised label is
surfaced as unknown rather than scored as zero stars.

## Limits

Stated plainly, because a tool that hides its limits is worse than one that has
them:

- **Not a clinical determination.** 17 of 28 criteria need a human.
- Frequency thresholds are **defaults**, not gene-specific. Real curation uses
  ClinGen expert-panel thresholds per gene.
- **Germline only.** Somatic variants are flagged, not classified; they belong
  under the AMP/ASCO/CAP tier system.
- Protein-notation **indels are not parsed** (`CFTR F508del` fails; its
  coordinates or rsID work).
- PS1 and PM5 are marked for human review rather than computed, though the
  machinery to look them up exists.
- PVS1 runs the first branch of the Abou Tayoun 2018 decision tree only.
  Initiation-codon loss is downgraded to Moderate; the last-exon and
  NMD-escape rules need exon structure this engine does not fetch, and are
  named in the criterion rather than guessed at.
- ClinVar RCV records are variant-plus-condition aggregates, not submissions.
  Genuine submitter-level discordance lives in SCV records, which are not in
  this payload — so the per-variant panel reports record-level disagreement,
  while the database-wide base rates are computed from the SCV release.
- BA1 has no exemption list. ClinGen publishes a short list of variants that
  should not receive BA1 despite high frequency; it is not checked.
- SpliceAI is absent for many variants in the hg38 index, which leaves BP7
  unresolved more often than it should be.

## Built with

Python 3.13, FastAPI, httpx, Pydantic, pytest. No front-end framework and no
build step. Charts are hand-written SVG.

Data: ClinVar, gnomAD v4, dbNSFP, dbSNP, Ensembl VEP, snpEff, Open Targets,
NCBI LitVar2, NCBI E-utilities. All public, none requiring a key.

Developed with AI assistance (Claude), as the GenHacks rules permit and require
to be disclosed. Project begun 19 November 2026, after the submission window
opened on 1 November.

## Licence

MIT.
