# Dissent

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

- ClinVar submitters **contradict each other**, and the aggregate label buries it.
- Many records carry `no assertion criteria provided` — an opinion with no
  stated method — rendered with the same visual weight as an expert review.
- Records go **stale**. BRAF V600E carries an assertion last evaluated in 2009.
- Predictors **split**. HFE C282Y: 14 call it damaging, 5 call it benign.
- **Germline criteria get applied to somatic variants**, which is a category
  error, and nothing warns you.

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

**2. A conflict panel.** Eleven detectors, ranked by severity:

| Detector | What it catches |
|---|---|
| Submitter discordance | Spread of classifications across ClinVar records |
| Weak evidence base | Share of records stating no assertion criteria |
| Stale assertion | Years since the most recent re-evaluation |
| Frequency vs claim | Pathogenic assertions above the BA1/BS1 thresholds |
| Predictor split | Disagreement across ~20 dbNSFP predictors |
| Computed vs reported | This worksheet landing opposite ClinVar's majority |
| Somatic/germline misuse | Germline ACMG applied to a somatic variant |
| Internal contradiction | Pathogenic and benign criteria both met |
| Ambiguous input | An rsID covering several alleles |
| Annotator disagreement | snpEff and VEP differing on the consequence |
| Unreviewed literature | Substantial literature, weak curation |

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

Tests — these never touch the network, so they run anywhere:

```bash
python -m pytest
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
No I/O, no state. That is what lets the 151 tests run offline against recorded
fixtures, and what lets a reviewer audit a rule by reading one function.

**PP3/BP4 use one calibrated predictor at its calibrated strength**
([Pejaver 2022](https://pubmed.ncbi.nlm.nih.gov/36413997/)), not a tally of
tools that happen to agree. Predictors trained on overlapping data are not
independent votes, which is why the conflict panel reports the split separately
rather than scoring it.

**PM2 is applied at Supporting**, per ClinGen SVI, not Moderate.

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

## Limits

Stated plainly, because a tool that hides its limits is worse than one that has
them:

- **Not a clinical determination.** 16 of 28 criteria need a human.
- Frequency thresholds are **defaults**, not gene-specific. Real curation uses
  ClinGen expert-panel thresholds per gene.
- **Germline only.** Somatic variants are flagged, not classified; they belong
  under the AMP/ASCO/CAP tier system.
- Protein-notation **indels are not parsed** (`CFTR F508del` fails; its
  coordinates or rsID work).
- PS1 and PM5 are marked for human review rather than computed, though the
  machinery to look them up exists.
- SpliceAI is absent for many variants in the hg38 index, which leaves BP7
  unresolved more often than it should be.

## Built with

Python 3.13, FastAPI, httpx, Pydantic, pytest. No front-end framework and no
build step. Charts are hand-written SVG.

Developed with AI assistance (Claude), as the GenHacks rules permit and require
to be disclosed. Project begun 19 November 2026, after the submission window
opened on 1 November.

## Licence

MIT.
