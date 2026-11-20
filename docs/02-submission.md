# GenHacks submission materials

Everything that needs to go on the Devpost page, ready to paste. Written to be
read by a judge who does not work in genetics.

---

## Project name

**Dissent**

## Tagline (one line)

Shows you the argument behind a genetic variant, not just the answer.

## Elevator pitch (Devpost's short field)

When laboratories disagree about what a genetic variant means, every public
database shows you one merged label and throws the disagreement away. We
measured how often that happens: across all 6.9 million submissions in ClinVar,
**34.5% of variants that more than one laboratory has classified are classified
differently by them** — and the rate is highest in the genes where the answer
matters most. Dissent keeps the disagreement and puts it first.

---

## What it does

Give it a genetic variant — `rs1800562`, `HFE C282Y`, `chr6:26092913G>A` — and
it produces three things in about three seconds. Give it a whole gene
panel and it orders the lot.

**A draft worksheet.** Clinical geneticists classify variants against 28 formal
criteria. Dissent evaluates all 28, with every judgement sourced and linked, and
states plainly which ones a computer cannot decide and what evidence would
settle them. Nothing is omitted, because a missing row reads as absence of
evidence.

**A conflict panel.** Twelve detectors, ranked by how much they should change
your reading, each carrying a measured base rate so a finding can be placed
rather than merely stated. They catch laboratories contradicting each other,
records with no stated method, conclusions a decade old sitting under a recent
one, prediction tools splitting, a frequency that contradicts a disease claim,
and criteria for inherited disease being applied to a tumour variant.

**A plain sentence.** Every report opens with one line that needs no
vocabulary. For HFE C282Y: *"7 laboratories call this variant disease-causing.
It is carried by about 1 in 8 people of Non-Finnish European descent. Those two
facts are hard to reconcile, and nothing in either database mentions the
other."*

**A triaged worklist.** Paste a gene panel at `/panel` and every variant is
interpreted against the same six sources, then sorted into four tiers:
contradicted, submitters split, thin evidence, no disagreement found. There is
no score — a single number invites trust in a precision that is not there, and
every row carries the reasons that put it in its tier, computed from ClinVar's
own structured fields rather than from the wording of a message. Export is CSV,
so the worklist can be annotated.

One thing is deliberately not a tier. "This tool says uncertain, the laboratory
says pathogenic" was one in the first version, and it captured 21 of 30 rows in
a test panel. That number is correct and it makes the tier useless: a
laboratory classifies from family segregation, functional assays and internal
case series, none of which is in a public API. So it is a flag on the row and a
count in the summary, labelled as expected rather than as a finding.

There is also a database-wide view at `/scale`: which genes disagreement
concentrates in, and where any of 15,326 genes places against the rest.

## Why it matters

Interpreting one variant means opening five databases and reconciling them by
hand — 20 to 60 minutes. But the time is not the problem. **The reconciliation
step is where the errors happen, and the tools that automate it hide exactly
that step behind a single label.**

The existing tools — VarSome, Franklin, OpenCRAVAT, wANNOVAR — aggregate well
and output a verdict. None of them treats disagreement as the output.

Our measurements, from ClinVar's own submission release:

| | |
|---|---|
| Variants with more than one record, classified differently | **34.5%** |
| Not re-evaluated in five years | 696,673 variants |
| No record stating any assertion criteria | 144,092 variants |
| Outright pathogenic-versus-benign conflict | 3,349 variants |

And the genes where it is worst are the ones with the most at stake: BRCA2
33.4%, TP53 33.1%, LDLR 29.8%, BRCA1 29.6%, MSH2 29.6%.

## How we built it

Python, FastAPI and six public APIs, queried concurrently in two phases. No
front-end framework, no build step, and **no API keys for anything** — a judge
can clone the repository and run it without signing up for a single service.

The architecture point that matters: every criterion and every detector is a
**pure function** over an evidence bundle. No I/O, no state. That is what makes
each rule auditable by reading one function, and what lets 411 tests run with
the socket constructor replaced — CI enforces that the suite never touches the
network, rather than a README asserting it.

Charts are hand-written SVG, so every font size is set deliberately rather than
inherited from a library's defaults.

## Challenges we ran into

**Four bugs that would each have produced a confident wrong answer.** Querying
the annotation API by rsID returned the wrong amino acid change on
multi-allelic sites — BRAF V600E came back as V600**G**. Matching by protein
position resolved BRCA1 C61G to an entirely different, benign variant, because
the position matched on a non-canonical transcript. A structural annotation was
selected as a variant's molecular consequence. And a phrase search for a protein
change ranked the gene we asked about below 25 unrelated ones, so the most
common TP53 mutation in human cancer silently failed to resolve.

**One domain error the tests caught.** We gated the strongest pathogenic
criterion on a gene constraint metric — which would have discarded every
recessive disease gene. Those metrics measure selection against *heterozygous*
loss of function, and carriers of recessive conditions are healthy, so CFTR
scores as unconstrained while loss of function is unquestionably the mechanism
of cystic fibrosis.

**One criterion that inverted its own answer.** We scored a benign criterion
from a count of people carrying two copies of a variant — until we noticed that
for HFE C282Y, those 3,248 people are not healthy controls. Carrying two copies
of that variant *is* what causes haemochromatosis. It is now reported as
context and never scored.

**Two ways to fail silently.** GraphQL endpoints answer schema errors with HTTP
200 and an `errors` array; without checking it, a renamed field would read as
"this variant is not in the database" while the status panel still said `ok`.
And an unrecognised ClinVar review-status label scored zero stars, so a single
vocabulary change would have turned every well-curated record in the database
into a high-severity evidence warning.

**A ranking that captured everything.** The panel's second tier was
originally "the computed call differs from the reported one". On a thirty-variant
test panel it took twenty-one rows. A tier that holds seventy per cent of the
work is not a priority list, it is the work. The fix was not to tune a weight
but to notice the difference was not evidence of anything: it is what a tool
built on public data should produce where a laboratory has private evidence.
The tier was replaced with one that measures disagreement *inside* ClinVar, and
the old signal became a flag labelled "expected".

**A detector that committed the error it looks for.** The same tier counted
every ClinVar record as a pathogenicity claim, so MTHFR A222V was reported as
"8 records disagree: benign vs uncertain" — when several of those records were
about warfarin dosing, which is not an answer to the same question. Reporting
non-disagreement as disagreement is precisely the failure this project was
built to catch.

**One class of input that cannot be rescued.** Much of the clinical literature
numbers protein residues from the mature protein; HGVS counts from the
initiator methionine. Factor V Leiden is R506Q in print and p.Arg534Gln in
HGVS. The offset is the signal peptide, it differs per gene, and no public
endpoint here exposes a mapping — searching dbNSFP by amino-acid position
returns the residue that sits at the HGVS position, which is a different
residue entirely. Guessing an offset would produce a confident answer about the
wrong variant. So these inputs fail, the error explains why, and one of them is
left in the demo panel on purpose.

## What we learned

That a tool in this domain is only as trustworthy as what it admits it cannot
do. Almost every design decision came down to the same choice: produce a
confident answer, or show the argument and say where it runs out. The published
classification framework itself has a gap we found and report rather than
paper over.

## What's next

Batch mode for a gene panel — "31 of the 400 genes on this cardiac panel have
variants your laboratories disagree about" — and computing the two criteria we
currently leave to a human but have the machinery to look up.

---

## Required disclosures

Per the GenHacks rules, open-source libraries, APIs and AI tools must be listed
in the project description.

**Built with:** Python 3.13, FastAPI, Uvicorn, httpx, Pydantic, pytest, respx.
No front-end framework; charts are hand-written SVG.

**Data sources, all public and keyless:** ClinVar, gnomAD v4, dbNSFP, dbSNP,
Ensembl VEP, snpEff, Open Targets, NCBI LitVar2, NCBI E-utilities.

**AI assistance:** built with Claude (Anthropic).

**Project start:** 19 November 2026, after the submission window opened on
1 November. The git history begins that day.

---

## Submission checklist

- [ ] Project page filled in on Devpost, including the disclosures above
- [ ] Track selected
- [ ] Demo video, **3:00 maximum**, on YouTube or Vimeo, **embedding enabled**
- [ ] Public repository link, with a real author identity on the commits
- [ ] All team members added on Devpost
- [ ] Live URL in the project page
- [ ] Opt in to the Deployxa prize if deployed there
- [ ] Submitted before **24 December** (official close is 00:00 ET on 26
      December, which is the end of 25 December — do not use that as the plan)

---

## Demo video, shot by shot

Three minutes exactly, with no slack. No intro card, no logo animation, no
"hi, I'm". The first frame is evidence. Captions throughout; several judges
will watch muted.

| Time | What is on screen | Narration |
|---|---|---|
| 0:00–0:12 | The real ClinVar page for HFE C282Y, cursor on "Pathogenic". Hard cut to the real gnomAD page, highlight 7% and the homozygote count. | "The official database says this variant causes disease. The population database says one in eight Northern Europeans carries it. Both pages are correct. Neither mentions the other." |
| 0:12–0:30 | Full screen: **34.5%**. Below it, small: 6,895,437 submissions · 4,470,690 variants. | "We measured every submission in ClinVar. When more than one laboratory weighs in, they disagree a third of the time. Every tool that merges this data gives you one label and throws the disagreement away." |
| 0:30–0:38 | Title card, static. *Dissent — shows you the argument, not just the answer.* | "Dissent keeps the argument." |
| 0:38–1:08 | Type `HFE C282Y`, hit Interpret. **Hold six full seconds** on the plain summary and the two verdict cards — green *Benign*, red *pathogenic*, and the line between them. Do not scroll. | "One query, six databases. Our worksheet says benign. The laboratories say pathogenic. It doesn't pick a winner — it shows both readings and why they differ." |
| 1:08–1:32 | Scroll to Conflicts. Hold on the first card and its measured base rate. | "Four conflicts, ranked. The top one: seven laboratories call it disease-causing, but it is above the frequency threshold where the guidelines say a variant is benign on its own. That combination occurs in less than one variant in a thousand." |
| 1:32–1:50 | The population chart, bars crossing the dashed 5% line. Cut to the 14-damaging / 5-benign split bar. | "The measurement behind it. And nineteen prediction tools — fourteen say damaging, five say benign. Most tools show you the fourteen." |
| 1:50–2:06 | The worksheet grid, each cell in plain English. Click one "needs a person" cell to open its detail. | "All twenty-eight criteria, in words — including the seventeen a computer can't decide. Those aren't hidden. Each says what evidence would settle it." |
| 2:06–2:30 | `/panel`. Click the carrier-screening panel, hit Triage. **Let the tiles fill on camera** — rows land every few seconds. Settle on the Contradicted group, three rows, then scroll once past Thin evidence. | "One variant saves you fifteen minutes. A panel changes what you read at all. Ten variants, same six sources, sorted by what a person still has to decide. Three are contradicted outright. One doesn't resolve — and the error tells you why." |
| 2:30–2:40 | `/scale`. The gene leaderboard. Pause on BRCA2, TP53, BRCA1. | "Across the whole database: the genes where laboratories disagree most are the ones where being wrong costs the most." |
| 2:40–2:46 | Fast cut: BRCA1 C61G, conflict panel quiet. | "When the evidence agrees, it says so. A tool that cries wolf teaches you to ignore it." |
| 2:46–2:54 | Split screen: `python -m dissent --panel panel.txt` printing the same worklist, then `pytest` → **411 passed**. | "Same engine in the terminal. Six public APIs, no keys. Four hundred and eleven tests, none touching the network." |
| 2:54–3:00 | Static end card: live URL (large), repository, track. "Built with Claude — disclosed per the rules." | "Dissent. The disagreement is the finding." |

**Before you hit record**, run this once. It fetches everything the demo
touches so nothing on camera waits on someone else's server, and exits non-zero
if any demo query would fail live:

```bash
python analysis/warm_cache.py
```

It is not a fixture dump — the cache holds the same responses the live code
fetches, keyed the same way — so nothing shown is a recording of something
other than what the tool does. Run it shortly before filming; entries expire
after a week.

**The panel shot is the one to get right.** It is the only moment that shows
the tool doing something a reader cannot do by hand, and the tiles filling in
real time is the proof it is live rather than a mock-up. Record it with a warm
cache so rows land quickly, but do not cut the waiting entirely — a progress
bar that visibly moves reads as real.

**If it runs long:** cut the BRCA1 shot to 4 seconds, or drop the terminal half
of the last demo shot. Protect the panel shot. Going over 3:00 risks the
submission.

**Avoid:** explaining what ACMG stands for (no seconds for it, nobody needs
it); narrating the architecture over the demo (put it on the project page);
apologising for an "Uncertain significance" result on camera.
