"""Plain-English readings of the ACMG/AMP criterion codes.

The codes are the field's vocabulary and the report keeps them, but `PVS1` and
`BS2` carry no meaning to anyone outside clinical genetics, and a worksheet of
28 such codes reads like a table of error numbers. Each line below says what
the criterion actually asks, in words that assume no background.

These are deliberately plain rather than precise: the exact wording of each
criterion is already in `CriterionResult.name`, and the evidence behind it is
in the rationale. This layer exists so a reader can follow the argument before
learning the vocabulary.
"""

from __future__ import annotations

PLAIN: dict[str, str] = {
    # Pathogenic
    "PVS1": "Switches the gene off, in a gene where that is known to cause disease",
    "PS1": "The same protein change is already known to cause disease",
    "PS2": "Arose new in the patient; neither parent carries it",
    "PS3": "Laboratory experiments show the protein stops working",
    "PS4": "Far more common in patients than in healthy people",
    "PM1": "Sits in a region of the protein where damage is known to matter",
    "PM2": "Absent, or nearly absent, from large healthy populations",
    "PM3": "In a recessive disease, paired with a known harmful variant",
    "PM4": "Makes the protein longer or shorter",
    "PM5": "A different change at the very same position causes disease",
    "PM6": "Appears to be new in the patient, but parentage was not confirmed",
    "PP1": "Tracks with the disease through the family",
    "PP2": "Harmless changes of this kind are rare in this gene",
    "PP3": "Prediction software expects the change to be damaging",
    "PP4": "The patient's symptoms point specifically to this gene",
    "PP5": "Another laboratory called it disease-causing (criterion retired)",
    # Benign
    "BA1": "Too common in healthy people to cause a rare disease",
    "BS1": "More common than the disease it is supposed to cause",
    "BS2": "Found in healthy people who should have been affected",
    "BS3": "Laboratory experiments show the protein still works",
    "BS4": "Does not track with the disease in affected families",
    "BP1": "In this gene, only gene-disabling changes cause disease",
    "BP2": "Appears alongside a variant that already explains the disease",
    "BP3": "A small insertion or deletion in a repetitive, unimportant stretch",
    "BP4": "Prediction software expects the change to be harmless",
    "BP5": "The patient's disease is already explained by something else",
    "BP6": "Another laboratory called it harmless (criterion retired)",
    "BP7": "A silent change that does not alter the protein or how it is read",
}


STATUS_PLAIN: dict[str, str] = {
    "met": "applies",
    "not_met": "does not apply",
    "insufficient": "data missing",
    "manual": "needs a person",
}


def plain_name(code: str) -> str:
    return PLAIN.get(code, code)
