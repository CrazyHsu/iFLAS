"""Isoelectric point exactly as Biopython 1.76 computed it (``Bio.SeqUtils.IsoelectricPoint``).

CPC2 feeds the peptide pI to the PU classifier. Biopython 1.78 rewrote the pK tables and
the root search, so the feature drifted from the published Python-2 run (BUG-145). This is
a port of the 1.76 bracket-and-bisection algorithm (Bjellqvist pK values), so the feature
stays value-identical regardless of the installed Biopython.
"""

_POSITIVE_PKS = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
_NEGATIVE_PKS = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
_PK_CTERMINAL = {"D": 4.55, "E": 4.75}
_PK_NTERMINAL = {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.7}
_CHARGED_AAS = ("K", "R", "H", "D", "E", "C", "Y")


def isoelectric_point(sequence):
    """Return the pI of a one-letter protein ``sequence`` (case-insensitive, non-empty)."""
    seq = str(sequence).upper()
    counts = {aa: float(seq.count(aa)) for aa in _CHARGED_AAS}
    counts["Nterm"] = counts["Cterm"] = 1.0
    pos_pks = dict(_POSITIVE_PKS, Nterm=_PK_NTERMINAL.get(seq[0], _POSITIVE_PKS["Nterm"]))
    neg_pks = dict(_NEGATIVE_PKS, Cterm=_PK_CTERMINAL.get(seq[-1], _NEGATIVE_PKS["Cterm"]))

    def charge(ph):
        positive = 0.0
        for aa, pk in pos_pks.items():
            cr = 10 ** (pk - ph)
            positive += counts[aa] * (cr / (cr + 1.0))
        negative = 0.0
        for aa, pk in neg_pks.items():
            cr = 10 ** (ph - pk)
            negative += counts[aa] * (cr / (cr + 1.0))
        return positive - negative

    # Bracket the root in 1-pH steps from pH 7, then bisect to 1e-4 (or an exact zero).
    ph = 7
    c = charge(ph)
    if c > 0.0:
        ph1 = ph
        while c > 0.0:
            ph = ph1 + 1.0
            c = charge(ph)
            if c > 0.0:
                ph1 = ph
        ph2 = ph
    else:
        ph2 = ph
        while c < 0.0:
            ph = ph2 - 1.0
            c = charge(ph)
            if c < 0.0:
                ph2 = ph
        ph1 = ph
    while ph2 - ph1 > 0.0001 and c != 0.0:
        ph = (ph1 + ph2) / 2.0
        c = charge(ph)
        if c > 0.0:
            ph1 = ph
        else:
            ph2 = ph
    return ph
