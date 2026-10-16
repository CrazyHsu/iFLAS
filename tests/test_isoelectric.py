"""BUG-145: CPC2's pI feature must reproduce the published Python-2 run (Biopython 1.76).

Expected values were computed in the ``iflas-legacy-py2`` env with Biopython 1.76
(``ProtParam.ProteinAnalysis(p).isoelectric_point()``). Biopython 1.87 returns different
numbers for several of them (e.g. DDDDEEEE 4.0500, KKKRRRHH 11.99997), which is the drift
the vendored algorithm removes.
"""

import pytest

from iflas.pu import cpc2
from iflas.pu.isoelectric import isoelectric_point

_BIOPYTHON_176 = {
    "INGAR": 9.75006103515625,
    "PETER": 4.53204345703125,
    "ACDEFGHIKLMNPQRSTVWY": 6.78448486328125,
    "DDDDEEEE": 3.35101318359375,
    "KKKRRRHH": 12.31060791015625,
    "MSTNPKPQRKTKRNTNRRPQDVKFPGG": 12.01654052734375,
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGL"
    "YTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRE"
    "RAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWH"
    "QALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL": 5.45269775390625,
    "GGGGG": 5.52496337890625,
    "CYCY": 5.51007080078125,
    "AAAAAAAAAAK": 8.79498291015625,
}


@pytest.mark.parametrize("peptide,expected", sorted(_BIOPYTHON_176.items()))
def test_isoelectric_point_matches_biopython_176_exactly(peptide, expected):
    assert isoelectric_point(peptide) == expected


def test_isoelectric_point_is_case_insensitive():
    assert isoelectric_point("ingar") == _BIOPYTHON_176["INGAR"]


def test_cpc2_protein_param_uses_the_vendored_algorithm():
    assert cpc2.protein_param("DDDDEEEE") == _BIOPYTHON_176["DDDDEEEE"]
