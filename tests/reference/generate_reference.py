"""Generate independent reference m/z values for mzPAF annotations.

This script does not import paftacular, tacular or peptacular. Every value is
built from the formulas in the mzPAF 1.0.1 specification (section 4.4.3 table,
4.4.4 to 4.4.7, 4.5 to 4.8) using pyteomics atomic masses (NIST) and pyteomics
amino acid compositions. Modification formulas are the Unimod compositions.

Regenerate with (from the repository root):

    uv run --no-project --with pyteomics==5.0.1 python tests/reference/generate_reference.py

Produced with pyteomics 5.0.1 on Python 3.13. The output is
tests/reference/mzpaf_reference.json. reference_molecules.json is an unmodified
copy of the official mzPAF Appendix B list (HUPO-PSI/mzPAF on GitHub). pyteomics
is not a dependency of paftacular. Only the frozen JSON is read by
tests/test_spec_reference.py.
"""

import json
import re
from pathlib import Path

from pyteomics import mass

OUT = Path(__file__).with_name("mzpaf_reference.json")

NIST = mass.nist_mass
PROTON = NIST["H+"][0][0]
ELECTRON = NIST["e*"][0][0]


def comp(formula: str) -> mass.Composition:
    """Composition from a formula that may contain ProForma isotope tokens such as [13C1]."""
    counts: dict[str, int] = {}
    for isotope, element, count in _tokens(formula):
        key = f"{element}[{isotope}]" if isotope else element
        counts[key] = counts.get(key, 0) + count
    return mass.Composition(counts)


def _tokens(formula: str):
    for match in re.finditer(r"\[(\d+)([A-Z][a-z]?)(-?\d*)\]|([A-Z][a-z]?)(-?\d*)", formula):
        if match.group(2):
            yield match.group(1), match.group(2), int(match.group(3) or 1)
        else:
            yield None, match.group(4), int(match.group(5) or 1)


def m(composition: mass.Composition) -> float:
    return mass.calculate_mass(composition=composition)


# Unimod monoisotopic compositions of the modifications used below.
MODS = {
    "Oxidation": "O",
    "Phospho": "HO3P",
    "Acetyl": "C2H2O",
    "Carbamidomethyl": "C2H3NO",
    "Amidated": "HNO-1",
    "Hex": "C6H10O5",
    "HexNAc(2)": "C16H26N2O10",
    "TMT6plex": "C8H20[13C4]N[15N]O2",
}

# Residue (-NH-CHR-CO-) compositions from pyteomics.
AA = {aa: mass.Composition(mass.std_aa_comp[aa]) for aa in "ACDEFGHIKLMNPQRSTVWY"}
BACKBONE = comp("C2H2NO")  # residue minus side chain R


def residues(proforma: str) -> tuple[list[tuple[str, mass.Composition]], mass.Composition, mass.Composition]:
    """Parse the small ProForma subset used here: [N]-SEQ[mod]...-[C]."""
    nterm = mass.Composition()
    cterm = mass.Composition()
    match = re.match(r"^\[([^\]]+)\]-", proforma)
    if match:
        nterm = comp(MODS[match.group(1)])
        proforma = proforma[match.end() :]
    match = re.search(r"-\[([^\]]+)\]$", proforma)
    if match:
        cterm = comp(MODS[match.group(1)])
        proforma = proforma[: match.start()]
    out = []
    for aa, mod in re.findall(r"([A-Z])(?:\[([^\]]+)\])?", proforma):
        c = mass.Composition(AA[aa])
        if mod:
            c += comp(MODS[mod])
        out.append((aa, c))
    return out, nterm, cterm


def total(parts) -> mass.Composition:
    c = mass.Composition()
    for part in parts:
        c += part
    return c


def mz_protonated(neutral: mass.Composition, charge: int) -> float:
    return (m(neutral) + charge * PROTON) / charge


# mzPAF 1.0.1 section 4.4.3 table, relative to the sum of residues.
SERIES = {
    "a": comp("C-1O-1"),
    "b": mass.Composition(),
    "c": comp("NH3"),
    "x": comp("CO2"),
    "y": comp("H2O"),
    "z": comp("H2O") - comp("NH2"),
}
PYTEOMICS_ION = {"a": "a", "b": "b", "c": "c", "x": "x", "y": "y", "z": "z-dot"}

cases: list[dict] = []


def add(annotation: str, value: float, charge: int, note: str, analyte: str | None = None) -> None:
    # Tolerance in Da on m/z. The default is 1e-6. Two known data differences widen it:
    # pyteomics 5.0.1 rounds 18O to 17.999161 (AME2020 17.99915961), and tacular keeps
    # Unimod modification masses to 6 decimals, so each named modification may add 5e-7.
    tolerance = 1e-6
    text = annotation + (analyte or "")
    if "18O" in text:
        tolerance = 5e-6
    elif any(f"[{name}]" in text for name in MODS):
        tolerance = 3e-6
    cases.append({"annotation": annotation, "analyte": analyte, "charge": charge, "mz": round(value, 9), "tolerance": tolerance, "note": note})


# Backbone series, resolved from analytes (termini and modifications included).
ANALYTES = ["PEPTIDEK", "MYPEPTIDEK", "[Acetyl]-M[Oxidation]YPEPS[Phospho]TIDEK", "PEPTC[Carbamidomethyl]DEK-[Amidated]"]
for analyte in ANALYTES:
    res, nterm, cterm = residues(analyte)
    n = len(res)
    for series, group in SERIES.items():
        for position in sorted({1, 2, n // 2, n - 1}):
            forward = series in "abc"
            fragment = res[:position] if forward else res[n - position :]
            terminal = nterm if forward else cterm
            neutral = total(c for _, c in fragment) + terminal + group
            for charge in (1, 2, 3):
                value = mz_protonated(neutral, charge)
                if not (nterm or cterm) and all(c == AA[a] for a, c in fragment):
                    seq = "".join(a for a, _ in fragment)
                    check = mass.fast_mass(seq, ion_type=PYTEOMICS_ION[series], charge=charge)
                    assert abs(check - value) < 1e-9, (series, seq, check, value)
                suffix = f"^{charge}" if charge > 1 else ""
                add(f"{series}{position}{suffix}", value, charge, "4.4.3 backbone series", analyte)

# Embedded sequences from section 4.4.3.
for text, seq, series, extra, adduct_formula, charge in [
    ("0@b2{LL}", "LL", "b", "", None, 1),
    ("0@y1{K}", "K", "y", "", None, 1),
    ("0@y1{K}-NH3", "K", "y", "-NH3", None, 1),
    ("0@b2{LC[Carbamidomethyl]}", "LC[Carbamidomethyl]", "b", "", None, 1),
    ("0@b1{[Acetyl]-M}", "[Acetyl]-M", "b", "", None, 1),
    ("0@y4{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", "M[Oxidation]ACK", "y", "-CH4OS", "HNa", 2),
]:
    res, nterm, cterm = residues(seq)
    neutral = total(c for _, c in res) + nterm + cterm + SERIES[series]
    if extra:
        neutral -= comp(extra[1:])
    if adduct_formula:
        value = (m(neutral + comp(adduct_formula)) - charge * ELECTRON) / charge
    else:
        value = mz_protonated(neutral, charge)
    add(text, value, charge, "4.4.3 embedded sequence example")

# Side-chain ions (d, v, w) from their parent ions and the spec 4.4.3 remarks.
# d = a minus the side-chain group beyond the beta carbon, w = z minus the same group,
# v = y minus the whole side chain R plus one hydrogen (C2H3NO2 on the n-1 residue sum).
LOST_GENERIC = {aa: AA[aa] - BACKBONE - comp("CH2") for aa in "CDEFHKLMNQRSWY"}
LOST = {("d", aa): group for aa, group in LOST_GENERIC.items()}
LOST.update({("w", aa): group for aa, group in LOST_GENERIC.items()})
LOST.update(
    {
        ("d", "V"): comp("CH3"),
        ("w", "V"): comp("CH3"),
        ("da", "T"): comp("CH3"),  # threonine da keeps OH
        ("db", "T"): comp("OH"),  # threonine db keeps CH3
        ("da", "I"): comp("CH3"),  # isoleucine da keeps C2H5
        ("db", "I"): comp("C2H5"),  # isoleucine db keeps CH3
        ("wa", "T"): comp("CH3"),
        ("wb", "T"): comp("OH"),
        ("wa", "I"): comp("CH3"),
        ("wb", "I"): comp("C2H5"),
    }
)
for series, sequence in [
    ("d", "PEL"),
    ("d", "PES"),
    ("d", "PEK"),
    ("d", "GV"),
    ("da", "PET"),
    ("db", "PET"),
    ("da", "PEI"),
    ("db", "PEI"),
    ("d", "[Acetyl]-L"),
    ("w", "LEK"),
    ("w", "SEK"),
    ("w", "VEK"),
    ("wa", "TEK"),
    ("wb", "TEK"),
    ("wa", "IEK"),
    ("wb", "IEK"),
    ("w", "LEK-[Amidated]"),
]:
    res, nterm, cterm = residues(sequence)
    residue = res[-1][0] if series.startswith("d") else res[0][0]
    parent = SERIES["a"] if series.startswith("d") else SERIES["z"]
    neutral = total(c for _, c in res) + nterm + cterm + parent - LOST[(series, residue)]
    for charge in (1, 2):
        suffix = f"^{charge}" if charge > 1 else ""
        add(f"{series}{len(res)}{{{sequence}}}{suffix}", mz_protonated(neutral, charge), charge, "4.4.3 side-chain ion")
for sequence in ["LEK", "GEK", "S[Phospho]EK", "PEK", "LEK-[Amidated]"]:
    res, nterm, cterm = residues(sequence)
    neutral = total(c for _, c in res) + cterm + SERIES["y"] - (AA[res[0][0]] - BACKBONE) - comp("H")
    if res[0][1] != AA[res[0][0]]:
        neutral -= res[0][1] - AA[res[0][0]]  # a side-chain modification leaves with the side chain
    add(f"v{len(res)}{{{sequence}}}", mz_protonated(neutral, 1), 1, "4.4.3 side-chain ion")
# Direct check of the spec table formulas for side-chain ions on a generic residue.
res, _, _ = residues("PEL")
assert abs(m(total(c for _, c in res[:2]) + comp("C2H4N")) - m(total(c for _, c in res) + SERIES["a"] - LOST[("d", "L")])) < 1e-9
res, _, _ = residues("LEK")
assert abs(m(total(c for _, c in res[1:]) + comp("C3H4O2")) - m(total(c for _, c in res) + SERIES["z"] - LOST[("w", "L")])) < 1e-9
assert abs(m(total(c for _, c in res[1:]) + comp("C2H3NO2")) - m(total(c for _, c in res) + SERIES["y"] - (AA["L"] - BACKBONE) - comp("H"))) < 1e-9

# Internal fragments (4.4.4) for MYPEPTIDEK: residues plus proton, yb by default.
res, _, _ = residues("MYPEPTIDEK")
for text, start, end, delta, charge in [
    ("m3:6", 3, 6, "", 1),
    ("m3:6-CO", 3, 6, "CO", 1),
    ("m3:6-CO-H2O^2", 3, 6, "COH2O", 2),
    ("m3:4", 3, 4, "", 1),
    ("m4:5", 4, 5, "", 1),
    ("m5:8-H2O", 5, 8, "H2O", 1),
]:
    neutral = total(c for _, c in res[start - 1 : end]) - comp(delta or "")
    add(text, mz_protonated(neutral, charge), charge, "4.4.4 internal fragment", "MYPEPTIDEK")

# Immonium ions (4.4.5): residue minus CO plus proton.
for aa in "ACDEFGHIKLMNPQRSTVWY":
    add(f"I{aa}", mz_protonated(AA[aa] - comp("CO"), 1), 1, "4.4.5 immonium")
for aa, mod in [("C", "Carbamidomethyl"), ("Y", "Phospho"), ("M", "Oxidation")]:
    add(f"I{aa}[{mod}]", mz_protonated(AA[aa] + comp(MODS[mod]) - comp("CO"), 1), 1, "4.4.5 modified immonium")
add("IL-CH2", mz_protonated(AA["L"] - comp("CO") - comp("CH2"), 1), 1, "4.4.5 immonium loss")
add("IC[+58.005]", mz_protonated(AA["C"] - comp("CO"), 1) + 58.005, 1, "4.4.5 immonium mass modification")

# Precursor ions (4.4.6) for a 4+ PEPTIDEK, using the adduct type column.
res, _, _ = residues("PEPTIDEK")
M = m(total(c for _, c in res) + comp("H2O"))
for text, hydrogens, charge in [
    ("p^4", 4, 4),
    ("p+H^3", 4, 3),
    ("p^3", 3, 3),
    ("p+2H^2", 4, 2),
    ("p^2", 2, 2),
    ("p+H^2", 3, 2),
    ("p+3H", 4, 1),
    ("p+2H", 3, 1),
    ("p+H", 2, 1),
    ("p", 1, 1),
]:
    # [M+kH]c+ carries k hydrogen nuclei and k - c electrons.
    add(text, (M + hydrogens * NIST["H"][1][0] - charge * ELECTRON) / charge, charge, "4.4.6 precursor table", "PEPTIDEK")
add("p-H3PO4^2", (M - m(comp("H3PO4")) + 2 * PROTON) / 2, 2, "4.4.6 precursor loss", "PEPTIDEK")
add("p-[Hex]", (M - m(comp(MODS["Hex"])) + PROTON), 1, "4.5 reference group loss (Unimod)", "PEPTIDEK")
add("p-[TMT6plex]-2H2O-HPO3", M - m(comp(MODS["TMT6plex"])) - 2 * m(comp("H2O")) - m(comp("HPO3")) + PROTON, 1, "4.5 loss chain", "PEPTIDEK")

# Reference ions (4.4.7, Appendix B): neutral formula plus proton.
refmol = json.loads(Path(__file__).with_name("reference_molecules.json").read_text())
for name, info in refmol.items():
    add(f"r[{name}]", mz_protonated(comp(info["chemical_formula"]), 1), 1, "Appendix B formula")
for name in ("Hex", "HexNAc(2)", "Phospho"):
    add(f"r[{name}]", mz_protonated(comp(MODS[name]), 1), 1, "4.4.7 Unimod reference name")
add("r[TMT127N][M+Na]", m(comp(refmol["TMT127N"]["chemical_formula"]) + comp("Na")) - ELECTRON, 1, "4.4.7 reference adduct")

# Neutral loss table (4.5): the tabulated exact mass must agree with pyteomics, then y1{K} minus each.
LOSS_TABLE = {
    "H": 1.007825,
    "NH3": 17.026549,
    "H2O": 18.010565,
    "CO": 27.994915,
    "CO2": 43.989829,
    "HCONH2": 45.021464,
    "HCOOH": 46.005479,
    "CH4OS": 63.998301,
    "SO3": 79.956818,
    "HPO3": 79.966331,
    "C2H5NOS": 91.009195,
    "C2H4O2S": 91.993211,
    "H3PO4": 97.976896,
}
y1k = total([AA["K"], comp("H2O")])
for formula, tabulated in LOSS_TABLE.items():
    # The printed table is rounded loosely (CH4OS is printed 63.998301, the formula gives 63.998285).
    assert abs(m(comp(formula)) - tabulated) < 2e-5, (formula, m(comp(formula)), tabulated)
    if abs(m(comp(formula)) - tabulated) > 1e-6:
        print(f"spec table differs: {formula} printed {tabulated} computed {m(comp(formula)):.6f}")
    add(f"y1{{K}}-{formula}", mz_protonated(y1k - comp(formula), 1), 1, "4.5 neutral loss table")
    add(f"y1{{K}}+{formula}", mz_protonated(y1k + comp(formula), 1), 1, "4.5 neutral gain")
add("y1{K}-2H2O", mz_protonated(y1k - comp("H2O") - comp("H2O"), 1), 1, "4.5 multiplied loss")
add("y2{EK}+CO-H2O", mz_protonated(total([AA["E"], AA["K"]]) + comp("H2O") + comp("CO") - comp("H2O"), 1), 1, "4.5 loss chain")
add("y2{EK}-[2H1]-NH3", mz_protonated(total([AA["E"], AA["K"]]) + comp("H2O") - comp("[2H1]") - comp("NH3"), 1), 1, "4.5 isotope loss")
add("y5{TIDEK}-H2[18O1][M+Na]", m(total(c for _, c in residues("TIDEK")[0]) + comp("H2O") - comp("H2[18O1]") + comp("Na")) - ELECTRON, 1, "4.5 isotope loss")
add("b3{PEP}-[Carbamidomethyl]", mz_protonated(total(c for _, c in residues("PEP")[0]) - comp(MODS["Carbamidomethyl"]), 1), 1, "4.4.4 Unimod loss")

# Isotopes (4.6).
C13 = NIST["C"][13][0] - NIST["C"][12][0]
N15 = NIST["N"][15][0] - NIST["N"][14][0]
H2 = NIST["H"][2][0] - NIST["H"][1][0]
O18 = NIST["O"][18][0] - NIST["O"][16][0]
y1 = mz_protonated(y1k, 1)
for text, shift in [
    ("+i", C13),
    ("+2i", 2 * C13),
    ("+3i", 3 * C13),
    ("-i", -C13),
    ("-2i", -2 * C13),
    ("+i13C", C13),
    ("+i15N", N15),
    ("+3i15N", 3 * N15),
    ("+i2H", H2),
    ("+2i18O", 2 * O18),
    ("+6i13C+2i15N", 6 * C13 + 2 * N15),
    ("+2i13C+i15N", 2 * C13 + N15),
]:
    add(f"y1{{K}}{text}", y1 + shift, 1, "4.6 isotope")
    add(f"y1{{K}}{text}^2", (y1 + PROTON + shift) / 2, 2, "4.6 isotope with charge")

# Adducts (4.7): neutral fragment plus adduct atoms minus one electron per charge.
y4 = total(c for _, c in residues("PEPT")[0]) + comp("H2O")
y5 = total(c for _, c in residues("TIDEK")[0]) + comp("H2O")
for text, neutral, adduct, charge in [
    ("y4{PEPT}[M+Na]", y4, "Na", 1),
    ("y4{PEPT}[M+NH4]", y4, "NH4", 1),
    ("y4{PEPT}[M+2Na]^2", y4, "Na2", 2),
    ("y4{PEPT}[M+2H+Na]^3", y4, "H2Na", 3),
    ("y4{PEPT}[M+H+Na]^2", y4, "HNa", 2),
    ("y4{PEPT}[M+H]", y4, "H", 1),
    ("y4{PEPT}[M+2H]^2", y4, "H2", 2),
    ("y5{TIDEK}-H2O[M+H+Na]^2", y5 - comp("H2O"), "HNa", 2),
    ("y4{PEPT}[M+[2H2]]^2", y4, "[2H2]", 2),
    ("y5{TIDEK}[M+[15N1]H4]", y5, "[15N1]H4", 1),
    ("y4{PEPT}-H2O+2i[M+H+Na]^2", y4 - comp("H2O"), "HNa", 2),
]:
    value = (m(neutral + comp(adduct)) - charge * ELECTRON) / charge
    if "+2i" in text:
        value += 2 * C13 / charge
    add(text, value, charge, "4.7 adduct")

# Chemical formulas (4.4.9): all nuclei of the charged species, minus electrons.
for formula, charge, iso in [
    ("C13H9", 1, 0),
    ("C12H9N", 1, 0),
    ("C13H9N", 1, 0),
    ("C13H10N", 1, 0),
    ("C13H11N", 1, 0),
    ("C13H12N", 1, 0),
    ("C14H10N", 1, 0),
    ("C14H11N", 1, 0),
    ("C14H10NO", 1, 0),
    ("C16H22O", 3, 1),
    ("C15[13C1]H22O", 3, 0),
]:
    text = f"f{{{formula}}}" + ("+i" if iso else "") + (f"^{charge}" if charge > 1 else "")
    add(text, (m(comp(formula)) + iso * C13 - charge * ELECTRON) / charge, charge, "4.4.9 formula")
add("f{C13H9}[M+H]", m(comp("C13H9")) - ELECTRON, 1, "4.4.9 formula adduct adds no atoms")

# SMILES (4.4.10): neutral molecule plus charge carriers.
for text, formula, adduct, charge in [
    ("s{CN=C=O}[M+H]", "C2H3NO", "H", 1),
    ("s{COc(c1)cccc1C#N}[M+H+Na]^2", "C8H7NO", "HNa", 2),
]:
    add(text, (m(comp(formula) + comp(adduct)) - charge * ELECTRON) / charge, charge, "4.4.10 smiles")

OUT.write_text(json.dumps(cases, indent=1) + "\n")
print(f"wrote {len(cases)} cases to {OUT}")
