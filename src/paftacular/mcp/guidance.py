"""Small, bundled references for hosts and models. No remote content is loaded."""

CONVENTIONS = [
    "Complete peptide, internal, and precursor calculations require sequence context. Offsets must be requested explicitly.",
    "mass_da is the charged species mass in Da. mz_th is mass divided by positive charge, in Th.",
    "All MCP calculations are monoisotopic. Average isotopomers are unsupported.",
    "Composition counts nuclei. Explicit adduct mass includes electron correction. Formula ions already contain the charged species atoms.",
    "Resolved context survives dictionary interchange but is absent from mzPAF text.",
    "Internal positions are one-based and inclusive. Resolved internal fragments exclude both analyte termini.",
    "The mzPAF internal-cleavage correction table differs from physical ion definitions. Do not infer cleavage identity from a correction.",
    "A mass match alone does not establish a fragment identity. Confidence values are supplied metadata, not computed probabilities.",
]

GUIDE = """paftacular MCP workflow guide

Start with get_capabilities to discover installed integrations and limits.
Use parse_annotations to inspect, validate, and normalize one or more mzPAF records.
Each record may contain comma-separated annotations. Failures retain their position.
Parsing validates notation and structural constraints, not all chemical references.

Use calculate_ion with annotation and a full ProForma analyte for mass or m/z.
For example y3^2 with PEPTIDE selects IDE. Do not pass IDE as the full analyte
unless IDE actually is the analyte. The annotation charge controls the calculation.
Use analytes with explicit reference and sequence fields for mixtures.
An omitted annotation reference selects analyte 1. Supplied and embedded context
must agree. Use resolve_annotation to retain this context without calculating.
Calculation and matching requests also accept a complete annotation dictionary
returned by another tool. This preserves resolved context across tool calls.

Use build_annotation to attach typed charge, losses, isotopes, adducts, mass error,
and confidence to a bare ion. Use serialize_annotation with a complete dictionary
returned by another tool. Its schema_version is distinct from the MCP response version.

Use generate_fragments to enumerate a/b/c/x/y/z fragments for a ProForma analyte.
It generates proper terminal fragments, excluding the full-length peptide.
Use positions and charges to keep the result below the 100-fragment limit.
Use calculate_ions for custom lists, including internal ions and neutral losses.

Use match_mz to compare an observed m/z against explicit candidate annotations.
It retains all candidates and errors, and ranks matching indices by absolute ppm error.
It does not search a sequence database, identify a peptide, or estimate confidence.

Request formula or composition only when needed. A numeric mass modification may
allow mass and m/z while preventing composition. Inspect status and property_errors.
Missing information is reported as an error, never as an invented value.
Resources contain scientific conventions and runnable tool-request examples.
"""

EXAMPLES = """Example MCP tool arguments

calculate_ion:
{"request":{"annotation":"y3^2","analyte":"PEPTIDE","properties":["mass","mz","composition"]}}

parse_annotations:
{"request":{"records":["y3^2,b2-H2O","invalid"]}}

build_annotation:
{"request":{"ion":"y3","charge":2,"neutral_losses":["-H2O"],"isotopes":["+i"]}}

resolve_annotation:
{"request":{"annotation":"2@y3","analytes":[{"reference":2,"sequence":"PEPTIDE"}]}}

generate_fragments:
{"request":{"analyte":"PEPTIDE","series":["b","y"],"charges":[1,2]}}

match_mz:
{"request":{"observed_mz":188.589359,"tolerance":10.0,"tolerance_unit":"ppm","candidates":[{"annotation":"y3^2","analyte":"PEPTIDE"}]}}
"""

RESOURCES = {
    "paftacular://guide": GUIDE,
    "paftacular://conventions": "\n".join(CONVENTIONS),
    "paftacular://examples": EXAMPLES,
}
