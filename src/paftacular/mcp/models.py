"""JSON contracts for AI clients, separate from the core interchange schema."""

from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from paftacular import __version__

MAX_RECORDS = 100
MAX_ANNOTATIONS = 1000
MAX_TEXT_BYTES = 16 * 1024
MAX_REQUEST_BYTES = 256 * 1024
MAX_RESULT_BYTES = 512 * 1024

Text = Annotated[str, Field(min_length=1, max_length=MAX_TEXT_BYTES)]
PositiveInt = Annotated[int, Field(ge=1, le=10000)]
Reference = Annotated[int, Field(ge=0, le=10000)]
Property = Literal["mass", "mz", "formula", "composition"]
BackboneSeries = Literal["a", "b", "c", "x", "y", "z"]


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True, allow_inf_nan=False)


class Analyte(Model):
    reference: Reference
    sequence: Text


class ContextRequest(Model):
    annotation: Text | dict[str, Any] = Field(description="mzPAF text or a complete schema_version 1 dictionary returned by another tool.")
    analyte: Text | None = Field(default=None, description="Full ProForma analyte, not the fragment sequence.")
    analytes: Annotated[list[Analyte], Field(min_length=1, max_length=MAX_RECORDS)] | None = None

    @model_validator(mode="after")
    def check_context(self):
        if self.analyte is not None and self.analytes is not None:
            raise ValueError("Supply analyte or analytes, not both")
        if self.analytes is not None and len({item.reference for item in self.analytes}) != len(self.analytes):
            raise ValueError("Analyte references must be unique")
        return self


class CalculationRequest(ContextRequest):
    properties: Annotated[list[Property], Field(min_length=1, max_length=4)] = Field(default_factory=lambda: ["mass", "mz"])
    mode: Literal["complete", "offsets"] = "complete"

    @model_validator(mode="after")
    def check_mode(self):
        if self.mode == "offsets" and (self.analyte is not None or self.analytes is not None):
            raise ValueError("Offset calculations cannot take analyte context")
        return self


class ParseRequest(Model):
    records: Annotated[list[Annotated[str, Field(max_length=MAX_TEXT_BYTES)]], Field(min_length=1, max_length=MAX_RECORDS)]


class BatchRequest(Model):
    requests: Annotated[list[CalculationRequest], Field(min_length=1, max_length=MAX_RECORDS)]


class SerializeRequest(Model):
    annotation: dict[str, Any] = Field(description="Complete schema_version 1 dictionary returned by parse or resolve tools.")


class BuildRequest(Model):
    ion: Text = Field(description="Bare mzPAF ion such as y3, m2:4, p, IM[Oxidation], or f{C2H4}.")
    charge: PositiveInt = 1
    analyte_reference: Reference | None = None
    is_auxiliary: bool = False
    neutral_losses: Annotated[list[Text], Field(max_length=32)] = Field(default_factory=list)
    isotopes: Annotated[list[Text], Field(max_length=32)] = Field(default_factory=list)
    adducts: Annotated[list[Text], Field(max_length=32)] = Field(default_factory=list)
    confidence: Annotated[float, Field(ge=0, le=1)] | None = None
    mass_error: float | None = None
    mass_error_unit: Literal["da", "ppm"] = "da"


class FragmentRequest(Model):
    analyte: Text
    series: Annotated[list[BackboneSeries], Field(min_length=1, max_length=6)] = Field(default_factory=lambda: ["b", "y"])
    charges: Annotated[list[PositiveInt], Field(min_length=1, max_length=10)] = Field(default_factory=lambda: [1])
    positions: Annotated[list[PositiveInt], Field(min_length=1, max_length=MAX_RECORDS)] | None = None
    properties: Annotated[list[Property], Field(min_length=1, max_length=4)] = Field(default_factory=lambda: ["mass", "mz"])


class MatchRequest(Model):
    observed_mz: Annotated[float, Field(gt=0)]
    candidates: Annotated[list[ContextRequest], Field(min_length=1, max_length=MAX_RECORDS)]
    tolerance: Annotated[float, Field(ge=0)] = 10.0
    tolerance_unit: Literal["ppm", "Th"] = "ppm"


class Error(Model):
    code: str
    message: str
    position: int | None = None
    annotation_index: int | None = None


class Envelope[T](Model):
    response_schema_version: Literal[1] = 1
    package_version: str = __version__
    data: T | None = None
    error: Error | None = None


class AnnotationView(Model):
    canonical: str
    annotation: dict[str, Any]
    sequence: str | None = None
    notices: list[str] = Field(default_factory=list)


class ParsedRecord(Model):
    index: int
    text: str
    annotations: list[AnnotationView] = Field(default_factory=list)
    error: Error | None = None


class ParsedBatch(Model):
    records: list[ParsedRecord]


class Calculation(Model):
    status: Literal["success", "partial", "error"]
    ion: AnnotationView
    mode: Literal["complete", "offsets"]
    context_source: Literal["analyte", "embedded", "resolved", "none", "ignored"]
    charge: int
    mass_type: Literal["monoisotopic"] = "monoisotopic"
    mass_basis: Literal["charged_species", "offsets_and_modifiers"]
    mass_da: float | None = None
    mz_th: float | None = None
    formula: str | None = None
    composition: dict[str, int] | None = None
    property_errors: dict[str, Error] = Field(default_factory=dict)


class CalculatedRecord(Model):
    index: int
    result: Calculation | None = None
    error: Error | None = None


class CalculatedBatch(Model):
    records: list[CalculatedRecord]


class Candidate(Model):
    index: int
    result: Calculation | None = None
    error: Error | None = None
    delta_th: float | None = None
    delta_ppm: float | None = None
    matched: bool = False


class Matches(Model):
    observed_mz: float
    tolerance: float
    tolerance_unit: Literal["ppm", "Th"]
    candidates: list[Candidate]
    matching_indices: list[int]
    interpretation: str = "Mass agreement is a candidate filter, not proof of fragment identity. Deltas are observed minus theoretical."


class Capabilities(Model):
    dependencies: dict[str, str | None]
    integrations: dict[str, bool]
    ion_series: list[str]
    resolvable_series: list[str]
    limits: dict[str, int]
    mass_types: list[str] = Field(default_factory=lambda: ["monoisotopic"])
    transport: str = "stdio"
    guide_resources: list[str]
    conventions: list[str]
