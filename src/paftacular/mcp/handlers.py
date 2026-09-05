"""Stateless adapters around the public paftacular API."""

import json
import logging
import math
from collections.abc import Callable
from importlib.metadata import PackageNotFoundError, version
from importlib.util import find_spec
from typing import Any

from pydantic import BaseModel

import paftacular as pft
from paftacular.constants import IonSeries

from .guidance import CONVENTIONS, RESOURCES
from .models import (
    MAX_ANNOTATIONS,
    MAX_RECORDS,
    MAX_REQUEST_BYTES,
    MAX_RESULT_BYTES,
    MAX_TEXT_BYTES,
    AnnotationView,
    BatchRequest,
    BuildRequest,
    CalculatedBatch,
    CalculatedRecord,
    Calculation,
    CalculationRequest,
    Candidate,
    Capabilities,
    ContextRequest,
    Error,
    FragmentRequest,
    Matches,
    MatchRequest,
    ParsedBatch,
    ParsedRecord,
    ParseRequest,
    SerializeRequest,
)

logger = logging.getLogger(__name__)
SEQUENCE_IONS = (pft.PeptideIon, pft.InternalFragment, pft.PrecursorIon)


class RequestError(ValueError):
    def __init__(self, code: str, message: str):
        self.code = code
        super().__init__(message)


def error_info(error: Exception) -> Error:
    if isinstance(error, RequestError):
        return Error(code=error.code, message=str(error))
    if isinstance(error, pft.PafParseError):
        return Error(code="parse_error", message=error.reason, position=error.position, annotation_index=error.annotation_index)
    if isinstance(error, ImportError):
        return Error(code="missing_dependency", message=f"{error}. For MCP with SMILES, install 'paftacular[mcp,smiles]'.")
    if isinstance(error, NotImplementedError):
        return Error(code="unsupported_calculation", message=str(error) or "This calculation is not supported")
    if isinstance(error, KeyError):
        return Error(code="unknown_reference", message=f"Unknown chemistry reference: {error}")
    if isinstance(error, ValueError | OverflowError):
        return Error(code="invalid_input", message=str(error))
    logger.error("Unexpected MCP calculation failure", exc_info=(type(error), error, error.__traceback__))
    return Error(code="internal_error", message="The calculation failed unexpectedly. Check the server diagnostics.")


def check_size(value: Any) -> None:
    """Bound decoded application values, including nested interchange dictionaries."""
    data = value.model_dump(mode="json") if isinstance(value, BaseModel) else value
    if len(json.dumps(data, ensure_ascii=True, allow_nan=False).encode()) > MAX_REQUEST_BYTES:
        raise RequestError("request_limit", "Request exceeds 256 KiB. Split the batch or reduce the input.")
    pending = [(data, 0)]
    while pending:
        item, depth = pending.pop()
        if depth > 32:
            raise RequestError("request_limit", "Request nesting exceeds 32 levels")
        if isinstance(item, str) and len(item.encode("utf-8")) > MAX_TEXT_BYTES:
            raise RequestError("request_limit", "Each input text must be at most 16 KiB")
        if isinstance(item, dict):
            pending.extend((part, depth + 1) for pair in item.items() for part in pair)
        elif isinstance(item, list):
            pending.extend((part, depth + 1) for part in item)


def view(annotation: pft.PafAnnotation) -> AnnotationView:
    notices = []
    if annotation.resolved_sequence is not None:
        notices.append("Canonical mzPAF text omits resolved context. Retain the annotation dictionary to preserve it.")
    return AnnotationView(canonical=annotation.serialize(), annotation=annotation.to_dict(), sequence=annotation.sequence, notices=notices)


def capabilities() -> Capabilities:
    dependencies = {}
    for name in ("mcp", "peptacular", "pysmiles", "tacular"):
        try:
            dependencies[name] = version(name)
        except PackageNotFoundError:
            dependencies[name] = None
    return Capabilities(
        dependencies=dependencies,
        integrations={name: find_spec(name) is not None for name in ("peptacular", "pysmiles")},
        ion_series=[series.value for series in IonSeries],
        resolvable_series=["a", "b", "c", "x", "y", "z", "internal", "precursor"],
        limits={
            "records": MAX_RECORDS,
            "annotations": MAX_ANNOTATIONS,
            "text_bytes": MAX_TEXT_BYTES,
            "request_bytes": MAX_REQUEST_BYTES,
            "result_bytes": MAX_RESULT_BYTES,
        },
        guide_resources=list(RESOURCES),
        conventions=CONVENTIONS,
    )


def parse_annotations(request: ParseRequest) -> ParsedBatch:
    records = []
    count = 0
    for result in pft.iter_parse(request.records):
        count += len(result.annotations)
        if count > MAX_ANNOTATIONS:
            raise RequestError("result_limit", "Parsing exceeds 1000 annotations. Split the records into smaller calls.")
        records.append(
            ParsedRecord(
                index=result.index,
                text=result.text,
                annotations=[view(item) for item in result.annotations],
                error=error_info(result.error) if result.error is not None else None,
            )
        )
    return ParsedBatch(records=records)


def contextualize(request: ContextRequest) -> pft.PafAnnotation:
    annotation = pft.parse_single(request.annotation) if isinstance(request.annotation, str) else pft.PafAnnotation.from_dict(request.annotation)
    if request.analyte is not None:
        annotation = annotation.resolve(request.analyte)
    elif request.analytes is not None:
        annotation = annotation.resolve({item.reference: item.sequence for item in request.analytes})
    return annotation


def resolve_annotation(request: ContextRequest) -> AnnotationView:
    if request.analyte is None and request.analytes is None:
        raise RequestError("missing_context", "Supply analyte or analytes to resolve and validate the selected sequence")
    return view(contextualize(request))


def serialize_annotation(request: SerializeRequest) -> AnnotationView:
    return view(pft.PafAnnotation.from_dict(request.annotation))


def build_annotation(request: BuildRequest) -> AnnotationView:
    bare = pft.parse_single(request.ion)
    if bare != pft.PafAnnotation(bare.ion_type):
        raise RequestError("invalid_input", "ion must be a bare ion. Supply modifiers through their separate fields.")
    annotation = pft.PafAnnotation(
        ion_type=bare.ion_type,
        charge=request.charge,
        analyte_reference=request.analyte_reference,
        is_auxiliary=request.is_auxiliary,
        neutral_losses=tuple(pft.NeutralLoss.parse(value) for value in request.neutral_losses),
        isotopes=tuple(pft.IsotopeSpecification.parse(value) for value in request.isotopes),
        adducts=tuple(pft.Adduct.parse(value) for value in request.adducts),
        confidence=request.confidence,
        mass_error=pft.MassError(request.mass_error, request.mass_error_unit) if request.mass_error is not None else None,
    )
    return view(pft.PafAnnotation.from_dict(annotation.to_dict()))


def calculate_ion(request: CalculationRequest) -> Calculation:
    annotation = contextualize(request)
    sequence_needed = isinstance(annotation.ion_type, SEQUENCE_IONS)
    complete = request.mode == "complete"
    if complete and sequence_needed and annotation.sequence is None:
        raise RequestError("missing_context", "Supply the full ProForma analyte for a complete fragment mass, or explicitly request mode='offsets'.")
    if complete and isinstance(annotation.ion_type, pft.PeptideIon) and annotation.ion_type.series not in ("a", "b", "c", "x", "y", "z"):
        raise RequestError("unsupported_calculation", "Complete peptide calculations support a/b/c/x/y/z series")
    if not complete and not sequence_needed:
        raise RequestError("invalid_input", "Offset mode applies only to peptide, internal, and precursor ions")
    context = "none"
    if request.analyte is not None or request.analytes is not None:
        context = "analyte"
    elif annotation.sequence is not None:
        context = ("resolved" if annotation.resolved_sequence is not None else "embedded") if complete else "ignored"
    result = Calculation(
        status="success",
        ion=view(annotation),
        mode=request.mode,
        context_source=context,
        charge=annotation.charge,
        mass_basis="charged_species" if complete else "offsets_and_modifiers",
    )
    operations: dict[str, tuple[str, Callable[[], Any]]] = {
        "mass": ("mass_da", lambda: annotation.mass(calculate_sequence=complete)),
        "mz": ("mz_th", lambda: annotation.mz(calculate_sequence=complete)),
        "formula": ("formula", lambda: annotation.formula(calculate_sequence=complete)),
        "composition": ("composition", lambda: annotation.dict_composition(calculate_sequence=complete)),
    }
    succeeded = 0
    for prop in dict.fromkeys(request.properties):
        field, operation = operations[prop]
        try:
            value = operation()
            if isinstance(value, float) and not math.isfinite(value):
                raise RequestError("invalid_result", "Calculation did not produce a finite number")
            setattr(result, field, value)
            succeeded += 1
        except Exception as error:
            result.property_errors[prop] = error_info(error)
    if result.property_errors:
        result.status = "partial" if succeeded else "error"
    return result


def calculate_ions(request: BatchRequest) -> CalculatedBatch:
    records = []
    for index, item in enumerate(request.requests):
        try:
            records.append(CalculatedRecord(index=index, result=calculate_ion(item)))
        except Exception as error:
            records.append(CalculatedRecord(index=index, error=error_info(error)))
    return CalculatedBatch(records=records)


def generate_fragments(request: FragmentRequest) -> CalculatedBatch:
    import peptacular as pt

    length = len(pt.parse(request.analyte).sequence)
    positions = list(dict.fromkeys(request.positions)) if request.positions is not None else range(1, length)
    series = list(dict.fromkeys(request.series))
    charges = list(dict.fromkeys(request.charges))
    if any(position >= length for position in positions):
        raise RequestError("invalid_input", "Generated terminal fragments must be shorter than the analyte")
    if not positions:
        raise RequestError("invalid_input", "The analyte has no proper terminal fragments")
    if len(positions) * len(series) * len(charges) > MAX_RECORDS:
        raise RequestError("request_limit", "Generation exceeds 100 fragments. Select fewer positions, series, or charges.")
    requests = [
        CalculationRequest(
            annotation=pft.PafAnnotation.make_peptide(ion, position, charge=charge).serialize(), analyte=request.analyte, properties=request.properties
        )
        for ion in series
        for position in positions
        for charge in charges
    ]
    return calculate_ions(BatchRequest(requests=requests))


def match_mz(request: MatchRequest) -> Matches:
    candidates = []
    for index, item in enumerate(request.candidates):
        candidate = Candidate(index=index)
        try:
            result = calculate_ion(CalculationRequest(**item.model_dump(), properties=["mz"]))
            candidate.result = result
            if result.mz_th is None:
                candidate.error = result.property_errors["mz"]
            elif result.mz_th <= 0:
                raise RequestError("invalid_result", "Theoretical m/z must be positive for matching")
            else:
                delta = request.observed_mz - result.mz_th
                ppm = delta / result.mz_th * 1e6
                if not math.isfinite(ppm):
                    raise RequestError("invalid_result", "Mass error is outside the finite numeric range")
                candidate.delta_th = delta
                candidate.delta_ppm = ppm
                candidate.matched = abs(ppm if request.tolerance_unit == "ppm" else delta) <= request.tolerance
        except Exception as error:
            candidate.error = error_info(error)
        candidates.append(candidate)
    matched = sorted((item for item in candidates if item.matched), key=lambda item: abs(item.delta_ppm or 0))
    return Matches(
        observed_mz=request.observed_mz,
        tolerance=request.tolerance,
        tolerance_unit=request.tolerance_unit,
        candidates=candidates,
        matching_indices=[item.index for item in matched],
    )
