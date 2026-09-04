import copy
import json
import pickle

import pytest

import paftacular as p


@pytest.mark.parametrize(
    "text",
    [
        "&2@y2{PE}-2H2O+i13C[M+H+Na]^2/0.000001ppm*0.123456789",
        "m2:4",
        "IK[Acetyl]",
        "r[Ref,X]",
        "_{Urocanic Acid}",
        "s{[13CH4]}",
        "f{C2H4}",
        "?2",
        "p",
    ],
)
def test_structured_json_roundtrip(text):
    annotation = p.parse_single(text)
    restored = p.PafAnnotation.from_dict(json.loads(json.dumps(annotation.to_dict())))
    assert restored == annotation
    assert restored.serialize() == annotation.serialize()


def test_structured_internal_preserves_backbone():
    annotation = p.PafAnnotation(p.InternalFragment(2, 3, None, p.IonSeries.B, p.IonSeries.X))
    assert p.PafAnnotation.from_dict(annotation.to_dict()) == annotation


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", True),
        ("schema_version", 2),
        ("charge", True),
        ("charge", 1.2),
        ("charge", 0),
        ("confidence", float("nan")),
        ("confidence", float("inf")),
        ("is_auxiliary", "false"),
        ("analyte_reference", -1),
        ("resolved_sequence", 12),
        ("neutral_losses", "-H2O"),
        ("mass_error", {"value": 1, "unit": "nonsense"}),
    ],
)
def test_reject_invalid_interchange_fields(field, value):
    data = p.parse_single("y2").to_dict()
    data[field] = value
    with pytest.raises(ValueError):
        p.PafAnnotation.from_dict(data)


@pytest.mark.parametrize("change", ["unknown", "missing", "ion_type", "ion_extra", "ion_position", "modifier"])
def test_reject_invalid_interchange_structure(change):
    data = p.parse_single("y2-H2O").to_dict()
    if change == "unknown":
        data["surprise"] = True
    elif change == "missing":
        del data["confidence"]
    elif change == "ion_type":
        data["ion"]["type"] = "NotAnIon"
    elif change == "ion_extra":
        data["ion"]["surprise"] = True
    elif change == "ion_position":
        data["ion"]["position"] = "2"
    else:
        data["neutral_losses"][0]["count"] = True
    with pytest.raises(ValueError):
        p.PafAnnotation.from_dict(data)


def test_structured_export_is_independent():
    annotation = p.parse_single("y2-H2O")
    before = annotation.to_dict()
    changed = annotation.to_dict()
    changed["neutral_losses"][0]["count"] = -2
    assert annotation.to_dict() == before


def test_legacy_export_stays_unversioned():
    data = p.parse_single("y2").as_dict()
    assert data["ion"] == "y2"
    assert "schema_version" not in data


def test_batch_retains_errors_and_record_indices():
    results = p.parse_batch(["y2,b3", "y2,b3!", "p", ""])
    assert [result.ok for result in results] == [True, False, True, True]
    assert [result.index for result in results] == list(range(4))
    assert len(results[0].annotations) == 2
    error = results[1].error
    assert isinstance(error, ValueError)
    assert error.position == 5
    assert error.annotation_index == 1
    assert error.text == "y2,b3!"
    assert results[1].annotations == ()
    assert results[3].annotations == ()
    assert pickle.loads(pickle.dumps(error)).reason == error.reason


def test_iter_parse_is_lazy():
    def records():
        yield "y2"
        raise RuntimeError("Read too far")

    iterator = p.iter_parse(records())
    assert next(iterator).ok
    with pytest.raises(RuntimeError, match="Read too far"):
        next(iterator)


@pytest.mark.parametrize("text", ["y2-H2O", "r[TMT126]", "IK", "_{sample}"])
def test_copy_components(text):
    annotation = p.parse_single(text)
    assert copy.copy(annotation) == annotation
    assert copy.deepcopy(annotation) == annotation
