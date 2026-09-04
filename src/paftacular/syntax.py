from collections.abc import Iterator

from .errors import PafParseError


def annotation_spans(text: str) -> Iterator[tuple[int, int]]:
    """Locate separators outside labels and balanced embedded ProForma."""
    start = 0
    index = 0
    brackets = 0
    braces = 0
    sequence = False
    for position, character in enumerate(text):
        if brackets:
            if character == "[":
                brackets += 1
            elif character == "]":
                brackets -= 1
        elif braces:
            if character == "[" and sequence:
                brackets = 1
            elif character == "{" and sequence:
                braces += 1
            elif character == "}":
                braces -= 1
        elif character == "[":
            brackets = 1
        elif character == "{":
            braces = 1
            sequence = position > start and text[position - 1].isdigit()
        elif character == ",":
            if not text[start:position].strip():
                raise PafParseError(text, position, index, "Empty annotation")
            yield start, position
            start = position + 1
            index += 1
    if braces or brackets:
        raise PafParseError(text, len(text), index, "Unclosed annotation delimiter")
    if not text[start:].strip():
        if index:
            raise PafParseError(text, len(text), index, "Trailing separator")
        return
    yield start, len(text)
