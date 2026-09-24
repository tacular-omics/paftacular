"""Exceptions raised by paftacular.

Every error caused by user input is a :class:`PaftacularError`, which is a ``ValueError``.
"""


class PaftacularError(ValueError):
    """Base class for paftacular errors: invalid annotations, components or calculations."""


class PafParseError(PaftacularError):
    """An mzPAF error with zero-based character and annotation positions."""

    def __init__(self, text: str, position: int, annotation_index: int, reason: str):
        self.text = text
        self.position = position
        self.annotation_index = annotation_index
        self.reason = reason
        super().__init__(f"Invalid mzPAF annotation at position {position} (annotation {annotation_index}): {reason}")

    def __reduce__(self):
        return type(self), (self.text, self.position, self.annotation_index, self.reason)


class PafUnknownReferenceError(PaftacularError, KeyError):
    """A reference molecule name that is not in the mzPAF reference list or Unimod.

    Raised when calculating an ``r[...]`` ion or a ``-[...]`` loss. It is a
    :class:`PaftacularError` (so a ``ValueError``) and also a ``KeyError``.
    """

    def __init__(self, name: str, reason: str | None = None):
        self.name = name
        self.reason = reason or f"Unknown reference molecule '{name}': not in the mzPAF reference list or Unimod"
        super().__init__(self.reason)

    def __str__(self) -> str:
        # KeyError.__str__ would wrap the message in quotes.
        return self.reason

    def __reduce__(self):
        return type(self), (self.name, self.reason)
