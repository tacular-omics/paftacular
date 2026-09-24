"""Exceptions raised by paftacular.

Every error caused by user input is a :class:`PaftacularError`, which is a ``ValueError``.
"""

from collections.abc import Callable
from functools import wraps


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


class PafUnsupportedCalculationError(PaftacularError):
    """A mass or composition that the annotation does not define.

    Raised for ``?`` (unannotated) and ``_{...}`` (named compound) ions, which carry no chemistry.
    """


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


def reraise_as_paftacular[**P, R](func: Callable[P, R]) -> Callable[P, R]:
    """Re-raise a ValueError from tacular or peptacular as a PaftacularError.

    peptacular parses sequences and modifications lazily, so its errors fire inside mass and
    composition calls, not where the text is first read. Wrap those calls with this.
    """

    @wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        try:
            return func(*args, **kwargs)
        except PaftacularError:
            raise
        except ValueError as error:
            raise PaftacularError(f"{type(error).__name__}: {error}") from error

    return wrapper
