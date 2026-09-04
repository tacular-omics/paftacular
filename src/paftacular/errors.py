class PafParseError(ValueError):
    """An mzPAF error with zero-based character and annotation positions."""

    def __init__(self, text: str, position: int, annotation_index: int, reason: str):
        self.text = text
        self.position = position
        self.annotation_index = annotation_index
        self.reason = reason
        super().__init__(f"Invalid mzPAF annotation at position {position} (annotation {annotation_index}): {reason}")

    def __reduce__(self):
        return type(self), (self.text, self.position, self.annotation_index, self.reason)
