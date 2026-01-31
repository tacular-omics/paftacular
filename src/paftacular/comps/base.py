"""Base classes and protocols for mzPAF components"""

from abc import ABC, abstractmethod
from collections import Counter

from tacular import ElementInfo


class Serializable(ABC):
    """Base class for serializable objects"""

    @abstractmethod
    def serialize(self) -> str:
        pass

    def __str__(self) -> str:
        return self.serialize()


class MassProvider(ABC):
    """Base class for objects that can provide mass"""

    @abstractmethod
    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass"""
        pass

    @property
    def monoisotopic_mass(self) -> float:
        return self.mass(monoisotopic=True)

    @property
    def average_mass(self) -> float:
        return self.mass(monoisotopic=False)


class CompositionProvider(ABC):
    """Base class for objects that can provide composition"""

    @property
    @abstractmethod
    def composition(self) -> Counter["ElementInfo"]:
        """Get elemental composition"""
        pass

    @property
    def dict_composition(self) -> dict[str, int]:
        """Get composition as a dictionary with element symbols as keys"""
        return {str(elem): count for elem, count in self.composition.items()}

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass from composition"""
        m = 0.0
        for elem, count in self.composition.items():
            m += elem.get_mass(monoisotopic) * count
        return m


class ScalableComposition(CompositionProvider):
    """Mixin for compositions that scale by count and sign"""

    count: int

    @property
    @abstractmethod
    def _single_composition(self) -> Counter["ElementInfo"]:
        """Get composition for single instance (before scaling)"""
        pass

    @property
    def composition(self) -> Counter["ElementInfo"]:
        """Get scaled composition"""
        return Counter({elem: count * self.count for elem, count in self._single_composition.items()})

    @property
    def _sign_prefix(self) -> str:
        """Get sign prefix for serialization"""
        sign_str = "+" if self.count > 0 else "-"
        count_str = "" if abs(self.count) == 1 else str(abs(self.count))
        return f"{sign_str}{count_str}"
