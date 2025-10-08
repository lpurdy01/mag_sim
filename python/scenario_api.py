"""Lightweight helpers for authoring simulation scenarios.

The API is intentionally tiny so agents and developers can build scenarios in
Python, validate them locally, and emit the JSON that the C++ ingestor
understands. The schema mirrors `ScenarioSpec` in the C++ layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Dict, List, Optional, Set, Union


@dataclass
class Domain:
    """Rectangular simulation window aligned with the solver grid."""

    Lx: float
    Ly: float
    nx: int
    ny: int

    def to_dict(self) -> Dict[str, object]:
        return {"Lx": float(self.Lx), "Ly": float(self.Ly), "nx": int(self.nx), "ny": int(self.ny)}


@dataclass
class Material:
    name: str
    mu_r: float

    def to_dict(self) -> Dict[str, object]:
        return {"name": self.name, "mu_r": float(self.mu_r)}


@dataclass
class UniformRegion:
    material: str

    def to_dict(self) -> Dict[str, object]:
        return {"type": "uniform", "material": self.material}


@dataclass
class Wire:
    x: float
    y: float
    radius: float
    I: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "type": "wire",
            "x": float(self.x),
            "y": float(self.y),
            "radius": float(self.radius),
            "I": float(self.I),
        }


@dataclass
class FieldMapOutput:
    """Request to write the full-field map to disk."""

    id: str
    path: Optional[str] = None
    quantity: str = "B"
    format: str = "csv"

    def validate(self) -> None:
        if not self.id:
            raise ValueError("FieldMapOutput id must be non-empty")
        if self.quantity != "B":
            raise ValueError("FieldMapOutput currently only supports quantity 'B'")
        if self.format != "csv":
            raise ValueError("FieldMapOutput currently only supports CSV format")
        if self.path is not None and not self.path:
            raise ValueError("FieldMapOutput path must be a non-empty string when provided")

    def to_dict(self) -> Dict[str, object]:
        data: Dict[str, object] = {
            "type": "field_map",
            "id": self.id,
            "quantity": self.quantity,
        }
        if self.path:
            data["path"] = self.path
        if self.format != "csv":
            data["format"] = self.format
        return data


@dataclass
class LineProbeOutput:
    """Sample the field along a horizontal or vertical line."""

    id: str
    axis: str
    value: float
    path: Optional[str] = None
    quantity: str = "Bmag"
    format: str = "csv"

    def validate(self) -> None:
        if not self.id:
            raise ValueError("LineProbeOutput id must be non-empty")
        if self.axis not in {"x", "y"}:
            raise ValueError("LineProbeOutput axis must be 'x' or 'y'")
        if self.quantity not in {"Bmag", "Bx", "By"}:
            raise ValueError("LineProbeOutput quantity must be 'Bmag', 'Bx', or 'By'")
        if self.format != "csv":
            raise ValueError("LineProbeOutput currently only supports CSV format")
        if self.path is not None and not self.path:
            raise ValueError("LineProbeOutput path must be a non-empty string when provided")

    def to_dict(self) -> Dict[str, object]:
        data: Dict[str, object] = {
            "type": "line_probe",
            "id": self.id,
            "axis": self.axis,
            "value": float(self.value),
            "quantity": self.quantity,
        }
        if self.path:
            data["path"] = self.path
        if self.format != "csv":
            data["format"] = self.format
        return data


ScenarioOutput = Union[FieldMapOutput, LineProbeOutput]


@dataclass
class Scenario:
    domain: Domain
    materials: List[Material] = field(default_factory=list)
    regions: List[UniformRegion] = field(default_factory=list)
    sources: List[Wire] = field(default_factory=list)
    outputs: List["ScenarioOutput"] = field(default_factory=list)
    units: str = "SI"
    version: str = "0.1"

    def _validate(self) -> None:
        if self.version != "0.1":
            raise ValueError("Only scenario version '0.1' is supported")
        if self.units != "SI":
            raise ValueError("Only SI units are currently supported")
        if not self.materials:
            raise ValueError("At least one material is required")
        names = {mat.name for mat in self.materials}
        if len(names) != len(self.materials):
            raise ValueError("Material names must be unique")
        if not self.regions:
            raise ValueError("At least one region must be provided")
        for region in self.regions:
            if region.material not in names:
                raise ValueError(f"Region references unknown material '{region.material}'")
            if not isinstance(region, UniformRegion):
                raise ValueError("Only UniformRegion is supported in v0.1")
        for wire in self.sources:
            if wire.radius <= 0.0:
                raise ValueError("Wire radius must be positive")
        seen_ids: Set[str] = set()
        for output in self.outputs:
            if not isinstance(output, (FieldMapOutput, LineProbeOutput)):
                raise ValueError(f"Unsupported output type: {type(output)}")
            output.validate()
            if output.id in seen_ids:
                raise ValueError(f"Duplicate output id '{output.id}'")
            seen_ids.add(output.id)

    def to_dict(self) -> Dict[str, object]:
        self._validate()
        data: Dict[str, object] = {
            "version": self.version,
            "units": self.units,
            "domain": self.domain.to_dict(),
            "materials": [mat.to_dict() for mat in self.materials],
            "regions": [region.to_dict() for region in self.regions],
            "sources": [wire.to_dict() for wire in self.sources],
        }
        if self.outputs:
            data["outputs"] = [output.to_dict() for output in self.outputs]
        return data

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent) + "\n"

    def save_json(self, path: Path | str) -> Path:
        data = self.to_json()
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(data, encoding="utf-8")
        return output_path


__all__ = [
    "Domain",
    "Material",
    "UniformRegion",
    "Wire",
    "FieldMapOutput",
    "LineProbeOutput",
    "Scenario",
]
