"""Lightweight helpers for authoring simulation scenarios.

The API is intentionally tiny so agents and developers can build scenarios in
Python, validate them locally, and emit the JSON that the C++ ingestor
understands. The schema mirrors `ScenarioSpec` in the C++ layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union


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
class HalfspaceRegion:
    normal_x: float
    normal_y: float
    offset: float
    material: str

    def to_dict(self) -> Dict[str, object]:
        return {
            "type": "halfspace",
            "normal": [float(self.normal_x), float(self.normal_y)],
            "offset": float(self.offset),
            "material": self.material,
        }


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
class Rotor:
    name: Optional[str] = None
    pivot: Tuple[float, float] = (0.0, 0.0)
    polygon_indices: List[int] = field(default_factory=list)
    magnet_indices: List[int] = field(default_factory=list)
    wire_indices: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        data: Dict[str, object] = {}
        if self.name:
            data["name"] = self.name
        px, py = self.pivot
        if px != 0.0 or py != 0.0:
            data["pivot"] = [float(px), float(py)]
        if self.polygon_indices:
            data["polygon_indices"] = [int(index) for index in self.polygon_indices]
        if self.magnet_indices:
            data["magnet_indices"] = [int(index) for index in self.magnet_indices]
        if self.wire_indices:
            data["wire_indices"] = [int(index) for index in self.wire_indices]
        return data


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
        if self.quantity not in {"B", "H", "BH", "energy_density"}:
            raise ValueError(
                "FieldMapOutput quantity must be 'B', 'H', 'BH', or 'energy_density'"
            )
        if self.format not in {"csv", "vti"}:
            raise ValueError("FieldMapOutput format must be 'csv' or 'vti'")
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
        if self.quantity not in {"Bmag", "Bx", "By", "Hx", "Hy", "Hmag", "energy_density"}:
            raise ValueError(
                "LineProbeOutput quantity must be one of 'Bmag', 'Bx', 'By', 'Hx', 'Hy', 'Hmag', or 'energy_density'"
            )
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


@dataclass
class StressTensorProbeOutput:
    """Closed contour probe evaluated via the Maxwell stress tensor."""

    id: str
    vertices: List[Tuple[float, float]]
    path: Optional[str] = None
    quantity: str = "torque"
    method: str = "stress_tensor"

    def validate(self) -> None:
        if not self.id:
            raise ValueError("StressTensorProbeOutput id must be non-empty")
        if self.method != "stress_tensor":
            raise ValueError("StressTensorProbeOutput only supports method='stress_tensor'")
        if self.quantity not in {"force", "torque", "force_and_torque"}:
            raise ValueError(
                "StressTensorProbeOutput quantity must be 'force', 'torque', or 'force_and_torque'"
            )
        if len(self.vertices) < 3:
            raise ValueError("StressTensorProbeOutput requires at least three vertices")
        for vertex in self.vertices:
            if len(vertex) != 2:
                raise ValueError("StressTensorProbeOutput vertices must be (x, y) tuples")
        if self.path is not None and not self.path:
            raise ValueError("StressTensorProbeOutput path must be a non-empty string when provided")

    def to_dict(self) -> Dict[str, object]:
        vertices = [[float(x), float(y)] for x, y in self.vertices]
        data: Dict[str, object] = {
            "type": "probe",
            "id": self.id,
            "probe_type": self.quantity,
            "method": self.method,
            "loop": {"type": "polygon", "vertices": vertices},
        }
        if self.path:
            data["path"] = self.path
        return data


@dataclass
class BackEmfProbeOutput:
    """Flux-change probe producing a back-EMF CSV across timeline frames."""

    id: str
    component: str = "Bmag"
    vertices: Optional[List[Tuple[float, float]]] = None
    rect: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None
    frames: Optional[List[int]] = None
    path: Optional[str] = None

    def validate(self) -> None:
        if not self.id:
            raise ValueError("BackEmfProbeOutput id must be non-empty")
        normalized = self.component.lower()
        if normalized not in {"bx", "by", "bmag", "b"}:
            raise ValueError("BackEmfProbeOutput component must be 'Bx', 'By', or 'Bmag'")
        has_vertices = self.vertices is not None and len(self.vertices) >= 3
        has_rect = self.rect is not None and len(self.rect) == 2
        if has_vertices and has_rect:
            raise ValueError("BackEmfProbeOutput cannot define both vertices and rect")
        if not has_vertices and not has_rect:
            raise ValueError("BackEmfProbeOutput requires either vertices or rect definition")
        if has_vertices:
            for vertex in self.vertices or []:
                if len(vertex) != 2:
                    raise ValueError("BackEmfProbeOutput vertices must be (x, y) tuples")
        if has_rect:
            x_range, y_range = self.rect  # type: ignore[misc]
            if len(x_range) != 2 or len(y_range) != 2:
                raise ValueError("BackEmfProbeOutput rect must provide (x_min, x_max) and (y_min, y_max)")
        if self.frames is not None:
            if len(self.frames) < 2:
                raise ValueError("BackEmfProbeOutput frames must contain at least two indices")
            for value in self.frames:
                if value < 0:
                    raise ValueError("BackEmfProbeOutput frame indices must be non-negative")
        if self.path is not None and not self.path:
            raise ValueError("BackEmfProbeOutput path must be a non-empty string when provided")

    def to_dict(self) -> Dict[str, object]:
        component_map = {"bx": "Bx", "by": "By", "bmag": "Bmag", "b": "Bmag"}
        component_key = component_map[self.component.lower()]
        data: Dict[str, object] = {"type": "back_emf_probe", "id": self.id, "component": component_key}
        if self.vertices is not None:
            vertices = [[float(x), float(y)] for x, y in self.vertices]
            data["region"] = {"type": "polygon", "vertices": vertices}
        elif self.rect is not None:
            (xmin, xmax), (ymin, ymax) = self.rect  # type: ignore[misc]
            data["region"] = {
                "type": "rect",
                "x_range": [float(xmin), float(xmax)],
                "y_range": [float(ymin), float(ymax)],
            }
        if self.frames is not None:
            data["frames"] = [int(index) for index in self.frames]
        if self.path:
            data["path"] = self.path
        return data


ScenarioOutput = Union[FieldMapOutput, LineProbeOutput, StressTensorProbeOutput, BackEmfProbeOutput]


@dataclass
class MagnetRectRegion:
    x_range: Tuple[float, float]
    y_range: Tuple[float, float]
    magnetization: Tuple[float, float]

    def to_dict(self) -> Dict[str, object]:
        return {
            "type": "rect",
            "x_range": [float(self.x_range[0]), float(self.x_range[1])],
            "y_range": [float(self.y_range[0]), float(self.y_range[1])],
            "magnetization": [float(self.magnetization[0]), float(self.magnetization[1])],
        }


@dataclass
class TimelineFrame:
    time: Optional[float] = None
    wire_currents: Optional[List[float]] = None
    wires: List[Dict[str, Union[int, float]]] = field(default_factory=list)
    rotor_angle_deg: Optional[float] = None
    rotor_angles: List[Dict[str, Union[int, float, str]]] = field(default_factory=list)
    magnets: List[Dict[str, Union[int, float, Tuple[float, float]]]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        data: Dict[str, object] = {}
        if self.time is not None:
            data["t"] = float(self.time)
        if self.wire_currents is not None:
            data["wire_currents"] = [float(value) for value in self.wire_currents]
        if self.wires:
            serialised = []
            for override in self.wires:
                if "index" not in override:
                    raise ValueError("Timeline wire override requires an 'index' key")
                entry: Dict[str, object] = {"index": int(override["index"])}
                if "current" in override:
                    entry["current"] = float(override["current"])
                elif "I" in override:
                    entry["I"] = float(override["I"])
                elif "scale" in override:
                    entry["scale"] = float(override["scale"])
                else:
                    raise ValueError("Timeline wire override requires current/I/scale value")
                serialised.append(entry)
            data["wires"] = serialised
        if self.rotor_angle_deg is not None:
            data["rotor_angle"] = float(self.rotor_angle_deg)
        if self.rotor_angles:
            serialised_angles: List[Dict[str, object]] = []
            for override in self.rotor_angles:
                entry: Dict[str, object] = {}
                if "index" in override:
                    entry["index"] = int(override["index"])
                elif "name" in override:
                    entry["name"] = str(override["name"])
                else:
                    raise ValueError("Timeline rotor angle override requires an 'index' or 'name'")
                if "angle_deg" in override:
                    entry["angle_deg"] = float(override["angle_deg"])
                elif "angle" in override:
                    entry["angle"] = float(override["angle"])
                else:
                    raise ValueError("Timeline rotor angle override requires 'angle_deg' or 'angle'")
                serialised_angles.append(entry)
            data["rotor_angles"] = serialised_angles
        if self.magnets:
            serialised = []
            for override in self.magnets:
                if "index" not in override:
                    raise ValueError("Timeline magnet override requires an 'index'")
                entry: Dict[str, object] = {"index": int(override["index"])}
                if "angle_deg" in override:
                    entry["angle_deg"] = float(override["angle_deg"])
                if "magnetization" in override:
                    mx, my = override["magnetization"]  # type: ignore[index]
                    entry["magnetization"] = [float(mx), float(my)]
                serialised.append(entry)
            data["magnets"] = serialised
        return data


@dataclass
class PolygonRegion:
    vertices: List[Tuple[float, float]]
    material: str

    def to_dict(self) -> Dict[str, object]:
        serialized = [[float(x), float(y)] for x, y in self.vertices]
        return {"type": "polygon", "vertices": serialized, "material": self.material}


RegionSpec = Union[UniformRegion, HalfspaceRegion, PolygonRegion]


@dataclass
class Scenario:
    domain: Domain
    materials: List[Material] = field(default_factory=list)
    regions: List[RegionSpec] = field(default_factory=list)
    sources: List[Wire] = field(default_factory=list)
    rotors: List[Rotor] = field(default_factory=list)
    magnet_regions: List[MagnetRectRegion] = field(default_factory=list)
    outputs: List["ScenarioOutput"] = field(default_factory=list)
    timeline: List[TimelineFrame] = field(default_factory=list)
    units: str = "SI"
    version: str = "0.2"

    def _validate(self) -> None:
        if self.version not in {"0.1", "0.2"}:
            raise ValueError("Scenario version must be '0.1' or '0.2'")
        if self.units != "SI":
            raise ValueError("Only SI units are currently supported")
        if not self.materials:
            raise ValueError("At least one material is required")
        names = {mat.name for mat in self.materials}
        if len(names) != len(self.materials):
            raise ValueError("Material names must be unique")
        if not self.regions:
            raise ValueError("At least one region must be provided")
        allow_halfspace = self.version == "0.2"
        uniform_seen = False
        for region in self.regions:
            if not isinstance(region, (UniformRegion, HalfspaceRegion, PolygonRegion)):
                raise ValueError(f"Unsupported region type: {type(region)}")
            if region.material not in names:
                raise ValueError(f"Region references unknown material '{region.material}'")
            if isinstance(region, UniformRegion):
                uniform_seen = True
            elif isinstance(region, HalfspaceRegion):
                if not allow_halfspace:
                    raise ValueError("HalfspaceRegion requires scenario version '0.2'")
                if region.normal_x == 0.0 and region.normal_y == 0.0:
                    raise ValueError("HalfspaceRegion normal must be non-zero")
            elif isinstance(region, PolygonRegion):
                if not allow_halfspace:
                    raise ValueError("PolygonRegion requires scenario version '0.2'")
                if len(region.vertices) < 3:
                    raise ValueError("PolygonRegion requires at least three vertices")
                for vertex in region.vertices:
                    if len(vertex) != 2:
                        raise ValueError("PolygonRegion vertices must be (x, y) pairs")
        if self.version == "0.1" and not uniform_seen:
            raise ValueError("Scenario v0.1 requires at least one UniformRegion")
        for wire in self.sources:
            if wire.radius <= 0.0:
                raise ValueError("Wire radius must be positive")
        polygon_count = sum(isinstance(region, PolygonRegion) for region in self.regions)
        rotor_name_set: Set[str] = set()
        for rotor in self.rotors:
            if not isinstance(rotor, Rotor):
                raise ValueError(f"Unsupported rotor type: {type(rotor)}")
            if rotor.name:
                if rotor.name in rotor_name_set:
                    raise ValueError(f"Duplicate rotor name '{rotor.name}'")
                rotor_name_set.add(rotor.name)
            for index in rotor.polygon_indices:
                if index < 0 or index >= polygon_count:
                    raise ValueError("Rotor polygon index out of range")
            for index in rotor.magnet_indices:
                if index < 0 or index >= len(self.magnet_regions):
                    raise ValueError("Rotor magnet index out of range")
            for index in rotor.wire_indices:
                if index < 0 or index >= len(self.sources):
                    raise ValueError("Rotor wire index out of range")
        seen_ids: Set[str] = set()
        for output in self.outputs:
            if not isinstance(output, (FieldMapOutput, LineProbeOutput, StressTensorProbeOutput, BackEmfProbeOutput)):
                raise ValueError(f"Unsupported output type: {type(output)}")
            output.validate()
            if output.id in seen_ids:
                raise ValueError(f"Duplicate output id '{output.id}'")
            seen_ids.add(output.id)
        for region in self.magnet_regions:
            if not isinstance(region, MagnetRectRegion):
                raise ValueError(f"Unsupported magnet region type: {type(region)}")
            if region.x_range[0] >= region.x_range[1] or region.y_range[0] >= region.y_range[1]:
                raise ValueError("MagnetRectRegion requires strictly increasing ranges")
        for frame in self.timeline:
            if not isinstance(frame, TimelineFrame):
                raise ValueError(f"Unsupported timeline frame type: {type(frame)}")

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
        if self.magnet_regions:
            data["magnet_regions"] = [region.to_dict() for region in self.magnet_regions]
        if self.rotors:
            data["rotors"] = [rotor.to_dict() for rotor in self.rotors]
        if self.outputs:
            data["outputs"] = [output.to_dict() for output in self.outputs]
        if self.timeline:
            data["timeline"] = [frame.to_dict() for frame in self.timeline]
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
    "HalfspaceRegion",
    "PolygonRegion",
    "Wire",
    "Rotor",
    "FieldMapOutput",
    "LineProbeOutput",
    "StressTensorProbeOutput",
    "BackEmfProbeOutput",
    "MagnetRectRegion",
    "TimelineFrame",
    "Scenario",
]
