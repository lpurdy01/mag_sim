from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import ezdxf
import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402


class DxfError(RuntimeError):
    """Raised when DXF artefacts cannot be processed."""


@dataclass
class DxfLayerInfo:
    name: str
    entity_count: int
    types: Sequence[str]


@dataclass
class DxfPrimitive:
    kind: str
    points: Sequence[Tuple[float, float]] | None = None
    center: Tuple[float, float] | None = None
    radius: float | None = None
    layer: str = "0"
    category: Optional[str] = None


_CATEGORY_COLOURS = {
    "domain": "#0d6efd",
    "material": "#198754",
    "magnet": "#d63384",
    "wire": "#fd7e14",
    "structural": "#20c997",
    "unassigned": "#6c757d",
}


def summarise_dxf(path: Path) -> List[DxfLayerInfo]:
    """Return basic layer statistics for the supplied DXF file."""

    try:
        doc = ezdxf.readfile(str(path))
    except (OSError, ezdxf.DXFStructureError, ezdxf.DXFTableEntryError) as exc:
        raise DxfError(f"Failed to read DXF '{path.name}': {exc}") from exc

    counts: Dict[str, int] = {}
    types: Dict[str, set[str]] = {}

    for entity in doc.modelspace():
        layer = entity.dxf.layer or "0"
        counts[layer] = counts.get(layer, 0) + 1
        types.setdefault(layer, set()).add(entity.dxftype())

    summary: List[DxfLayerInfo] = []
    for layer, count in sorted(counts.items()):
        summary.append(DxfLayerInfo(name=layer, entity_count=count, types=sorted(types[layer])))

    return summary


def _polyline_points(entity) -> List[Tuple[float, float]]:
    if entity.dxftype() == "LWPOLYLINE":
        return [(pt[0], pt[1]) for pt in entity.get_points("xy")]
    if entity.dxftype() == "POLYLINE":
        return [(v.dxf.location.x, v.dxf.location.y) for v in entity.vertices()]
    if entity.dxftype() == "LINE":
        return [
            (float(entity.dxf.start.x), float(entity.dxf.start.y)),
            (float(entity.dxf.end.x), float(entity.dxf.end.y)),
        ]
    return []


def load_primitives(path: Path, layers: Iterable[str], categories: Dict[str, Optional[str]]) -> List[DxfPrimitive]:
    """Extract drawable primitives for the requested layers."""

    selected = {layer.lower() for layer in layers}
    if not selected:
        return []

    try:
        doc = ezdxf.readfile(str(path))
    except (OSError, ezdxf.DXFStructureError, ezdxf.DXFTableEntryError) as exc:
        raise DxfError(f"Failed to read DXF '{path.name}': {exc}") from exc

    primitives: List[DxfPrimitive] = []
    for entity in doc.modelspace():
        layer = (entity.dxf.layer or "0")
        if layer.lower() not in selected:
            continue

        category = categories.get(layer)
        if entity.dxftype() in {"LWPOLYLINE", "POLYLINE", "LINE"}:
            points = _polyline_points(entity)
            if not points:
                continue
            primitives.append(DxfPrimitive(kind="polyline", points=points, layer=layer, category=category))
        elif entity.dxftype() == "CIRCLE":
            center = (float(entity.dxf.center.x), float(entity.dxf.center.y))
            primitives.append(
                DxfPrimitive(
                    kind="circle",
                    center=center,
                    radius=float(entity.dxf.radius),
                    layer=layer,
                    category=category,
                )
            )

    return primitives


def render_preview(primitives: Sequence[DxfPrimitive], output_path: Path) -> Path:
    """Render DXF primitives into an image."""

    if not primitives:
        raise DxfError("No DXF layers selected for preview.")

    xmin = math.inf
    xmax = -math.inf
    ymin = math.inf
    ymax = -math.inf

    for primitive in primitives:
        if primitive.kind == "polyline" and primitive.points:
            xs, ys = zip(*primitive.points)
            xmin = min(xmin, min(xs))
            xmax = max(xmax, max(xs))
            ymin = min(ymin, min(ys))
            ymax = max(ymax, max(ys))
        elif primitive.kind == "circle" and primitive.center and primitive.radius is not None:
            cx, cy = primitive.center
            r = primitive.radius
            xmin = min(xmin, cx - r)
            xmax = max(xmax, cx + r)
            ymin = min(ymin, cy - r)
            ymax = max(ymax, cy + r)

    if not math.isfinite(xmin) or not math.isfinite(xmax) or xmin == xmax or ymin == ymax:
        xmin, xmax, ymin, ymax = -0.5, 0.5, -0.5, 0.5

    padding_x = max(0.05 * (xmax - xmin), 1e-3)
    padding_y = max(0.05 * (ymax - ymin), 1e-3)
    xmin -= padding_x
    xmax += padding_x
    ymin -= padding_y
    ymax += padding_y

    fig, ax = plt.subplots(figsize=(6, 6))
    colour_cache: Dict[str, str] = {}

    for primitive in primitives:
        category = primitive.category or "unassigned"
        colour = colour_cache.get(category)
        if colour is None:
            colour = _CATEGORY_COLOURS.get(category, _CATEGORY_COLOURS["unassigned"])
            colour_cache[category] = colour

        if primitive.kind == "polyline" and primitive.points:
            xs, ys = zip(*primitive.points)
            ax.plot(xs, ys, color=colour, linewidth=1.2)
        elif primitive.kind == "circle" and primitive.center and primitive.radius is not None:
            circle = Circle(primitive.center, primitive.radius, fill=False, linewidth=1.2, color=colour)
            ax.add_patch(circle)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return output_path


def export_scenario_to_dxf(scenario_path: Path, output_dir: Path) -> List[Path]:
    """Convert a scenario JSON into a set of DXF artefacts."""

    from python.visualize_scenario_field import compute_domain_bounds, gather_boundaries, load_scenario

    spec, wires = load_scenario(scenario_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    exported: List[Path] = []

    # Domain boundary
    try:
        xmin, xmax, ymin, ymax = compute_domain_bounds(spec)
    except ValueError as exc:
        raise DxfError(str(exc)) from exc

    domain_doc = ezdxf.new("R2010")
    msp = domain_doc.modelspace()
    msp.add_lwpolyline(
        [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)],
        dxfattribs={"layer": "domain"},
    )
    domain_path = output_dir / "domain.dxf"
    domain_doc.saveas(domain_path)
    exported.append(domain_path)

    # Material regions
    boundaries = gather_boundaries(spec)
    region_doc = ezdxf.new("R2010")
    region_msp = region_doc.modelspace()
    for index, (xs, ys, material) in enumerate(boundaries.get("polygons", [])):
        if len(xs) < 2:
            continue
        coords = list(zip(xs, ys))
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        layer_name = f"material_{material or index}"
        region_msp.add_lwpolyline(coords, dxfattribs={"layer": layer_name})
    if len(region_msp):
        regions_path = output_dir / "materials.dxf"
        region_doc.saveas(regions_path)
        exported.append(regions_path)

    magnet_doc = ezdxf.new("R2010")
    magnet_msp = magnet_doc.modelspace()
    for idx, (xs, ys) in enumerate(boundaries.get("magnet_polygons", [])):
        if len(xs) < 2:
            continue
        coords = list(zip(xs, ys))
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        magnet_msp.add_lwpolyline(coords, dxfattribs={"layer": f"magnet_{idx}"})
    for idx, rect in enumerate(boundaries.get("magnet_rects", [])):
        xmin_r, xmax_r, ymin_r, ymax_r = rect
        coords = [
            (xmin_r, ymin_r),
            (xmax_r, ymin_r),
            (xmax_r, ymax_r),
            (xmin_r, ymax_r),
            (xmin_r, ymin_r),
        ]
        magnet_msp.add_lwpolyline(coords, dxfattribs={"layer": f"magnet_rect_{idx}"})
    if len(magnet_msp):
        magnets_path = output_dir / "magnets.dxf"
        magnet_doc.saveas(magnets_path)
        exported.append(magnets_path)

    if wires:
        wire_doc = ezdxf.new("R2010")
        wire_msp = wire_doc.modelspace()
        for idx, wire in enumerate(wires):
            wire_msp.add_circle(
                center=(wire.x, wire.y),
                radius=wire.radius,
                dxfattribs={"layer": f"wire_{idx}"},
            )
        wires_path = output_dir / "wires.dxf"
        wire_doc.saveas(wires_path)
        exported.append(wires_path)

    return exported
