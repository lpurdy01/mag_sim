#!/usr/bin/env python3
"""Lightweight sanity check for VTK ImageData exports."""

import argparse
import sys

try:
    import vtk  # type: ignore
    from vtk.util.numpy_support import vtk_to_numpy  # type: ignore
except ImportError as exc:  # pragma: no cover - optional dependency
    print("vtk module is required to run this script. Install it via 'pip install vtk'.")
    sys.exit(1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "path",
        nargs="?",
        default="outputs/field_frame_000.vti",
        help="VTK ImageData file to inspect (default: outputs/field_frame_000.vti)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(args.path)
    reader.Update()
    data = reader.GetOutput()

    cell_data = data.GetCellData()
    bx = cell_data.GetArray("Bx")
    by = cell_data.GetArray("By")

    if bx is None or by is None:
        print("VTK verification: FAILED (Bx/By arrays missing)")
        return 1

    bx_np = vtk_to_numpy(bx)
    by_np = vtk_to_numpy(by)

    if bx_np.shape != by_np.shape or bx_np.size == 0:
        print("VTK verification: FAILED (array shape mismatch)")
        return 1

    print("VTK verification: PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())

