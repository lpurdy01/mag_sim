#!/usr/bin/env python3
"""Lightweight sanity check for VTK ImageData exports."""

import argparse
import sys

try:
    import numpy as np
    import vtk  # type: ignore
    from vtk.util.numpy_support import vtk_to_numpy  # type: ignore
except ImportError as exc:  # pragma: no cover - optional dependency
    print("vtk and numpy modules are required. Install them via 'pip install vtk numpy'.")
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
    bvec = cell_data.GetArray("B")

    if bx is None or by is None or bvec is None:
        print("VTK verification: FAILED (B arrays missing)")
        return 1

    if bvec.GetNumberOfComponents() != 3:
        print("VTK verification: FAILED (B vector components != 3)")
        return 1

    bx_np = vtk_to_numpy(bx)
    by_np = vtk_to_numpy(by)
    bvec_np = vtk_to_numpy(bvec).reshape((-1, 3))

    if bx_np.shape != by_np.shape or bx_np.size == 0:
        print("VTK verification: FAILED (array shape mismatch)")
        return 1

    if bvec_np.shape[0] != bx_np.shape[0]:
        print("VTK verification: FAILED (B vector size mismatch)")
        return 1

    if not np.allclose(bvec_np[:, 0], bx_np):
        print("VTK verification: FAILED (B vector X component mismatch)")
        return 1
    if not np.allclose(bvec_np[:, 1], by_np):
        print("VTK verification: FAILED (B vector Y component mismatch)")
        return 1
    if not np.allclose(bvec_np[:, 2], 0.0):
        print("VTK verification: FAILED (B vector Z component expected to be zero)")
        return 1

    print("VTK verification: PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())

