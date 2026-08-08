#!/usr/bin/env python3
"""Extrude a one-cell-thick unit-voxel Medit .mesh topology."""

from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path
from typing import Iterable


AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
AXIS_NAMES = ("x", "y", "z")
CELL_CORNERS = tuple(itertools.product((0, 1), repeat=3))
CornerPattern = tuple[tuple[int, int, int], ...]
Point = tuple[float, float, float]
CellMap = dict[Point, int]


def _clean_lines(path: Path) -> list[str]:
    lines: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            lines.append(line)
    return lines


def _section(lines: list[str], name: str) -> tuple[int, list[str]]:
    try:
        start = lines.index(name)
    except ValueError as exc:
        raise ValueError(f"Missing {name} section") from exc
    count = int(lines[start + 1])
    records = lines[start + 2 : start + 2 + count]
    if len(records) != count:
        raise ValueError(f"Incomplete {name} section")
    return count, records


def read_mesh(path: Path) -> tuple[list[tuple[float, float, float, int]], list[tuple[list[int], int]]]:
    lines = _clean_lines(path)
    vertex_count, vertex_lines = _section(lines, "Vertices")
    hex_count, hex_lines = _section(lines, "Hexahedra")

    vertices: list[tuple[float, float, float, int]] = []
    for line in vertex_lines:
        values = line.split()
        if len(values) < 4:
            raise ValueError(f"Invalid vertex record: {line}")
        vertices.append((float(values[0]), float(values[1]), float(values[2]), int(values[3])))

    hexahedra: list[tuple[list[int], int]] = []
    for line in hex_lines:
        values = line.split()
        if len(values) < 9:
            raise ValueError(f"Invalid hexahedron record: {line}")
        node_ids = [int(value) for value in values[:8]]
        if any(node_id < 1 or node_id > vertex_count for node_id in node_ids):
            raise ValueError(f"Hexahedron contains an invalid vertex id: {line}")
        hexahedra.append((node_ids, int(values[8])))

    if len(vertices) != vertex_count or len(hexahedra) != hex_count:
        raise ValueError("Mesh section count mismatch")
    return vertices, hexahedra


def _key(values: Iterable[float]) -> tuple[float, ...]:
    return tuple(round(value, 12) for value in values)


def validate_unit_voxel_mesh(
    vertices: list[tuple[float, float, float, int]],
    hexahedra: list[tuple[list[int], int]],
) -> tuple[dict[Point, int], CellMap, CornerPattern]:
    if not vertices or not hexahedra:
        raise ValueError("Mesh must contain vertices and hexahedra")

    vertex_refs: dict[Point, int] = {}
    for x, y, z, reference in vertices:
        point_key = _key((x, y, z))
        if point_key in vertex_refs:
            raise ValueError(f"Duplicate vertex at {point_key}")
        vertex_refs[point_key] = reference

    cells: CellMap = {}
    corner_pattern: CornerPattern | None = None
    expected_corners = set(CELL_CORNERS)
    for node_ids, reference in hexahedra:
        points = [vertices[node_id - 1][:3] for node_id in node_ids]
        axes = [sorted({point[axis] for point in points}) for axis in range(3)]
        if any(len(values) != 2 for values in axes):
            raise ValueError("Every cell must be an axis-aligned hexahedron")
        for axis, values in enumerate(axes):
            if not math.isclose(values[1] - values[0], 1.0, rel_tol=0.0, abs_tol=1e-9):
                raise ValueError(
                    f"Every cell edge must have length 1; found {AXIS_NAMES[axis]} span "
                    f"{values[1] - values[0]:.17g}"
                )

        current_pattern: CornerPattern = tuple(
            tuple(0 if math.isclose(point[axis], axes[axis][0], abs_tol=1e-9) else 1
                  for axis in range(3))
            for point in points
        )
        if set(current_pattern) != expected_corners:
            raise ValueError("A hexahedron does not contain all eight unit-cell corners")
        if corner_pattern is None:
            corner_pattern = current_pattern
        elif current_pattern != corner_pattern:
            raise ValueError("Hexahedra do not use one consistent local vertex ordering")

        origin = _key((axes[0][0], axes[1][0], axes[2][0]))
        if origin in cells:
            raise ValueError(f"Duplicate hexahedron at {origin}")
        cells[origin] = reference

    if corner_pattern is None:
        raise ValueError("Mesh contains no hexahedra")
    return vertex_refs, cells, corner_pattern


def thin_axes(cells: CellMap) -> list[int]:
    return [axis for axis in range(3) if len({origin[axis] for origin in cells}) == 1]


def extrude_voxel_topology(
    vertex_refs: dict[Point, int],
    cells: CellMap,
    corner_pattern: CornerPattern,
    axis: int,
    layers: int,
    side: str,
) -> tuple[list[tuple[float, float, float, int]], list[tuple[list[int], int]]]:
    if layers < 1:
        raise ValueError("Layer count must be at least 1")
    if side not in ("positive", "negative"):
        raise ValueError("Side must be positive or negative")

    axis_origins = sorted({origin[axis] for origin in cells})
    if len(axis_origins) != 1:
        raise ValueError(
            f"Selected {AXIS_NAMES[axis]} axis has {len(axis_origins)} cell layers; "
            "the extrusion axis must have exactly one layer"
        )
    base = axis_origins[0]
    direction = 1.0 if side == "positive" else -1.0

    output_cells: CellMap = {}
    for origin, reference in cells.items():
        for layer in range(layers + 1):
            shifted = list(origin)
            shifted[axis] += direction * layer
            shifted_key = _key(shifted)
            if shifted_key in output_cells:
                raise ValueError(f"Extrusion produced a duplicate cell at {shifted_key}")
            output_cells[shifted_key] = reference

    sorted_origins = sorted(output_cells, key=lambda point: (point[2], point[1], point[0]))
    vertex_ids: dict[Point, int] = {}
    output_vertices: list[tuple[float, float, float, int]] = []
    for origin in sorted_origins:
        for corner in CELL_CORNERS:
            point = _key(origin[current_axis] + corner[current_axis] for current_axis in range(3))
            if point in vertex_ids:
                continue
            source_point = list(point)
            source_point[axis] = min(base + 1.0, max(base, source_point[axis]))
            source_key = _key(source_point)
            if source_key not in vertex_refs:
                raise ValueError(f"Cannot find source vertex reference for {source_key}")
            output_vertices.append((point[0], point[1], point[2], vertex_refs[source_key]))
            vertex_ids[point] = len(output_vertices)

    output_hexahedra: list[tuple[list[int], int]] = []
    for origin in sorted_origins:
        node_ids = []
        for corner in corner_pattern:
            point = _key(origin[current_axis] + corner[current_axis] for current_axis in range(3))
            node_ids.append(vertex_ids[point])
        output_hexahedra.append((node_ids, output_cells[origin]))

    return output_vertices, output_hexahedra


def _format_number(value: float) -> str:
    if abs(value) < 1e-14:
        value = 0.0
    return format(value, ".17g")


def write_mesh(
    path: Path,
    vertices: list[tuple[float, float, float, int]],
    hexahedra: list[tuple[list[int], int]],
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as output:
        output.write("MeshVersionFormatted 1\n")
        output.write("Dimension 3\n")
        output.write("Vertices\n")
        output.write(f"{len(vertices)}\n")
        for x, y, z, reference in vertices:
            output.write(
                f"{_format_number(x)} {_format_number(y)} {_format_number(z)} {reference}\n"
            )
        output.write("Hexahedra\n")
        output.write(f"{len(hexahedra)}\n")
        for node_ids, reference in hexahedra:
            output.write(" ".join(str(node_id) for node_id in node_ids))
            output.write(f" {reference}\n")
        output.write("End\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extrude an arbitrary one-cell-thick unit-voxel Medit topology."
    )
    parser.add_argument("input", type=Path, help="Input .mesh file")
    parser.add_argument("-o", "--output", type=Path, help="Output .mesh file")
    parser.add_argument("--axis", choices=AXIS_INDEX, required=True, help="One-cell-thick axis")
    parser.add_argument("--layers", type=int, required=True, help="Additional unit-cell layers")
    parser.add_argument(
        "--side",
        choices=("negative", "positive"),
        default="positive",
        help="Extrusion direction (default: positive)",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.layers < 1:
        raise SystemExit("--layers must be at least 1")
    if not args.input.is_file():
        raise SystemExit(f"Input file does not exist: {args.input}")

    output = args.output or args.input.with_name(
        f"{args.input.stem}_thick_{args.axis}_{args.side}_{args.layers}.mesh"
    )
    if output.exists() and not args.force:
        raise SystemExit(f"Output already exists: {output} (use --force to overwrite)")

    try:
        vertices, hexahedra = read_mesh(args.input)
        vertex_refs, cells, corner_pattern = validate_unit_voxel_mesh(vertices, hexahedra)
        new_vertices, new_hexahedra = extrude_voxel_topology(
            vertex_refs,
            cells,
            corner_pattern,
            AXIS_INDEX[args.axis],
            args.layers,
            args.side,
        )
    except (ValueError, KeyError) as exc:
        raise SystemExit(f"Cannot extrude mesh: {exc}") from exc

    output.parent.mkdir(parents=True, exist_ok=True)
    write_mesh(output, new_vertices, new_hexahedra)
    validate_unit_voxel_mesh(*read_mesh(output))

    print(f"Input : {args.input}")
    print(f"Output: {output}")
    print(f"Axis  : {args.axis}, thickness 1 -> {args.layers + 1}, side {args.side}")
    print(f"Cells : {len(hexahedra)} -> {len(new_hexahedra)}")
    print(f"Points: {len(vertices)} -> {len(new_vertices)}")


if __name__ == "__main__":
    main()
