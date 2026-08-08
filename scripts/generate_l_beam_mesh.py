#!/usr/bin/env python3
"""Generate a one-cell-thick, unit-hexahedron L-beam Medit mesh."""

from __future__ import annotations

import argparse
from pathlib import Path


LatticePoint = tuple[int, int, int]
Cell = tuple[int, int]


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate an L-shaped x-z design domain with one unit cell in the y direction."
        )
    )
    parser.add_argument("--length", "--L", type=positive_int, default=64, help="Outer x length")
    parser.add_argument("--height", "--H", type=positive_int, default=64, help="Outer z height")
    parser.add_argument(
        "--leg-width",
        "--Ld",
        type=positive_int,
        default=24,
        help="Width Ld of the left vertical leg",
    )
    parser.add_argument(
        "--leg-height",
        "--Hd",
        type=positive_int,
        default=24,
        help="Height Hd of the bottom horizontal leg",
    )
    parser.add_argument("-o", "--output", type=Path, help="Output Medit .mesh path")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists",
    )
    return parser.parse_args()


def validate_dimensions(length: int, height: int, leg_width: int, leg_height: int) -> None:
    if leg_width > length:
        raise ValueError(f"Ld ({leg_width}) cannot exceed L ({length})")
    if leg_height > height:
        raise ValueError(f"Hd ({leg_height}) cannot exceed H ({height})")


def build_cells(length: int, height: int, leg_width: int, leg_height: int) -> list[Cell]:
    """Return cells as (i, k), where k=0 is the top cell row."""
    bottom_start = height - leg_height
    return [
        (i, k)
        for k in range(height)
        for i in range(length)
        if i < leg_width or k >= bottom_start
    ]


def build_vertex_ids(cells: list[Cell]) -> dict[LatticePoint, int]:
    points: set[LatticePoint] = set()
    for i, k in cells:
        for dk in (0, 1):
            for j in (0, 1):
                for di in (0, 1):
                    points.add((i + di, j, k + dk))

    # Match the repository's regular mesh ordering: z rows, then y, then x.
    return {
        point: vertex_id
        for vertex_id, point in enumerate(
            sorted(points, key=lambda point: (point[2], point[1], point[0])),
            start=1,
        )
    }


def default_output_path(
    length: int,
    height: int,
    leg_width: int,
    leg_height: int,
) -> Path:
    repository_root = Path(__file__).resolve().parents[1]
    directory = (
        repository_root
        / "input"
        / "topopt"
        / "data"
        / "mmc"
        / f"l_beam_{length}_1_{height}"
    )
    return directory / f"l_beam_{length}_1_{height}_Ld{leg_width}_Hd{leg_height}.mesh"


def format_coordinate(value: float) -> str:
    if value.is_integer():
        return str(int(value))
    return format(value, ".17g")


def write_mesh(
    path: Path,
    length: int,
    height: int,
    cells: list[Cell],
    vertex_ids: dict[LatticePoint, int],
) -> None:
    half_length = length / 2.0
    half_height = height / 2.0

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as output:
        output.write("MeshVersionFormatted 1\n")
        output.write("Dimension 3\n")
        output.write("Vertices\n")
        output.write(f"{len(vertex_ids)}\n")

        for i, j, k in sorted(vertex_ids, key=vertex_ids.get):
            x = i - half_length
            y = j - 0.5
            z = half_height - k
            output.write(
                f"{format_coordinate(x)} {format_coordinate(y)} "
                f"{format_coordinate(z)} -1\n"
            )

        output.write("Hexahedra\n")
        output.write(f"{len(cells)}\n")
        for i, k in cells:
            corners = (
                (i, 0, k),
                (i + 1, 0, k),
                (i + 1, 1, k),
                (i, 1, k),
                (i, 0, k + 1),
                (i + 1, 0, k + 1),
                (i + 1, 1, k + 1),
                (i, 1, k + 1),
            )
            output.write(" ".join(str(vertex_ids[corner]) for corner in corners))
            output.write(" 1\n")

        output.write("End\n")


def main() -> int:
    args = parse_args()
    validate_dimensions(args.length, args.height, args.leg_width, args.leg_height)

    output_path = args.output or default_output_path(
        args.length,
        args.height,
        args.leg_width,
        args.leg_height,
    )
    output_path = output_path.resolve()
    if output_path.exists() and not args.force:
        raise FileExistsError(
            f"Output already exists: {output_path}\nUse --force to overwrite it."
        )

    cells = build_cells(args.length, args.height, args.leg_width, args.leg_height)
    vertex_ids = build_vertex_ids(cells)
    write_mesh(output_path, args.length, args.height, cells, vertex_ids)

    expected_cells = (
        args.leg_width * args.height
        + (args.length - args.leg_width) * args.leg_height
    )
    if len(cells) != expected_cells:
        raise RuntimeError(
            f"Internal cell-count mismatch: generated {len(cells)}, expected {expected_cells}"
        )

    print(f"Generated: {output_path}")
    print(f"Dimensions: L={args.length}, y=1, H={args.height}")
    print(f"Legs: Ld={args.leg_width}, Hd={args.leg_height}")
    print(f"Vertices: {len(vertex_ids)}")
    print(f"Hexahedra: {len(cells)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
