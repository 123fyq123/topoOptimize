#!/usr/bin/env python3
"""Generate a structured hexahedral chair-shaped Medit mesh."""

from __future__ import annotations

import argparse
from pathlib import Path


LatticePoint = tuple[int, int, int]
Cell = tuple[int, int, int]


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a chair design domain consisting of a solid seat/base "
            "and a straight vertical backrest."
        )
    )
    parser.add_argument("--width", type=positive_int, default=24)
    parser.add_argument("--depth", type=positive_int, default=24)
    parser.add_argument("--seat-height", type=positive_int, default=16)
    parser.add_argument("--back-height", type=positive_int, default=24)
    parser.add_argument("--back-thickness", type=positive_int, default=6)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists",
    )
    return parser.parse_args()


def validate_dimensions(
    width: int,
    depth: int,
    seat_height: int,
    back_height: int,
    back_thickness: int,
) -> None:
    if back_thickness > depth:
        raise ValueError(
            f"back thickness ({back_thickness}) cannot exceed depth ({depth})"
        )
    if min(width, depth, seat_height, back_height) <= 0:
        raise ValueError("all chair dimensions must be positive")


def build_cells(
    width: int,
    depth: int,
    seat_height: int,
    back_height: int,
    back_thickness: int,
) -> list[Cell]:
    """Return the union of the lower seat/base and the rear vertical backrest."""
    cells: list[Cell] = []

    # 下部完整长方体是座椅初始设计域，也为优化后的支腿结构提供材料空间。
    for k in range(seat_height):
        for j in range(depth):
            for i in range(width):
                cells.append((i, j, k))

    # 椅背只占后侧若干层，并沿 z 方向竖直延伸，避免初始域出现倾斜背板。
    rear_start = depth - back_thickness
    for k in range(seat_height, seat_height + back_height):
        for j in range(rear_start, depth):
            for i in range(width):
                cells.append((i, j, k))

    return cells


def build_vertex_ids(cells: list[Cell]) -> dict[LatticePoint, int]:
    points: set[LatticePoint] = set()
    for i, j, k in cells:
        for dk in (0, 1):
            for dj in (0, 1):
                for di in (0, 1):
                    points.add((i + di, j + dj, k + dk))

    return {
        point: vertex_id
        for vertex_id, point in enumerate(
            sorted(points, key=lambda point: (point[2], point[1], point[0])),
            start=1,
        )
    }


def default_output_path(
    width: int,
    depth: int,
    seat_height: int,
    back_height: int,
    back_thickness: int,
) -> Path:
    total_height = seat_height + back_height
    repository_root = Path(__file__).resolve().parents[1]
    directory = (
        repository_root
        / "input"
        / "topopt"
        / "data"
        / "mmc"
        / f"chair_{width}_{depth}_{total_height}"
    )
    return directory / (
        f"chair_{width}_{depth}_{total_height}"
        f"_seat{seat_height}_back{back_thickness}x{back_height}.mesh"
    )


def format_coordinate(value: float) -> str:
    if value.is_integer():
        return str(int(value))
    return format(value, ".17g")


def write_mesh(
    path: Path,
    width: int,
    depth: int,
    total_height: int,
    cells: list[Cell],
    vertex_ids: dict[LatticePoint, int],
) -> None:
    half_width = width / 2.0
    half_depth = depth / 2.0
    half_height = total_height / 2.0

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as output:
        output.write("MeshVersionFormatted 1\n")
        output.write("Dimension 3\n")
        output.write("Vertices\n")
        output.write(f"{len(vertex_ids)}\n")

        for i, j, k in sorted(vertex_ids, key=vertex_ids.get):
            x = i - half_width
            y = j - half_depth
            z = k - half_height
            output.write(
                f"{format_coordinate(x)} {format_coordinate(y)} "
                f"{format_coordinate(z)} -1\n"
            )

        output.write("Hexahedra\n")
        output.write(f"{len(cells)}\n")
        for i, j, k in cells:
            corners = (
                (i, j, k),
                (i + 1, j, k),
                (i + 1, j + 1, k),
                (i, j + 1, k),
                (i, j, k + 1),
                (i + 1, j, k + 1),
                (i + 1, j + 1, k + 1),
                (i, j + 1, k + 1),
            )
            output.write(" ".join(str(vertex_ids[corner]) for corner in corners))
            output.write(" 1\n")

        output.write("End\n")


def main() -> int:
    args = parse_args()
    validate_dimensions(
        args.width,
        args.depth,
        args.seat_height,
        args.back_height,
        args.back_thickness,
    )

    total_height = args.seat_height + args.back_height
    output_path = args.output or default_output_path(
        args.width,
        args.depth,
        args.seat_height,
        args.back_height,
        args.back_thickness,
    )
    output_path = output_path.resolve()
    if output_path.exists() and not args.force:
        raise FileExistsError(
            f"Output already exists: {output_path}\nUse --force to overwrite it."
        )

    cells = build_cells(
        args.width,
        args.depth,
        args.seat_height,
        args.back_height,
        args.back_thickness,
    )
    vertex_ids = build_vertex_ids(cells)
    write_mesh(
        output_path,
        args.width,
        args.depth,
        total_height,
        cells,
        vertex_ids,
    )

    expected_cells = (
        args.width * args.depth * args.seat_height
        + args.width * args.back_thickness * args.back_height
    )
    if len(cells) != expected_cells:
        raise RuntimeError(
            f"Internal cell-count mismatch: generated {len(cells)}, "
            f"expected {expected_cells}"
        )

    print(f"Generated: {output_path}")
    print(
        f"Dimensions: width={args.width}, depth={args.depth}, "
        f"total height={total_height}"
    )
    print(
        f"Seat/base height={args.seat_height}, backrest="
        f"{args.back_thickness} x {args.back_height}"
    )
    print(f"Vertices: {len(vertex_ids)}")
    print(f"Hexahedra: {len(cells)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
