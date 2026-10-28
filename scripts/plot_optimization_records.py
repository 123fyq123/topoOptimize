#!/usr/bin/env python3
"""Generate compact comparison figures from MMC-SIMP paper TSV records."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SYNC_COLOR = "#0000FF"
SEQUENTIAL_COLOR = "#FFFF00"
TARGET_COLOR = "#555555"
MAX_COMPARISON_FE_BUDGET = 100.0


def configure_chinese_font() -> None:
	# 优先使用 Windows 常见中文字体，避免图题和坐标轴中的中文显示为方框。
	matplotlib.rcParams["font.sans-serif"] = [
		"Microsoft YaHei",
		"SimHei",
		"Noto Sans CJK SC",
		"DejaVu Sans",
	]
	matplotlib.rcParams["axes.unicode_minus"] = False


@dataclass
class RunRecord:
	config_path: Path
	iteration_path: Path
	config: dict[str, str]
	rows: list[dict[str, str]]

	@property
	def method(self) -> str:
		return self.config["method"]

	@property
	def run_id(self) -> str:
		return self.config["run_id"]

	@property
	def final_row(self) -> dict[str, str]:
		coupled_rows = self.coupled_rows
		if not coupled_rows:
			raise ValueError(f"Run has no coupled-compliance rows: {self.run_id}")
		return coupled_rows[-1]

	@property
	def coupled_rows(self) -> list[dict[str, str]]:
		# 顺序版只有材料阶段的柔度才包含 E(eta)；同步版排除纯 MMC 暖启动，
		# 保证两条曲线比较的都是同一个耦合柔度定义。
		coupled_stage = "simp_material" if self.method == "sequential" else "simultaneous"
		return [row for row in self.rows if row["stage"] == coupled_stage]


def finite_number(value: str | None) -> float | None:
	try:
		number = float(value)  # type: ignore[arg-type]
	except (TypeError, ValueError):
		return None
	return number if math.isfinite(number) else None


def read_config(path: Path) -> dict[str, str]:
	with path.open("r", encoding="utf-8-sig", newline="") as stream:
		return {
			row["parameter"]: row["value"]
			for row in csv.DictReader(stream, delimiter="\t")
		}


def read_iterations(path: Path) -> list[dict[str, str]]:
	with path.open("r", encoding="utf-8-sig", newline="") as stream:
		return list(csv.DictReader(stream, delimiter="\t"))


def discover_runs(search_root: Path) -> list[RunRecord]:
	runs: list[RunRecord] = []
	for config_path in search_root.rglob("paper_run_config.tsv"):
		# 新记录均位于 records/<run_id>；忽略旧的手工归档和 latest 文件。
		if "records" not in config_path.parts:
			continue
		iteration_path = config_path.with_name("paper_iterations.tsv")
		if not iteration_path.exists():
			continue
		config = read_config(config_path)
		if "run_id" not in config or "method" not in config:
			continue
		runs.append(
			RunRecord(
				config_path=config_path,
				iteration_path=iteration_path,
				config=config,
				rows=read_iterations(iteration_path),
			)
		)
	return runs


def parameter_key(run: RunRecord) -> tuple[str, ...]:
	return (
		run.config.get("model_signature", ""),
		run.config.get("boundary_signature", ""),
		run.config.get("material_E1", ""),
		run.config.get("material_E2", ""),
		run.config.get("material_ratio_f", ""),
		run.config.get("mmc_target_rate", ""),
	)


def select_pair(
	runs: Iterable[RunRecord],
	sequential_run_id: str | None,
	simultaneous_run_id: str | None,
) -> tuple[RunRecord, RunRecord]:
	run_list = list(runs)
	by_id = {run.run_id: run for run in run_list}
	if sequential_run_id:
		sequential = by_id.get(sequential_run_id)
		if sequential is None or sequential.method != "sequential":
			raise ValueError(f"Unknown sequential run id: {sequential_run_id}")
	else:
		sequential_candidates = [run for run in run_list if run.method == "sequential"]
		if not sequential_candidates:
			raise ValueError("No sequential record was found.")
		sequential = max(
			sequential_candidates,
			key=lambda run: run.config_path.stat().st_mtime,
		)

	if simultaneous_run_id:
		simultaneous = by_id.get(simultaneous_run_id)
		if simultaneous is None or simultaneous.method != "simultaneous":
			raise ValueError(f"Unknown simultaneous run id: {simultaneous_run_id}")
	else:
		matching = [
			run
			for run in run_list
			if run.method == "simultaneous" and parameter_key(run) == parameter_key(sequential)
		]
		if not matching:
			raise ValueError(
				"No simultaneous run matches the latest sequential model, boundary and material parameters."
			)
		simultaneous = max(matching, key=lambda run: run.config_path.stat().st_mtime)

	if parameter_key(sequential) != parameter_key(simultaneous):
		raise ValueError("The selected runs do not use the same model, boundary and material parameters.")
	return sequential, simultaneous


def series(
	run: RunRecord,
	field: str,
	stage: str | None = None,
	max_fe_budget: float | None = None,
) -> tuple[list[float], list[float]]:
	x_values: list[float] = []
	y_values: list[float] = []
	for row in run.rows:
		if stage is not None and row["stage"] != stage:
			continue
		x = finite_number(row.get("cumulative_fe_solves"))
		y = finite_number(row.get(field))
		if (
			x is not None
			and y is not None
			and (max_fe_budget is None or x <= max_fe_budget)
		):
			x_values.append(x)
			y_values.append(y)
	return x_values, y_values


def row_at_or_before_budget(run: RunRecord, budget: float) -> dict[str, str]:
	candidates: list[dict[str, str]] = []
	for row in run.coupled_rows:
		fe_solves = finite_number(row.get("cumulative_fe_solves"))
		compliance = finite_number(row.get("compliance"))
		if fe_solves is not None and fe_solves <= budget and compliance is not None:
			candidates.append(row)
	if not candidates:
		raise ValueError(
			f"Run has no coupled state at or before FE budget {budget:g}: {run.run_id}"
		)
	return candidates[-1]


def common_fe_budget(sequential: RunRecord, simultaneous: RunRecord) -> float:
	# 论文图统一比较前 100 次有限元求解；即使顺序版继续跑到 160 次，
	# 也不使用其额外计算预算下的终值参与柔度对比。
	return min(
		MAX_COMPARISON_FE_BUDGET,
		final_value(sequential, "cumulative_fe_solves"),
		final_value(simultaneous, "cumulative_fe_solves"),
	)


def style_axis(axis: plt.Axes) -> None:
	axis.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
	axis.set_axisbelow(True)
	axis.spines["top"].set_visible(False)
	axis.spines["right"].set_visible(False)


def save_figure(figure: plt.Figure, output_path: Path) -> None:
	figure.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
	plt.close(figure)


def experiment_caption(run: RunRecord) -> str:
	load_x = float(run.config["neumann_centroid_x"])
	load_z = float(run.config["neumann_centroid_z"])
	return (
		f"模型：{run.config['model_id']}；边界条件：{run.config['boundary_signature'][:8]}；"
		f"载荷中心：({load_x:.3f}, {load_z:.3f})；"
		f"材料弹性模量：{run.config['material_E1']}/{run.config['material_E2']}；"
		f"目标材料体积比：{run.config['material_ratio_f']}"
	)


def threshold_note(sequential: RunRecord, simultaneous: RunRecord) -> str:
	sequential_threshold = sequential.config.get("material_active_threshold", "unknown")
	simultaneous_threshold = simultaneous.config.get("material_active_threshold", "unknown")
	if sequential_threshold == simultaneous_threshold:
		return f"材料活跃阈值：{sequential_threshold}"
	return (
		"注意：材料活跃阈值不同，"
		f"顺序优化={sequential_threshold}，同步优化={simultaneous_threshold}"
	)


def plot_compliance(
	sequential: RunRecord,
	simultaneous: RunRecord,
	output_dir: Path,
) -> None:
	figure, axis = plt.subplots(figsize=(7.2, 4.6))
	budget = common_fe_budget(sequential, simultaneous)
	sx, sy = series(
		sequential,
		"compliance",
		stage="simp_material",
		max_fe_budget=budget,
	)
	mx, my = series(
		simultaneous,
		"compliance",
		stage="simultaneous",
		max_fe_budget=budget,
	)
	axis.plot(mx, my, color=SYNC_COLOR, linewidth=2.2, label="同步优化")
	axis.plot(
		sx,
		sy,
		color=SEQUENTIAL_COLOR,
		linewidth=2.2,
		marker="o",
		markersize=2.8,
		markevery=max(1, len(sx) // 8),
		label="顺序优化（材料阶段）",
	)
	seq_budget_row = row_at_or_before_budget(sequential, budget)
	sim_budget_row = row_at_or_before_budget(simultaneous, budget)
	seq_budget_x = final_number(seq_budget_row, "cumulative_fe_solves")
	sim_budget_x = final_number(sim_budget_row, "cumulative_fe_solves")
	seq_budget_y = final_number(seq_budget_row, "compliance")
	sim_budget_y = final_number(sim_budget_row, "compliance")
	axis.axvline(
		budget,
		color=TARGET_COLOR,
		linestyle=":",
		linewidth=1.2,
		label=f"统一比较位置：第{budget:g}次有限元求解",
	)
	axis.scatter([seq_budget_x], [seq_budget_y], color=SEQUENTIAL_COLOR, s=38, zorder=4)
	axis.scatter([sim_budget_x], [sim_budget_y], color=SYNC_COLOR, s=38, zorder=4)
	axis.set_xlim(left=0.0, right=budget * 1.02)
	axis.set_xlabel("累计有限元求解次数")
	axis.set_ylabel("耦合柔度")
	axis.legend(frameon=False, fontsize=9)
	style_axis(axis)
	figure.suptitle("同步与顺序优化的柔度变化", fontsize=12, weight="bold", y=0.98)
	figure.text(0.5, 0.925, experiment_caption(sequential), ha="center", fontsize=8)
	figure.text(
		0.5,
		0.02,
		(
			f"第{budget:g}次有限元求解：顺序优化={seq_budget_y:.6g}，"
			f"同步优化={sim_budget_y:.6g}；"
			+ threshold_note(sequential, simultaneous)
		),
		ha="center",
		fontsize=8,
		color="#8B0000",
	)
	figure.tight_layout(rect=(0.0, 0.06, 1.0, 0.89))
	save_figure(figure, output_dir / "01_compliance_vs_cumulative_fe_solves.png")


def final_value(run: RunRecord, field: str) -> float:
	value = finite_number(run.final_row.get(field))
	if value is None:
		raise ValueError(f"Missing final {field} for {run.run_id}")
	return value


def final_number(row: dict[str, str], field: str) -> float:
	value = finite_number(row.get(field))
	if value is None:
		raise ValueError(f"Missing finite {field} in a comparison row.")
	return value


def plot_equal_budget_compliance(
	sequential: RunRecord,
	simultaneous: RunRecord,
	output_dir: Path,
) -> None:
	budget = common_fe_budget(sequential, simultaneous)
	seq_row = row_at_or_before_budget(sequential, budget)
	sim_row = row_at_or_before_budget(simultaneous, budget)
	seq_value = final_number(seq_row, "compliance")
	sim_value = final_number(sim_row, "compliance")
	values = [seq_value, sim_value]
	reduction = (seq_value - sim_value) / seq_value * 100.0

	figure, axis = plt.subplots(figsize=(7.2, 4.8))
	bars = axis.bar(
		["顺序优化", "同步优化"],
		values,
		color=[SEQUENTIAL_COLOR, SYNC_COLOR],
		width=0.58,
	)
	axis.set_ylabel("耦合柔度")
	axis.set_ylim(0.0, max(values) * 1.28)
	for bar, value in zip(bars, values):
		axis.text(
			bar.get_x() + bar.get_width() / 2.0,
			bar.get_height(),
			f"{value:.6g}",
			ha="center",
			va="bottom",
			fontsize=10,
		)
	style_axis(axis)
	figure.suptitle(
		f"相同计算预算下的柔度（第{budget:g}次有限元求解）",
		fontsize=12,
		weight="bold",
	)
	figure.text(0.5, 0.105, experiment_caption(sequential), ha="center", fontsize=8)
	figure.text(
		0.5,
		0.065,
		f"同步优化相对顺序优化的柔度降低幅度：{reduction:.2f}%",
		ha="center",
		fontsize=8,
		color=TARGET_COLOR,
	)
	figure.text(
		0.5,
		0.025,
		threshold_note(sequential, simultaneous),
		ha="center",
		fontsize=8,
		color="#8B0000",
	)
	figure.subplots_adjust(left=0.12, right=0.97, bottom=0.25, top=0.82)
	save_figure(figure, output_dir / "02_compliance_at_equal_fe_budget.png")


def plot_constraints(
	sequential: RunRecord,
	simultaneous: RunRecord,
	output_dir: Path,
) -> None:
	figure, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
	volume_target = float(sequential.config["mmc_target_rate"])
	ratio_target = float(sequential.config["material_ratio_f"])

	for run, color, label in (
		(sequential, SEQUENTIAL_COLOR, "顺序优化"),
		(simultaneous, SYNC_COLOR, "同步优化"),
	):
		budget = common_fe_budget(sequential, simultaneous)
		x, y = series(run, "mmc_volume", max_fe_budget=budget)
		axes[0].plot(x, y, color=color, linewidth=2.0, label=label)
		x, y = series(run, "material_ratio", max_fe_budget=budget)
		axes[1].plot(x, y, color=color, linewidth=2.0, label=label)

	axes[0].axhline(volume_target, color=TARGET_COLOR, linestyle="--", linewidth=1.2)
	axes[1].axhline(ratio_target, color=TARGET_COLOR, linestyle="--", linewidth=1.2)
	axes[0].set_title("MMC体积分数", fontsize=10, weight="bold")
	axes[1].set_title("材料体积比", fontsize=10, weight="bold")
	for axis in axes:
		axis.set_xlim(left=0.0, right=budget * 1.02)
		axis.set_xlabel("累计有限元求解次数")
		axis.legend(frameon=False, fontsize=8)
		style_axis(axis)
	figure.suptitle("体积约束与材料比例约束变化", fontsize=12, weight="bold")
	figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
	save_figure(figure, output_dir / "03_volume_and_material_ratio_constraints.png")


def write_summary(
	sequential: RunRecord,
	simultaneous: RunRecord,
	output_dir: Path,
) -> None:
	fields = [
		("compliance", "compliance"),
		("mmc_volume", "mmc_volume"),
		("material_ratio", "material_ratio"),
		("material_grayness", "material_grayness"),
		("geometry_change", "geometry_change"),
		("material_change", "material_change"),
		("cumulative_fe_solves", "cumulative_fe_solves"),
		("elapsed_seconds", "elapsed_seconds"),
		("active_cell_count", "active_cell_count"),
	]
	with (output_dir / "comparison_summary.tsv").open(
		"w", encoding="utf-8", newline=""
	) as stream:
		writer = csv.writer(stream, delimiter="\t")
		writer.writerow(
			[
				"method",
				"run_id",
				"comparison_point",
				"model_id",
				"boundary_signature",
				"material_active_threshold",
				*[name for name, _ in fields],
			]
		)
		budget = common_fe_budget(sequential, simultaneous)
		for run in (sequential, simultaneous):
			# 汇总表只记录统一的前 100 次求解比较点，不再混入顺序版
			# 第 160 次求解的终值，防止后续论文统计误用不同计算预算。
			row = row_at_or_before_budget(run, budget)
			writer.writerow(
				[
					run.method,
					run.run_id,
					f"equal_fe_budget_{budget:g}",
					run.config["model_id"],
					run.config["boundary_signature"],
					run.config.get("material_active_threshold", ""),
					*[row.get(field, "") for _, field in fields],
				]
			)


def default_output_directory(sequential: RunRecord) -> Path:
	records_root = sequential.config_path.parent.parent
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	boundary = sequential.config["boundary_signature"][:8]
	return records_root / (
		f"plots_compare_{sequential.config['model_id']}_bc{boundary}_equal_budget_{timestamp}"
	)


def parse_arguments() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Plot the latest compatible sequential and simultaneous MMC-SIMP records."
	)
	parser.add_argument(
		"--search-root",
		type=Path,
		default=Path("out/build/bin"),
		help="Directory containing the optimization result trees.",
	)
	parser.add_argument("--sequential-run-id", default=None)
	parser.add_argument("--simultaneous-run-id", default=None)
	parser.add_argument("--output-dir", type=Path, default=None)
	return parser.parse_args()


def main() -> None:
	configure_chinese_font()
	arguments = parse_arguments()
	runs = discover_runs(arguments.search_root)
	sequential, simultaneous = select_pair(
		runs,
		arguments.sequential_run_id,
		arguments.simultaneous_run_id,
	)
	output_dir = arguments.output_dir or default_output_directory(sequential)
	output_dir.mkdir(parents=True, exist_ok=True)

	plot_compliance(sequential, simultaneous, output_dir)
	plot_equal_budget_compliance(sequential, simultaneous, output_dir)
	plot_constraints(sequential, simultaneous, output_dir)
	write_summary(sequential, simultaneous, output_dir)

	print(f"Sequential run: {sequential.run_id}")
	print(f"Simultaneous run: {simultaneous.run_id}")
	print(f"Output directory: {output_dir.resolve()}")
	budget = common_fe_budget(sequential, simultaneous)
	seq_row = row_at_or_before_budget(sequential, budget)
	sim_row = row_at_or_before_budget(simultaneous, budget)
	seq_compliance = final_number(seq_row, "compliance")
	sim_compliance = final_number(sim_row, "compliance")
	print(
		f"Equal FE budget: {budget:g}; sequential compliance={seq_compliance:.12g}; "
		f"simultaneous compliance={sim_compliance:.12g}; "
		f"reduction={(seq_compliance - sim_compliance) / seq_compliance * 100.0:.4f}%"
	)
	if sequential.config.get("material_active_threshold") != simultaneous.config.get(
		"material_active_threshold"
	):
		print(
			"WARNING: material_active_threshold differs: "
			f"sequential={sequential.config.get('material_active_threshold')}, "
			f"simultaneous={simultaneous.config.get('material_active_threshold')}"
		)


if __name__ == "__main__":
	main()
