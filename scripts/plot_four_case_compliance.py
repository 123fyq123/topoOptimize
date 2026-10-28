#!/usr/bin/env python3
"""绘制四种 MMC-SIMP 工况下顺序与同步优化的柔度曲线。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_optimization_records import (
	MAX_COMPARISON_FE_BUDGET,
	SEQUENTIAL_COLOR,
	SYNC_COLOR,
	RunRecord,
	common_fe_budget,
	configure_chinese_font,
	discover_runs,
	final_number,
	parameter_key,
	row_at_or_before_budget,
	series,
	style_axis,
)


CASE_DEFINITIONS = (
	("hex_mesh_64_1_32", "9ab66eb2", "右下角竖向受力"),
	("hex_mesh_64_1_32", "a733a7eb", "右端中点竖向受力"),
	("hex_mesh_64_1_32", "df28892e", "右端斜向受力"),
	("hex_mesh_64_1_64", "16a0b7c1", "64×1×64右端水平受力"),
)


def latest_pair(
	runs: list[RunRecord],
	model_id: str,
	boundary_prefix: str,
) -> tuple[RunRecord, RunRecord]:
	"""选择同一模型和边界条件下最新且材料参数一致的一对记录。"""
	candidates = [
		run
		for run in runs
		if run.config.get("model_id") == model_id
		and run.config.get("boundary_signature", "").startswith(boundary_prefix)
	]
	sequential_runs = [run for run in candidates if run.method == "sequential"]
	simultaneous_runs = [run for run in candidates if run.method == "simultaneous"]
	compatible_pairs = [
		(sequential, simultaneous)
		for sequential in sequential_runs
		for simultaneous in simultaneous_runs
		if parameter_key(sequential) == parameter_key(simultaneous)
	]
	if not compatible_pairs:
		raise ValueError(
			f"没有找到完整配对：model={model_id}, boundary={boundary_prefix}"
		)
	return max(
		compatible_pairs,
		key=lambda pair: (
			pair[0].config_path.stat().st_mtime,
			pair[1].config_path.stat().st_mtime,
		),
	)


def plot_case(
	axis: plt.Axes,
	sequential: RunRecord,
	simultaneous: RunRecord,
	title: str,
) -> tuple[float, float, float, float]:
	"""在一个子图中绘制同一工况的两种耦合柔度，并返回比较指标。"""
	budget = min(
		MAX_COMPARISON_FE_BUDGET,
		common_fe_budget(sequential, simultaneous),
	)
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
	axis.plot(sx, sy, color=SEQUENTIAL_COLOR, linewidth=2.2, label="顺序优化")
	axis.plot(mx, my, color=SYNC_COLOR, linewidth=2.2, label="同步优化")

	seq_row = row_at_or_before_budget(sequential, budget)
	sim_row = row_at_or_before_budget(simultaneous, budget)
	seq_value = final_number(seq_row, "compliance")
	sim_value = final_number(sim_row, "compliance")
	reduction = (seq_value - sim_value) / seq_value * 100.0
	axis.scatter([budget], [seq_value], color=SEQUENTIAL_COLOR, s=24, zorder=4)
	axis.scatter([budget], [sim_value], color=SYNC_COLOR, s=24, zorder=4)
	axis.set_title(title, fontsize=11, weight="bold")
	axis.set_xlim(0.0, budget * 1.02)
	axis.set_xlabel("累计有限元求解次数")
	axis.set_ylabel("归一化柔度")
	axis.set_facecolor("#F7F7F7")
	axis.text(
		0.98,
		0.95,
		f"第{budget:g}次：同步比顺序低 {reduction:.2f}%",
		transform=axis.transAxes,
		ha="right",
		va="top",
		fontsize=8.5,
		bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82},
	)
	style_axis(axis)
	return budget, seq_value, sim_value, reduction


def parse_arguments() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="绘制四工况柔度对比图。")
	parser.add_argument("--search-root", type=Path, default=Path("out/build/bin"))
	parser.add_argument(
		"--output",
		type=Path,
		default=Path(
			"out/build/bin/mmc_simp_result/records/plots_20260928_four_cases/"
			"00_four_case_compliance_summary.png"
		),
	)
	return parser.parse_args()


def main() -> None:
	arguments = parse_arguments()
	configure_chinese_font()
	runs = discover_runs(arguments.search_root)
	figure, axes = plt.subplots(2, 2, figsize=(12.2, 8.0))

	for axis, (model_id, boundary_prefix, title) in zip(
		axes.flat, CASE_DEFINITIONS
	):
		sequential, simultaneous = latest_pair(runs, model_id, boundary_prefix)
		budget, seq_value, sim_value, reduction = plot_case(
			axis, sequential, simultaneous, title
		)
		print(
			f"{title}: budget={budget:g}, sequential={seq_value:.12g}, "
			f"simultaneous={sim_value:.12g}, reduction={reduction:.4f}%"
		)

	handles, labels = axes.flat[0].get_legend_handles_labels()
	figure.legend(
		handles,
		labels,
		loc="upper center",
		ncol=2,
		frameon=False,
		bbox_to_anchor=(0.5, 0.955),
	)
	figure.suptitle("四种工况下顺序与同步优化的柔度变化", fontsize=14, weight="bold")
	figure.text(
		0.5,
		0.018,
		"黄色：顺序优化；蓝色：同步优化。右下角工况的顺序版活跃阈值仍为0.3，其余工况均为0.01。",
		ha="center",
		fontsize=9,
		color="#8B0000",
	)
	figure.tight_layout(rect=(0.02, 0.05, 0.98, 0.90), h_pad=2.0, w_pad=1.5)
	arguments.output.parent.mkdir(parents=True, exist_ok=True)
	figure.savefig(arguments.output, dpi=300, bbox_inches="tight", facecolor="white")
	plt.close(figure)
	print(f"输出：{arguments.output.resolve()}")


if __name__ == "__main__":
	main()
