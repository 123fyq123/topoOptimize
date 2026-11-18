"""绘制 MMC-SIMP 论文研究路线与推荐验证算例示意图。"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from matplotlib import font_manager


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "docs" / "figures"


def configure_chinese_font():
    """优先使用 Windows 微软雅黑，避免论文示意图中文乱码。"""
    font_path = Path(r"C:\Windows\Fonts\msyh.ttc")
    if font_path.exists():
        font_manager.fontManager.addfont(str(font_path))
        plt.rcParams["font.family"] = font_manager.FontProperties(
            fname=str(font_path)
        ).get_name()
    plt.rcParams["axes.unicode_minus"] = False


def box(ax, xy, width, height, text, color):
    patch = Rectangle(
        xy,
        width,
        height,
        facecolor=color,
        edgecolor="#263238",
        linewidth=1.3,
        joinstyle="round",
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=10,
    )
    return patch


def arrow(ax, start, end, color="#455a64", style="-|>"):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            mutation_scale=13,
            linewidth=1.4,
            color=color,
            connectionstyle="arc3,rad=0.0",
        )
    )


def draw_research_route():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    fig.suptitle("体细分 MMC-SIMP 多材料拓扑优化研究路线", fontsize=16, weight="bold")

    ax = axes[0]
    ax.set_title("创新点一：顺序优化（基线方法）", fontsize=13, weight="bold")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    box(ax, (0.6, 7.4), 2.2, 1.2, "体细分/等几何\n结构分析", "#d9edf7")
    box(ax, (3.8, 7.4), 2.2, 1.2, "单材料 MMC\n优化几何", "#ffe082")
    box(ax, (7.0, 7.4), 2.2, 1.2, "固定 MMC 域\n多材料 SIMP", "#90caf9")
    arrow(ax, (2.8, 8.0), (3.8, 8.0))
    arrow(ax, (6.0, 8.0), (7.0, 8.0))
    box(ax, (3.0, 4.3), 4.0, 1.3, "显式边界 + 域内材料布局", "#c8e6c9")
    arrow(ax, (8.1, 7.4), (6.2, 5.6))
    ax.text(
        5,
        2.5,
        "特点：实现简单、收敛稳定\n局限：材料性能不反作用于 MMC 受力路径",
        ha="center",
        va="center",
        fontsize=11,
        color="#37474f",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#fafafa", edgecolor="#b0bec5"),
    )

    ax = axes[1]
    ax.set_title("创新点二：体细分下的耦合同步协同优化", fontsize=13, weight="bold")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    box(ax, (0.5, 7.5), 2.2, 1.2, "MMC 几何变量\n$d$", "#ffe082")
    box(ax, (7.3, 7.5), 2.2, 1.2, "SIMP 材料变量\n$\\eta$", "#90caf9")
    box(
        ax,
        (3.0, 7.2),
        4.0,
        1.8,
        "统一耦合场\n$a_e=m_e(d)E_e(\\eta_e)$",
        "#c8e6c9",
    )
    arrow(ax, (2.7, 8.1), (3.0, 8.1))
    arrow(ax, (7.3, 8.1), (7.0, 8.1))
    box(ax, (3.4, 4.8), 3.2, 1.2, "同一体细分分析\n$K(a)U=F$", "#d9edf7")
    arrow(ax, (5.0, 7.2), (5.0, 6.0))
    box(ax, (0.7, 2.2), 3.6, 1.5, "一致 MMC 灵敏度\n$\\partial C/\\partial d$\n$\\rightarrow$ MMA", "#fff3cd")
    box(ax, (5.7, 2.2), 3.6, 1.5, "一致材料灵敏度\n$\\partial C/\\partial\\eta$\n$\\rightarrow$ OC/MMA", "#dbeafe")
    arrow(ax, (4.3, 4.8), (2.7, 3.7))
    arrow(ax, (5.7, 4.8), (7.3, 3.7))
    arrow(ax, (2.0, 2.2), (1.6, 7.5), color="#ef6c00")
    arrow(ax, (8.0, 2.2), (8.4, 7.5), color="#1565c0")
    ax.text(
        5,
        0.8,
        "同一位移场、同一目标函数、同一约束状态；材料与几何双向反馈",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#37474f",
    )

    path = OUTPUT_DIR / "mmc_simp_research_route.png"
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def fixed_edge(ax, x, y0, y1):
    ax.plot([x, x], [y0, y1], color="#d32f2f", linewidth=4)
    for y in [y0 + i * (y1 - y0) / 7 for i in range(8)]:
        ax.plot([x - 0.18, x], [y - 0.12, y + 0.12], color="#d32f2f", linewidth=1)


def load_arrow(ax, start, end, label="$F$"):
    arrow(ax, start, end, color="#1565c0")
    ax.text(start[0] + 0.15, start[1] + 0.12, label, color="#1565c0", fontsize=10)


def support(ax, x, y):
    ax.add_patch(
        Polygon([[x, y], [x - 0.22, y - 0.32], [x + 0.22, y - 0.32]],
                closed=True, facecolor="#9e9e9e", edgecolor="#424242")
    )


def setup_benchmark_ax(ax, title):
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_aspect("equal")
    ax.axis("off")


def draw_extruded_prism(ax, front, depth=(0.35, 0.28), color="#b2dfdb"):
    back = [(x + depth[0], y + depth[1]) for x, y in front]
    ax.add_patch(Polygon(back, closed=True, facecolor="#80cbc4", edgecolor="#37474f"))
    ax.add_patch(Polygon(front, closed=True, facecolor=color, edgecolor="#37474f", linewidth=1.2))
    for p, q in zip(front, back):
        ax.plot([p[0], q[0]], [p[1], q[1]], color="#37474f", linewidth=1)


def draw_benchmarks():
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6), constrained_layout=True)
    fig.suptitle("建议的二维/三维多材料拓扑优化验证模型", fontsize=16, weight="bold")

    ax = axes[0, 0]
    setup_benchmark_ax(ax, "二维悬臂梁：基础正确性")
    ax.add_patch(Rectangle((0, 0), 4.5, 2, facecolor="#e0f2f1", edgecolor="#37474f"))
    fixed_edge(ax, 0, 0, 2)
    load_arrow(ax, (4.2, 1.0), (4.2, 0.15))
    ax.set_xlim(-0.4, 5.0); ax.set_ylim(-0.5, 2.5)

    ax = axes[0, 1]
    setup_benchmark_ax(ax, "二维 MBB 梁：标准对比")
    ax.add_patch(Rectangle((0, 0), 5, 1.8, facecolor="#e0f2f1", edgecolor="#37474f"))
    support(ax, 0.3, 0); support(ax, 4.7, 0)
    load_arrow(ax, (2.5, 2.5), (2.5, 1.7))
    ax.set_xlim(-0.3, 5.3); ax.set_ylim(-0.6, 2.8)

    ax = axes[0, 2]
    setup_benchmark_ax(ax, "二维 L 型梁：非凸域与应力集中")
    l_shape = [(0, 0), (4.2, 0), (4.2, 1.4), (1.5, 1.4), (1.5, 4.3), (0, 4.3)]
    ax.add_patch(Polygon(l_shape, closed=True, facecolor="#e0f2f1", edgecolor="#37474f"))
    ax.plot([0, 1.5], [4.3, 4.3], color="#d32f2f", linewidth=4)
    load_arrow(ax, (4.0, 1.1), (4.0, 0.25))
    ax.set_xlim(-0.4, 4.8); ax.set_ylim(-0.4, 4.9)

    ax = axes[1, 0]
    setup_benchmark_ax(ax, "三维悬臂梁：三维耦合能力")
    front = [(0, 0), (4.5, 0), (4.5, 1.7), (0, 1.7)]
    draw_extruded_prism(ax, front)
    fixed_edge(ax, 0, 0, 1.7)
    load_arrow(ax, (4.65, 1.05), (4.65, 0.15))
    ax.set_xlim(-0.4, 5.4); ax.set_ylim(-0.5, 2.5)

    ax = axes[1, 1]
    setup_benchmark_ax(ax, "三维 L 型梁：复杂设计域")
    front = [(0, 0), (4.0, 0), (4.0, 1.3), (1.4, 1.3), (1.4, 4.1), (0, 4.1)]
    draw_extruded_prism(ax, front, depth=(0.32, 0.25))
    ax.plot([0, 1.4], [4.1, 4.1], color="#d32f2f", linewidth=4)
    load_arrow(ax, (4.15, 1.0), (4.15, 0.15))
    ax.set_xlim(-0.4, 4.8); ax.set_ylim(-0.4, 4.8)

    ax = axes[1, 2]
    setup_benchmark_ax(ax, "三维扭转梁：空间载荷路径")
    front = [(0, 0), (4.5, 0), (4.5, 1.7), (0, 1.7)]
    draw_extruded_prism(ax, front, depth=(0.42, 0.34))
    fixed_edge(ax, 0, 0, 1.7)
    load_arrow(ax, (4.75, 1.9), (4.75, 2.65), "$+F$")
    load_arrow(ax, (4.45, 0.05), (4.45, -0.7), "$-F$")
    ax.text(3.0, 2.35, "力偶形成扭矩", fontsize=9, color="#455a64")
    ax.set_xlim(-0.4, 5.5); ax.set_ylim(-1.0, 3.0)

    fig.text(
        0.5,
        -0.01,
        "红色：固定边界；蓝色：外载荷。二维用于公式与收敛验证，三维用于体现体细分分析和工程扩展价值。",
        ha="center",
        fontsize=10,
        color="#455a64",
    )
    path = OUTPUT_DIR / "mmc_simp_recommended_benchmarks.png"
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def main():
    configure_chinese_font()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for output in (draw_research_route(), draw_benchmarks()):
        print(output)


if __name__ == "__main__":
    main()
