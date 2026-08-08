#!/usr/bin/env python3
"""Tkinter front end for thicken_hex_mesh.py."""

from __future__ import annotations

import ctypes
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from thicken_hex_mesh import (
    AXIS_INDEX,
    AXIS_NAMES,
    extrude_voxel_topology,
    read_mesh,
    thin_axes,
    validate_unit_voxel_mesh,
    write_mesh,
)


SIDE_LABELS = {
    "正方向": "positive",
    "负方向": "negative",
}


class MeshThickenerApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("六面体网格加厚工具")
        self.root.geometry("720x510")
        self.root.minsize(660, 470)

        self.input_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.axis_var = tk.StringVar(value="x")
        self.layers_var = tk.IntVar(value=1)
        self.side_label_var = tk.StringVar(value=next(iter(SIDE_LABELS)))
        self.overwrite_var = tk.BooleanVar(value=False)
        self.preview_var = tk.StringVar(value="请选择一个规则六面体 .mesh 文件")
        self.status_var = tk.StringVar(value="就绪")
        self._auto_output = True
        self._mesh_cache: tuple | None = None
        self._cached_input: Path | None = None

        self._configure_style()
        self._build_ui()
        for variable in (self.axis_var, self.layers_var, self.side_label_var):
            variable.trace_add("write", self._settings_changed)

    def _configure_style(self) -> None:
        style = ttk.Style()
        if "vista" in style.theme_names():
            style.theme_use("vista")
        style.configure("Title.TLabel", font=("Microsoft YaHei UI", 16, "bold"))
        style.configure("Section.TLabelframe.Label", font=("Microsoft YaHei UI", 10, "bold"))
        style.configure("Primary.TButton", font=("Microsoft YaHei UI", 10, "bold"), padding=(18, 8))
        style.configure("Status.TLabel", foreground="#2457a6")

    def _build_ui(self) -> None:
        outer = ttk.Frame(self.root, padding=20)
        outer.pack(fill="both", expand=True)

        ttk.Label(outer, text="六面体网格加厚", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            outer,
            text="适用于任意拓扑形状：单元边长为 1，所选厚度轴只有 1 层",
        ).pack(anchor="w", pady=(3, 16))

        files = ttk.LabelFrame(outer, text="文件", padding=12, style="Section.TLabelframe")
        files.pack(fill="x")
        files.columnconfigure(1, weight=1)
        ttk.Label(files, text="输入网格").grid(row=0, column=0, sticky="w", padx=(0, 10), pady=5)
        input_entry = ttk.Entry(files, textvariable=self.input_var)
        input_entry.grid(row=0, column=1, sticky="ew", pady=5)
        input_entry.bind("<Return>", self._load_input)
        input_entry.bind("<FocusOut>", self._load_input)
        ttk.Button(files, text="浏览...", command=self._browse_input).grid(row=0, column=2, padx=(8, 0), pady=5)
        ttk.Label(files, text="输出网格").grid(row=1, column=0, sticky="w", padx=(0, 10), pady=5)
        output_entry = ttk.Entry(files, textvariable=self.output_var)
        output_entry.grid(row=1, column=1, sticky="ew", pady=5)
        output_entry.bind("<Key>", self._output_edited)
        ttk.Button(files, text="浏览...", command=self._browse_output).grid(row=1, column=2, padx=(8, 0), pady=5)

        settings = ttk.LabelFrame(outer, text="加厚设置", padding=12, style="Section.TLabelframe")
        settings.pack(fill="x", pady=(14, 0))
        settings.columnconfigure(1, weight=1)
        settings.columnconfigure(3, weight=1)

        ttk.Label(settings, text="加厚轴").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=6)
        self.axis_combo = ttk.Combobox(
            settings,
            textvariable=self.axis_var,
            values=("x", "y", "z"),
            state="readonly",
            width=12,
        )
        self.axis_combo.grid(row=0, column=1, sticky="ew", padx=(0, 18), pady=6)

        ttk.Label(settings, text="增加层数").grid(row=0, column=2, sticky="w", padx=(0, 8), pady=6)
        ttk.Spinbox(settings, from_=1, to=1000, textvariable=self.layers_var, width=12).grid(
            row=0, column=3, sticky="ew", pady=6
        )

        ttk.Label(settings, text="添加方向").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=6)
        ttk.Combobox(
            settings,
            textvariable=self.side_label_var,
            values=tuple(SIDE_LABELS),
            state="readonly",
        ).grid(row=1, column=1, sticky="ew", padx=(0, 18), pady=6)

        ttk.Checkbutton(settings, text="允许覆盖已有输出文件", variable=self.overwrite_var).grid(
            row=2, column=0, columnspan=4, sticky="w", pady=(6, 0)
        )

        preview = ttk.LabelFrame(outer, text="网格预览", padding=12, style="Section.TLabelframe")
        preview.pack(fill="x", pady=(14, 0))
        ttk.Label(preview, textvariable=self.preview_var).pack(anchor="w")

        footer = ttk.Frame(outer)
        footer.pack(fill="x", side="bottom", pady=(18, 0))
        ttk.Label(footer, textvariable=self.status_var, style="Status.TLabel").pack(side="left")
        ttk.Button(footer, text="关闭", command=self.root.destroy).pack(side="right")
        ttk.Button(footer, text="生成网格", command=self._generate, style="Primary.TButton").pack(
            side="right", padx=(0, 10)
        )

    def _browse_input(self) -> None:
        filename = filedialog.askopenfilename(
            title="选择输入网格",
            initialdir=str(Path(__file__).resolve().parents[1] / "input" / "topopt" / "data"),
            filetypes=(("Medit mesh", "*.mesh"), ("所有文件", "*.*")),
        )
        if not filename:
            return
        self.input_var.set(filename)
        self._auto_output = True
        self._load_input()

    def _browse_output(self) -> None:
        initial = self.output_var.get().strip()
        filename = filedialog.asksaveasfilename(
            title="选择输出位置",
            initialdir=str(Path(initial).parent) if initial else None,
            initialfile=Path(initial).name if initial else None,
            defaultextension=".mesh",
            filetypes=(("Medit mesh", "*.mesh"), ("所有文件", "*.*")),
        )
        if filename:
            self.output_var.set(filename)
            self._auto_output = False

    def _output_edited(self, _event: tk.Event) -> None:
        self._auto_output = False

    def _settings_changed(self, *_args: object) -> None:
        if self._auto_output:
            self._set_default_output()
        self._refresh_preview()

    def _set_default_output(self) -> None:
        raw_input = self.input_var.get().strip()
        if not raw_input:
            return
        input_path = Path(raw_input)
        side = SIDE_LABELS.get(self.side_label_var.get(), "positive")
        try:
            layers = max(1, int(self.layers_var.get()))
        except (tk.TclError, ValueError):
            layers = 1
        output = input_path.with_name(
            f"{input_path.stem}_thick_{self.axis_var.get()}_{side}_{layers}.mesh"
        )
        self.output_var.set(str(output))

    def _load_input(self, _event: tk.Event | None = None) -> None:
        input_path = Path(self.input_var.get().strip())
        try:
            vertices, hexahedra = read_mesh(input_path)
            validated = validate_unit_voxel_mesh(vertices, hexahedra)
            available_axes = thin_axes(validated[1])
            if not available_axes:
                raise ValueError("没有找到只有 1 个单元厚度的轴")
        except (OSError, ValueError, KeyError) as exc:
            self._mesh_cache = None
            self._cached_input = None
            self.preview_var.set(f"无法读取网格：{exc}")
            self.status_var.set("输入网格无效")
            return
        self._mesh_cache = (vertices, hexahedra, validated)
        self._cached_input = input_path
        available_axis_names = tuple(AXIS_NAMES[axis] for axis in available_axes)
        self.axis_combo.configure(values=available_axis_names)
        if self.axis_var.get() not in available_axis_names:
            self.axis_var.set(available_axis_names[0])
        self._set_default_output()
        self._refresh_preview()
        self.status_var.set("输入网格检查通过")

    def _refresh_preview(self) -> None:
        if self._mesh_cache is None:
            return
        try:
            layers = int(self.layers_var.get())
        except (tk.TclError, ValueError):
            self.preview_var.set("层数必须是正整数")
            return
        if layers < 1:
            self.preview_var.set("层数必须至少为 1")
            return

        cells = self._mesh_cache[2][1]
        axis_name = self.axis_var.get()
        axis = AXIS_INDEX[axis_name]
        current_layers = len({origin[axis] for origin in cells})
        if current_layers != 1:
            self.preview_var.set(f"{axis_name} 轴有 {current_layers} 层，不能作为 2.5D 厚度轴")
            return
        self.preview_var.set(
            f"厚度轴 {axis_name}：1 → {layers + 1} 层    "
            f"单元数量：{len(cells)} → {len(cells) * (layers + 1)}"
        )

    def _generate(self) -> None:
        input_path = Path(self.input_var.get().strip())
        output_text = self.output_var.get().strip()
        output_path = Path(output_text) if output_text else None
        if not input_path.is_file():
            messagebox.showerror("输入错误", "请选择有效的输入 .mesh 文件。")
            return
        if output_path is None:
            messagebox.showerror("输出错误", "请选择输出文件。")
            return
        if self._cached_input != input_path or self._mesh_cache is None:
            self._load_input()
            if self._mesh_cache is None:
                messagebox.showerror("网格错误", self.preview_var.get())
                return
        if output_path.exists() and not self.overwrite_var.get():
            messagebox.showerror("文件已存在", "输出文件已经存在。勾选允许覆盖或选择其他文件名。")
            return

        try:
            layers = int(self.layers_var.get())
            if layers < 1:
                raise ValueError("层数必须至少为 1")
            vertices, hexahedra, validated = self._mesh_cache
            vertex_refs, cells, corner_pattern = validated
            side = SIDE_LABELS[self.side_label_var.get()]
            new_vertices, new_hexahedra = extrude_voxel_topology(
                vertex_refs,
                cells,
                corner_pattern,
                AXIS_INDEX[self.axis_var.get()],
                layers,
                side,
            )
            output_path.parent.mkdir(parents=True, exist_ok=True)
            write_mesh(output_path, new_vertices, new_hexahedra)
            validate_unit_voxel_mesh(*read_mesh(output_path))
        except (OSError, ValueError, KeyError) as exc:
            self.status_var.set("生成失败")
            messagebox.showerror("生成失败", str(exc))
            return

        self.status_var.set(f"生成成功：{output_path.name}")
        messagebox.showinfo(
            "生成成功",
            f"输出文件：\n{output_path}\n\n"
            f"厚度：1 → {layers + 1} 层\n"
            f"单元数量：{len(hexahedra)} → {len(new_hexahedra)}",
        )


def main() -> None:
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except (AttributeError, OSError):
        pass
    root = tk.Tk()
    MeshThickenerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
