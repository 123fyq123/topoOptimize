#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace mmc_layout {

// Shared chair layout for pure MMC and simultaneous MMC-SIMP. The target
// domain is a 24 x 24 x 40 chair with a 16-cell base and a 6 x 24 backrest.
template <typename TComponent>
std::vector<TComponent> MakeChairLayout(const std::vector<double>& bounds) {
	if (bounds.size() != 6) return {};
	const double dx = bounds[1] - bounds[0];
	const double dy = bounds[3] - bounds[2];
	const double dz = bounds[5] - bounds[4];
	if (!(dx > 0.0 && dy > 0.0 && dz > 0.0)) return {};

	// These fractions reproduce the generated 24 x 24 x 40 chair domain.
	const double base_height = 0.40 * dz;
	const double back_height = dz - base_height;
	const double back_thickness = 0.25 * dy;
	const double seat_z = bounds[4] + base_height;
	const double back_front_y = bounds[3] - back_thickness;

	constexpr int base_columns = 4;
	constexpr int base_rows = 4;
	constexpr int base_layers = 2;
	constexpr int back_columns = 4;
	constexpr int back_rows = 4;
	const double base_cell_x = dx / static_cast<double>(base_columns);
	const double base_cell_y = dy / static_cast<double>(base_rows);
	const double base_cell_z = base_height / static_cast<double>(base_layers);
	const double back_cell_x = dx / static_cast<double>(back_columns);
	const double back_cell_z = back_height / static_cast<double>(back_rows);
	const double radius = 0.07 * std::min({
		base_cell_x,
		base_cell_y,
		base_cell_z,
		back_thickness,
		back_cell_x,
		back_cell_z });

	std::vector<TComponent> components;
	components.reserve(437);

	const auto add_component = [&](double center_x,
		double center_y,
		double center_z,
		double vector_x,
		double vector_y,
		double vector_z) {
		const double planar_length = std::sqrt(
			vector_x * vector_x + vector_y * vector_y);
		const double half_length = std::sqrt(
			planar_length * planar_length + vector_z * vector_z);
		if (!(half_length > 0.0)) return;
		components.emplace_back(
			center_x,
			center_y,
			center_z,
			half_length,
			radius,
			radius,
			0.0,
			-std::atan2(vector_z, planar_length),
			std::atan2(vector_y, vector_x));
	};

	const auto add_crossed_cell = [&](double center_x,
		double center_y,
		double center_z,
		double cell_x,
		double cell_y,
		double cell_z) {
		const double half_x = 0.5 * cell_x;
		const double half_y = 0.5 * cell_y;
		const double half_z = 0.5 * cell_z;

		// 三个主平面都布置双对角杆；从 X/Y/Z 任一方向观察都能看到清晰 X 形。
		add_component(center_x, center_y, center_z, half_x, half_y, 0.0);
		add_component(center_x, center_y, center_z, half_x, -half_y, 0.0);
		add_component(center_x, center_y, center_z, half_x, 0.0, half_z);
		add_component(center_x, center_y, center_z, half_x, 0.0, -half_z);
		add_component(center_x, center_y, center_z, 0.0, half_y, half_z);
		add_component(center_x, center_y, center_z, 0.0, half_y, -half_z);

		// 三根轴向杆与六根斜杆在中心相交，组成三个相互垂直的局部“米”字。
		add_component(center_x, center_y, center_z, half_x, 0.0, 0.0);
		add_component(center_x, center_y, center_z, 0.0, half_y, 0.0);
		add_component(center_x, center_y, center_z, 0.0, 0.0, half_z);
	};

	for (int layer = 0; layer < base_layers; ++layer) {
		for (int row = 0; row < base_rows; ++row) {
			for (int column = 0; column < base_columns; ++column) {
				add_crossed_cell(
					bounds[0] + (static_cast<double>(column) + 0.5) * base_cell_x,
					bounds[2] + (static_cast<double>(row) + 0.5) * base_cell_y,
					bounds[4] + (static_cast<double>(layer) + 0.5) * base_cell_z,
					base_cell_x,
					base_cell_y,
					base_cell_z);
			}
		}
	}

	const double back_center_y = 0.5 * (back_front_y + bounds[3]);
	for (int row = 0; row < back_rows; ++row) {
		for (int column = 0; column < back_columns; ++column) {
			add_crossed_cell(
				bounds[0] + (static_cast<double>(column) + 0.5) * back_cell_x,
				back_center_y,
				seat_z + (static_cast<double>(row) + 0.5) * back_cell_z,
				back_cell_x,
				back_thickness,
				back_cell_z);
		}
	}

	// Local collectors cover the current central seat load and backrest load.
	const double center_x = 0.5 * (bounds[0] + bounds[1]);
	const double seat_load_y = bounds[2] + 0.375 * dy;
	add_component(center_x, seat_load_y, seat_z, 0.5 * base_cell_x, 0.0, 0.0);
	add_component(center_x, seat_load_y, seat_z, 0.0, 0.5 * base_cell_y, 0.0);

	const double back_load_z = seat_z + 0.5 * back_height;
	add_component(
		center_x, back_front_y, back_load_z,
		0.5 * back_cell_x, 0.0, 0.0);
	add_component(
		center_x, back_front_y, back_load_z,
		0.0, 0.5 * back_thickness, 0.0);
	add_component(
		center_x, back_front_y, back_load_z,
		0.0, 0.0, 0.5 * back_cell_z);

	return components;
}

} // namespace mmc_layout
