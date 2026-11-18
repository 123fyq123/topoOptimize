#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace mmc_layout {

// Shared by the pure-MMC and simultaneous solvers. The target domain is the
// 48 x 10 x 48 L beam whose two in-plane arms are both 16 cells wide.
template <typename TComponent>
std::vector<TComponent> MakeThreeDimensionalLBeamLayout(
	const std::vector<double>& bounds) {
	if (bounds.size() != 6) return {};
	const double dx = bounds[1] - bounds[0];
	const double dy = bounds[3] - bounds[2];
	const double dz = bounds[5] - bounds[4];
	if (!(dx > 0.0 && dy > 0.0 && dz > 0.0)) return {};

	const double arm_width = std::min(dx, dz) / 3.0;
	const int cells_across_arm = 2;
	const int extension_cell_count = 4;
	const double arm_cell_size =
		arm_width / static_cast<double>(cells_across_arm);
	const double horizontal_cell_size =
		(dx - arm_width) / static_cast<double>(extension_cell_count);
	const double vertical_cell_size =
		(dz - arm_width) / static_cast<double>(extension_cell_count);
	if (!(horizontal_cell_size > 0.0 && vertical_cell_size > 0.0)) return {};

	const double center_y = 0.5 * (bounds[2] + bounds[3]);
	const double layer_offset = 0.25 * dy;
	const double radius = 0.08 * std::min(
		std::min(arm_cell_size, horizontal_cell_size),
		std::min(vertical_cell_size, dy));
	std::vector<TComponent> components;
	components.reserve(106);

	const auto add_component = [&](double center_x,
		double component_center_y,
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
			component_center_y,
			center_z,
			half_length,
			radius,
			radius,
			0.0,
			-std::atan2(vector_z, planar_length),
			std::atan2(vector_y, vector_x));
	};

	const auto add_three_dimensional_cell = [&](double cell_center_x,
		double cell_center_z,
		double cell_size_x,
		double cell_size_z) {
		const double half_x = 0.5 * cell_size_x;
		const double half_z = 0.5 * cell_size_z;
		// Four body diagonals connect all opposite corners of the local 3D cell.
		for (const double direction_y : {-0.5 * dy, 0.5 * dy}) {
			add_component(
				cell_center_x, center_y, cell_center_z,
				half_x, direction_y, half_z);
			add_component(
				cell_center_x, center_y, cell_center_z,
				half_x, direction_y, -half_z);
		}
		add_component(
			cell_center_x, center_y, cell_center_z,
			0.0, 0.5 * dy, 0.0);
	};

	const auto add_grid = [&](double origin_x,
		double origin_z,
		double cell_size_x,
		double cell_size_z,
		int column_count,
		int row_count) {
		for (int row = 0; row < row_count; ++row) {
			for (int column = 0; column < column_count; ++column) {
				add_three_dimensional_cell(
					origin_x + (static_cast<double>(column) + 0.5) * cell_size_x,
					origin_z + (static_cast<double>(row) + 0.5) * cell_size_z,
					cell_size_x,
					cell_size_z);
			}
		}
	};

	add_grid(
		bounds[0], bounds[4],
		arm_cell_size, arm_cell_size,
		cells_across_arm, cells_across_arm);
	add_grid(
		bounds[0], bounds[4] + arm_width,
		arm_cell_size, vertical_cell_size,
		cells_across_arm, extension_cell_count);
	add_grid(
		bounds[0] + arm_width, bounds[4],
		horizontal_cell_size, arm_cell_size,
		extension_cell_count, cells_across_arm);

	const double top_center_x = bounds[0] + 0.5 * arm_width;
	const double right_center_z = bounds[4] + 0.5 * arm_width;
	for (const double offset_y : {-layer_offset, layer_offset}) {
		add_component(
			top_center_x, center_y + offset_y, bounds[5],
			0.5 * arm_width, 0.0, 0.0);
		add_component(
			bounds[1], center_y + offset_y, right_center_z,
			0.0, 0.0, 0.5 * arm_width);
	}
	add_component(
		top_center_x, center_y, bounds[5],
		0.0, 0.5 * dy, 0.0);
	add_component(
		bounds[1], center_y, right_center_z,
		0.0, 0.5 * dy, 0.0);

	return components;
}

} // namespace mmc_layout
