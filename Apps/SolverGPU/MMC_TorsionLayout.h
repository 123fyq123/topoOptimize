#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace mmc_layout {

// Shared by the pure-MMC and simultaneous solvers, whose component types
// expose the same nine-argument constructor.
template <typename TComponent>
std::vector<TComponent> MakeTorsionBeamLayout(const std::vector<double>& bounds) {
	if (bounds.size() != 6) return {};
	const double dx = bounds[1] - bounds[0];
	const double dy = bounds[3] - bounds[2];
	const double dz = bounds[5] - bounds[4];
	if (!(dx > 0 && dy > 0 && dz > 0)) return {};

	// Three X sections and two Y rows. For 48 x 16 x 16, the half spans
	// are (8, 4, 4) and the initial cross-section half width is 1.2.
	const double hx = dx / 6, hy = dy / 4, hz = dz / 4;
	const double mid_z = bounds[4] + dz / 2;
	const double radius = 0.075 * std::min(dx, std::min(dy, dz));
	std::vector<TComponent> components;
	components.reserve(56);
	const auto add = [&](double x, double y, double z, double vx, double vy, double vz) {
		// v points from the center to one end of the component.
		const double planar_length = std::sqrt(vx * vx + vy * vy);
		const double half_length = std::sqrt(planar_length * planar_length + vz * vz);
		components.emplace_back(x, y, z, half_length, radius, radius,
			0.0, -std::atan2(vz, planar_length), std::atan2(vy, vx));
	};
	const double layer_offsets[2] = {-hz, hz};

	for (int ix = 0; ix < 3; ++ix) {
		const double x = bounds[0] + (2 * ix + 1) * hx;
		for (int iy = 0; iy < 2; ++iy) {
			const double y = bounds[2] + (2 * iy + 1) * hy;
			for (double offset : layer_offsets) {
				add(x, y, mid_z + offset, hx, 0, 0);
				add(x, y, mid_z + offset, 0, hy, 0);
				add(x, y, mid_z, hx, hy, offset);
				add(x, y, mid_z, hx, -hy, offset);
			}
			add(x, y, mid_z, 0, 0, hz);
		}
	}
	// Two symmetric Z-edge collectors connect all four loaded corners.
	const double load_edge_y[2] = {bounds[2], bounds[3]};
	for (double y : load_edge_y) {
		add(bounds[1], y, mid_z, 0, 0, dz / 2);
	}
	return components;
}

} // namespace mmc_layout
