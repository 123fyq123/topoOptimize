#pragma once

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace optimization_paper_record {

inline double NaN() {
	return std::numeric_limits<double>::quiet_NaN();
}

inline double ElapsedSeconds(const std::chrono::steady_clock::time_point& begin) {
	return std::chrono::duration<double>(
		std::chrono::steady_clock::now() - begin).count();
}

inline double L2Norm(const std::vector<double>& values) {
	double squared_norm = 0.0;
	for (double value : values) {
		squared_norm += value * value;
	}
	return std::sqrt(squared_norm);
}

struct IterationRecord {
	std::string method;
	std::string stage;
	int global_iteration = 0;
	int stage_iteration = 0;
	double compliance = NaN();
	double shape_compliance = NaN();
	double mmc_volume = NaN();
	double material_v1 = NaN();
	double material_v2 = NaN();
	double material_ratio = NaN();
	double ratio_residual = NaN();
	double mean_material_variable = NaN();
	double material_grayness = NaN();
	double geometry_change = NaN();
	double material_change = NaN();
	double geometry_gradient_norm = NaN();
	double material_gradient_norm = NaN();
	int active_cell_count = 0;
	int cumulative_fe_solves = 0;
	double iteration_seconds = NaN();
	double elapsed_seconds = NaN();
};

class IterationWriter {
public:
	IterationWriter() = default;

	IterationWriter(const std::string& file_path, bool append) {
		Open(file_path, append);
	}

	bool Open(const std::string& file_path, bool append) {
		if (file_path.empty()) {
			return false;
		}
		const std::filesystem::path path(file_path);
		if (path.has_parent_path()) {
			std::filesystem::create_directories(path.parent_path());
		}
		const bool write_header = !append || !std::filesystem::exists(path)
			|| std::filesystem::file_size(path) == 0;
		output_.open(
			path,
			std::ios::out | (append ? std::ios::app : std::ios::trunc));
		if (!output_.good()) {
			return false;
		}
		output_ << std::setprecision(17);
		if (write_header) {
			output_
				<< "method\tstage\tglobal_iteration\tstage_iteration"
				<< "\tcompliance\tshape_compliance\tmmc_volume"
				<< "\tmaterial_v1\tmaterial_v2\tmaterial_ratio\tratio_residual"
				<< "\tmean_material_variable\tmaterial_grayness"
				<< "\tgeometry_change\tmaterial_change"
				<< "\tgeometry_gradient_norm\tmaterial_gradient_norm"
				<< "\tactive_cell_count\tcumulative_fe_solves"
				<< "\titeration_seconds\telapsed_seconds\n";
		}
		return true;
	}

	bool Good() const {
		return output_.good();
	}

	void Write(const IterationRecord& record) {
		if (!output_.is_open()) {
			return;
		}
		output_
			<< record.method << '\t'
			<< record.stage << '\t'
			<< record.global_iteration << '\t'
			<< record.stage_iteration << '\t'
			<< record.compliance << '\t'
			<< record.shape_compliance << '\t'
			<< record.mmc_volume << '\t'
			<< record.material_v1 << '\t'
			<< record.material_v2 << '\t'
			<< record.material_ratio << '\t'
			<< record.ratio_residual << '\t'
			<< record.mean_material_variable << '\t'
			<< record.material_grayness << '\t'
			<< record.geometry_change << '\t'
			<< record.material_change << '\t'
			<< record.geometry_gradient_norm << '\t'
			<< record.material_gradient_norm << '\t'
			<< record.active_cell_count << '\t'
			<< record.cumulative_fe_solves << '\t'
			<< record.iteration_seconds << '\t'
			<< record.elapsed_seconds << '\n';
		output_.flush();
	}

private:
	std::ofstream output_;
};

inline bool WriteConfig(
	const std::string& file_path,
	const std::vector<std::pair<std::string, std::string>>& entries) {
	const std::filesystem::path path(file_path);
	if (path.has_parent_path()) {
		std::filesystem::create_directories(path.parent_path());
	}
	std::ofstream output(path, std::ios::out | std::ios::trunc);
	if (!output.good()) {
		return false;
	}
	output << "parameter\tvalue\n";
	for (const auto& entry : entries) {
		output << entry.first << '\t' << entry.second << '\n';
	}
	return output.good();
}

template <typename T>
inline std::string ToString(const T& value) {
	std::ostringstream output;
	output << std::setprecision(17) << value;
	return output.str();
}

} // namespace optimization_paper_record
