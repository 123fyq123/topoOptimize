#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace mmc_trace {

struct VectorSummary {
	double minimum = std::numeric_limits<double>::quiet_NaN();
	double maximum = std::numeric_limits<double>::quiet_NaN();
	double mean = std::numeric_limits<double>::quiet_NaN();
	double l2_norm = std::numeric_limits<double>::quiet_NaN();
};

inline VectorSummary Summarize(const std::vector<double>* values) {
	VectorSummary summary;
	if (values == nullptr || values->empty()) return summary;
	summary.minimum = (*values)[0];
	summary.maximum = (*values)[0];
	double sum = 0.0;
	double squared_sum = 0.0;
	for (double value : *values) {
		summary.minimum = std::min(summary.minimum, value);
		summary.maximum = std::max(summary.maximum, value);
		sum += value;
		squared_sum += value * value;
	}
	summary.mean = sum / static_cast<double>(values->size());
	summary.l2_norm = std::sqrt(squared_sum);
	return summary;
}

inline double Dot(
	const std::vector<double>& first,
	const std::vector<double>& second) {
	const size_t count = std::min(first.size(), second.size());
	double value = 0.0;
	for (size_t index = 0; index < count; ++index) {
		value += first[index] * second[index];
	}
	return value;
}

inline double TangentGradientNorm(
	const std::vector<double>* objective_gradient,
	const std::vector<double>* equality_gradient) {
	if (objective_gradient == nullptr || equality_gradient == nullptr ||
		objective_gradient->size() != equality_gradient->size() ||
		objective_gradient->empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	const double denominator = Dot(*equality_gradient, *equality_gradient);
	if (denominator <= 1e-30) {
		return Summarize(objective_gradient).l2_norm;
	}
	const double scale = Dot(*objective_gradient, *equality_gradient) / denominator;
	double squared_norm = 0.0;
	for (size_t index = 0; index < objective_gradient->size(); ++index) {
		const double projected =
			(*objective_gradient)[index] - scale * (*equality_gradient)[index];
		squared_norm += projected * projected;
	}
	return std::sqrt(squared_norm);
}

inline std::string SampleValues(
	const std::vector<double>* values,
	int requested_count = 8) {
	if (values == nullptr || values->empty() || requested_count <= 0) return "";
	const int value_count = static_cast<int>(values->size());
	const int sample_count = std::min(value_count, requested_count);
	std::ostringstream samples;
	samples << std::setprecision(17);
	int previous_index = -1;
	for (int sample = 0; sample < sample_count; ++sample) {
		const int index = sample_count == 1
			? 0
			: static_cast<int>(std::llround(
				static_cast<double>(sample) * (value_count - 1) / (sample_count - 1)));
		if (index == previous_index) continue;
		if (previous_index >= 0) samples << ';';
		samples << index << ':' << (*values)[index];
		previous_index = index;
	}
	return samples.str();
}

struct Record {
	std::string mode;
	std::string action;
	int iteration = 0;
	int component_count = 0;
	int vertex_count = 0;
	int cell_count = 0;
	double objective = std::numeric_limits<double>::quiet_NaN();
	double volume = std::numeric_limits<double>::quiet_NaN();
	double volume_constraint = std::numeric_limits<double>::quiet_NaN();
	double material_v1 = std::numeric_limits<double>::quiet_NaN();
	double material_v2 = std::numeric_limits<double>::quiet_NaN();
	double ratio_residual = std::numeric_limits<double>::quiet_NaN();
	double positive_ratio_constraint = std::numeric_limits<double>::quiet_NaN();
	double negative_ratio_constraint = std::numeric_limits<double>::quiet_NaN();
	double active_constraint_residual = std::numeric_limits<double>::quiet_NaN();
	double geometry_change_max = std::numeric_limits<double>::quiet_NaN();
	double material_change_max = std::numeric_limits<double>::quiet_NaN();
	double geometry_change_tolerance = std::numeric_limits<double>::quiet_NaN();
	double material_change_tolerance = std::numeric_limits<double>::quiet_NaN();
	double constraint_tolerance = std::numeric_limits<double>::quiet_NaN();
	int geometry_change_ok = -1;
	int material_change_ok = -1;
	int constraint_ok = -1;
	int convergence_triggered = 0;
	double previous_objective_linear_geometry = std::numeric_limits<double>::quiet_NaN();
	double previous_objective_linear_material = std::numeric_limits<double>::quiet_NaN();
	double previous_geometry_constraint_linear_geometry = std::numeric_limits<double>::quiet_NaN();
	double previous_ratio_constraint_linear_geometry = std::numeric_limits<double>::quiet_NaN();
	double previous_ratio_constraint_linear_material = std::numeric_limits<double>::quiet_NaN();
	const std::vector<double>* global_phi = nullptr;
	const std::vector<double>* heaviside = nullptr;
	const std::vector<double>* heaviside_derivative = nullptr;
	const std::vector<double>* occupancy = nullptr;
	const std::vector<double>* physical_variables = nullptr;
	const std::vector<double>* mmc_objective_gradient = nullptr;
	const std::vector<double>* geometry_constraint_gradient = nullptr;
	const std::vector<double>* material_eta = nullptr;
	const std::vector<double>* material_modulus = nullptr;
	const std::vector<double>* material_objective_gradient = nullptr;
	const std::vector<double>* ratio_material_gradient = nullptr;
	const std::vector<double>* geometry_step = nullptr;
	const std::vector<double>* material_step = nullptr;
};

class Writer {
public:
	explicit Writer(const std::string& file_path) {
		const std::filesystem::path path(file_path);
		if (!path.parent_path().empty()) {
			std::error_code error;
			std::filesystem::create_directories(path.parent_path(), error);
			if (error) return;
		}
		output_.open(file_path, std::ios::out | std::ios::trunc);
		if (output_.good()) {
			output_ << std::setprecision(17);
			WriteHeader();
		}
	}

	bool Good() const { return output_.good(); }

	void Write(const Record& record) {
		if (!output_.good()) return;
		const VectorSummary global_phi = Summarize(record.global_phi);
		const VectorSummary heaviside = Summarize(record.heaviside);
		const VectorSummary heaviside_derivative = Summarize(record.heaviside_derivative);
		const VectorSummary occupancy = Summarize(record.occupancy);
		const VectorSummary physical_variables = Summarize(record.physical_variables);
		const VectorSummary mmc_gradient = Summarize(record.mmc_objective_gradient);
		const VectorSummary geometry_constraint_gradient =
			Summarize(record.geometry_constraint_gradient);
		const VectorSummary eta = Summarize(record.material_eta);
		const VectorSummary modulus = Summarize(record.material_modulus);
		const VectorSummary material_gradient =
			Summarize(record.material_objective_gradient);
		const VectorSummary ratio_material_gradient =
			Summarize(record.ratio_material_gradient);
		const VectorSummary geometry_step = Summarize(record.geometry_step);
		const VectorSummary material_step = Summarize(record.material_step);

		output_ << record.mode << '\t' << record.iteration << '\t' << record.action
			<< '\t' << record.component_count << '\t' << record.vertex_count
			<< '\t' << record.cell_count << '\t' << record.objective
			<< '\t' << record.volume << '\t' << record.volume_constraint
			<< '\t' << record.material_v1 << '\t' << record.material_v2
			<< '\t' << record.ratio_residual
			<< '\t' << record.positive_ratio_constraint
			<< '\t' << record.negative_ratio_constraint
			<< '\t' << record.active_constraint_residual
			<< '\t' << record.geometry_change_max << '\t' << record.material_change_max
			<< '\t' << record.geometry_change_tolerance
			<< '\t' << record.material_change_tolerance
			<< '\t' << record.constraint_tolerance
			<< '\t' << record.geometry_change_ok << '\t' << record.material_change_ok
			<< '\t' << record.constraint_ok << '\t' << record.convergence_triggered
			<< '\t' << record.previous_objective_linear_geometry
			<< '\t' << record.previous_objective_linear_material
			<< '\t' << record.previous_geometry_constraint_linear_geometry
			<< '\t' << record.previous_ratio_constraint_linear_geometry
			<< '\t' << record.previous_ratio_constraint_linear_material;
		WriteSummary(global_phi);
		WriteSummary(heaviside);
		WriteSummary(heaviside_derivative);
		WriteSummary(occupancy);
		WriteSummary(physical_variables);
		WriteSummary(mmc_gradient);
		WriteSummary(geometry_constraint_gradient);
		WriteSummary(eta);
		WriteSummary(modulus);
		WriteSummary(material_gradient);
		WriteSummary(ratio_material_gradient);
		WriteSummary(geometry_step);
		WriteSummary(material_step);
		output_ << '\t' << TangentGradientNorm(
			record.material_objective_gradient,
			record.ratio_material_gradient)
			<< '\t' << SampleValues(record.physical_variables)
			<< '\t' << SampleValues(record.global_phi)
			<< '\t' << SampleValues(record.heaviside)
			<< '\t' << SampleValues(record.heaviside_derivative)
			<< '\t' << SampleValues(record.occupancy)
			<< '\t' << SampleValues(record.mmc_objective_gradient)
			<< '\t' << SampleValues(record.geometry_constraint_gradient)
			<< '\t' << SampleValues(record.material_eta)
			<< '\t' << SampleValues(record.material_objective_gradient)
			<< '\n';
		output_.flush();
	}

private:
	void WriteHeader() {
		output_
			<< "mode\titeration\taction\tcomponent_count\tvertex_count\tcell_count"
			<< "\tobjective\tvolume\tvolume_constraint\tmaterial_v1\tmaterial_v2"
			<< "\tratio_residual\tpositive_ratio_constraint\tnegative_ratio_constraint"
			<< "\tactive_constraint_residual\tgeometry_change_max\tmaterial_change_max"
			<< "\tgeometry_change_tolerance\tmaterial_change_tolerance\tconstraint_tolerance"
			<< "\tgeometry_change_ok\tmaterial_change_ok\tconstraint_ok\tconvergence_triggered"
			<< "\tprevious_objective_linear_geometry\tprevious_objective_linear_material"
			<< "\tprevious_geometry_constraint_linear_geometry"
			<< "\tprevious_ratio_constraint_linear_geometry"
			<< "\tprevious_ratio_constraint_linear_material";
		WriteSummaryHeader("global_phi");
		WriteSummaryHeader("heaviside");
		WriteSummaryHeader("heaviside_derivative");
		WriteSummaryHeader("occupancy");
		WriteSummaryHeader("physical_variables");
		WriteSummaryHeader("mmc_objective_gradient");
		WriteSummaryHeader("geometry_constraint_gradient");
		WriteSummaryHeader("material_eta");
		WriteSummaryHeader("material_modulus");
		WriteSummaryHeader("material_objective_gradient");
		WriteSummaryHeader("ratio_material_gradient");
		WriteSummaryHeader("geometry_step");
		WriteSummaryHeader("material_step");
		output_
			<< "\tmaterial_tangent_gradient_l2"
			<< "\tphysical_variables_samples\tglobal_phi_samples\theaviside_samples"
			<< "\theaviside_derivative_samples\toccupancy_samples"
			<< "\tmmc_objective_gradient_samples\tgeometry_constraint_gradient_samples"
			<< "\tmaterial_eta_samples\tmaterial_objective_gradient_samples\n";
	}

	void WriteSummaryHeader(const char* prefix) {
		output_ << '\t' << prefix << "_min\t" << prefix << "_max\t"
			<< prefix << "_mean\t" << prefix << "_l2";
	}

	void WriteSummary(const VectorSummary& summary) {
		output_ << '\t' << summary.minimum << '\t' << summary.maximum
			<< '\t' << summary.mean << '\t' << summary.l2_norm;
	}

	std::ofstream output_;
};

} // namespace mmc_trace
