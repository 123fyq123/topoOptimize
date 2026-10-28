#pragma once

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace optimization_paper_record {

inline double NaN();

template <typename T>
inline std::string ToString(const T& value);

struct BoundarySummary {
	bool readable = false;
	int point_count = 0;
	int dirichlet_count = 0;
	int neumann_count = 0;
	double dirichlet_centroid_x = NaN();
	double dirichlet_centroid_y = NaN();
	double dirichlet_centroid_z = NaN();
	double neumann_centroid_x = NaN();
	double neumann_centroid_y = NaN();
	double neumann_centroid_z = NaN();
	double resultant_x = 0.0;
	double resultant_y = 0.0;
	double resultant_z = 0.0;
};

struct RunMetadata {
	std::string schema_version = "2";
	std::string run_id;
	std::string started_at_utc;
	std::string method;
	std::string model_id;
	std::string model_signature;
	std::string boundary_file;
	std::string boundary_signature;
	std::string record_directory;
	BoundarySummary boundary;
};

inline double NaN() {
	return std::numeric_limits<double>::quiet_NaN();
}

inline std::string SanitizePathToken(std::string value) {
	for (char& character : value) {
		const unsigned char code = static_cast<unsigned char>(character);
		if (!std::isalnum(code) && character != '-' && character != '_') {
			character = '_';
		}
	}
	return value.empty() ? "unknown" : value;
}

inline std::string FileSignature(const std::string& file_path) {
	std::ifstream input(file_path, std::ios::binary);
	if (!input.good()) {
		return "unreadable";
	}
	// FNV-1a provides a stable experiment identity; it is not used for security.
	std::uint64_t signature_value = 14695981039346656037ull;
	char buffer[4096];
	while (input.read(buffer, sizeof(buffer)) || input.gcount() > 0) {
		for (std::streamsize index = 0; index < input.gcount(); ++index) {
			signature_value ^= static_cast<unsigned char>(buffer[index]);
			signature_value *= 1099511628211ull;
		}
	}
	std::ostringstream output;
	output << std::hex << std::setfill('0') << std::setw(16) << signature_value;
	return output.str();
}

inline BoundarySummary ReadBoundarySummary(const std::string& boundary_file) {
	BoundarySummary summary;
	std::ifstream input(boundary_file);
	int declared_count = 0;
	if (!(input >> declared_count) || declared_count < 0) {
		return summary;
	}
	double dirichlet_x = 0.0, dirichlet_y = 0.0, dirichlet_z = 0.0;
	double neumann_x = 0.0, neumann_y = 0.0, neumann_z = 0.0;
	for (int index = 0; index < declared_count; ++index) {
		double x = 0.0, y = 0.0, z = 0.0;
		double vx = 0.0, vy = 0.0, vz = 0.0;
		int type = 0;
		if (!(input >> x >> y >> z >> type >> vx >> vy >> vz)) {
			return BoundarySummary{};
		}
		if (type == 1) {
			dirichlet_x += x;
			dirichlet_y += y;
			dirichlet_z += z;
			++summary.dirichlet_count;
		}
		else if (type == 2) {
			neumann_x += x;
			neumann_y += y;
			neumann_z += z;
			summary.resultant_x += vx;
			summary.resultant_y += vy;
			summary.resultant_z += vz;
			++summary.neumann_count;
		}
	}
	summary.readable = true;
	summary.point_count = declared_count;
	if (summary.dirichlet_count > 0) {
		const double count = static_cast<double>(summary.dirichlet_count);
		summary.dirichlet_centroid_x = dirichlet_x / count;
		summary.dirichlet_centroid_y = dirichlet_y / count;
		summary.dirichlet_centroid_z = dirichlet_z / count;
	}
	if (summary.neumann_count > 0) {
		const double count = static_cast<double>(summary.neumann_count);
		summary.neumann_centroid_x = neumann_x / count;
		summary.neumann_centroid_y = neumann_y / count;
		summary.neumann_centroid_z = neumann_z / count;
	}
	return summary;
}

inline std::tm UtcTime(std::time_t value) {
	std::tm result{};
#ifdef _WIN32
	gmtime_s(&result, &value);
#else
	gmtime_r(&value, &result);
#endif
	return result;
}

inline RunMetadata CreateRunMetadata(
	const std::string& base_output_directory,
	const std::string& method,
	const std::string& input_model_file,
	const std::string& boundary_file) {
	const auto now = std::chrono::system_clock::now();
	const std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	const std::tm utc = UtcTime(now_time);
	const auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
		now.time_since_epoch()).count() % 1000;
	std::ostringstream timestamp;
	timestamp << std::put_time(&utc, "%Y%m%dT%H%M%S")
		<< std::setfill('0') << std::setw(3) << milliseconds << 'Z';
	std::ostringstream iso_time;
	iso_time << std::put_time(&utc, "%Y-%m-%dT%H:%M:%S")
		<< '.' << std::setfill('0') << std::setw(3) << milliseconds << 'Z';

	RunMetadata metadata;
	metadata.started_at_utc = iso_time.str();
	metadata.method = method;
	metadata.model_id = std::filesystem::path(input_model_file).stem().string();
	metadata.model_signature = FileSignature(input_model_file);
	metadata.boundary_file = boundary_file;
	metadata.boundary_signature = FileSignature(boundary_file);
	metadata.boundary = ReadBoundarySummary(boundary_file);
	const std::string boundary_token = metadata.boundary_signature.substr(
		0,
		std::min<std::size_t>(8, metadata.boundary_signature.size()));
	metadata.run_id = timestamp.str() + "_" + SanitizePathToken(method)
		+ "_" + SanitizePathToken(metadata.model_id) + "_bc" + boundary_token;

	const std::filesystem::path records_root =
		std::filesystem::path(base_output_directory) / "records";
	std::filesystem::create_directories(records_root);
	std::filesystem::path candidate = records_root / metadata.run_id;
	for (int suffix = 2; std::filesystem::exists(candidate); ++suffix) {
		candidate = records_root / (metadata.run_id + "_" + std::to_string(suffix));
	}
	metadata.record_directory = candidate.string();
	std::filesystem::create_directories(candidate);
	return metadata;
}

inline std::vector<std::pair<std::string, std::string>> MetadataEntries(
	const RunMetadata& metadata) {
	return {
		{ "record_schema_version", metadata.schema_version },
		{ "run_id", metadata.run_id },
		{ "started_at_utc", metadata.started_at_utc },
		{ "method", metadata.method },
		{ "model_id", metadata.model_id },
		{ "model_signature", metadata.model_signature },
		{ "boundary_file", metadata.boundary_file },
		{ "boundary_signature", metadata.boundary_signature },
		{ "boundary_point_count", ToString(metadata.boundary.point_count) },
		{ "dirichlet_point_count", ToString(metadata.boundary.dirichlet_count) },
		{ "dirichlet_centroid_x", ToString(metadata.boundary.dirichlet_centroid_x) },
		{ "dirichlet_centroid_y", ToString(metadata.boundary.dirichlet_centroid_y) },
		{ "dirichlet_centroid_z", ToString(metadata.boundary.dirichlet_centroid_z) },
		{ "neumann_point_count", ToString(metadata.boundary.neumann_count) },
		{ "neumann_centroid_x", ToString(metadata.boundary.neumann_centroid_x) },
		{ "neumann_centroid_y", ToString(metadata.boundary.neumann_centroid_y) },
		{ "neumann_centroid_z", ToString(metadata.boundary.neumann_centroid_z) },
		{ "load_resultant_x", ToString(metadata.boundary.resultant_x) },
		{ "load_resultant_y", ToString(metadata.boundary.resultant_y) },
		{ "load_resultant_z", ToString(metadata.boundary.resultant_z) },
		{ "record_directory", metadata.record_directory }
	};
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
