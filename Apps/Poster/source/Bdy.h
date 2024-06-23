#pragma once
#include <vector>
#include <string>
namespace IGA {
	class BdyBase {
	public:
		enum BdyType {
			Dirichlet = 1,
			Neumann,
			Robin
		};
		typedef std::vector<double> BdyValue;
		BdyType bdy_type;
		BdyValue bdy_value;
	};
	struct HexBdyBase : public BdyBase {
		enum ElementType {
			Point = 1,
			Line,
			Face,
			Cell
		};

		typedef std::vector<double> VertexPos;
		typedef std::vector<VertexPos> ElementPoints;
		ElementPoints element_points;
		ElementType element_type;

	};
	typedef std::vector<HexBdyBase> HexBdy;

	static std::vector<int> protected_cell;

	struct CCBezierBdyBase :public BdyBase {
		size_t cell_id;
		size_t local_base_id;
		size_t global_base_id;
	};
	typedef std::vector<CCBezierBdyBase> CCBezierBdy;

	void BdyConvert(const std::string& filename, HexBdy& hexbdy, CCBezierBdy& ccbdy);

	void WriteBdy(const std::string& filename, HexBdy& hexbdy, std::vector<int>& pro);
	//void WriteBdy(const std::string& filename, CCBezierBdy& ccbdy);
	void ReadBdy(const std::string& filename, HexBdy& hexbdy, std::vector<int>& pro);
	//void ReadBdy(const std::string& filename, CCBezierBdy& ccbdy);
}