#ifndef HEX_TOPLOPY_KERNEL_H
#define HEX_TOPLOPY_KERNEL_H
#include "ComHeader.h"

//规则的边界边度数
static const unsigned char REGULAR_EDGE_BDY_N = 3;
//规则的内部边度数
static const unsigned char REGULAR_EDGE_INNER_N = 4;
//
static const unsigned char REGULAR_VERTEX_INNER_N = 6;
static const unsigned char REGULAR_VERTEX_BDY_N = 5;

//没有意义 只是用来标示中间迭代过程中引入的奇异边
static const unsigned char UNREGULAR_EDGE_ADDIN_N = 10;
const uint16_t FACE_LINK_V_INDEX[6 * 4] = {
	0,4,5,1,
	2,6,7,3,
	1,5,6,2,
	3,7,4,0,
	4,7,6,5,
	1,0,3,2
};
const uint16_t EDGE_LINK_V_INDEX[12*2] = {
	0,1, 1,5, 5,4, 4,0,
	2,3, 3,7, 7,6, 6,2,
	1,2, 6,5, 3,0, 4,7
};
	
const uint16_t FACE_LINK_E_INDEX[6 * 4] = {
	0,1,2,3,
	4,5,6,7,
	8,7,9,1,
	10,3,11,5,
	2,9,6,11,
	0,10,4,8
};

class MeshHandle {
public:
	// Default constructor
	explicit MeshHandle(int _idx) : idx_(_idx) {};
		
	MeshHandle& operator=(int _idx) {
		idx_ = _idx;
		return *this;
	}

	MeshHandle& operator=(const MeshHandle& _idx) {
		idx_ = _idx.idx_;
		return *this;
	}

	inline bool is_valid() const { return idx_ != -1; }

	inline bool operator<(const MeshHandle& _idx) const { return (this->idx_ < _idx.idx_); }

	inline bool operator<(int _idx) const { return idx_ < _idx; }

	inline bool operator>(const MeshHandle& _idx) const { return (this->idx_ > _idx.idx_); }

	inline bool operator>(int _idx) const { return idx_ > _idx; }

	inline bool operator==(const MeshHandle& _h) const { return _h.idx_ == this->idx_; }

	inline bool operator!=(const MeshHandle& _h) const { return _h.idx_ != this->idx_; }

	inline const int& idx() const { return idx_; }

	void idx(const int& _idx) { idx_ = _idx; }

	inline operator int() const { return idx_; }//隐式转换

	void reset() { idx_ = -1; }

private:
	int idx_;
};

	
struct compare_OVM
{
	bool operator()(const MeshHandle &a, const MeshHandle &b) const
	{
		return a.idx()<b.idx();
	}
};
//其他重载

//bool operator==(const int& _lhs, const MeshHandle& _rhs) {
//	return _lhs == _rhs.idx();
//}
//
//bool operator==(const unsigned int& _lhs, const MeshHandle& _rhs) {
//
//	return _lhs == (unsigned int)_rhs.idx();
//}
//
//bool operator!=(const int& _lhs, const MeshHandle& _rhs) {
//
//	return !(_lhs == _rhs);
//}
//
//bool operator!=(const unsigned int& _lhs, const MeshHandle& _rhs) {
//
//	return !(_lhs == _rhs);
//}
//
//std::ostream& operator<<(std::ostream& _ostr, const MeshHandle& _handle) {
//	_ostr << _handle.idx();
//	return _ostr;
//}
//
//std::istream& operator >> (std::istream& _istr, MeshHandle& _handle) {
//	int val = 0;
//	_istr >> val;
//	_handle.idx(val);
//	return _istr;
//}


class VertexHandle : public MeshHandle { public: explicit VertexHandle(int _idx = -1) : MeshHandle(_idx) {} };
class EdgeHandle : public MeshHandle { public: explicit EdgeHandle(int _idx = -1) : MeshHandle(_idx) {} };
class FaceHandle : public MeshHandle { public: explicit FaceHandle(int _idx = -1) : MeshHandle(_idx) {} };
class CellHandle : public MeshHandle { public: explicit CellHandle(int _idx = -1) : MeshHandle(_idx) {} };
	

#endif // !HEX_TOPLOPY_KERNEL_H
