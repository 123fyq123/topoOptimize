#pragma once
#include"SimulationHead.h"
#include"SparseSystem.h"
class assembler {
	
protected:
	sparseSystem m_system;
public:

	virtual ~assembler()
	{}

	virtual void refresh();

	virtual void assemble();

	sparseSystem& system() { return m_system; }
};