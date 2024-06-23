#pragma once
#include"head.h"
#include"SparseSystem.h"
class assembler {
	
protected:
	SparseSystem m_system;
public:

	virtual void refresh();

	virtual void assemble();

	SparseSystem& system() { return m_system; }
};