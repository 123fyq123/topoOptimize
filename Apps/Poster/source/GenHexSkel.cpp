#include "GenHexSkel.h"
using namespace std;

void GenHexSkel(std::string InputFile) {
	std::string filename(InputFile);
	BMesh::Unitest_common::main_genNodeHex(filename);
}