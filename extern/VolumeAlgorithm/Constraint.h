#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include "MeshKernel/Mesh.h"

using namespace std;

void read_constraint_data(string filename, unordered_map<int, bool>& res);

void write_constraint_data(string filename, unordered_map<int, bool>& res);

void updateAllHandles_kepp_constraint(MeshKernel::VolumeMesh& mesh, unordered_map<int, bool>& force_mp, unordered_map<int, bool>& spc_mp);
