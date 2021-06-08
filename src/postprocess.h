#pragma once
#include "mesh.h"
#include "nr3.h"
#include "shapequad.h"
class postprocess
{
public:
	postprocess();
	~postprocess();

	void PostProcess(mesh &inmesh, const MatDoub & nodalsol, std::vector<std::vector<double>> &sol);
};

