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
    //void PostProcess(mesh * inmesh, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly);
        /**
	 * @brief Computes the solution in a specific element and xi,eta coord
     * @param [in] allcoords 
     * @param [in] meshtopology   
     * @param [in] el element id
     * @param [in] nodalsol mesh solution
	 * @param [in] xi   integration point coord
     * @param [in] eta integration point coord
     * @param [out] sol solution in the specific xi,eta coord
	 */
	//void SolPt(const std::vector<std::vector< std::vector<Doub > > > &allcoords,const MatInt &meshtopology, const Int &el, const  NRmatrix<Doub>  &nodalsol, const Doub //&xi, const Doub &eta, NRmatrix<Doub>  &xycoords, MatDoub &sol);
    void SolPt(mesh * inmesh, const Int &el, const  NRmatrix<Doub>  &nodalsol, const Doub &xi, const Doub &eta, NRmatrix<Doub>  &xycoords, MatDoub &sol);
};

