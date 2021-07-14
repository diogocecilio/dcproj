#pragma once
//#include "nr3.h"
#include <math.h>
#include <cmath>
//#include "elastoplasticbase.h"
//#include "shapequad.h"
#include <stdio.h>
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include "material.h"
#include "shapequad.h"

class mesh
{
public:

	mesh(std::vector<std::vector< std::vector<Doub > > > &allcoords, NRmatrix<Doub> &meshnodes, NRmatrix<Int> &meshtopology);

	mesh(material * mat , std::vector<std::vector< std::vector<Doub > > >& allcoords, NRmatrix<Doub>& meshnodes, NRmatrix<Int>& meshtopology);

	mesh(material* mat, std::vector<std::vector< std::vector<Doub > > >& allcoords, NRmatrix<Doub>& meshnodes, NRmatrix<Int>& meshtopology,MatDoub  HHAT );

	mesh(mesh &copy);

	mesh();

	~mesh();

	static void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords);
	MatDoub FindSolution(VecDoub coord, MatDoub datatosearch);
	MatDoub  TransferSolution( mesh &out, MatDoub datatosearch);
	std::vector<std::vector< std::vector<Doub > > >  GetAllCoords();
	MatDoub GetMeshNodes();
	MatInt GetMeshTopology();
	void Assemble(MatDoub& KG, MatDoub& Fint, MatDoub& Fbody);
    void Assemble(SparseMatrix<double>  &KG, VectorXd &Fint, VectorXd &Fbody);

	inline void SetHhat(NRmatrix<Doub> HHAT)
	{
		fHHAT = HHAT;
	}

	inline void GetHhat(NRmatrix<Doub> &HHAT)
	{
		 HHAT=fHHAT;
	}
	//material GetMaterial();


	material* fmaterial;
private:

	 
	 std::vector<std::vector< std::vector<Doub > > > fallcoords;
	 MatDoub fmeshnodes;
	 MatInt fmeshtopology;
	 NRmatrix<Doub> fHHAT;
	 NRvector<MatDoub> fhhatvel;
	
};

// Generic functor

struct FuncdSearch {
	VecDoub co;
	MatDoub psis, GradPsi, elcoords, psist,xycoords;
	shapequad shape = shapequad(2, 1);
	Doub operator() (VecDoub_I &x)
	{
		Doub xi = x[0];
		Doub eta = x[1];
		shape.shapes(psis, GradPsi, xi, eta);
		psis.Transpose(psist);
		psist.Mult(elcoords, xycoords);
		Doub dx = xycoords[0][0] - co[0];
		Doub dy = xycoords[0][1] - co[1];
		Doub dist = sqrt(dx*dx + dy*dy);
		return dist;
	}

	void set(VecDoub tosearch,MatDoub elcoordsout)
	{
		co = tosearch;
		elcoords = elcoordsout;
	}
};



