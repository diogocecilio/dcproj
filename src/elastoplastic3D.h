#pragma once

//#include "elastmat2D.h"

#include "material.h"
#include "shapehexahedron.h"
#include "error.h"
using namespace std;

template <class YC>
class elastoplastic3D : public material{
public:
    elastoplastic3D(NRmatrix<Doub> bodyforce,Int order);
    elastoplastic3D( NRmatrix<Doub>  bodyforce, Int order, NRmatrix<Doub>  HHAT);
	elastoplastic3D(elastoplastic3D & copy);
	elastoplastic3D();
	~elastoplastic3D();
//
    void DirichletBC(SparseMatrix<double>  &KG, VectorXd &Fint,std::vector<int> ids, Int  dir, Int val);

	void Contribute(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, NRvector<Doub> intptsw, NRmatrix<Doub>  elcoords, NRmatrix<Doub>  eldisplace);

	void CalcStiff(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, const NRmatrix<Doub>   &elcoords, NRmatrix<Doub>  eldisplace);

	void assembleBandN(NRmatrix<Doub>  &B, NRmatrix<Doub>  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi);
	void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords);
	void DirichletBC(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val);
	void ContributeLineNewan(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val);
	void ContributeCurvedLine(NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology, Doub force);
    
    void SolPt(mesh * inmesh, const Int &el, const  NRmatrix<Doub>  &nodalsol, const Doub &xi, const Doub &eta, NRmatrix<Doub>  &xycoords, MatDoub &sol);
	void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly);
    void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, Int var,const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol);
    void PostProcessStrain(mesh * inmesh, NRvector<NRvector<NRtensor<Doub>>> &sol);
	void PostProcessIntegrationPointVar(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol);

	void SetRandomField(NRmatrix<Doub>  hhat);
	void SetRandomFieldLocal(NRvector<NRmatrix<Doub> >  hhatvel);
	void SetMemory(Int ngloblapoints, Int systemsize);
	void UpdateDisplacement(NRmatrix<Doub>  displace);
    void UpdateDisplacement(VectorXd displace);
	void UpdatePlasticStrain();
	void UpdateBodyForce(NRmatrix<Doub>  newbodyforce);
	void ResetPlasticStrain();
	void ResetDisplacement();
	void ResetMat();
	void ResetMemory();
	NRmatrix<Doub> GetSolution();
    void ComputeSolAndDSol(mesh * inmesh,NRmatrix<Doub>&sol,NRmatrix<Doub>&dsol);
    void ComputeSolAndDSol(mesh * inmesh,NRvector<NRmatrix<Doub>>&sol,NRvector<NRmatrix<Doub>>&dsol);
    void ComputeSolution(mesh * inmesh,NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace,NRvector<Doub> ptsw,NRmatrix<Doub> &sol,NRmatrix<Doub> &dsol);
	void ResetCounter()
	{
		fglobalcounter = 0;
	}


	void SetMatConstants(NRvector<Doub>& consts)
	{
		fYC.SetMatConstants(consts);
	}

	void GetMatConstants(NRvector<Doub>& consts)
	{
		fYC.GetMatConstants(consts);
	}
	

    void SetTangentMatrixType(bool type)
    {
         fYC.SetTangentMatrixType(type);
    }
     
	virtual void GetMaterialData()
	{

	}

	NRmatrix<Doub>  GetBodyForce()
	{
		return fbodyforce;
	}
	virtual NRvector<NRtensor<Doub> > GetPlasticStrain();

    NRvector<Doub> ComputePhi(NRtensor<Doub> eps)
    {
        NRvector<Doub> valphi = fYC.phi(eps);
        return valphi;
    }
//
//
private:

	NRmatrix<Doub>  fbodyforce;
	Int fOrder;
	//mesh fmesh;



public:
	NRmatrix<Doub>  fHHAT;
	NRvector<NRmatrix<Doub> > fhhatvel;
	YC fYC;

	NRmatrix<Doub>  fdisplace;
	NRvector<NRtensor<Doub> > fepspvec, fepspsolitern;
	Int fglobalcounter;

 


};



