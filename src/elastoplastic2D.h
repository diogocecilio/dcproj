#pragma once

//#include "elastmat2D.h"

#include "material.h"
#include "shapequad.h"

using namespace std;

template <class YC>
class elastoplastic2D : public material{
public:
	elastoplastic2D(Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int order, NRmatrix<Doub>   HHAT);
	elastoplastic2D(Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int order);
	elastoplastic2D(Doub young, Doub nu, Doub sigy, Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int order, NRmatrix<Doub>   HHAT);

	elastoplastic2D(elastoplastic2D & copy);
	elastoplastic2D();
	~elastoplastic2D();

//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order);
//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT);
//	elastmat2D();
//	~elastmat2D();
//
	void Contribute(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody,NRvector<Doub> ptsw, NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace);

    void ContributeEig(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody,NRvector<Doub> ptsw, NRmatrix<Doub>  elcoords, NRmatrix<Doub>  eldisplace);

	void CacStiff(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, const NRmatrix<Doub>   &elcoords, NRmatrix<Doub>  eldisplace);
	void Assemble(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, NRmatrix<Doub>  &KG, NRmatrix<Doub>  &Fint,NRmatrix<Doub>  &Fbody);

    void Assemble(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, SparseMatrix<double>  &KG, VectorXd &Fint, VectorXd &Fbody);

	//void Assemble(MatDoub &KG, MatDoub &Fint, MatDoub &Fbody);
	void assembleBandN(NRmatrix<Doub>  &B, NRmatrix<Doub>  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi);
    //	B.assign(3, psis.nrows() * 2, 0.);
	//N.assign(2, psis.nrows() * 2, 0.);
    void assembleBandN(MatrixXd&B, MatrixXd  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi);
	void assembleConstitutiveMatrix(MatDoub &C, Doub mult);
	void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords);
	void DirichletBC(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val);
    void DirichletBC(SparseMatrix<double> & KG, VectorXd& FG, std::vector<int> ids, Int  dir, Int val);
	void ContributeLineNewan(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val);
	void ContributeCurvedLine(NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology, Doub force);
	void SolPt(const std::vector<std::vector< std::vector<Doub > > > &allcoords,const MatInt &meshtopology, const Int &el, const  NRmatrix<Doub>  &solG, const Doub &xi, const Doub &eta, NRmatrix<Doub>  &xycoords, MatDoub &sol);
//
	void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly);
    void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, Int var,const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol);

	void PostProcessIntegrationPointVar(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol);

	//void SetMemory(MatDoub displace, NRvector<TensorDoub> epspvec, NRvector<TensorDoub>  epspsolitern, Int globalcounter);
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
    NRvector<NRtensor<Doub> > GetPlasticStrain();

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

	virtual void GetMaterialData()
	{

	}

	NRmatrix<Doub>  GetBodyForce()
	{
		return fbodyforce;
	}
//
//
private:

	NRmatrix<Doub>  fbodyforce;
	Int fplanestress;
	Doub fthickness;
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



