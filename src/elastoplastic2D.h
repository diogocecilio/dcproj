#pragma once

//#include "elastmat2D.h"

#include "material.h"
#include "shapequad.h"

using namespace std;

template <class YC>
class elastoplastic2D : public material{
public:
	elastoplastic2D(Doub thickness, MatDoub bodyforce, Int planestress, Int order, MatDoub  HHAT);
	elastoplastic2D(Doub thickness, MatDoub bodyforce, Int planestress, Int order);
	elastoplastic2D(Doub young, Doub nu, Doub sigy, Doub thickness, MatDoub bodyforce, Int planestress, Int order, MatDoub  HHAT);

	elastoplastic2D(elastoplastic2D & copy);
	elastoplastic2D();
	~elastoplastic2D();

//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order);
//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT);
//	elastmat2D();
//	~elastmat2D();
//
	void Contribute(MatDoub &ek, MatDoub &efint, MatDoub &efbody, Doub xi, Doub eta, Doub w, MatDoub elcoords, MatDoub eldisplace);
	void CacStiff(MatDoub &ek, MatDoub &efint, MatDoub &efbody, const MatDoub  &elcoords, MatDoub eldisplace);
	void Assemble(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub &KG, MatDoub &Fint,MatDoub &Fbody);
	//void Assemble(MatDoub &KG, MatDoub &Fint, MatDoub &Fbody);
	void assembleBandN(MatDoub &B, MatDoub &N, const MatDoub &psis, const MatDoub &GradPhi);
	void assembleConstitutiveMatrix(MatDoub &C, Doub mult);
	void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords);
	void DirichletBC(MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val);
	void ContributeLineNewan(MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val);
	void ContributeCurvedLine(MatDoub &KG, MatDoub &FG, MatDoub meshnodes, MatInt linetopology, Doub force);
	void SolPt(const std::vector<std::vector< std::vector<Doub > > > &allcoords,const MatInt &meshtopology, const Int &el, const  MatDoub &solG, const Doub &xi, const Doub &eta, MatDoub &xycoords, MatDoub &sol);
//
	void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, const MatDoub & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly);
    void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, Int var,const MatDoub & nodalsol, std::vector<std::vector<double>> &sol);

	void PostProcessIntegrationPointVar(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, const MatDoub & nodalsol, std::vector<std::vector<double>> &sol);

	//void SetMemory(MatDoub displace, NRvector<TensorDoub> epspvec, NRvector<TensorDoub>  epspsolitern, Int globalcounter);
	void SetRandomField(MatDoub hhat);
	void SetRandomFieldLocal(NRvector<MatDoub>  hhatvel);
	void SetMemory(Int ngloblapoints, Int systemsize);
	void UpdateDisplacement(MatDoub displace);
	void UpdatePlasticStrain();
	void UpdateBodyForce(MatDoub newbodyforce);
	void ResetPlasticStrain();
	void ResetDisplacement();
	void ResetMat();
	void ResetMemory();
	NRmatrix<Doub> GetSolution();

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

	MatDoub GetBodyForce()
	{
		return fbodyforce;
	}
//
//
private:

	MatDoub fbodyforce;
	Int fplanestress;
	Doub fthickness;
	Int fOrder;
	//mesh fmesh;



public:
	MatDoub fHHAT;
	NRvector<MatDoub> fhhatvel;
	YC fYC;

	MatDoub fdisplace;
	NRvector<TensorDoub> fepspvec, fepspsolitern;
	Int fglobalcounter;

 


};



