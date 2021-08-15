#pragma once
#include "nr3.h"
#include "vonmises.h"
#include "druckerprager.h"
#include "elastoplasticbase.h"
//#include "mesh.h"


class material
{
public:
	material();
	~material();

	//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order);
	//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT);
	//	elastmat2D();
	//	~elastmat2D();
	//
	virtual void Contribute(MatDoub& ek, MatDoub& efint, MatDoub& efbody,NRvector<Doub> ptsw, MatDoub elcoords, MatDoub eldisplace)=0;
	virtual void CacStiff(MatDoub& ek, MatDoub& efint, MatDoub& efbody, const MatDoub& elcoords, MatDoub eldisplace)=0;
	virtual void Assemble(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub& KG, MatDoub& Fint, MatDoub& Fbody)=0;
	//virtual void Assemble(MatDoub &KG, MatDoub &Fint, MatDoub &Fbody)=0;
	virtual void assembleBandN(MatDoub& B, MatDoub& N, const MatDoub& psis, const MatDoub& GradPhi)=0;
	//virtual void assembleConstitutiveMatrix(MatDoub& C, Doub mult)=0;
	virtual void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub& elcoords)=0;
	virtual void DirichletBC(MatDoub& KG, MatDoub& FG, std::vector<int> ids, Int  dir, Int val)=0;
    virtual void DirichletBC(SparseMatrix<double>  &KG, VectorXd &Fint,std::vector<int> ids, Int  dir, Int val)=0;
	virtual void ContributeLineNewan(MatDoub& KG, MatDoub& FG, std::vector<int> ids, Int  dir, Int val)=0;
	virtual void ContributeCurvedLine(MatDoub& KG, MatDoub& FG, MatDoub meshnodes, MatInt linetopology, Doub force)=0;
	virtual void SolPt(const std::vector<std::vector< std::vector<Doub > > >& allcoords, const MatInt& meshtopology, const Int& el, const  MatDoub& solG, const Doub& xi, const Doub& eta, MatDoub& xycoords, MatDoub& sol)=0;
	//
	virtual void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, const MatDoub& nodalsol, std::vector<std::vector<double>>& solx, std::vector<std::vector<double>>& soly)=0;
	virtual void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, Int var, const MatDoub& nodalsol, std::vector<std::vector<double>>& sol)=0;

	virtual void PostProcessIntegrationPointVar(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, const MatDoub& nodalsol, std::vector<std::vector<double>>& sol)=0;

	//void SetMemory(MatDoub displace, NRvector<TensorDoub> epspvec, NRvector<TensorDoub>  epspsolitern, Int globalcounter);
	virtual void SetRandomField(MatDoub hhat)=0;
	virtual void SetMemory(Int ngloblapoints, Int systemsize)=0;
	virtual void UpdateDisplacement(MatDoub displace)=0;
    virtual void UpdateDisplacement(VectorXd displace)=0;
	virtual void UpdatePlasticStrain()=0;
	virtual void UpdateBodyForce(MatDoub newbodyforce)=0;
	virtual void ResetPlasticStrain()=0;
	virtual void ResetDisplacement()=0;
	virtual void ResetMat()=0;
	virtual void ResetMemory()=0;
	virtual NRmatrix<Doub> GetSolution() = 0;
    virtual NRvector<NRtensor<Doub> > GetPlasticStrain()=0;
	virtual void SetMatConstants(NRvector<Doub>& consts)=0;
	virtual void GetMatConstants(NRvector<Doub>& consts) = 0;
	virtual void ResetCounter() = 0;
	virtual void SetRandomFieldLocal(NRvector<MatDoub>  hhatvel)=0;
	//virtual elastoplasticbase GetYieldCriterion()=0;
	virtual MatDoub GetBodyForce() =0;
	//virtual void GetMaterialRandomField()=0;

	//MatDoub fHHAT;
	//NRvector<MatDoub> fhhatvel;
	//Int fglobalcounter;
	// 
//	void Contribute(MatDoub& ek, MatDoub& ef, Doub xi, Doub eta, Doub w, MatDoub elcoords);
//	void CacStiff(MatDoub& ek, MatDoub& ef, const MatDoub& elcoords);
//	void Assemble(MatDoub& KG, MatDoub& FG, const std::vector<std::vector< std::vector<Doub > > >& allcoords, const MatDoub& meshnodes, const MatInt meshtopology);
//	void assembleBandN(MatDoub& B, MatDoub& N, const MatDoub& psis, const MatDoub& GradPhi);
//	void assembleConstitutiveMatrix(MatDoub& C, Doub mult);
//void GetElCoords(std::vector<std::vector< std::vector<Doub > > >& allcoords, Int el, MatDoub& elcoords);
//	void DirichletBC(MatDoub& KG, MatDoub& FG, std::vector<int> ids, Int  dir, Int val);
	//void ContributeLineNewan(MatDoub& KG, MatDoub& FG, std::vector<int> ids, Int  dir, Int val);
	//virtual void SolPt(const Int& el, const  MatDoub& solG, const Doub& xi, const Doub& eta, MatDoub& xycoords, MatDoub& sol)=0;
//
//	void PostProcess(const MatDoub& nodalsol, std::vector<std::vector<double>>& solx, std::vector<std::vector<double>>& soly);
	//virtual void PostProcess(const MatDoub& nodalsol, std::vector<std::vector<double>>& sol)=0;

};

