#pragma once

#include "material.h"
#include "mesh.h"
#include "shape.h"
#include "shapequad.h"
#include "shapetri.h"
class elastmat2D : public material
{
public:
    elastmat2D (Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order );
    elastmat2D ( Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT );
    elastmat2D();
    ~elastmat2D();

    void Contribute ( MatDoub &ek, MatDoub &ef, Doub xi, Doub eta, Doub w, MatDoub elcoords );
	//void CalcStiff ( MatDoub& ek, MatDoub& efint, MatDoub& efbody, const MatDoub& elcoords, MatDoub eldisplace ) =0;
    void CalcStiff ( MatDoub &ek, MatDoub &ef, const MatDoub  &elcoords );
	void CalcStiff ( MatDoub& ek, MatDoub& efint, MatDoub& efbody, const MatDoub& elcoords, MatDoub eldisplace ){DebugStop();}
    //void Assemble(MatDoub &KG, MatDoub &FG, const std::vector<std::vector< std::vector<Doub > > > &allcoords, const MatDoub &meshnodes, const MatInt meshtopology);
    void assembleBandN ( MatDoub &B, MatDoub &N, const MatDoub &psis, const MatDoub &GradPhi );
    void assembleConstitutiveMatrix ( MatDoub &ce );

    void DirichletBC ( MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val );
    void ContributeLineNewan ( MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val );

	
	 virtual void Contribute ( MatDoub& ek, MatDoub& efint, MatDoub& efbody,NRvector<Doub> ptsw, MatDoub elcoords, MatDoub eldisplace ){DebugStop();};
    
	virtual void SetTangentMatrixType ( bool type ){DebugStop();}


void GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords )
{
    elcoords.assign ( allcoords[el].size(), 2, 0. );
    for ( Int j = 0; j < allcoords[el].size(); j++ ) {
        Doub x = allcoords[el][j][0];
        Doub y = allcoords[el][j][1];
        elcoords[j][0] = x;
        elcoords[j][1] = y;
    }
}
// 	void  ComputeLoadVector( MatDoub &KG, MatDoub & FG, std::vector<int> ids,MatDoub & force );
     //virtual void GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub& elcoords ){DebugStop();}

    virtual void DirichletBC ( SparseMatrix<double>  &KG, VectorXd &Fint,std::vector<int> ids, Int  dir, Int val ){DebugStop();}

    void AssembleLoadvector ( NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology, Doub  force );
	
	void ContributeCurvedLine ( NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology,Doub force ){DebugStop();}
    //virtual void SolPt(const std::vector<std::vector< std::vector<Doub > > >& allcoords, const MatInt& meshtopology, const Int& el, const  MatDoub& solG, const Doub& xi, const Doub& eta, MatDoub& xycoords, MatDoub& sol)=0;
   virtual  void SolPt ( mesh * inmesh, const Int &el, const  NRmatrix<Doub>  &nodalsol, const Doub &xi, const Doub &eta, NRmatrix<Doub>  &xycoords, MatDoub &sol ){ DebugStop();}
    //
    virtual void PostProcess ( std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, const MatDoub& nodalsol,std::vector<std::vector<double>>& solx, std::vector<std::vector<double>>& soly ) {DebugStop();}

    //virtual void PostProcess(mesh * inmesh, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly)=0;
   virtual  void PostProcess ( std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, Int var, const MatDoub& nodalsol, std::vector<std::vector<double>>& sol ) {DebugStop();}
    // virtual void PostProcessStrain(mesh * inmesh, NRvector<NRvector<NRtensor<Doub>>> &sol)=0;
   virtual  void PostProcessIntegrationPointVar ( std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, const MatDoub& nodalsol, std::vector<std::vector<double>>& sol ){DebugStop();}

    virtual void SetRandomField ( MatDoub hhat ) {DebugStop();}
    virtual void SetMemory ( Int ngloblapoints, Int systemsize ){DebugStop();}
    virtual void UpdateDisplacement ( VectorXd displace ) {DebugStop();}
    virtual void UpdatePlasticStrain() {DebugStop();}
    virtual void UpdateBodyForce ( MatDoub newbodyforce ) {DebugStop();}
    virtual void ResetPlasticStrain() {DebugStop();}
    virtual void ResetDisplacement() {DebugStop();}
    virtual void ResetMat() {DebugStop();}
    virtual void ResetMemory() {DebugStop();}


   virtual  NRvector<NRtensor<Doub> > GetPlasticStrain(){DebugStop();}
    virtual void SetMatConstants ( NRvector<Doub>& consts ){DebugStop();}
    virtual void GetMatConstants ( NRvector<Doub>& consts ) {DebugStop();}
    virtual void ResetCounter() {DebugStop();}
    virtual void SetRandomFieldLocal ( NRvector<MatDoub>  hhatvel ) {DebugStop();}
    //virtual elastoplasticbase GetYieldCriterion()=0;
    virtual MatDoub GetBodyForce() {DebugStop();}
    
   virtual  void ComputeSolAndDSol ( mesh * inmesh,NRvector<NRmatrix<Doub>>&sol,NRvector<NRmatrix<Doub>>&dsol ) {DebugStop();}
   virtual  NRvector<Doub> ComputePhi ( NRtensor<Doub> eps ){
	   NRvector<Doub>  dumb(3,0.);
	   return dumb;
}
   virtual  void ComputeSolution ( mesh * inmesh,NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace,NRvector<Doub> ptsw,NRmatrix<Doub> &sol,NRmatrix<Doub> &dsol ) {DebugStop();}
   
    void ComputeSolAndDSol ( mesh * inmesh,NRmatrix<Doub>&sol,NRmatrix<Doub>&dsol );
	void UpdateDisplacement ( NRmatrix<Doub>  displace )
	{
		fdisplace = displace;
	}

	NRmatrix<Doub> GetSolution()
	{
		return fdisplace;	
	}
	
private:
    Doub fyoung;
    Doub fnu;
    Doub fbodyforce;
    Int fplanestress;
    Doub fthickness;
    Int fOrder;
   // mesh fmesh;
	NRmatrix<Doub>  fdisplace;
	shape *fshape;
public:
    MatDoub fHHAT;
    MatDoub fhhatvel;
	
};

