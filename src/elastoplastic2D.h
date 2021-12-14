#pragma once
#include "material.h"
#include "shape.h"
#include "mesh.h"
#include "error.h"

    /** @brief This class implements an elastoplastic material and store the hardening variables and strain state.
    Detailed description follows here.
    @author Diogo Cecílio
    @date dez 2021
    @brief Reference: <a href="https://ascelibrary.org/doi/abs/10.1061/%28ASCE%29EM.1943-7889.0001737">link text</a>
    */

using namespace std;

template <class YC>
class elastoplastic2D : public material{
public:
    /**
	 * @brief Class constructor
	 * @param [in] thickness 
     * @param [in] bodyforce 
     * @param [in] planestress  plane stress = 1 and plane strain = 0 
     * @param [in] order interpolation order
     * @param [in] HHAT random field
	 */
	elastoplastic2D(Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int order, NRmatrix<Doub>   HHAT);
        /**
	 * @brief Class constructor
	 * @param [in] thickness 
     * @param [in] bodyforce 
     * @param [in] planestress  plane stress = 1 and plane strain = 0 
     * @param [in] order interpolation order
	 */
	elastoplastic2D(Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int order);

    /**
     * @brief Copy constructor
     */
	elastoplastic2D(elastoplastic2D & copy);
    
       /**
     * @brief Default constructor
     */
	elastoplastic2D();
       /**
     * @brief Default destructor
     */
	~elastoplastic2D();

    /**
	 * @brief    It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param [in] ptsw integration point and weigth
     * @param [in] elcoords element nodal coordinates
     * @param [in] eldisplace  element nodal displacement
	 */
	void Contribute(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody,NRvector<Doub> ptsw, NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace);

/**
	 * @brief Computes the element stiffness matrix and right hand side
     * @param [in] elcoords element nodal coordinates
     * @param [in] eldisplace  element nodal displacement
     * @param [out] ek element matrix
	 * @param [out] efint element right hand side internal forces
     * @param [out] efbody element right hand side body forces
	 */
	void CalcStiff(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, const NRmatrix<Doub>   &elcoords, NRmatrix<Doub>  eldisplace);
    
    /**
	 * @brief Assemble the displacement matrix B and the shape matrix N
	 */
	void assembleBandN(NRmatrix<Doub>  &B, NRmatrix<Doub>  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi);

    /**
	 * @brief Assemble the displacement matrix B and the shape matrix N
	 */
    void assembleBandN(MatrixXd&B, MatrixXd  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi);

    /**
	 * @brief Get element coordinates
	 */
	void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords);
    
    /**
	 * @brief Imposes Dirichlet boundary conditions
	 */
	void DirichletBC(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val);
    
    /**
	 * @brief Imposes Dirichlet boundary conditions
	 */
    void DirichletBC(SparseMatrix<double> & KG, VectorXd& FG, std::vector<int> ids, Int  dir, Int val);
    
    /**
	* @brief Uniforme load right hand side
	 */
	void ContributeLineNewan(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val);
    
    /**
	* @brief Uniforme load right hand side over a curced boundary
	 */
	void ContributeCurvedLine(NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology, Doub force);
    
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
    
    /**
	 * @brief Post process the solution with a specific refinement
     * @param [in] allcoords 
     * @param [in] meshnodes 
     * @param [in] meshtopology   
     * @param [in] nodalsol element id
     * @param [out] solx x and y coordinate and post processed displacement in x direction [[xcoord,ycoord,ux],...[]]
	 * @param [out] soly  x and y coordinate and post processed displacement in y direction [[xcoord,ycoord,uy],...,[]]
	 */
	void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly);
    
    /**
	 * @brief Post process the solution with a specific refinement
     * @param [in] allcoords 
     * @param [in] meshnodes 
     * @param [in] meshtopology
     * @param [in] var what to post process: ux:0 and uy:1     
     * @param [in] nodalsol element id
     * @param [out] sol x and y coordinate and post processed displacement [[xcoord,ycoord,ux or uy],...[]]
	 */
    void PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, Int var,const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol);

    
    /**
	 * @brief Post process the solution with a specific refinement
     * @param [in] allcoords 
     * @param [in] meshnodes 
     * @param [in] meshtopology
     * @param [in] var what to post process: ux:0 and uy:1     
     * @param [in] nodalsol element id
     * @param [out] sol x and y coordinate and post processed displacement [[xcoord,ycoord,ux or uy],...[]]
	 */
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
    NRvector<NRtensor<Doub> > GetPlasticStrain();
    void ComputeSolution(mesh * inmesh,NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace,NRvector<Doub> ptsw,NRmatrix<Doub> &sol,NRmatrix<Doub> &dsol);

    void ComputeSolAndDSol(mesh * inmesh,NRmatrix<Doub>&sol,NRmatrix<Doub>&dsol);
    void ComputeSolAndDSol(mesh * inmesh,NRvector<NRmatrix<Doub>>&sol,NRvector<NRmatrix<Doub>>&dsol);

    Doub ComputePhi(NRtensor<Doub> eps)
    {
        Doub valphi=fYC.phi(eps);
        return valphi;
    }

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
Int getIndex(vector<int> indexvec, int inode)
{
    Int id;
    std::vector<int>::iterator it = std::find(indexvec.begin(), indexvec.end(), inode);
    if (it != indexvec.end() || indexvec.size()!=0)
    {
        id = std::distance(indexvec.begin(), it);
        //std::cout << "Element Found" << std::endl;
    }else
    {
        id=Int(10e12);
    }
    return id;
}

void deleteduplicates(std::vector<int> &v)
{
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it) {
        end = std::remove(it + 1, end, *it);
    }

    v.erase(end, v.end());
}

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



