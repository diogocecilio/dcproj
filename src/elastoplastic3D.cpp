#include "elastoplastic3D.h"
#include "vonmises.h"
#include "druckerprager.h"

template <class YC>
elastoplastic3D<YC>::elastoplastic3D(NRmatrix<Doub> bodyforce,Int order)
{
	fbodyforce = bodyforce;
	fOrder = order;
//	fHHAT = HHAT;
	//fmesh = inmesh;

}

template <class YC>
elastoplastic3D<YC>::elastoplastic3D( NRmatrix<Doub>  bodyforce, Int order, NRmatrix<Doub>  HHAT)
{
	fbodyforce = bodyforce;
	fOrder = order;
	fHHAT = HHAT;

}

template <class YC>
elastoplastic3D<YC>::elastoplastic3D(elastoplastic3D & copy)
{
}

template <class YC>
elastoplastic3D<YC>::elastoplastic3D()
{
}



template <class YC>
elastoplastic3D<YC>::~elastoplastic3D()
{
}

template <class YC>
void elastoplastic3D<YC>::SetMemory(Int nglobalpts, Int systemsize)
{
	fdisplace.assign(systemsize,1, 0.);
 
	fepspvec.assign(nglobalpts, 0.);
	fepspsolitern.assign(nglobalpts, 0.);
	fglobalcounter = 0;

}

template <class YC>
void elastoplastic3D<YC>::ResetMemory()
{
	Int systemsize = fdisplace.nrows();
	fdisplace.assign(systemsize, 1, 0.);

	Int nglobalpts = fepspvec.size();
	fepspvec.assign(nglobalpts, 0.);
	fepspsolitern.assign(nglobalpts, 0.);
	fglobalcounter = 0;
}

template <class YC>
void elastoplastic3D<YC>::SetRandomField(NRmatrix<Doub>  hhat)
{
	fHHAT = hhat;
}
template <class YC>
void elastoplastic3D<YC>::SetRandomFieldLocal(NRvector<NRmatrix<Doub> >  hhatvel)
{
	fhhatvel = hhatvel;
}


template <class YC>
void elastoplastic3D<YC>::UpdateDisplacement(NRmatrix<Doub>  displace)
{
	fdisplace = displace;
}

template <class YC>
void elastoplastic3D<YC>::UpdateDisplacement(VectorXd displace)
{
    double sz =  displace.size();
    fdisplace.resize(sz,1);
	for(int i=0;i<sz;i++)fdisplace[i][0] = displace(i);
}

template <class YC>
void elastoplastic3D<YC>::UpdatePlasticStrain()
{
	fepspsolitern = fepspvec;
	fglobalcounter = 0;

}

template <class YC>
void elastoplastic3D<YC>::UpdateBodyForce(NRmatrix<Doub>  newbodyforce)
{
	fbodyforce = newbodyforce;
}

template <class YC>
void elastoplastic3D<YC>::ResetPlasticStrain()
{

	
	//NRvector<TensorDoub> fepspvec, fepspsolitern;
	Int rows = fepspvec.size();
	for (Int i = 0;i < rows;i++)
	{
		for (Int j = 0;j < 6;j++)
		{
			fepspvec[i].XX() = 0. ;
			fepspvec[i].YY() = 0. ;
			fepspvec[i].ZZ() = 0. ;
			fepspvec[i].XZ() = 0. ;
			fepspvec[i].YZ() = 0. ;
			fepspvec[i].XY() = 0. ;

			fepspsolitern[i].XX() = 0.;
			fepspsolitern[i].YY() = 0.;
			fepspsolitern[i].ZZ() = 0.;
			fepspsolitern[i].XZ() = 0.;
			fepspsolitern[i].YZ() = 0.;
			fepspsolitern[i].XY() = 0.;
		}
	}

}

template <class YC>
NRvector<NRtensor<Doub> >   elastoplastic3D<YC>::GetPlasticStrain()
{
	return fepspvec;
}

template <class YC>
void elastoplastic3D<YC>::ResetMat()
{
	fglobalcounter = 0;
	ResetDisplacement();
	ResetPlasticStrain();
	ResetMemory();

}

template <class YC>
void elastoplastic3D<YC>::ResetDisplacement()
{
	Int rows = fdisplace.nrows();
	fdisplace.assign(rows, 1,0.);
}

template <class YC>
void elastoplastic3D<YC>::Contribute(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, NRvector<Doub>  intptsw, NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace)
{
	int type = 1;
	shapehexahedron objshapes(fOrder, type);

    Doub xi,eta,zeta,w;
    xi = intptsw[0];
    eta = intptsw[1];
    zeta = intptsw[2];
    w = intptsw[3];
	NRmatrix<Doub>  psis, GradPsi, elcoordst, xycoords, Jac, InvJac(3, 3);
	NRmatrix<Doub > GradPhi, B, BT, N, NT, psist, C, BC, BCS, stress(6, 1, 0.), temp, CS, KSt;
	objshapes.shapes(psis, GradPsi, xi, eta,zeta);
	psis.Transpose(psist);
   // elcoords.Print();
	psist.Mult(elcoords, xycoords);
	NRvector<NRmatrix<Doub> > hhat(fhhatvel.size());

	if (fhhatvel.size() != 0)
	{
		for (Int ivar = 0;ivar < fhhatvel.size();ivar++)
		{
			psist.Mult(fhhatvel[ivar], hhat[ivar]);
		}
		
	}
	GradPsi.Mult(elcoords, Jac);
    Doub DetJ = -Jac[0][2]*Jac[1][1]*Jac[2][0] + Jac[0][1]*Jac[1][2]*Jac[2][0] + Jac[0][2]*Jac[1][0]*Jac[2][1] - Jac[0][0]*Jac[1][2]*Jac[2][1] - Jac[0][1]* Jac[1][0]* Jac[2][2] + Jac[0][0]* Jac[1][1]* Jac[2][2];
   // std::cout << "\n DetJ  =  "<< DetJ << endl;
	if (DetJ <= 0)
    //if (true)
	{

		std::cout << "\n DetJ < 0 " << endl;
		std::cout << "\n xi " << xi << endl;
		std::cout << "\n eta " << eta << endl;
        std::cout << "\n zeta " << zeta << endl;
         std::cout << "\n w " << w << endl;
        Jac.Print();
		GradPsi.Print();
		elcoords.Print();
		GradPsi.Mult(elcoords, Jac);
		xycoords.Print();

		psis.Print();

		objshapes.shapes(psis, GradPsi, xi, eta,zeta);
		std::cout << "\n Det < 0 " <<std::endl;
        DebugStop();
	}
	InvJac[0][0]= (-Jac[1][2]* Jac[2][1] + Jac[1][1]* Jac[2][2])/DetJ;
    InvJac[0][1]= (Jac[0][2]* Jac[2][1] - Jac[0][1]*Jac[2][2])/DetJ;
    InvJac[0][2]= (-Jac[0][2]*Jac[1][1] + Jac[0][1]*Jac[1][2])/DetJ;

    InvJac[1][0]= (Jac[1][2]*Jac[2][0] - Jac[1][0]*Jac[2][2])/DetJ;
    InvJac[1][1]= (-Jac[0][2]* Jac[2][0] + Jac[0][0]* Jac[2][2])/DetJ;
    InvJac[1][2]= (Jac[0][2]* Jac[1][0] - Jac[0][0]* Jac[1][2])/DetJ;

    InvJac[2][0]= (-Jac[1][1]* Jac[2][0]+ Jac[1][0]* Jac[2][1])/DetJ;
    InvJac[2][1]= (Jac[0][1]* Jac[2][0]- Jac[0][0]* Jac[2][1])/DetJ;
    InvJac[2][2]= (-Jac[0][1]* Jac[1][0] + Jac[0][0]* Jac[1][1])/DetJ;

	InvJac.Mult(GradPsi, GradPhi);
	NRmatrix<Doub>  gradu;
	GradPhi.Mult(eldisplace, gradu);

	//{ {dudx, dudy, dudz}, { dvdx, dvdy, dvdz }, {dwdx, dwdy, dwdz}} = gradprevsol;
	Doub ex = gradu[0][0];// dudx;
	Doub ey = gradu[1][1];
    Doub ez = gradu[2][2];
	Doub exy = (gradu[0][1] + gradu[1][0]);
    Doub exz = (gradu[0][2] + gradu[2][0]);
    Doub eyz = (gradu[1][2] + gradu[2][1]);

	NRtensor<Doub>   epst(0.), epsp(0.), projstress(0.), projstrain(0.), epspeint(0.);
    
	NRmatrix<Doub>  Dep;
	Doub  projgamma=0.;
	epst.XX() = ex;epst.YY() = ey;epst.ZZ()=ez;epst.XY() = exy;epst.XZ()=exz;epst.YZ()=eyz;
	//epsp = epspsoliternGLOBAL[globalcounter];
	epsp = fepspsolitern[fglobalcounter];

	//std::cout << "epst = " << std::endl;
	//epst.Print();
	//std::cout << "epsp = " << std::endl;
	//epsp.Print();
	//cout << "\n antes c  = " << fYC.GetCoes() << endl;
	if (fhhatvel.size() != 0)
	{
		for (Int ivar = 0;ivar < fhhatvel.size();ivar++)
		{
			psist.Mult(fhhatvel[ivar], hhat[ivar]);
		}
		fYC.updateatributes(hhat);
	}
	//cout << "\n c  = " << fYC.GetCoes() <<endl;

	fYC.closestpointproj(epst,epsp,projstress,projstrain,Dep,projgamma);
	if (fhhatvel.size() != 0)
	{
		fYC.restoreoriginalatributes();
	}
	//epspeint = epst - (epse);
	epspeint = epst;
	epspeint = epspeint - projstrain;

	fepspvec[fglobalcounter] = epspeint;
	//epspvecGLOBAL[globalcounter] = epspeint;

	fglobalcounter++;
	assembleBandN(B, N, psis, GradPhi);
	N.Transpose(NT);
	B.Transpose(BT);
	
	//Dep.Print();

    //B.Print();
   // N.Print();
   // BC=(BT*Dep)*B;
	BT.Mult(Dep, BC);
	BC.Mult(B, ek);

    //#define _XX_ 0#define _YY_ 1#define _ZZ_ 2#define _XZ_ 3#define _YZ_ 4#define _XY_ 5
	stress[0][0] = projstress.XX();stress[1][0] = projstress.YY();stress[2][0] = projstress.ZZ();
    stress[5][0] = projstress.XY();stress[3][0] = projstress.XZ();stress[3][0] = projstress.YZ();
	ek *= w*DetJ;
	//std::cout << " stress "<< std::endl;
	//stress.Print();
	//ek.Print();
	//ef = (Transpose[BB].stress) weight DetJ;
	//ef2 = (Transpose[NShapes].{0, -bodyforce}) weight DetJ;
	BT.Mult(stress, efint);
    //fbodyforce.Print();
	NT.Mult(fbodyforce, efbody);
	efint *= w*DetJ;
	efbody *= w*DetJ;
}
template <class YC>
void elastoplastic3D<YC>::CacStiff(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, const NRmatrix<Doub>   &elcoords,NRmatrix<Doub>  eldisplace)
{
	MatDoub ptsweigths, ekt, eftint,eftbody;
	Doub xi, eta, w;
	Int nnodes = elcoords.nrows();
	ek.assign(nnodes * 3, nnodes * 3, 0.);
	efint.assign(nnodes * 3, 1, 0.);
	efbody.assign(nnodes * 3, 1, 0.);
	shapehexahedron shape = shapehexahedron(fOrder, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();

    NRvector<Doub> ptsw(4);
	for (Int ipt = 0;ipt < npts;ipt++)
	{

		ptsw[0]=ptsweigths[ipt][0];
		ptsw[1]= ptsweigths[ipt][1];
		ptsw[2]= ptsweigths[ipt][2];
        ptsw[3]=ptsweigths[ipt][3];
		//std::cout << "integration point :" << " xi = " << xi << "  eta = "<<  eta << std::endl;
		Contribute(ekt, eftint, eftbody, ptsw, elcoords, eldisplace);
		ek += ekt;
		efint += eftint;
		efbody += eftbody;
	}
}

template <class YC>
void elastoplastic3D<YC>::assembleBandN(NRmatrix<Doub>  &B, NRmatrix<Doub>  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi)
{
	B.assign(6, psis.nrows() * 3, 0.);
	N.assign(3, psis.nrows() * 3, 0.);

	Int j = 0, k = 0,w=0,indexaux=1;
	for (Int i = 0;i < psis.nrows() * 3;i++)
	{
		if ((i) % 3 == 0)// 0 3 6 9 ...
		{

            B[0][i] = GradPhi[0][j];
			B[1][i] = 0;
            B[2][i] = 0;
            B[3][i] = 0;
            B[4][i] = GradPhi[2][j];
			B[5][i] = GradPhi[1][j];

			N[0][i] = psis[j][0];
			N[1][i] = 0;
            N[2][i] = 0;
			j++;
		}
		else
        {
            if(i == indexaux) // 1 4 7 9 ...
            {
                B[0][i] = 0;
                B[1][i] = GradPhi[1][k];
                B[2][i] = 0;
                B[3][i] = GradPhi[2][k];
                B[4][i] = 0;
                B[5][i] = GradPhi[0][k];

                N[0][i] = 0;
                N[1][i] = psis[k][0];
                N[2][i] = 0;

                k++;
                indexaux+=3;
            }
            else
            {
                B[0][i] = 0;
                B[1][i] = 0;
                B[2][i] = GradPhi[2][w];
                B[3][i] = GradPhi[1][w];
                B[4][i] = GradPhi[0][w];
                B[5][i] = 0;


                N[0][i] = 0;
                N[1][i] = 0;
                N[2][i] = psis[w][0];
                w++;

            }
        }

    }
}

template <class YC>
void elastoplastic3D<YC>::GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords)
{
	elcoords.assign(allcoords[el].size(), 3, 0.);
	for (Int j = 0; j < allcoords[el].size(); j++)
	{
		Doub x = allcoords[el][j][0];
		Doub y = allcoords[el][j][1];
        Doub z = allcoords[el][j][2];
		elcoords[j][0] = x;
		elcoords[j][1] = y;
        elcoords[j][2] = z;
	}
}
template <class YC>
void elastoplastic3D<YC>::DirichletBC(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val)
{
	Int nodes = ids.size();
	Int sz = KG.nrows();
	Int fu = 0, cu = 0;
	for (Int i = 0;i < nodes;i++)
	{
		Int pso = ids[i];
		if (dir == 0)
		{
			for (Int j = 0; j < sz;j++)
			{
				KG[3 * pso][j] = 0;
				KG[j][3 * pso] = 0;
			}
			KG[3 * pso][3 * pso] = 1;
			FG[3 * pso][0] = val;
		}
		if (dir == 1)
		{
			for (Int j = 0; j < sz;j++)
			{
				KG[3 * pso + 1][j] = 0;
				KG[j][3 * pso + 1] = 0;
			}
			KG[3 * pso + 1][3 * pso + 1] = 1;
			FG[3 * pso + 1][0] = val;
		}
		if (dir == 2)
        {
            for (Int j = 0; j < sz;j++)
			{
				KG[3 * pso + 2][j] = 0;
				KG[j][3 * pso + 2] = 0;
			}
			KG[3 * pso + 2][3 * pso + 2] = 1;
			FG[3 * pso + 2][0] = val;
        }
	}
}

template <class YC>
void elastoplastic3D<YC>::DirichletBC(SparseMatrix<double>  &KG, VectorXd &Fint,std::vector<int> ids, Int  dir, Int val)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
void elastoplastic3D<YC>::ContributeLineNewan(NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
void elastoplastic3D<YC>::ContributeCurvedLine(NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology,Doub force)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
void elastoplastic3D<YC>::SolPt(const std::vector<std::vector< std::vector<Doub > > > &allcoords, const MatInt &meshtopology, const Int &el, const  NRmatrix<Doub>  &solG, const Doub &xi, const Doub &eta, NRmatrix<Doub>  &xycoords, NRmatrix<Doub>  &sol)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
void elastoplastic3D<YC>::PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnode, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}
template <class YC>
//void elastoplastic2D<YC>::PostProcess(mesh & inmesh,Int var,const MatDoub & nodalsol, std::vector<std::vector<double>> &sol)
void elastoplastic3D<YC>::PostProcess(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnode, MatInt meshtopology, Int var, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
void elastoplastic3D<YC>::PostProcessIntegrationPointVar(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnode, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
void elastoplastic3D<YC>::Assemble(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, NRmatrix<Doub>  &KG, NRmatrix<Doub>  &Fint,NRmatrix<Doub>  &Fbody)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
void elastoplastic3D<YC>::Assemble(std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnodes, MatInt meshtopology, SparseMatrix<double>  &KG, VectorXd &Fint, VectorXd &Fbody)
{
    std::cout<< "Not Implemented." << std::endl;
    DebugStop();
}

template <class YC>
NRmatrix<Doub>  elastoplastic3D<YC>::GetSolution()
{
	return fdisplace;
}


template class elastoplastic3D<vonmises>;
template class elastoplastic3D<druckerprager>;
