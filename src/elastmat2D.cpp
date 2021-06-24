#include "elastmat2D.h"
#include "gridmesh.h"

using namespace std;
elastmat2D::elastmat2D(mesh &inmesh, Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order) 
{
	fyoung = young;
	fnu = nu;
	fbodyforce = bodyforce;
	fplanestress = planestress;
	fthickness = thickness;
	fOrder = order;
	fmesh = inmesh;
}

elastmat2D::elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT)
{
	fyoung = young;
	fnu = nu;
	fbodyforce = bodyforce;
	fplanestress = planestress;
	fthickness = thickness;
	fOrder = order;
	fHHAT = HHAT;
}

elastmat2D::elastmat2D()
{

}

elastmat2D::~elastmat2D()
{
    

}


void elastmat2D::Assemble(MatDoub &KG, MatDoub &FG, const std::vector<std::vector< std::vector<Doub > > > &allcoords, const MatDoub &meshnodes, const MatInt meshtopology)
{
	MatDoub ek, ef, elcoords, eltopology;
	std::vector<std::vector< std::vector<Doub > > > alco = fmesh.GetAllCoords();
	GetElCoords(alco, 0, elcoords);
	Int rows = elcoords.nrows();
	Int sz = 2 * meshnodes.nrows();
	Int cols = rows;
	//KG.assign(sz, sz, 0.);
	//FG.assign(sz, 1, 0.);
	Int nels = allcoords.size();

	Int fu = 0;

	for (Int iel = 0;iel < nels;iel++)
	{
		if (fHHAT.nrows() != 0)
		{
			fhhatvel.assign(rows, 1, 0.);
			for (Int inode = 0;inode < rows;inode++)
			{
				fhhatvel[inode][0] = fHHAT[meshtopology[iel][inode]][0];
			}
		}

		GetElCoords(alco, iel, elcoords);
		CacStiff(ek, ef, elcoords);
		for (Int irow = 0;irow < rows;irow++)
		{
			Int rowglob = meshtopology[iel][irow];
			for (Int icol = 0;icol < cols;icol++)
			{
				Int colglob = meshtopology[iel][icol];
				KG[2 * rowglob + fu][2 * colglob + fu] += ek[2 * irow + fu][2 * icol + fu];
				KG[2 * rowglob + fu][2 * colglob + 1 + fu] += ek[2 * irow + fu][2 * icol + 1 + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + fu] += ek[2 * irow + 1 + fu][2 * icol + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + 1 + fu] += ek[2 * irow + 1 + fu][2 * icol + 1 + fu];

			}
			FG[2 * rowglob + fu][0] += ef[2 * irow + fu][0];
			FG[2 * rowglob + 1 + fu][0] += ef[2 * irow + 1 + fu][0];
		}

	}

}

void elastmat2D::CacStiff(MatDoub &ek, MatDoub &ef, const MatDoub  &elcoords)
{

	MatDoub ptsweigths, ekt, eft;
	Doub xi, eta, w;
	Int nnodes = elcoords.nrows();
	ek.assign(nnodes * 2, nnodes * 2, 0.);
	ef.assign(nnodes * 2, 1, 0.);
	shapequad shape = shapequad(fOrder, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();

	for (Int ipt = 0;ipt < npts;ipt++)
	{
		xi = ptsweigths[ipt][0];
		eta = ptsweigths[ipt][1];
		w = ptsweigths[ipt][2];
		Contribute(ekt, eft, xi, eta, w, elcoords);
		ek += ekt;
		ef += eft;
	}

}

void elastmat2D::Contribute(MatDoub &ek, MatDoub &ef, Doub xi, Doub eta, Doub w, MatDoub elcoords)
{
	MatDoub psis, GradPsi, elcoordst, xycoords, Jac, InvJac(2, 2), GradPhi, B, BT, N, NT, psist, C, BC, BCS, stress(3, 1, 0.), bodyforce(2, 1), temp, CS, KSt;
	bodyforce[0][0] = 0;
	bodyforce[1][0] = -fbodyforce;

	int type = 1;
	shapequad objshapes(fOrder, type);

	objshapes.shapes(psis, GradPsi, xi, eta);
	psis.Transpose(psist);
	psist.Mult(elcoords, xycoords);
	MatDoub hhat;
	if (fhhatvel.nrows() != 0)
	{
		psist.Mult(fhhatvel, hhat);
	}

	GradPsi.Mult(elcoords, Jac);
	Int nnodes = psis.nrows();
	Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
	if (DetJ <= 0)
	{

		std::cout << "\n DetJ < 0 " << endl;
		std::cout << "\n xi " << xi << endl;
		std::cout << "\n eta " << eta << endl;
		GradPsi.Print();
		elcoords.Print();
		GradPsi.Mult(elcoords, Jac);
		xycoords.Print();

		psis.Print();

		objshapes.shapes(psis, GradPsi, xi, eta);
		return;
	}
	InvJac[0][0] = Jac[1][1] / DetJ;   InvJac[0][1] = -Jac[0][1] / DetJ;
	InvJac[1][0] = -Jac[1][0] / DetJ;	InvJac[1][1] = Jac[0][0] / DetJ;
	InvJac.Mult(GradPsi, GradPhi);
	assembleBandN(B, N, psis, GradPhi);
	N.Transpose(NT);
	B.Transpose(BT);
	assembleConstitutiveMatrix(C, 1.);
	BT.Mult(C, BC);
	BC.Mult(B, ek);

	if (fhhatvel.nrows() != 0)
	{
		Doub mult = hhat[0][0];
		assembleConstitutiveMatrix(CS, mult);
		BT.Mult(CS, BCS);
		BCS.Mult(B, KSt);
		ek += KSt;
	}
	else {

	}

	ek *= w*DetJ*fthickness;
	BT.Mult(stress, ef);
	NT.Mult(bodyforce, temp);
	ef -= temp;
	ef *= w*DetJ;
}

void elastmat2D::assembleBandN(MatDoub &B, MatDoub &N, const MatDoub &psis, const MatDoub &GradPhi)
{
	B.assign(3, psis.nrows() * 2, 0.);
	N.assign(2, psis.nrows() * 2, 0.);

	Int j = 0, k = 0;
	for (Int i = 0;i < psis.nrows() * 2;i++)
	{
		if (i % 2 == 0)
		{
			B[0][i] = GradPhi[0][j];
			B[1][i] = 0;
			B[2][i] = GradPhi[1][j];

			N[0][i] = psis[j][0];
			N[1][i] = 0;
			j++;
		}
		else {
			B[0][i] = 0;
			B[1][i] = GradPhi[1][k];
			B[2][i] = GradPhi[0][k];

			N[0][i] = 0;
			N[1][i] = psis[k][0];

			k++;
		}

	}
}
void elastmat2D::assembleConstitutiveMatrix(MatDoub &C, Doub mult)
{
	Doub nusqr = fnu*fnu;
	Doub young = mult*fyoung, nu = fnu;
	C.assign(3, 3, 0.);
	C[0][0] = young / (1 - nusqr);   C[0][1] = nu*young / (1 - nusqr);C[0][2] = 0.;
	C[1][0] = nu*young / (1 - nusqr);C[1][1] = young / (1 - nusqr);   C[1][2] = 0.;
	C[2][0] = 0.;                    C[2][1] = 0.;                    C[2][2] = young / (2 * (1 + nu));
}

void elastmat2D::GetElCoords(std::vector<std::vector< std::vector<Doub > > > &allcoords, Int el, MatDoub & elcoords)
{

	elcoords.assign(allcoords[el].size(), 2, 0.);
	for (Int j = 0; j < allcoords[el].size(); j++)
	{
		Doub x = allcoords[el][j][0];
		Doub y = allcoords[el][j][1];
		elcoords[j][0] = x;
		elcoords[j][1] = y;
	}
}

void elastmat2D::DirichletBC(MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val)
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
				KG[2 * pso][j] = 0;
				KG[j][2 * pso] = 0;
			}
			KG[2 * pso][2 * pso] = 1;
			FG[2 * pso][0] = val;
		}
		else
		{
			for (Int j = 0; j < sz;j++)
			{
				KG[2 * pso + 1][j] = 0;
				KG[j][2 * pso + 1] = 0;
			}
			KG[2 * pso + 1][2 * pso + 1] = 1;
			FG[2 * pso + 1][0] = val;
		}


	}

}
void elastmat2D::ContributeLineNewan(MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val)
{
	MatDoub psis, gradpsis, ptsws;

	int type = 1;
	shapequad objshapes(fOrder, type);

	objshapes.pointsandweigths1D(ptsws);
	Int npts = ptsws.nrows();

	Doub xi, w;
	MatDoub DetJ;
	for (Int ipt = 0;ipt < npts;ipt++)
	{
		xi = ptsws[ipt][0];
		w = ptsws[ipt][2];
		objshapes.shapes1D(psis, gradpsis, xi);
		psis.Mult(gradpsis, DetJ);

	}

}

void elastmat2D::SolPt( const Int &el, const  MatDoub &solG, const Doub &xi, const Doub &eta, MatDoub &xycoords, MatDoub &sol)
{
	int type = 1;
	shapequad objshapes(fOrder, type);
	MatDoub psis, GradPsi, elcoords, psist, solel;
		std::vector<std::vector< std::vector<Doub > > > alco = fmesh.GetAllCoords();
	GetElCoords(alco, el, elcoords);
	objshapes.shapes(psis, GradPsi, xi, eta);
	Int nodes = psis.nrows();
	Int nstatevars = 1;
	solel.assign(nodes, 1, 0);
	psis.Transpose(psist);
	psist.Mult(elcoords, xycoords);

	for (Int inode = 0; inode < nodes;inode++)
	{
		solel[inode][0] = solG[fmesh.GetMeshTopology()[el][inode]][0];
	}
	psist.Mult(solel, sol);

}



void elastmat2D::PostProcess( const MatDoub & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly)
{
	int type = 1;
	shapequad objshapes(fOrder, type);
	std::vector<std::vector< std::vector<Doub > > > alco = fmesh.GetAllCoords();
	MatDoub elcoords, eltopology, psis, gradpsis, xycoords, psist;
	GetElCoords(alco, 0, elcoords);
	Int rows = elcoords.nrows();
	Int cols = rows;
	Int nels = fmesh.GetAllCoords().size();
	Doub refine = 0.1;

	for (Int iel = 0;iel < nels;iel++)
	{
		GetElCoords(alco, iel, elcoords);
		for (Doub xi = -1.;xi < 1 - refine;xi += refine)
		{
			std::vector<double> sol(3);
			for (Doub eta = -1.; eta < 1 - refine;eta += refine)
			{
				Doub approx = 0., approy = 0.;
				objshapes.shapes(psis, gradpsis, xi, eta);
				psis.Transpose(psist);
				psist.Mult(elcoords, xycoords);
				sol[0] = xycoords[0][0];
				sol[1] = xycoords[0][1];
				for (Int inode = 0;inode < elcoords.nrows();inode++)
				{
					approx += psis[inode][0] * nodalsol[fmesh.GetMeshTopology()[iel][inode] * 2][0];
					approy += psis[inode][0] * nodalsol[fmesh.GetMeshTopology()[iel][inode] * 2 + 1][0];
				}
				sol[2] = approx;
				solx.push_back(sol);
				sol[2] = approy;
				soly.push_back(sol);

			}
		}
	}
}

void elastmat2D::PostProcess( const MatDoub & nodalsol, std::vector<std::vector<double>> &sol)
{
	int type = 1;
	shapequad objshapes(fOrder, type);
std::vector<std::vector< std::vector<Doub > > > alco = fmesh.GetAllCoords();
	MatDoub elcoords, eltopology, psis, gradpsis, xycoords, psist;
	GetElCoords(alco, 0, elcoords);
	Int rows = elcoords.nrows();
	Int cols = rows;
	Int nels = fmesh.GetAllCoords().size();
	Doub refine = 0.05;

	for (Int iel = 0;iel < nels;iel++)
	{
		GetElCoords(alco, iel, elcoords);
		for (Doub xi = -1.;xi < 1 - refine;xi += refine)
		{
			std::vector<double> soli(3);
			for (Doub eta = -1.; eta < 1 - refine;eta += refine)
			{
				Doub approx = 0., approy = 0.;
				objshapes.shapes(psis, gradpsis, xi, eta);
				psis.Transpose(psist);
				psist.Mult(elcoords, xycoords);
				soli[0] = xycoords[0][0];
				soli[1] = xycoords[0][1];
				for (Int inode = 0;inode < elcoords.nrows();inode++)
				{
					approx += psis[inode][0] * nodalsol[fmesh.GetMeshTopology()[iel][inode]][0];
				}
				soli[2] = approx;
				sol.push_back(soli);

			}
		}
	}
}
