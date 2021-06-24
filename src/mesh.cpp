
//#include "mins_ndim.h"
#include <stdio.h>
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include "mesh.h"

mesh::mesh(std::vector<std::vector< std::vector<Doub > > > &allcoords, NRmatrix<Doub> &meshnodes, NRmatrix<Int> &meshtopology)
{
	fallcoords = allcoords;
	fmeshnodes = meshnodes;
	fmeshtopology = meshtopology;
}

mesh::mesh(material *mat, std::vector<std::vector< std::vector<Doub > > >& allcoords, NRmatrix<Doub>& meshnodes, NRmatrix<Int>& meshtopology)
{
	fmaterial = mat;
	fallcoords = allcoords;
	fmeshnodes = meshnodes;
	fmeshtopology = meshtopology;


}


mesh::mesh(material* mat, std::vector<std::vector< std::vector<Doub > > >& allcoords, NRmatrix<Doub>& meshnodes, NRmatrix<Int>& meshtopology,MatDoub  HHAT)
{
	fmaterial = mat;
	fallcoords = allcoords;
	fmeshnodes = meshnodes;
	fmeshtopology = meshtopology;
	fHHAT = HHAT;

}

mesh::mesh()
{
}

mesh::mesh(mesh &copy)
{
	copy = mesh(fallcoords, fmeshnodes, fmeshtopology);
}

mesh::~mesh()
{
    delete fmaterial;
}

std::vector<std::vector< std::vector<Doub > > >  mesh::GetAllCoords()
{
	int szdebug = fallcoords.size();
	return fallcoords;
}

MatDoub mesh::GetMeshNodes()
{
	return fmeshnodes;
}

MatInt mesh::GetMeshTopology()
{
	return fmeshtopology;
}

//material mesh::GetMaterial()
//{
//	return *fmaterial;
//}


void mesh::GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords)
{

	elcoords.assign(allcoords[el].size(), 2, 0.);
	Int sz = allcoords[el].size();
	for (Int j = 0; j <sz; j++)
	{
		Doub x = allcoords[el][j][0];
		Doub y = allcoords[el][j][1];
		elcoords[j][0] = x;
		elcoords[j][1] = y;
	}
}

MatDoub mesh::FindSolution(VecDoub coord, MatDoub datatosearch)
{
	//MatDoub datatosearch [ data]

	MatDoub solout,returnsol(1,3,0.);
	Int el = 0;
	Doub xi = 0.;
	Doub eta = 0.;
	MatDoub psis, GradPsi, elcoords, psist, solel, xycoords, sol;
	shapequad shape = shapequad(2, 1);
	shape.shapes(psis, GradPsi, xi, eta);
	Doub l0 = 200;
	Doub l = 200;
	std::vector<int> possible;
	int counter = 0;
	bool  breaktrue = false;
	while (counter<40)
	{
		for (int iel = 0;iel < fallcoords.size();iel++)
		{
			GetElCoords(fallcoords, iel, elcoords);
			Int nodes = psis.nrows();
			Int nstatevars = 1;
			solel.assign(nodes, 1, 0);
			psis.Transpose(psist);
			psist.Mult(elcoords, xycoords);
			Doub dx = xycoords[0][0] - coord[0];
			Doub dy = xycoords[0][1] - coord[1];
			Doub dist = sqrt(dx*dx + dy*dy);
			if (dist<l0)
			{
				possible.push_back(iel);
			}
		}
		//cout << "l0 = " << l0 << endl;
		l0 *= 0.5;
		counter++;
	}
	std::cout << "in co  = " << endl;
	coord.Print();
	el = possible.size() - 1;
	std::cout << "el = " <<possible[el] << endl;
	GetElCoords(fallcoords, possible[el], elcoords);
	std::cout << "elcoords = " <<  endl;
	elcoords.Print();
	MatDoub elnodesol(elcoords.nrows(), 1, 0.);
	for (Int inode = 0;inode < elcoords.nrows();inode++)
	{
		elnodesol[inode][0] = datatosearch[ fmeshtopology[possible[el]][inode] ][0];
	}
	std::cout << "elnodesol = " << endl;
	elnodesol.Print();

	VecDoub p(2, -0.);
	FuncdSearch func2;
	func2.set(coord, elcoords);
	//Linemethod<FuncdSearch> line(func2);
	//Doub min = line.linmin();
	//std::cout << "min = " << min  << endl;


	//Powell<FuncdSearch> powell(func2);


	for (Doub xi = -1.;xi <= 1.;xi += 0.0001)
	{
		for (Doub eta = -1.;eta <= 1.;eta += 0.0001)
		{
			shape.shapes(psis, GradPsi, xi, eta);
			psis.Transpose(psist);
			psist.Mult(elcoords, xycoords);
			Doub dx = xycoords[0][0] - coord[0];
			Doub dy = xycoords[0][1] - coord[1];
			Doub dist = sqrt(dx*dx + dy*dy);
			if (dist < 0.1)
			{
				p[0] = xi;
				p[1] = eta;
				break;
			}

		}
	}


	//p = powell.minimize(p);

	std::cout << "solu = " << endl;
	//p.Print();

	shape.shapes(psis, GradPsi, p[0], p[1]);
	psis.Transpose(psist);
	psist.Mult(elcoords, xycoords);
	psist.Mult(elnodesol, solout);
	returnsol[0][0] = xycoords[0][0];
	returnsol[0][1] = xycoords[0][1];
	returnsol[0][2] = solout[0][0];
	returnsol.Print();
	return returnsol;


	
	
}

MatDoub mesh::TransferSolution( mesh &out, MatDoub datatosearch)
{
	Int outsize = out.fmeshnodes.nrows();//Um valor para cada nó
	MatDoub outsol(outsize,1,0.);// newval na malha coarse
	Int outnels = out.fallcoords.size();
	MatDoub outelcoords,returnsol;
	VecDoub co(2,0.);
	//preciso fornecer para a malha fina a coordenada do no da malha grossa para que a busca na malha fina seja realizada
	//portanto preciso iterar na malha grossa out
	for (int iel = 0;iel < outnels;iel++)
	{
		GetElCoords(fallcoords, iel, outelcoords);
		for (Int inode = 0;inode < outelcoords.nrows();inode++)
		{
			co[0] = outelcoords[inode][0];
			co[1] = outelcoords[inode][1];
			returnsol = FindSolution(co, datatosearch); //(acha coordena x y e val na malha fina a partir da coordenada 
			//fornecido pela malha grossa)
			outsol[out.fmeshtopology[iel][inode]][0] = returnsol[0][2];
		}
		
	}
	return outsol;
}


void mesh::Assemble(MatDoub& KG, MatDoub& Fint, MatDoub& Fbody)
{
	//std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
	//MatDoub meshnodes = fmesh.GetMeshNodes();
	//MatInt meshtopology = fmesh.GetMeshTopology();

	//cout << "all cc" << allcoords[0].size() <<endl;

	MatDoub ek, efint, efbody, elcoords, eltopology;
	GetElCoords(fallcoords, 0, elcoords);
	Int nnodes = fmeshnodes.nrows();
	Int rows = elcoords.nrows();
	Int sz = 2 * nnodes;
	Int cols = rows;
	KG.assign(sz, sz, 0.);
	Fint.assign(sz, 1, 0.);
	Fbody.assign(sz, 1, 0.);
	Int nels = fallcoords.size();

	//uglob = Table[
		//Table[{displacement[[2 topol[[k, j]] - 1]],
		//	displacement[[2 topol[[k, j]]]]}, { j, 1,
			//Length[topol[[k]]] }], { k, 1, nels }];
	NRmatrix<NRvector<Doub>> uglob;
	uglob.resize(nels, nnodes);
	for (int i = 0; i < nels; i++)
	{
		for (int j = 0; j < nnodes; j++) {
			uglob[i][j].assign(2, 0.);
		}
	}

	for (Int iel = 0; iel < nels; iel++)
	{
		for (Int node = 0; node < rows; node++)
		{
			uglob[iel][node][0] = fmaterial->GetSolution()[2 * fmeshtopology[iel][node]][0];
			uglob[iel][node][1] = fmaterial->GetSolution()[2 * fmeshtopology[iel][node] + 1][0];
		}
	}

	//	uglob.Print2();

	Int fu = 0;
	for (Int iel = 0; iel < nels; iel++)
	{
		if (fHHAT.nrows() != 0)
		{
			fhhatvel.resize(fHHAT.ncols());
			for (Int ivar = 0; ivar < fHHAT.ncols(); ivar++) {
				fhhatvel[ivar].assign(rows, 1, 0.);
				for (Int inode = 0; inode < rows; inode++)
				{
					fhhatvel[ivar][inode][0] = fHHAT[fmeshtopology[iel][inode]][ivar];
				}
			}
		}
		fmaterial->SetRandomField(fHHAT);
		fmaterial->SetRandomFieldLocal(fhhatvel);
		MatDoub elementdisplace(elcoords.nrows(), 2, 0.);
		for (Int i = 0; i < elcoords.nrows(); i++)for (Int j = 0; j < 2; j++)elementdisplace[i][j] = uglob[iel][i][j];
		GetElCoords(fallcoords, iel, elcoords);
		fmaterial->CacStiff(ek, efint, efbody, elcoords, elementdisplace);
		for (Int irow = 0; irow < rows; irow++)
		{
			Int rowglob = fmeshtopology[iel][irow];
			for (Int icol = 0; icol < cols; icol++)
			{
				Int colglob = fmeshtopology[iel][icol];
				KG[2 * rowglob + fu][2 * colglob + fu] += ek[2 * irow + fu][2 * icol + fu];
				KG[2 * rowglob + fu][2 * colglob + 1 + fu] += ek[2 * irow + fu][2 * icol + 1 + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + fu] += ek[2 * irow + 1 + fu][2 * icol + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + 1 + fu] += ek[2 * irow + 1 + fu][2 * icol + 1 + fu];

			}
			Fbody[2 * rowglob + fu][0] += efbody[2 * irow + fu][0];
			Fbody[2 * rowglob + 1 + fu][0] += efbody[2 * irow + 1 + fu][0];

			Fint[2 * rowglob + fu][0] += efint[2 * irow + fu][0];
			Fint[2 * rowglob + 1 + fu][0] += efint[2 * irow + 1 + fu][0];
		}

	}
	fmaterial->ResetCounter();
}

