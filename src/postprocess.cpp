#include "postprocess.h"



postprocess::postprocess()
{
}


postprocess::~postprocess()
{
}



void postprocess::PostProcess(mesh &inmesh,const MatDoub & nodalsol, std::vector<std::vector<double>> &sol)
{

	std::vector<std::vector< std::vector<Doub > > > allcoords = inmesh.GetAllCoords();
	MatDoub meshnodes = inmesh.GetMeshNodes();
	MatInt meshtopology = inmesh.GetMeshTopology();
	MatDoub elcoords, eltopology, psis, gradpsis, xycoords, psist;
	inmesh.GetElCoords(allcoords, 0, elcoords);
	Int rows = elcoords.nrows();
	Int cols = rows;
	Int nels = allcoords.size();
	Doub refine = 0.05;
	shapequad shape = shapequad(2, 1);
	for (Int iel = 0;iel < nels;iel++)
	{
		inmesh.GetElCoords(allcoords, iel, elcoords);
		for (Doub xi = -1.;xi < 1 - refine;xi += refine)
		{
			std::vector<double> soli(3);
			for (Doub eta = -1.; eta < 1 - refine;eta += refine)
			{
				Doub approx = 0., approy = 0.;
				shape.shapes(psis, gradpsis, xi, eta);
				psis.Transpose(psist);
				psist.Mult(elcoords, xycoords);
				soli[0] = xycoords[0][0];
				soli[1] = xycoords[0][1];
				for (Int inode = 0;inode < elcoords.nrows();inode++)
				{
					approx += psis[inode][0] * nodalsol[meshtopology[iel][inode]][0];
				}
				soli[2] = approx;
				sol.push_back(soli);

			}
		}
	}
}
