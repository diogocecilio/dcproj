#include "gridmesh.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector


gridmesh::gridmesh(Doub L, Doub h, Int nx, Int ny, Int order)
{
	fL = L;
	fh = h;
	fnx = nx;
	fny = ny;
	forder = order;
}


gridmesh::~gridmesh()
{
}


void gridmesh::CreateMesh(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology)
{

	if (forder == 2)
	{
		Doub dx = fL / (2 * fnx);
		Doub dy = fh / (2 * fny);
		std::vector< std::vector<Doub> > meshnodes;
		Doub x = 0, y = 0;
		if (forder == 2)
		{
			for (Int i = 0; i < 2 * fny + 1; i++)
			{
				if (i % 2 == 0)
				{
					for (Int j = 0; j < 2 * fnx + 1; j++)
					{
						std::vector<Doub> co(2);
						co[0] = x;
						co[1] = y;
						meshnodes.push_back(co);
						x += dx;
					}
				}
				else
				{
					for (Int k = 0; k < fnx + 1; k++)
					{
						std::vector<Doub> co(2);
						co[0] = x;
						co[1] = y;
						meshnodes.push_back(co);
						x += 2 * dx;
					}
				}
				x = 0;
				y += dy;

			}

		}
		int nnodes = meshnodes.size();
		meshcoords.resize(nnodes, 2);

		for (int i = 0; i < meshnodes.size(); i++)
		{
			meshcoords[i][0] = meshnodes[i][0];
			meshcoords[i][1] = meshnodes[i][1];
		}

		//meshcoords.Print();


		std::vector< std::vector<Int> > temp;
		//meshtopology = {};
		Int b = 0;
		Int a = 1;
		Int l = 0;
		Int c = 3 * fnx + 2;
		std::vector<Int> data(8);
		for (Int i = 0;i < fny;i++)
		{
			for (Int j = 0;j < fnx;j++)
			{
				data[0] = a;
				data[1] = a + 2;
				data[2] = 3 * fnx + 4 + a;
				data[3] = 3 * fnx + 3 + b;
				data[4] = a + 1;
				data[5] = 2 * fnx + 3 + l;
				data[6] = 3 * fnx + 4 + b;
				data[7] = 2 * fnx + 2 + l;
				temp.push_back(data);
				a += 2;
				b += 2;
				l += 1;
			}
			l = 3 * fnx + 2 + c*(i);
			a = 3 * fnx + 3 + c*(i);
			b = 3 * fnx + 2 + c*(i);
		}
		nnodes = temp.size();
		meshtopology.resize(nnodes, forder * 4);

		for (int i = 0; i < nnodes; i++)
		{
			for (Int j = 0;j < forder * 4;j++)
			{
				meshtopology[i][j] = temp[i][j] - 1;
			}
		}

		std::vector<Doub> temp3(2);
		for (Int i = 0; i < meshtopology.nrows();i++)
		{
			std::vector< std::vector<Doub> > temp2;
			for (Int j = 0; j < meshtopology.ncols(); j++)
			{
				Int top = meshtopology[i][j];
				temp3[0] = meshcoords[top][0];
				temp3[1] = meshcoords[top][1];
				temp2.push_back(temp3);
			}
			allcoords.push_back(temp2);
		}

	}
	else {
		Int topolsz = 4;
		Doub dx = fL / (fnx - 1.);
		Doub dy = fh / (fny - 1.);
		Doub x = 0, y = 0;
		std::vector< std::vector<Doub> > meshnodes;
		std::vector<Doub> co(2);
		//for (Int i = 0;i < fny;i++)
		//{
		//	//y = 0;
		//	for (Int j = 0;j < fnx;j++)
		//	{
		//		 co[0] = x;
		//		 co[1] = y;
		//		 meshnodes.push_back(co);
		//		 x += dx;
		//	}
		//	y += dy;
		//	dx *= -1;
		//	x += dx;
		//}

		for (Int i = 0;i < fnx;i++)
		{
			y = 0;
			for (Int j = 0;j < fny;j++)
			{
				co[0] = x;
				co[1] = y;
				meshnodes.push_back(co);
				y += dy;
			}
			x += dx;
		}

		int nnodes = meshnodes.size();
		meshcoords.resize(nnodes, 2);

		for (int i = 0; i < meshnodes.size(); i++)
		{
			meshcoords[i][0] = meshnodes[i][0];
			meshcoords[i][1] = meshnodes[i][1];
		}


		std::vector<Int> data(4);
		std::vector< std::vector<Int> > temp;
		for (Int j = 0;j <= fny - 2;j++)
		{
			for (Int i = 1;i < fnx*fny - fny;i += fny)
			{

				/*			data[0] = i+j;
				data[3] = i+j+fny;
				data[2] = i + j + fny + 1;
				data[1] = i + j + 1;*/
				data[0] = i + j;
				data[1] = i + j + fny;
				data[2] = i + j + fny + 1;
				data[3] = i + j + 1;
				temp.push_back(data);
			}
		}
		//if (fny == 2)
		//{
		//	Int i = 1, j = 0;
		//	//data[0] = i + j;
		//	//data[3] = i + j + fny;
		//	//data[2] = i + j + fny + 1;
		//	//data[1] = i + j + 1;
		//	data[0] = i + j;
		//	data[1] = i + j + fny;
		//	data[2] = i + j + fny + 1;
		//	data[3] = i + j + 1;
		//	temp.push_back(data);
		//}

		nnodes = temp.size();
		meshtopology.resize(nnodes, forder * 4);

		for (int i = 0; i < nnodes; i++)
		{
			for (Int j = 0;j < forder * 4;j++)
			{
				meshtopology[i][j] = temp[i][j] - 1;
			}
		}

		std::vector<Doub> temp3(2);
		for (Int i = 0; i < meshtopology.nrows();i++)
		{
			std::vector< std::vector<Doub> > temp2;
			for (Int j = 0; j < meshtopology.ncols(); j++)
			{
				Int top = meshtopology[i][j];
				temp3[0] = meshcoords[top][0];
				temp3[1] = meshcoords[top][1];
				temp2.push_back(temp3);
			}
			allcoords.push_back(temp2);
		}
		/*meshtopology =
		Flatten[Table[
		Table[{i + j, i + j + ny, i + j + ny + 1, i + j + 1}, { i, 1,
		nx ny - ny, ny }], { j, 0, ny - 2 }], 1];8*/
	}

	fallcoords = allcoords;
	fmeshcoords = meshcoords;
	fmeshtopology = meshtopology;


}

void gridmesh::PrintAllCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, MatInt  meshtopology)
{
	for (Int i = 0; i < meshtopology.nrows(); i++)
	{
		//std::cout << "\n" << std::endl;
		for (Int j = 0; j < meshtopology.ncols(); j++)
		{
			Doub x = allcoords[i][j][0];
			Doub y = allcoords[i][j][1];
			//std::cout << " [ " << x << " " << y << " ] , " << std::endl;
		}
	}
}

 void gridmesh::GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords)
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

 void gridmesh::FindIdsInPath(const MatDoub & path, std::vector<std::vector< std::vector<Doub > > > &allcoords, MatInt & meshtopology, std::vector<int> & idpath)
{
	MatDoub elcoords;
	Int nels = allcoords.size();
	GetElCoords(allcoords, 0, elcoords);
	Int nnodes = elcoords.nrows();
	for (Int iel = 0;iel < nels;iel++)
	{
		GetElCoords(allcoords, iel, elcoords);
		for (Int inode = 0;inode < nnodes;inode++)
		{
			Doub x = elcoords[inode][0];
			Doub y = elcoords[inode][1];

			for (Int ipath = 0;ipath < path.nrows();ipath++)
			{
				Doub copathx = path[ipath][0];
				Doub copathy = path[ipath][1];

				if (fabs(x - copathx) < 10.e-4 && fabs(y - copathy) < 10.e-4)
				{
					idpath.push_back(meshtopology[iel][inode]);
					ipath = path.nrows();
				}
			}
		}

	}

	sort(idpath.begin(), idpath.end());
	idpath.erase(unique(idpath.begin(), idpath.end()), idpath.end());

}

void gridmesh::Line(VecDoub a, VecDoub b, Int ndivs, MatDoub & path)
{
	Doub x0 = a[0];
	Doub xf = b[0];

	Doub y0 = a[1];
	Doub yf = b[1];

	Doub dx = (xf - x0) / ndivs;
	Doub dy = (yf - y0) / ndivs;

	path.assign(ndivs, 2, 0.);

	for (Int idiv = 0;idiv < ndivs;idiv++)
	{
		path[idiv][0] = x0 + idiv * dx;
		path[idiv][1] = y0 + idiv * dy;
	}


}

void gridmesh::PrintGMeshVTK(std::vector<std::vector< std::vector<Doub > > >  allcoords, MatInt meshtopology, std::ofstream &file)
{
	file.clear();
	int nelements = allcoords.size();

	MatDoub elcoords;

	std::stringstream node, connectivity, type;

	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;

	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";

	int actualNode = -1, size = 0, nVALIDelements = 0;

	for (int el = 0; el < nelements; el++)
	{
		GetElCoords(allcoords, el, elcoords);
		Int elNnodes = elcoords.nrows();
		size += (1 + elNnodes);
		connectivity << elNnodes;
		for (int t = 0; t < elNnodes; t++)
		{
			double x = elcoords[t][0];
			double y = elcoords[t][1];
			node << x << " " << y << " " << 0;
			node << std::endl;
			actualNode++;
			connectivity << " " << meshtopology[el][t];
		}
		connectivity << std::endl;
		Int elType = 23;

		type << elType << std::endl;
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();


	file << "CELLS " << nVALIDelements << " ";

	file << size << std::endl;
	file << connectivity.str() << std::endl;

	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();

	file.close();
}

