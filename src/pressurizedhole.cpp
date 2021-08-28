#include "pressurizedhole.h"
#include "vtkmesh.h"
//class VTKGraphMesh;
pressurizedhole::pressurizedhole()
{

}

pressurizedhole::pressurizedhole(mesh* inmesh)
{
    fmesh = inmesh;
}

pressurizedhole::~pressurizedhole()
{

}


void  pressurizedhole::IterativeProcess()
{


    MatDoub KG,Fint,Fbody,FG;
    elastoplastic2D< vonmises > mat0;
    mesh mesh0;
    CreateMatAndMesh(mesh0,mat0);

    mesh0.Assemble(KG,Fint, Fbody);
    int sz = KG.nrows();
    LoadBC(mesh0,  KG,  FG);

    NRmatrix<Doub> u;
    SolveEigen( KG, Fint, u);

    //u.Print();
    //std::cout << "Displace = " <<u[127*3-1][0] << std::endl;

    u.assign(sz,1,0.);
	Doub finalload = 0.19209;
	//Doub fac[] = { 0.1 / finalload, 0.14 / finalload, 0.18 / finalload, 0.19 / finalload, 1. };
    Doub fac[] = {0.019209/ finalload, 0.038418/ finalload, 0.057627/ finalload, 0.076836/ finalload, 0.096045/ finalload, 0.115254/ finalload, 0.134463/ finalload, 0.153672/ finalload, 0.172881/ finalload, 0.19209/ finalload};
	Int steps = 10;
	Int counterout = 1;
	MatDoub solpost(1000, 2, 0.);
	for (Int iload = 0; iload < steps; iload++)
	{
		std::cout << "load step = " << iload << std::endl;
		Int counter = 0, maxcount = 30;
		Doub err1 = 10., err2 = 10., tol = 10.e-5;
		MatDoub dw(sz, 1, 0.), res(sz, 1, 0.), FINT,FBODY, R;
		while (counter <  maxcount && err1 > tol)
		{
			MatDoub FGint = FG;
			mesh0.Assemble(KG,Fint, Fbody);

			//KG.Print();
			//FINT.Print();

			FGint *= fac[iload];
			FGint -= Fint;
			R = FGint;

            InsertBC(mesh0,  KG,  R);

            SolveEigen( KG, R, dw);

			u += dw;

			mesh0.fmaterial->UpdateDisplacement(u);

			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;

			rnorm = R.NRmatrixNorm();

			normdw = dw.NRmatrixNorm();

			normfg = FG.NRmatrixNorm();

			unorm = u.NRmatrixNorm();

			err1 = rnorm / (fac[iload] *normfg);

			err2 = normdw / unorm;

			std::cout << " Iteration number = " << counter << " |  |R|/|FE| = " << err1 << " | deltau/u " << err2 << std::endl;
			counter++;
		}
		mesh0.fmaterial->UpdatePlasticStrain();
		counterout++;
		solpost[iload][0] = fabs(u[0][0]);
		solpost[iload][1] = fabs(fac[iload] * finalload);
        std::vector<string> scalar_names;
        std::vector<string> vector_names;
        //  TPZStack<std::string> scalar_names,vector_names, tensor_names;
        vector_names.push_back("Displacement");
        vector_names.push_back("Strain");
        //vector_names.push_back("SqrtJ2(EPSP)");
        //vector_names.push_back("Stress");
        Int dim=2;
        string slopestr="ring";
        VTKGraphMesh vtkobj(&mesh0,dim,scalar_names,vector_names,slopestr);
        vtkobj.DrawSolution( counterout, counter);

	}

	std::ofstream file8("loadvsdisplacementlu.txt");
	OutPutFile(solpost, file8);
	std::vector<std::vector<double>> solx, soly;
    std::vector<std::vector<std::vector<Doub>>> allcoords = mesh0.GetAllCoords();
    MatInt meshtopology=mesh0.GetMeshTopology();
    MatDoub meshnodes=mesh0.GetMeshNodes();
	mesh0.fmaterial->PostProcess(allcoords,meshnodes,meshtopology, u, solx, soly);
    NRmatrix<Doub> sol,dsol;
    mesh0.fmaterial->ComputeSolAndDSol(&mesh0,sol,dsol);

	std::ofstream file("soly.txt");
	OutPutPost(soly, file);


}


void pressurizedhole::CreateMatAndMesh(mesh &getmesh, material &mat)
{
    string nodestr = "/home/diogocecilio/projects/dcproj/nodes-pressure-fino.txt";
	string elsstr = "/home/diogocecilio/projects/dcproj/elements-pressure-fino.txt";

	MatDoub hhatinho;
	MatDoub  meshcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

    std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);
	std::ofstream filemesh2("meshtopology.txt");
	OutPutPost(meshtopology, filemesh2);

    int dim=2;
    MatDoub KG,Fint,Fbody;
    NRmatrix<Doub> bodyforce(2,1);
    bodyforce.assign(2,1,0.);
    bodyforce[1][0]=0.;
    Int order=2;

    MatDoub ptsweigths;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows() * npts;
	Int sz = dim * meshcoords.nrows();
    Doub thickness=1.;
    Int planestress=0;

    elastoplastic2D< vonmises >* mat0 = new elastoplastic2D< vonmises >(thickness,bodyforce,planestress,order);

	mat0->fYC.setup(210.,0.3, 0.24);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);

    mesh* mesh0 = new mesh(dim,mat0, allcoords, meshcoords, meshtopology);
    getmesh=*mesh0;
    mat =*mat0;
}


void pressurizedhole::LoadBC(mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F)
{
    std::vector<std::vector<std::vector<Doub>>> allcoords = mesh0.GetAllCoords();
    MatInt meshtopology=mesh0.GetMeshTopology();
    MatInt  linetopology;
	std::vector<int> idpathcirc;
	MatDoub pathcirc;
	int ndivs = 1000;
	Doub delta;
	pathcirc.assign(ndivs + 1, 2, 0.);
	Int i = 0;
	delta = (M_PI / 2.) / (Doub(ndivs));
	for (Doub theta = 0;theta < M_PI / 2.; theta += delta) {
		pathcirc[i][0] = 100. * cos(theta);
		pathcirc[i][1] = 100. * sin(theta);
		i++;
	}

	gridmesh::FindIdsInPath(pathcirc, allcoords, meshtopology, idpathcirc);
    //for (int i = 0;i < idpathcirc.size();i++)std::cout << " ID  = " << idpathcirc[i] << endl;

	std::vector<std::vector<int>>  linetopol = LineTopology(idpathcirc, 2);
	ToMatInt(linetopol, linetopology);

    Doub pressure = 0.19209;
	mesh0.fmaterial->ContributeCurvedLine(K, F, mesh0.GetMeshNodes(), linetopology, pressure);
}


void pressurizedhole::InsertBC(mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F)
{
    NRvector<double> constcoorddata(3,0.);
    constcoorddata[0]=100.;
    constcoorddata[1]=0.;
    constcoorddata[2]=0.;
    int constcoord=0;
    std::vector<int> idsbottom,idsleft;

    std::vector<std::vector<std::vector<Doub>>> allcoords = mesh0.GetAllCoords();
    MatInt meshtopology=mesh0.GetMeshTopology();
    NRvector<int>constcoord2(3);
    constcoord2[0]=0;//livre x
    constcoord2[1]=1;//fixo y
    constcoord2[2]=1;//fixo z
    FindIds(constcoorddata, constcoord2, allcoords ,meshtopology , idsbottom);

    constcoorddata[0]=0.;
    constcoorddata[1]=100.;
    constcoorddata[2]=0.;
    constcoord2[0]=1;//livre x
    constcoord2[1]=0;//fixo y
    constcoord2[2]=1;//fixo z
    FindIds(constcoorddata, constcoord2, allcoords ,meshtopology , idsleft);
    //int nids=idsleft.size();
    //for (auto& x : foo().items()) { /* .. */ }
    //for( int i = 0;i< nids;i++ )std::cout<<idsface[i]<<std::endl;
    //for( int i = 0;i< nids;i++ )std::cout<<mesh0.GetMeshNodes()[idsleft[i]][0] << mesh0.GetMeshNodes()[idsleft[i]][1]<<std::endl;

    double val =0.;
    int dir =1;
    mesh0.fmaterial->DirichletBC(K,F, idsbottom, dir,val);

    dir =0;
    mesh0.fmaterial->DirichletBC(K,F, idsleft, dir,val);

}

void pressurizedhole::SolveElasticHole()
{
    MatDoub KG,Fint,Fbody;
    elastoplastic2D< vonmises > mat0;
    mesh mesh0;
    CreateMatAndMesh(mesh0,mat0);

    mesh0.Assemble(KG,Fint, Fbody);

    InsertBC(mesh0,  KG,  Fint);

    NRmatrix<Doub> u;
    SolveEigen( KG, Fint, u);
    u.Print();
    std::cout << "Displace = " <<u[127*3-1][0] << std::endl;
}
void pressurizedhole::SolveEigen(MatDoub A, MatDoub b, MatDoub& x)
{

	x.assign(A.nrows(), 1, 0.);
	MatrixXd AA(A.nrows(), A.nrows());
	VectorXd bbb(A.nrows());
	for (int i = 0; i < A.nrows(); i++)
	{
		for (int j = 0; j < A.ncols(); j++)
		{
			AA(i, j) = A[i][j];
		}
	}
	for (int i = 0; i < A.nrows(); i++)bbb(i) = b[i][0];
	VectorXd xxx = AA.llt().solve(bbb);//mais rápido
	for (int i = 0; i < A.nrows(); i++)x[i][0] = xxx(i);
}
void pressurizedhole::SolveEigenSparse(MatDoub A, MatDoub b, MatDoub& x)
{

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int sz=A.nrows();

    //30 = approx number of nonzero term in K matrix
    tripletList.reserve(sz*30);
   // tripletList.reserve(80000);

    x.assign(sz, 1, 0.);
	SparseMatrix<double> AA(sz, sz);
	VectorXd bbb(sz);
	for (int i = 0; i < sz; i++)
	{
		for (int j = 0; j < sz; j++)
		{
            if(fabs(A[i][j])>1.e-12)
            {
                tripletList.push_back(T(i,j,A[i][j]));
            }
		}
		bbb(i) = b[i][0];
	}

    AA.setFromTriplets(tripletList.begin(), tripletList.end());

    AA.makeCompressed();
    SimplicialLLT< SparseMatrix<double> > solver;
    VectorXd xx = solver.compute(AA).solve(bbb);
    for(int i=0;i<sz;i++)x[i][0]=xx(i);
}

void pressurizedhole::ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord)
{
	std::vector<std::vector<Int>> topol;
	string line, temp;

	ifstream myfile(filenameel);
	//
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::vector<string> tokens;
			istringstream iss(line);
			while (iss >> temp)
				tokens.push_back(temp);
			std::vector<Int> input_int = vecstr_to_vec<Int>(tokens);
			for (int k = 0; k < input_int.size(); k++)
			{
				input_int[k] = input_int[k] - 1;
			}
			topol.push_back(input_int);
		}
		myfile.close();
	}
	else std::cout << "Unable to open file";

	meshtopology.CopyFromVector(topol);


	std::vector<std::vector<Doub>> coords;
	string line2, temp2;
	ifstream myfile2(filenamecoord);
	if (myfile2.is_open())
	{
		while (getline(myfile2, line2))
		{
			std::vector<string> tokens;
			istringstream iss(line2);
			while (iss >> temp2)
				tokens.push_back(temp2);
			std::vector<Doub> input_doub = vecstr_to_vec<Doub>(tokens);

			//std::vector<Doub> input_doub2(input_doub.size() - 1);
			//for (int k = 1;k < input_doub.size();k++)
			//{
			//	input_doub2[k] = input_doub[k];
			//}

			coords.push_back(input_doub);
		}
		myfile2.close();
	}
	else std::cout << "Unable to open file";

	meshcoords.CopyFromVector(coords);



	std::vector<Doub> temp33(3);
	for (Int i = 0; i < meshtopology.nrows(); i++)
	{
		std::vector< std::vector<Doub> > temp22;
		for (Int j = 0; j < meshtopology.ncols(); j++)
		{
			Int top = meshtopology[i][j];
			temp33[0] = meshcoords[top][0];
			temp33[1] = meshcoords[top][1];
            temp33[2] = meshcoords[top][2];
			temp22.push_back(temp33);
		}
		allcoords.push_back(temp22);
	}



}


template <class T>
std::vector<T> pressurizedhole::vecstr_to_vec(std::vector<string> vs)
{
	std::vector<T> ret;
	for (std::vector<string>::iterator it = vs.begin() +1; it != vs.end() ; ++it)
	{
		istringstream iss(*it);
		T temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

void  pressurizedhole::FindIds(NRvector<double> constcoorddata,NRvector<int> constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& ids)
{
    //constcoorddata vector containig info about face to search id. It must contain any coodinate locate in the in the face
	MatDoub elcoords;
	int nels = allcoords.size();
	GetElCoords(allcoords, 0, elcoords);
	Int nnodes = elcoords.nrows();
    int sum=0;
    //constcoord.size() = 1 face
    //constcoord.size() = 2 linha
    //constcoord.size() = 3 pontos

    std::vector<int> dirs;
    for(int iconst=0;iconst<constcoord.size();iconst++)
    {
        sum+=constcoord[iconst];
        if(constcoord[iconst]==1)dirs.push_back(iconst);

    }
    for (Int iel = 0; iel < nels; iel++)
    {
        GetElCoords(allcoords, iel, elcoords);
        for (Int inode = 0; inode < nnodes; inode++)
        {

            if(sum==1)
            {
                if(fabs(elcoords[inode][dirs[0]] - constcoorddata[dirs[0]])<1.e-4  )
                {
                        ids.push_back(meshtopology[iel][inode]);
                }
            }else if (sum==2)
            {
                if(fabs(elcoords[inode][dirs[0]] - constcoorddata[dirs[0]])<1.e-4 && abs(elcoords[inode][dirs[1]] - constcoorddata[dirs[1]])<1.e-4 )
                {
                    ids.push_back(meshtopology[iel][inode]);
                }
            }
            else if (sum==3)
            {
                if(fabs(elcoords[inode][dirs[0]] - constcoorddata[dirs[0]])<1.e-4 && abs(elcoords[inode][dirs[1]] - constcoorddata[dirs[1]])<1.e-4 && abs(elcoords[inode][dirs[2]] - constcoorddata[dirs[2]])<1.e-4)
                {
                    ids.push_back(meshtopology[iel][inode]);
                }
            }


        }


    }

	sort(ids.begin(), ids.end());
	ids.erase(unique(ids.begin(), ids.end()), ids.end());
}

void pressurizedhole::GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords)
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


std::vector<std::vector<int>>  pressurizedhole::LineTopology(std::vector<int> ids, Int order)
{
	Int k = 0;

	std::vector<std::vector<int>> vg;
	for (int j = 0;j < ids.size() / order;j++)
	{
		std::vector<int> v;
		for (int i = 0;i < order + 1;i++)
		{
			v.push_back(ids[i + k]);
		}
		vg.push_back(v);
		k += order;
	}
	return vg;
}

void pressurizedhole::ToMatInt(std::vector<std::vector<int>> in, MatInt & out)
{
	Int  rows = in.size();
	Int cols = in[0].size();
	out.assign(rows, cols, 0.);
	for (Int i = 0;i < rows;i++)for (Int j = 0;j < cols;j++)out[i][j] = in[i][j];
}
