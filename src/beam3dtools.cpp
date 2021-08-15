#include "beam3dtools.h"

beam3dtools::beam3dtools()
{

}

beam3dtools::beam3dtools(mesh* inmesh)
{
    fmesh = inmesh;
}

beam3dtools::~beam3dtools()
{

}

Doub  beam3dtools::computelamda(MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l)
{
    Doub dlamb=0.;
	Int sz = dwb.nrows();
	Doub aa = 0.;
	for (Int i = 0; i < sz; i++)aa += dwb[i][0] * dwb[i][0];
	Doub bb = 0.;
	MatDoub dwcopy = dw;
	dwcopy += dws;
	for (Int i = 0; i < sz; i++)bb += dwb[i][0] * dwcopy[i][0];
	bb *= 2;
	Doub cc = 0.;
	for (Int i = 0; i < sz; i++)cc += dwcopy[i][0] * dwcopy[i][0];

	cc -= l * l;
	Doub delta = bb * bb - 4. * aa * cc;

    if(delta<0)
    {
        std::cout<< "a = "<< aa << std::endl;
        std::cout<< "b = "<< bb << std::endl;
        std::cout<< "c = "<< cc << std::endl;
        std::cout<< "delta = "<< delta << std::endl;
        std::cout << " deta negativo. "<<std::endl;
        DebugStop();
    }
	dlamb = (-bb + sqrt(delta)) / (2. * aa);
	return dlamb;


}

void beam3dtools::IterativeProcess( )
{
   MatDoub KG,Fint,Fbody,FG;
    elastoplastic3D< vonmises > mat0;
    mesh mesh0;
    CreateMatAndMesh(mesh0,mat0);

    mesh0.Assemble(KG,Fint, Fbody);
    int sz = KG.nrows();
    FG.assign(sz,1,0.);
    LoadBC(mesh0,  KG,  FG);
    InsertBC(mesh0,  KG,  FG);
    NRmatrix<Doub> u;
    SolveEigen( KG, FG, u);

    //u.Print();
    std::cout << "Displace = " <<u[127*3-1][0] << std::endl;

    u.assign(sz,1,0.);
	Doub finalload = 2.;
	Doub fac[] = { 0.2,0.4,0.6,0.8,1. };

	Int steps = 5;
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

            //dw.Print();
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

	}

	std::ofstream file8("loadvsdisplacementlu.txt");
	OutPutFile(solpost, file8);
	std::vector<std::vector<double>> solx, soly;
    std::vector<std::vector<std::vector<Doub>>> allcoords = mesh0.GetAllCoords();
    MatInt meshtopology=mesh0.GetMeshTopology();
    MatDoub meshnodes=mesh0.GetMeshNodes();
	mesh0.fmaterial->PostProcess(allcoords,meshnodes,meshtopology, u, solx, soly);
	std::ofstream file("soly.txt");
	OutPutPost(soly, file);

}
void beam3dtools::CreateMatAndMeshCube(mesh&getmesh, material &mat)
{
    string nodestr = "/home/diogo/projects/dcproj/one-el-p2-nodes-v2.txt";
	string elsstr = "/home/diogo/projects/dcproj/one-el-p2-els.txt";

	MatDoub hhatinho;
	MatDoub  meshcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

    int dim=3;
    MatDoub KG,Fint,Fbody;
    NRmatrix<Doub> bodyforce(3,1);
    bodyforce.assign(3,1,0.);
    bodyforce[1][0]=0.;
    Int order=2;

    MatDoub ptsweigths;
	shapehexahedron shape = shapehexahedron(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows() * npts;
	Int sz = dim * meshcoords.nrows();

    elastoplastic3D< vonmises >* mat0 = new elastoplastic3D< vonmises >(bodyforce,order);

	mat0->fYC.setup(20000.,0.0, 21000.);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);

    mesh* mesh0 = new mesh(dim,mat0, allcoords, meshcoords, meshtopology);

    getmesh=*mesh0;
    mat =*mat0;
}

void beam3dtools::CreateMatAndMesh(mesh &getmesh, material &mat)
{
    string nodestr = "/home/diogo/projects/dcproj/beam3D-nodes.txt";
	string elsstr = "/home/diogo/projects/dcproj/beam3D-elements.txt";

	MatDoub hhatinho;
	MatDoub  meshcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

    int dim=3;
    MatDoub KG,Fint,Fbody;
    NRmatrix<Doub> bodyforce(3,1);
    bodyforce.assign(3,1,0.);
    bodyforce[1][0]=0.;
    Int order=2;

    MatDoub ptsweigths;
	shapehexahedron shape = shapehexahedron(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows() * npts;
	Int sz = dim * meshcoords.nrows();

    elastoplastic3D< vonmises >* mat0 = new elastoplastic3D< vonmises >(bodyforce,order);

	mat0->fYC.setup(21000000.,0.3, 21000.);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);

    mesh* mesh0 = new mesh(dim,mat0, allcoords, meshcoords, meshtopology);
    getmesh=*mesh0;
    mat =*mat0;


}




void beam3dtools::InsertBC(mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F)
{
    NRvector<double> constcoorddata(3,0.);
    constcoorddata[0]=0.;
    constcoorddata[1]=0.;
    constcoorddata[2]=0.2;
    int constcoord=0;
    std::vector<int> idsface;

    std::vector<std::vector<std::vector<Doub>>> allcoords = mesh0.GetAllCoords();
    MatInt meshtopology=mesh0.GetMeshTopology();
    //FindIdsInFace(constcoorddata, constcoord, allcoords ,meshtopology , idsface);

    NRvector<int>constcoord2(3);
    constcoord2[0]=1;//fixo x
    constcoord2[1]=0;//fixo y
    constcoord2[2]=0;//livre z
    FindIds(constcoorddata, constcoord2, allcoords ,meshtopology , idsface);

    int nids=idsface.size();
    //for (auto& x : foo().items()) { /* .. */ }
   // for( int i = 0;i< nids;i++ )std::cout<<idsface[i]<<std::endl;

    int dir =0;
    double val =0.;
    mesh0.fmaterial->DirichletBC(K,F, idsface, dir,val);

    dir =1;
    mesh0.fmaterial->DirichletBC(K,F, idsface, dir,val);

    dir =2;
    mesh0.fmaterial->DirichletBC(K,F, idsface, dir,val);

}



void beam3dtools::LoadBC(mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F)
{

    F.assign(F.nrows(),1,0.);
    F[127*3-1][0]=-2.;
}

void beam3dtools::InsertCubeBC(mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F)
{
    NRvector<double> constcoorddata(3,0.);
    constcoorddata[0]=0.;
    constcoorddata[1]=0.;
    constcoorddata[2]=0.;
    int constcoord=0;
    std::vector<int> idsface;

    std::vector<std::vector<std::vector<Doub>>> allcoords = mesh0.GetAllCoords();
    MatInt meshtopology=mesh0.GetMeshTopology();
    //FindIdsInFace(constcoorddata, constcoord, allcoords ,meshtopology , idsface);

    NRvector<int>constcoord2(3);
    constcoord2[0]=0;//fixo x
    constcoord2[1]=0;//fixo y
    constcoord2[2]=1;//livre z
    FindIds(constcoorddata, constcoord2, allcoords ,meshtopology , idsface);

    int nids=idsface.size();
    //for (auto& x : foo().items()) { /* .. */ }
    for( int i = 0;i< nids;i++ )std::cout<<idsface[i]<<std::endl;

    int dir =0;
    double val =0.;
    mesh0.fmaterial->DirichletBC(K,F, idsface, dir,val);

    dir =1;
    mesh0.fmaterial->DirichletBC(K,F, idsface, dir,val);

    dir =2;
    mesh0.fmaterial->DirichletBC(K,F, idsface, dir,val);

    F.assign(F.nrows(),1,0.);
    //{3.33333, 3.33333, 3.33333, 3.33333, -13.3333, -13.3333, -13.3333, -13.3333}
    F[3*5-1][0]=3.333333;
    F[3*6-1][0]=3.333333;
    F[3*7-1][0]=3.333333;
    F[3*8-1][0]=3.333333;
    F[3*17-1][0]=-13.3333;
    F[3*18-1][0]=-13.3333;
    F[3*19-1][0]=-13.3333;
    F[3*20-1][0]=-13.3333;
}


void beam3dtools::SolveElasticCube()
{
    MatDoub KG,Fint,Fbody;
    elastoplastic3D< vonmises > mat0;
    mesh mesh0;
    CreateMatAndMeshCube(mesh0,mat0);

    mesh0.Assemble(KG,Fint, Fbody);

    InsertCubeBC(mesh0,  KG,  Fint);

    NRmatrix<Doub> u;
    SolveEigen( KG, Fint, u);
   // u.Print();
    std::cout << "Displace = " <<u[5*3-1][0] << std::endl;

   // Int dimension=3;
   // VTKGraphMesh vtkobj(&mesh0,dimension,u);

   // vtkobj.DrawSolution( 1,  0.);

}

void beam3dtools::SolveElasticBeam()
{
    MatDoub KG,Fint,Fbody;
    elastoplastic3D< vonmises > mat0;
    mesh mesh0;
    CreateMatAndMesh(mesh0,mat0);

    mesh0.Assemble(KG,Fint, Fbody);

    InsertBC(mesh0,  KG,  Fint);
    LoadBC(mesh0,  KG,  Fint);

    NRmatrix<Doub> u;
    SolveEigen( KG, Fint, u);

    //string beamstr="beam";
    //VTKGraphMesh vtkobj(&mesh0,3,u,beamstr);

    //vtkobj.DrawSolution( 1,  0.);

    u.Print();
    std::cout << "Displace = " <<u[127*3-1][0] << std::endl;
}
void beam3dtools::SolveEigen(MatDoub A, MatDoub b, MatDoub& x)
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
	VectorXd xxx = AA.lu().solve(bbb);//mais rápido
	//VectorXd xxx = AA.fullPivHouseholderQr().solve(bbb);
	//VectorXd xxx = AA.fullPivLu().solve(bbb);
	//VectorXd xxx = AA.ldlt().solve(bbb);
	//VectorXd xxx = AA.lu().solve(bbb);
	for (int i = 0; i < A.nrows(); i++)x[i][0] = xxx(i);
}
void beam3dtools::SolveEigenSparse(MatDoub A, MatDoub b, MatDoub& x)
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

void beam3dtools::ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord)
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
	//meshtopology.Print();

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
	//meshcoords.Print();


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
std::vector<T> beam3dtools::vecstr_to_vec(std::vector<string> vs)
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


void  beam3dtools::FindIds(NRvector<double> constcoorddata,NRvector<int> constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& ids)
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


void beam3dtools::GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords)
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



