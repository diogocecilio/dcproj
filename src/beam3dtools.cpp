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
        //DebugStop();
    }
	dlamb = (-bb + sqrt(delta)) / (2. * aa);
	return dlamb;


}

std::vector<std::vector<double>>   beam3dtools::IterativeProcess( int ndesi, Doub dlamb0,Doub alphatol, int niter)
{
    elastoplastic3D< vonmises > mat0;
    mesh mesh0;
    CreateMatAndMesh(mesh0,mat0);

    Int sz = 3 * mesh0.GetMeshNodes().nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.),FEXT(sz, 1, 0.), ptsweigths;



	mesh0.fmaterial->ResetPlasticStrain();
	mesh0.fmaterial->ResetDisplacement();
	mesh0.fmaterial->ResetMat();

	MatDoub displace, displace0;
	displace.assign(sz, 1, 0.);

	Doub l = 0, l0 = 0, lamb = 1., dlamb = 0., diff=0.,diff2 = 100;
	Int counterout = 0;

	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);


	MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;

    mesh0.Assemble(KG,FINT, FBODY);

    //FBODY.Print();
    R=FEXT;
	R += FBODY;
	R *= lamb;
	R -= FINT;

	bool check = false;
	InsertBC2(mesh0,  KG,  FBODY);
    InsertBC2(mesh0,  KG,  R);
	SolveEigen(KG, R, dws);
	SolveEigen(KG, FBODY, dwb);
    dws.Print();
    dwb.Print();
	MatDoub dwbt, mult;
	dwb.Transpose(dwbt);
	dwbt.Mult(dwb, mult);
	l0 = sqrt(dlamb0 * dlamb0 * mult[0][0]);
	l = l0;
	dlamb = 0.1;
	lamb = 0.1;
	Doub scalar = 1.5;
	Doub rtol =0.001;

	cout << " \n SYSTEM SIZE  = " << KG.nrows() << std::endl;
	Doub rnorm = 10., lambn0=0.;
    Int counter = 0, maxcount = 20;

    FEXT[127 *3 -1][0]=-2.0;
	do
	{
		Doub err1 = 10., err2 = 10., tol = 1.e-5;
        if(counter<maxcount)lambn0 = lamb;
		diff = 10;
		displace0 = displace;
		rnorm = 10.;
        Doub meantime=0.;
        counter = 0;
		do
		{

            chrono::steady_clock sc;
            auto start = sc.now();     // start timer
			mesh0.Assemble(KG, FINT, FBODY);
            FBODY+=FEXT;
            R = FBODY;
			R *= lamb;
			R -= FINT;

            InsertBC2(mesh0,  KG,  FBODY);
            InsertBC2(mesh0,  KG,  R);

			SolveEigen(KG, R, dws);
			SolveEigen(KG, FBODY, dwb);

			dlamb = computelamda(dwb, dws, dw, l);
			if (isnan(dlamb) == 1) {
				std::cout << " isnan(dlamb) detected. " << endl;
                DebugStop();
			}
			if(dlamb==0)
            {
             //lamb+=0.1;
            }
            else{
			lamb += dlamb;
            }
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;

			displace += dww;
			mesh0.fmaterial->UpdateDisplacement(displace);

			rnorm = R.NRmatrixNorm();
			Doub normdw = dww.NRmatrixNorm();
			Doub unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;

            auto end = sc.now();
            auto time_span = static_cast<chrono::duration<double>>(end - start);

			std::cout << " Time in newton step :  " << time_span.count() <<" Iteration number = " << counter << " |  |R| = " << rnorm <<" | lamb  = " << lamb <<std::endl;
			counter++;
			if (counter == 1)rnorm = 10;


		} while (counter < maxcount && rnorm > rtol);
        meantime/=counter;

        std::cout << " exter iter = " << counterout << "  | newton iters = " << counter  << " |  |R| = " << rnorm << " |  lamb  = " << lamb << std::endl;

		if (rnorm > 0.5)
		{
			counterout++;
			dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

			fmesh->fmaterial->UpdateDisplacement(displace0);
			displace = displace0;
				if (check==true)
				{
					dlamb0 *= 0.1/scalar;
					scalar *= 2.;

				}
				else
				{
					dlamb0 *= 1/scalar;
					scalar *= 1.5;
				}
			dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

			mesh0.Assemble(KG, FINT, FBODY);
			R=FEXT;
            R += FBODY;
			R *= lamb;
			R -= FINT;
            InsertBC2(mesh0,  KG,  FBODY);
            InsertBC2(mesh0,  KG,  R);
			SolveEigen(KG, R, dws);
			SolveEigen(KG, FBODY, dwb);
			dwb.Transpose(dwbt);
			dwbt.Mult(dwb, mult);
			l0 = sqrt(dlamb0 * dlamb0 * mult[0][0]);
			l = l0;
			dlamb = computelamda(dwb, dws, dw, l);
			lamb = 0;
		}
		else {
			scalar = 1.;
			l *= Doub(ndesi) / Doub(counter);
			if (l > 10.)l = 10.;
			diff2 = fabs(lambn0 - lamb);
			mesh0.fmaterial->UpdatePlasticStrain();
		}

		counterout++;
		dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

		solcount[0] = 0.;
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = rnorm;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = 0.;
		uvf[1] = lamb;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);

	} while (counterout <= niter && fabs(diff2) > alphatol );


	return solpost;

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

	mat0->fYC.setup(200000.,0., 2000);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);

    mesh* mesh0 = new mesh(dim,mat0, allcoords, meshcoords, meshtopology);
    getmesh=*mesh0;
    mat =*mat0;
}


void beam3dtools::InsertBC2(mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F)
{
    NRvector<double> constcoorddata(3,0.);
    constcoorddata[0]=0.;
    constcoorddata[1]=0.;
    constcoorddata[2]=0.2;
    int constcoord=0;
    std::vector<int> idsface;

    std::vector<std::vector<std::vector<Doub>>> allcoords = mesh0.GetAllCoords();
    MatInt meshtopology=mesh0.GetMeshTopology();
    FindIdsInFace(constcoorddata, constcoord, allcoords ,meshtopology , idsface);

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
    FindIdsInFace(constcoorddata, constcoord, allcoords ,meshtopology , idsface);

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

    //MatDoub F,u;
    F.assign(F.nrows(),1,0.);
    F[127*3-1][0]=-2.;
}

void beam3dtools::SolveElasticBeam()
{
    MatDoub KG,Fint,Fbody;
    elastoplastic3D< vonmises > mat0;
    mesh mesh0;
    CreateMatAndMesh(mesh0,mat0);

    mesh0.Assemble(KG,Fint, Fbody);

    InsertBC(mesh0,  KG,  Fint);

    NRmatrix<Doub> u;
    SolveEigen( KG, Fint, u);
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
	VectorXd xxx = AA.llt().solve(bbb);//mais rápido
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
	meshtopology.Print();

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
	meshcoords.Print();


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

void beam3dtools::FindIdsInFace(NRvector<double> constcoorddata, int constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& idsface)
{
    //constcoorddata vector containig info about face to search id. It must contain any coodinate locate in the in the face
	MatDoub elcoords;
	int nels = allcoords.size();
	GetElCoords(allcoords, 0, elcoords);
	Int nnodes = elcoords.nrows();
	for (Int iel = 0; iel < nels; iel++)
	{
		GetElCoords(allcoords, iel, elcoords);
		for (Int inode = 0; inode < nnodes; inode++)
		{
			Doub x = elcoords[inode][0];
			Doub y = elcoords[inode][1];
            Doub z = elcoords[inode][2];

            if(constcoord==0)//xconstant
            {
                if (fabs(x - constcoorddata[constcoord])<1.e-4  )
                {
					idsface.push_back(meshtopology[iel][inode]);
                }
            }
            else if(constcoord==1)//yconstant
            {
                if (fabs(y - constcoorddata[constcoord])<1.e-4  )
                {
					idsface.push_back(meshtopology[iel][inode]);
                }
            }
            else if(constcoord==2)//zconstant
            {
                if (fabs(z - constcoorddata[constcoord])<1.e-4 )
                {
					idsface.push_back(meshtopology[iel][inode]);
                }
            }
		}

	}

	sort(idsface.begin(), idsface.end());
	idsface.erase(unique(idsface.begin(), idsface.end()), idsface.end());

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

    /*cout << 1%1 << endl;
    NRmatrix<Doub> bodyforce(3,1),ek,efi,efb,elcoords,eldisplace;
    eldisplace.assign(8,3,0.);
    NRvector<Doub> ptsw(4);
    ptsw[0]= -1.;
    ptsw[1]=0.1;
    ptsw[2]=0.2;
    ptsw[3]=0.8;

    bodyforce[1][0]=-20.;
    Int order=2;
    //druckerprager mat = druckerprager(20000,0.3,10., 30*M_PI/180.);
    elastoplastic3D< druckerprager >* mat0 = new elastoplastic3D< druckerprager >(bodyforce,order);
   // void elastoplastic2D<YC>::Contribute(NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody,NRvector<Doub> ptsw, NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace)
    elcoords.assign(8,3,0.);
    elcoords[0][0]=0;elcoords[0][1]=0;elcoords[0][2]=0;
    elcoords[1][0]=1;elcoords[1][1]=0;elcoords[1][2]=0;
    elcoords[2][0]=1;elcoords[2][1]=1;elcoords[2][2]=0;
    elcoords[3][0]=0;elcoords[3][1]=1;elcoords[3][2]=0;
    elcoords[4][0]=0;elcoords[4][1]=0;elcoords[4][2]=1;
    elcoords[5][0]=1;elcoords[5][1]=0;elcoords[5][2]=1;
    elcoords[6][0]=1;elcoords[6][1]=1;elcoords[6][2]=1;
    elcoords[7][0]=0;elcoords[7][1]=1;elcoords[7][2]=1;

    mat0->fYC.setup(20000., 0.3, 10., 30.*M_PI/180.);
	mat0->SetMemory(1000, 1000);
	mat0->UpdateBodyForce(bodyforce);
    //mat0->Contribute(ek,efi,efb,ptsw,elcoords,eldisplace);
   // mat0->CacStiff(ek,efi,efb,elcoords,eldisplace);
    //elastoplastic3D<druckerprager> *materaial = new elastoplastic3D<druckerprager>(bodyforce,order);
    //ek.Print();
   // efi.Print();
   // efb.Print();
    //return 0;
    */
/*
    NRmatrix<Doub> bodyforce(3,1),ek,efi,efb,elcoords,eldisplace;
    eldisplace.assign(8,3,0.);
    NRvector<Doub> ptsw(4);
    ptsw[0]= -1.;
    ptsw[1]=0.1;
    ptsw[2]=0.2;
    ptsw[3]=0.8;

    bodyforce[1][0]=-20.;
    Int order=2;
    elastoplastic3D< druckerprager >* mat0 = new elastoplastic3D< druckerprager >(bodyforce,order);



    mat0->fYC.setup(20000., 0.3, 10., 30.*M_PI/180.);
	mat0->SetMemory(1000, 1000);
	mat0->UpdateBodyForce(bodyforce);



    double coords[20][3]= {
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, 1},
        {-1, 1, 1}, {0, -1, -1}, {1, 0, -1}, {0, 1, -1}, {-1, 0, -1}, {0, -1, 1}, {1, 0, 1},
        {0, 1, 1}, {-1, 0, 1}, {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}
    };
    elcoords.assign(20,3,0.);
    eldisplace.assign(20,3,0.);
    for(int i=0;i<20;i++)for(int j=0;j<3;j++)elcoords[i][j]=coords[i][j];
    //mat0->Contribute(ek,efi,efb,ptsw,elcoords,eldisplace);

    mat0->CacStiff(ek,efi,efb,elcoords,eldisplace);
   */
/*
    string nodestr = "/home/diogo/projects/dcproj/one-el-p2-nodes.txt";
	string elsstr = "/home/diogo/projects/dcproj/one-el-p2-els.txt";

	MatDoub hhatinho;
	MatDoub  meshcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("one-el-p2-nodes-debug.txt");
	OutPutPost(meshcoords, filemesh1);
	std::ofstream filemesh2("one-el-p2-els-debug.txt");
	OutPutPost(meshtopology, filemesh2);
    int dim=3;
    MatDoub KG,Fint,Fbody;
    NRmatrix<Doub> bodyforce(3,1);
    bodyforce[1][0]=-20.;
    Int order=2;

    MatDoub ptsweigths;
	shapehexahedron shape = shapehexahedron(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows() * npts;
	Int sz = dim * meshcoords.nrows();

    elastoplastic3D< druckerprager >* mat0 = new elastoplastic3D< druckerprager >(bodyforce,order);

	mat0->fYC.setup(20000.,0.3, 10., 30*M_PI/180.);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);

    mesh* mesh0 = new mesh(dim,mat0, allcoords, meshcoords, meshtopology);
    mesh0->Assemble(KG,Fint, Fbody);

    std::vector<int> ids(8);
    ids[0]=1;
    ids[1]=2;
    ids[2]=3;
    ids[3]=4;
    ids[4]=5;
    ids[5]=6;
    ids[6]=7;
    ids[7]=8;

    int dir =1;
    double val =0.;
    mesh0->fmaterial->DirichletBC(KG,Fint, ids, dir,val);

    dir =2;
    mesh0->fmaterial->DirichletBC(KG,Fint, ids, dir,val);

    dir =3;
    mesh0->fmaterial->DirichletBC(KG,Fint, ids, dir,val);

    KG.Print();

*/




