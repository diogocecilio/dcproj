


#include "slopeproject.h"
#include <chrono>
//class VTKGraphMesh;
#ifdef _WIN32    /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
#include <direct.h>
#endif
//#include "vtkmesh.h"
#include "shapetri.h"
slopeproject::slopeproject()
{

}

slopeproject::slopeproject ( mesh* inmesh, KLGalerkinRF* inklgalerking )
{
    fmesh = inmesh;
    fklgalerking = inklgalerking;
   // findbcidsfatsn ( fmesh, fidsvector );
	findbcids( fmesh, fidsvector );
}
slopeproject::slopeproject ( mesh* inmesh, KLGalerkinRF* inklgalerking, NRmatrix<MatDoub> randomfield )
{
    fmesh = inmesh;
    fklgalerking = inklgalerking;
    frandomfield = randomfield;
   // findbcidsfatsn ( fmesh, fidsvector );
	findbcids( fmesh, fidsvector );
}
slopeproject::slopeproject ( mesh inmesh, KLGalerkinRF inklgalerking, NRmatrix<MatDoub> randomfield )
{
    fmesh = &inmesh;
    fklgalerking = &inklgalerking;
    frandomfield = randomfield;
   // findbcidsfatsn ( fmesh, fidsvector );
	findbcids( fmesh, fidsvector );
}
slopeproject::~slopeproject()
{
    delete fmesh;
    delete fklgalerking;
}

slopeproject::slopeproject ( slopeproject& copy )
{

}

void slopeproject::CreateRandomField ( string namefolder )
{

    KLGalerkinRF* objKLGalerkinRF = fklgalerking;

    int check;

    char* cstr = new char[namefolder.length() + 1];
    strcpy ( cstr, namefolder.c_str() );
#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
    check = mkdir ( cstr,777 );
#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
    check = mkdir ( cstr );
#endif

    string datafile = namefolder;
    datafile += "/datarandom.txt";
    std::ofstream file ( datafile );
    file << " samples = " << objKLGalerkinRF->GetSamples() <<
         " | expansion order = " << objKLGalerkinRF->GetExpansionorder() << " | func type = " << objKLGalerkinRF->Getftype() << endl;
    file << "Lx = " << objKLGalerkinRF->GetLx() << " | Ly = " << objKLGalerkinRF->GetLy() << " variance = " << objKLGalerkinRF->GetSig() << endl;
    VecComplex  val;
    MatDoub  vec, HHAT;
    NRmatrix<MatDoub> randomfield;
    NRmatrix<Doub> cosionfield, frictionfield;

    std::vector<std::vector<double>> errpost;
    objKLGalerkinRF->SolveGenEigValProblem ( val, vec, randomfield, errpost );
    cosionfield = randomfield[0][0];
    frictionfield = randomfield[1][0];
    datafile = namefolder;
    datafile += "/coesionfield.txt";
    std::ofstream coesfile ( datafile );
    OutPutPost ( cosionfield, coesfile );

    datafile = namefolder;
    datafile += "/frictionfield.txt";
    std::ofstream frictionfile ( datafile );
    OutPutPost ( frictionfield, frictionfile );

    delete objKLGalerkinRF;
}

//void slopeproject::findbcids(mesh* gmesh, std::vector<std::vector<int>>& idsvector)
//{
//	Int ndivs = 100000;
//	MatDoub pathbottom, pathleft, pathright, pathdisplace;
//	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
//	VecDoub a(2), b(2);
//	a[0] = 0.; a[1] = 0.;
//	b[0] = 80.; b[1] = 0;
//	Line(a, b, ndivs, pathbottom);
//	FindIdsInPath(pathbottom, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsbottom);
//
//	idsvector.push_back(idsbottom);
//
//	a[0] = 0.; a[1] = 0.;
//	b[0] = 0.; b[1] = 40.;
//	Line(a, b, ndivs, pathleft);
//	FindIdsInPath(pathleft, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsleft);
//
//	idsvector.push_back(idsleft);
//	a[0] = 80.; a[1] = 0.;
//	b[0] = 80.; b[1] = 20.;
//	Line(a, b, ndivs, pathright);
//	FindIdsInPath(pathright, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsright);
//
//	idsvector.push_back(idsright);
//
//	a[0] = 29.99; a[1] = 39.99;
//	b[0] = 30.; b[1] = 40.;
//	Line(a, b, ndivs, pathdisplace);
//	FindIdsInPath(pathdisplace, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), iddisplace);
//
//	idsvector.push_back(iddisplace);
//
//}

void slopeproject::SolveEigen ( SparseMatrix<double> A, VectorXd b, VectorXd& x )
{
    //start = sc.now();
    SimplicialLLT< SparseMatrix<double> > solver;
    x = solver.compute ( A ).solve ( b );
    // end = sc.now();
    /// time_span = static_cast<chrono::duration<double>>(end - start);
    // cout << "Operation took: " << time_span.count() << " seconds.";
}

void slopeproject::findbcids45 ( mesh* gmesh, std::vector<std::vector<int>>& idsvector )
{
    //malha 606
    //50x20
    Int ndivs = 100000;
    MatDoub pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    VecDoub a ( 2 ), b ( 2 );
    a[0] = 0.;
    a[1] = 0.;
    b[0] = 50.;
    b[1] = 0;
    Line ( a, b, ndivs, pathbottom );

    std::vector<std::vector<std::vector<double> > > ggmesh = gmesh->GetAllCoords();

    NRmatrix<int> meshtop = gmesh->GetMeshTopology();

    FindIdsInPath ( pathbottom, ggmesh,meshtop, idsbottom );

    idsvector.push_back ( idsbottom );

    a[0] = 0.;
    a[1] = 0.;
    b[0] = 0.;
    b[1] = 20.;
    Line ( a, b, ndivs, pathleft );
    FindIdsInPath ( pathleft, ggmesh,meshtop, idsleft );

    idsvector.push_back ( idsleft );
    a[0] = 50.;
    a[1] = 0.;
    b[0] = 50.;
    b[1] = 10;
    Line ( a, b, ndivs, pathright );
    FindIdsInPath ( pathright, ggmesh,meshtop, idsright );

    idsvector.push_back ( idsright );

    a[0] = 19.99;
    a[1] = 19.99;
    //a[0] = 14.99; a[1] = 19.99;
    b[0] = 20.;
    b[1] = 20.;
    //b[0] = 15.; b[1] = 20.;
    Line ( a, b, ndivs, pathdisplace );
    FindIdsInPath ( pathdisplace,ggmesh,meshtop, iddisplace );

    //cout<< "ID displace = "<<iddisplace[0]<<endl;

    idsvector.push_back ( iddisplace );

}

void slopeproject::findbcids ( mesh* gmesh, std::vector<std::vector<int>>& idsvector )
{
    //30x10
    ////mesh 2x1h5
    Int ndivs = 100000;
    MatDoub pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    VecDoub a ( 2 ), b ( 2 );
    a[0] = 0.;
    a[1] = 0.;
    b[0] = 50.;
    b[1] = 0;
    Line ( a, b, ndivs, pathbottom );

    std::vector<std::vector<std::vector<double> > > ggmesh = gmesh->GetAllCoords();

    NRmatrix<int> meshtop = gmesh->GetMeshTopology();

    FindIdsInPath ( pathbottom, ggmesh,meshtop, idsbottom );

    idsvector.push_back ( idsbottom );

    a[0] = 0.;
    a[1] = 0.;
    b[0] = 0.;
    b[1] = 20.;
    Line ( a, b, ndivs, pathleft );
    FindIdsInPath ( pathleft, ggmesh,meshtop, idsleft );

    idsvector.push_back ( idsleft );
    a[0] = 50.;
    a[1] = 0.;
    b[0] = 50.;
    b[1] = 10.;
    Line ( a, b, ndivs, pathright );
    FindIdsInPath ( pathright, ggmesh,meshtop, idsright );

    idsvector.push_back ( idsright );

    a[0] = 19.99;
    a[1] = 19.99;
    b[0] = 20.;
    b[1] = 20.;
    Line ( a, b, ndivs, pathdisplace );
    FindIdsInPath ( pathdisplace,ggmesh,meshtop, iddisplace );

    // cout<< "ID displace = "<<iddisplace[0]<<endl;

    idsvector.push_back ( iddisplace );

}


void slopeproject::findbcidsfat5022 ( mesh* gmesh, std::vector<std::vector<int>>& idsvector )
{
    //30x10
    ////mesh 2x1h5
    Int ndivs = 100000;
    MatDoub pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    VecDoub a ( 2 ), b ( 2 );
    a[0] = 0.;
    a[1] = 0.;
    b[0] = 50.;
    b[1] = 0;
    Line ( a, b, ndivs, pathbottom );

    std::vector<std::vector<std::vector<double> > > ggmesh = gmesh->GetAllCoords();

    NRmatrix<int> meshtop = gmesh->GetMeshTopology();

    FindIdsInPath ( pathbottom, ggmesh,meshtop, idsbottom );

    idsvector.push_back ( idsbottom );

    a[0] = 0.;
    a[1] = 0.;
    b[0] = 0.;
    b[1] = 22.;
    Line ( a, b, ndivs, pathleft );
    FindIdsInPath ( pathleft, ggmesh,meshtop, idsleft );

    idsvector.push_back ( idsleft );
    a[0] = 50.;
    a[1] = 0.;
    b[0] = 50.;
    b[1] = 17.;
    Line ( a, b, ndivs, pathright );
    FindIdsInPath ( pathright, ggmesh,meshtop, idsright );

    idsvector.push_back ( idsright );

    a[0] = 19.99;
    a[1] = 21.99;
    b[0] = 20.;
    b[1] = 22.;
    Line ( a, b, ndivs, pathdisplace );
    FindIdsInPath ( pathdisplace,ggmesh,meshtop, iddisplace );

    // cout<< "ID displace = "<<iddisplace[0]<<endl;

    idsvector.push_back ( iddisplace );

}

void slopeproject::findbcidsfatsn ( mesh* gmesh, std::vector<std::vector<int>>& idsvector )
{
    //30x10
    ////mesh 2x1h5
    Int ndivs = 100000;
    MatDoub pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    VecDoub a ( 2 ), b ( 2 );
    a[0] = 0.;
    a[1] = 0.;
    b[0] = 75.;
    b[1] = 0;
    Line ( a, b, ndivs, pathbottom );

    std::vector<std::vector<std::vector<double> > > ggmesh = gmesh->GetAllCoords();

    NRmatrix<int> meshtop = gmesh->GetMeshTopology();

    FindIdsInPath ( pathbottom, ggmesh,meshtop, idsbottom );

    idsvector.push_back ( idsbottom );

    a[0] = 0.;
    a[1] = 0.;
    b[0] = 0.;
    b[1] = 40.;
    Line ( a, b, ndivs, pathleft );
    FindIdsInPath ( pathleft, ggmesh,meshtop, idsleft );

    idsvector.push_back ( idsleft );
    a[0] = 75.;
    a[1] = 0.;
    b[0] = 75.;
    b[1] = 30.;
    Line ( a, b, ndivs, pathright );
    FindIdsInPath ( pathright, ggmesh,meshtop, idsright );

    idsvector.push_back ( idsright );

    a[0] = 34.99;
    a[1] = 39.99;
    b[0] = 35.;
    b[1] = 40.;
    Line ( a, b, ndivs, pathdisplace );
    FindIdsInPath ( pathdisplace,ggmesh,meshtop, iddisplace );

    // cout<< "ID displace = "<<iddisplace[0]<<endl;

    idsvector.push_back ( iddisplace );

}

void slopeproject::findbcidsfat ( mesh* gmesh, std::vector<std::vector<int>>& idsvector )
{
    //30x10
    ////mesh 2x1h5
    Int ndivs = 100000;
    MatDoub pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    VecDoub a ( 2 ), b ( 2 );
    a[0] = 0.;
    a[1] = 0.;
    b[0] = 38.;
    b[1] = 0;
    Line ( a, b, ndivs, pathbottom );

    std::vector<std::vector<std::vector<double> > > ggmesh = gmesh->GetAllCoords();

    NRmatrix<int> meshtop = gmesh->GetMeshTopology();

    FindIdsInPath ( pathbottom, ggmesh,meshtop, idsbottom );

    idsvector.push_back ( idsbottom );

    a[0] = 0.;
    a[1] = 0.;
    b[0] = 0.;
    b[1] = 19.;
    Line ( a, b, ndivs, pathleft );
    FindIdsInPath ( pathleft, ggmesh,meshtop, idsleft );

    idsvector.push_back ( idsleft );
    a[0] = 38.;
    a[1] = 0.;
    b[0] = 38.;
    b[1] = 14.;
    Line ( a, b, ndivs, pathright );
    FindIdsInPath ( pathright, ggmesh,meshtop, idsright );

    idsvector.push_back ( idsright );

    a[0] = 15.99;
    a[1] = 18.99;
    b[0] = 16.;
    b[1] = 19.;
    Line ( a, b, ndivs, pathdisplace );
    FindIdsInPath ( pathdisplace,ggmesh,meshtop, iddisplace );

    //  cout<< "ID displace = "<<iddisplace[0]<<endl;

    idsvector.push_back ( iddisplace );

}


/*
std::vector<std::vector<double>>  slopeproject::IterativeProcessShearRed(Doub fac, Doub delta,Doub tol) {

	mesh* meshint = fmesh;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	std::vector<std::vector<int>> idsvector;

	findbcids(meshint, idsvector);

	material* mat = meshint->fmaterial;

	idsbottom = idsvector[0];
	idsleft = idsvector[1];
	idsright = idsvector[2];
	iddisplace = idsvector[3];

	mat->ResetMat();
	Int sz = 2 * meshint->GetMeshNodes().nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	MatDoub displace, displace0, R, sol;
	displace.assign(sz, 1, 0.);
	displace0.assign(sz, 1, 0.);

	Doub l = 0, l0 = 0, lamb = 1., lambn = 0, dlamb = 0., lamb3, diff = 100, diff2 = 100;
	Int counterout = 0, maxcountout = 80;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);

	Int niter = 500;
	Int postcounter = 0;

	NRvector<Doub> matconsts(4, 0.);
	mat->GetMatConstants(matconsts);
	Doub  young = matconsts[0];
	Doub nu = matconsts[1];
	Doub c = matconsts[2];
	Doub phi = matconsts[3];

	Doub res = 10, facn = 0, FS = fac, FSmin = 0, FSmax = 100000;
	Doub norm = 100000.;
	bool boll = false;

	Doub c0 = c;
	Doub phi0 = phi;

	NRmatrix<Doub> hhat0;
	meshint->GetHhat(hhat0);
	NRmatrix<Doub> hhatcopy  = hhat0;



	if (hhat0.nrows() != 0)
	{
		for (int irow = 0; irow < hhat0.nrows(); irow++)
		{
			hhatcopy[irow][0] = hhat0[irow][0] / FS;
			hhatcopy[irow][1] = atan(tan(hhat0[irow][1]) / FS);
		}
		meshint->SetHhat(hhatcopy);
		c = c0 / FS;
		phi = atan(tan(phi0) / FS);
		matconsts[2] = c;
		matconsts[3] = phi;
		mat->SetMatConstants(matconsts);
	}
	else {
		c = c0 / FS;
		phi = atan(tan(phi0) / FS);
		matconsts[2] = c;
		matconsts[3] = phi;
		mat->SetMatConstants(matconsts);
	}

	do
	{
		Int counter = 0;
		//displace.assign(sz, 1, 0.);
		//material->UpdateDisplacement(displace);
		sol.assign(sz, 1, 0.);
		norm = 1000.;

		mat->ResetMat();
        fmesh->fmaterial->SetTangentMatrixType(false);
        Doub tol2 = 0.01;
		do
		{
            chrono::steady_clock sc;
            auto start = sc.now();

			FINT.assign(sz, 1, 0.);
			meshint->Assemble(KG, FINT, FBODY);
			R = FBODY;
			//R *= fac;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);

		//	SolveEigen(KG, R, sol);
            SolveEigenSparse(KG, R, sol);

			displace += sol;
			Doub u = fabs(displace[2 * iddisplace[0] + 1][0]);
			mat->UpdateDisplacement(displace);
			norm = R.NRmatrixNorm();

            auto end = sc.now();
            auto time_span = static_cast<chrono::duration<double>>(end - start);
            cout << "Iter  " << counter ;
            cout << "| Time: " << time_span.count();
			std::cout << " R norm = " << norm << " | phi0/phi = " << tan(phi0) / tan(phi) << " | c0/c = " << c0 / c << " | c = " << c <<" | phi = " << phi << std::endl;
			counter++;
			postcounter++;

            if (isnan(norm) == 1) {
				std::cout << "NAN" << endl;
				return solpost;
			}
			//if(counter>20)fmesh->fmaterial->SetElastic();

			//cout<< "fbody norm = "<< FBODY.NRmatrixNorm() <<endl;
			if(counter>10)fmesh->fmaterial->SetTangentMatrixType(true);

		} while (norm > tol2 && counter <30 && norm < 1000.);

        if(norm<1){
		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = FS;
		solcount[2] = 0;
		solcount[3] = 0;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		uvf[1] = FS;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);
}

		res = (fac - facn) / fac;
		if (norm>= tol2) {
			displace = displace0;
			//R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);
			//facn = FS;
			FSmax = FS;
			FS = (FSmin + FSmax) / 2.;
			//fac = FS;
			//mat->ResetMat();
			//mat->UpdateDisplacement(displace0);
			boll = true;
		}
		else {
			if (boll == true)
			{
				displace0 = displace;
				facn = FS;
				FSmin = FS;
				FS = 1. / ((1. / FSmin + 1. / FSmax) / 2.);
				fac = FS;

			}
			else {
				displace0 = displace;
				facn = FS;
				FSmin = FS;
				FS += delta;
				fac = FS;
				//mat->UpdatePlasticStrain();
			}

		}



		//Caso tenha random field a reducao da resistencia tem que ser aplicada a todos os pontos
		if (hhat0.nrows() != 0)
		{
			for (int irow = 0; irow < hhat0.nrows(); irow++)
			{
				hhatcopy[irow][0] = hhat0[irow][0] / FS;
				hhatcopy[irow][1] = atan(tan(hhat0[irow][1]) / FS);
			}
			meshint->SetHhat(hhatcopy);
			c = c0 / FS;
			phi = atan(tan(phi0) / FS);
			matconsts[2] = c;
			matconsts[3] = phi;
			mat->SetMatConstants(matconsts);
		}
		else {
			c = c0 / FS;
			phi = atan(tan(phi0) / FS);
			matconsts[2] = c;
			matconsts[3] = phi;
			mat->SetMatConstants(matconsts);
		}


		counterout++;

		std::cout << " iter = " << counterout << " | FS= " << FS << " | FSmax " << FSmax << " | FSmin = " << FSmin << std::endl;



	}  while ((FSmax - FSmin) / FS > tol); //while ((FSmax - FSmin) / FS > tol && ( norm <1.e20));
	//std::cout << "FOS eeeeee = " << FS << std::endl;

	if (true)
	{
        std::cout << "FOS eeeeee = " << FS << std::endl;
        MatDoub solpost23;
		solpost23.CopyFromVector(solpost2);
		string names = "loadvsdisplacementSRM.dat";
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);

		std::vector<std::vector<double>> epsppost;
		mat->PostProcessIntegrationPointVar(meshint->GetAllCoords(), meshint->GetMeshNodes(), meshint->GetMeshTopology(), mat->GetSolution(), epsppost);
		//string name3 = "/home/diogocecilio/Dropbox/slope-reliability/results/mathematicas-new/sqrtJ2SRM-Cho-Determ";
		//string ext3 = ".dat";
        string name3 = "sqrtJ2SRM-Cho-Determ.dat";
		//string ext3 = ".dat";
		//name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		mat->PostProcess(meshint->GetAllCoords(), meshint->GetMeshNodes(), meshint->GetMeshTopology(), mat->GetSolution(), solx, soly);
		filename = "/home/diogocecilio/Dropbox/slope-reliability/results/mathematicas-new/solySRM.dat";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "/home/diogocecilio/Dropbox/slope-reliability/results/mathematicas-new/solxSRM.dat";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

	return solpost;
}*/


std::vector<std::vector<double>>  slopeproject::IterativeProcessShearRed ( Doub fac, Doub delta,Doub tol )
{

    std::ofstream file8 ( "uxf-srm-mc2x1.nb" );
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    std::vector<std::vector<int>> idsvector;
idsvector=fidsvector;
    //findbcids ( fmesh, idsvector );
  //  findbcidsfatsn ( fmesh, idsvector );
//findbcidsfat5022( fmesh, idsvector );
    idsbottom = idsvector[0];
    idsleft = idsvector[1];
    idsright = idsvector[2];
    iddisplace = idsvector[3];

    Int sz = 2 * fmesh->GetMeshNodes().nrows();
    MatDoub KG ( sz, sz, 0. ), FBODY ( sz, 1, 0. ), FINT ( sz, 1, 0. ), ptsweigths,displace,displace0,R;


    Int counterout=0;

    fmesh->fmaterial->ResetPlasticStrain();
    fmesh->fmaterial->ResetDisplacement();
    fmesh->fmaterial->ResetMat();


    Int type=1;



    MatDoub  sol;
    displace.assign ( sz, 1, 0. );
    displace0.assign ( sz, 1, 0. );


    std::vector<double> solcount ( 7, 0. ), uvf ( 2, 0. );
    std::vector<std::vector<double>> solpost, solpost2;
    solpost.push_back ( solcount );

    Int niter = 500;
    Int postcounter = 0;


    NRvector<Doub> matconsts ( 4, 0. );
    fmesh->fmaterial->GetMatConstants ( matconsts );
    Doub  young = matconsts[0];
    Doub nu = matconsts[1];
    Doub c = matconsts[2];
    Doub phi = matconsts[3];

    //       consts[0] = fyoung;
    //    consts[1] = fnu;
    //    consts[2] = fc ;
    //    consts[3] = fPhi;
    //    consts[4] = fPsi;

    Doub res = 10, facn = 0, FS = fac, FSmin = 0, FSmax =10000,diff=10,diff2=10;
    Doub norm = 100000.;
    bool boll = false;

    Doub c0 = c;
    Doub phi0 = phi;

    NRmatrix<Doub> hhat0;
    fmesh->GetHhat ( hhat0 );
    NRmatrix<Doub> hhatcopy  = hhat0;

//    NRvector<Doub> matconsts2 ( 3 );
//    bool vonmises=true;

//     if ( vonmises ) {
//
//         c = c0 / FS;
//         matconsts2[2] = c;
//         mat->SetMatConstants ( matconsts2 );
//     } else {

    if ( hhat0.nrows() != 0 ) {
        for ( int irow = 0; irow < hhat0.nrows(); irow++ ) {
            hhatcopy[irow][0] = hhat0[irow][0] / FS;
            hhatcopy[irow][1] = atan ( tan ( hhat0[irow][1] ) / FS );
        }
        fmesh->SetHhat ( hhatcopy );
        c = c0 / FS;
        phi = atan ( tan ( phi0 ) / FS );
        matconsts[2] = c;
        matconsts[3] = phi;
        fmesh->fmaterial->SetMatConstants ( matconsts );
    } else {
        c = c0 / FS;
        phi = atan ( tan ( phi0 ) / FS );
        matconsts[2] = c;
        matconsts[3] = phi;
        fmesh->fmaterial->SetMatConstants ( matconsts );
    }



//    }
   fmesh->fmaterial->SetTangentMatrixType ( false );
    do {

        Int counter = 0;
        //displace.assign(sz, 1, 0.);
        //material->UpdateDisplacement(displace);
        sol.assign ( sz, 1, 0. );
        norm = 1000.;

        fmesh->fmaterial->ResetMat();
   //      fmesh->fmaterial->SetTangentMatrixType ( false );
        Doub tol2 = 1;
        do {
            chrono::steady_clock sc;
            auto start = sc.now();

            FINT.assign ( sz, 1, 0. );
            fmesh->Assemble ( KG, FINT, FBODY );
            R = FBODY;
            //R *= fac;
            R -= FINT;

            InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );

            norm = R.NRmatrixNorm();

            //cout << norm << endl;
            SolveEigenSparse ( type,KG, R, sol );
            //  SolveEigen3(KG, R, sol);

            //MatDoub KGinv;
            //KG.ComputeInverse(KGinv);
            //KGinv.Mult(R,sol);
            displace += sol;
            Doub u = fabs ( displace[2 * iddisplace[0] + 1][0] );
            fmesh->fmaterial->UpdateDisplacement ( displace );
            norm = R.NRmatrixNorm();

            auto end = sc.now();
            auto time_span = static_cast<chrono::duration<double>> ( end - start );
          //  cout << "Iter  " << counter ;
          //  cout << "| Time: " << time_span.count();
           // std::cout << " R norm = " << norm << " | phi0/phi = " << tan ( phi0 ) / tan ( phi ) << " | c0/c = " << c0 / c << " | c = " << c <<" | phi = " << phi << //std::endl;
            counter++;
            postcounter++;

            if ( isnan ( norm ) == 1 ) {
                std::cout << "NAN" << endl;
                return solpost;
            }
            //if(counter>20)fmesh->fmaterial->SetElastic();

            //cout<< "fbody norm = "<< FBODY.NRmatrixNorm() <<endl;
            //if ( counter>10 ) fmesh->fmaterial->SetTangentMatrixType ( true );
            // if(norm > 4000.)

        } while ( norm > tol2 && counter <20 && norm < 200000. );

        if ( norm<1 ) {
            solcount[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
            solcount[1] = FS;
            solcount[2] = 0;
            solcount[3] = 0;
            solcount[4] = diff;
            solcount[5] = diff2;
            solcount[6] = counterout;

            uvf[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
            uvf[1] = FS;

            solpost2.push_back ( uvf );
            solpost.push_back ( solcount );
        }

        res = ( fac - facn ) / fac;
        if ( norm>= tol2 ) {
            displace = displace0;
            //R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);
            //facn = FS;
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
            //fac = FS;
            //mat->ResetMat();
            //mat->UpdateDisplacement(displace0);
            boll = true;
        } else {
            if ( boll == true ) {
                displace0 = displace;
                facn = FS;
                FSmin = FS;
                FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
                fac = FS;

            } else {
                displace0 = displace;
                facn = FS;
                FSmin = FS;
                FS += delta;
                fac = FS;
                //mat->UpdatePlasticStrain();
            }

        }
        //Caso tenha random field a reducao da resistencia tem que ser aplicada a todos os pontos
        if ( hhat0.nrows() != 0 ) {
            for ( int irow = 0; irow < hhat0.nrows(); irow++ ) {
                hhatcopy[irow][0] = hhat0[irow][0] / FS;
                hhatcopy[irow][1] = atan ( tan ( hhat0[irow][1] ) / FS );
            }
            fmesh->SetHhat ( hhatcopy );
            c = c0 / FS;
            phi = atan ( tan ( phi0 ) / FS );
            matconsts[2] = c;
            matconsts[3] = phi;
            fmesh->fmaterial->SetMatConstants ( matconsts );
        } else {
            c = c0 / FS;
            phi = atan ( tan ( phi0 ) / FS );
            matconsts[2] = c;
            matconsts[3] = phi;
            fmesh->fmaterial->SetMatConstants ( matconsts );
        }

        // fmesh->fmaterial->ResetMat();
        counterout++;

        std::cout << " iter = " << counterout << " | FS= " << FS << " | FSmax " << FSmax << " | FSmin = " << FSmin << std::endl;



    }  while ( ( FSmax - FSmin ) / FS > tol ); //while ((FSmax - FSmin) / FS > tol && ( norm <1.e20));
    std::cout << " iter = " << counterout << " | FS= " << FS << " | FSmax " << FSmax << " | FSmin = " << FSmin << std::endl;

//	fmesh->fmaterial->UpdateDisplacement ( displace );
    return solpost;
}


std::vector<std::vector<double>>  slopeproject::IterativeProcessShearRed2 ( Doub fac, Doub delta,Doub tol )
{

    mesh* meshint = fmesh;
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    std::vector<std::vector<int>> idsvector;

    findbcids ( meshint, idsvector );

    material* mat = meshint->fmaterial;

    idsbottom = idsvector[0];
    idsleft = idsvector[1];
    idsright = idsvector[2];
    iddisplace = idsvector[3];

    mat->ResetMat();
    Int sz = 2 * meshint->GetMeshNodes().nrows();
    SparseMatrix<double>  KG ( sz,sz );
    VectorXd FBODY ( sz ), FINT ( sz ),R ( sz ),sol ( sz ),displace ( sz ), displace0 ( sz );
    MatDoub ptsweigths;
    Doub rtemp = 10.;


    Doub l = 0, l0 = 0, lamb = 1., lambn = 0, dlamb = 0., lamb3, diff = 100, diff2 = 100;
    Int counterout = 0, maxcountout = 80;
    std::vector<double> solcount ( 7, 0. ), uvf ( 2, 0. );
    std::vector<std::vector<double>> solpost, solpost2;
    solpost.push_back ( solcount );

    Int niter = 500;
    Int postcounter = 0;

    NRvector<Doub> matconsts ( 4, 0. );
    mat->GetMatConstants ( matconsts );
    Doub  young = matconsts[0];
    Doub nu = matconsts[1];
    Doub c = matconsts[2];
    Doub phi = matconsts[3];

    Doub res = 10, facn = 0, FS = fac, FSmin = 0, FSmax = 100000;
    Doub norm = 100000.;
    bool boll = false;

    Doub c0 = c;
    Doub phi0 = phi;

    NRmatrix<Doub> hhat0;
    meshint->GetHhat ( hhat0 );
    NRmatrix<Doub> hhatcopy  = hhat0;



    if ( hhat0.nrows() != 0 ) {
        for ( int irow = 0; irow < hhat0.nrows(); irow++ ) {
            hhatcopy[irow][0] = hhat0[irow][0] / FS;
            hhatcopy[irow][1] = atan ( tan ( hhat0[irow][1] ) / FS );
        }
        meshint->SetHhat ( hhatcopy );
        c = c0 / FS;
        phi = atan ( tan ( phi0 ) / FS );
        matconsts[2] = c;
        matconsts[3] = phi;
        mat->SetMatConstants ( matconsts );
    } else {
        c = c0 / FS;
        phi = atan ( tan ( phi0 ) / FS );
        matconsts[2] = c;
        matconsts[3] = phi;
        mat->SetMatConstants ( matconsts );
    }



    do {
        Int counter = 0;
        //displace.assign(sz, 1, 0.);
        //material->UpdateDisplacement(displace);
        sol.resize ( sz );
        norm = 1000.;

        mat->ResetMat();
        do {
            chrono::steady_clock sc;
            auto start = sc.now();

            FINT.resize ( sz );
            meshint->Assemble ( KG, FINT, FBODY );
            R = FBODY;
            //R *= fac;
            R -= FINT;

            //InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);

            SolveEigen ( KG, R, sol );

            displace += sol;
            Doub u = fabs ( displace ( 2 * iddisplace[0] + 1 ) );
            mat->UpdateDisplacement ( displace );
            norm = R.norm();

            auto end = sc.now();
            auto time_span = static_cast<chrono::duration<double>> ( end - start );
            cout << "Operation took: " << time_span.count() << " seconds !!!";

            std::cout << " R norm = " << norm << " | phi0/phi = " << tan ( phi0 ) / tan ( phi ) << " | c0/c = " << c0 / c << " | c = " << c <<" | phi = " << phi << std::endl;
            counter++;
            postcounter++;

            if ( isnan ( norm ) == 1 ) {
                std::cout << "NAN" << endl;
                return solpost;
            }

        } while ( norm > 0.01 && counter <15 && norm < 5000. );

        if ( norm<1 ) {
            solcount[0] = fabs ( displace ( 2 * iddisplace[0] + 1 ) );
            solcount[1] = FS;
            solcount[2] = 0;
            solcount[3] = 0;
            solcount[4] = diff;
            solcount[5] = diff2;
            solcount[6] = counterout;

            uvf[0] = fabs ( displace ( 2 * iddisplace[0] + 1 ) );
            uvf[1] = FS;

            solpost2.push_back ( uvf );
            solpost.push_back ( solcount );
        }

        res = ( fac - facn ) / fac;
        if ( norm>= 1. ) {
            displace = displace0;
            //R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);
            //facn = FS;
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
            //fac = FS;
            //mat->ResetMat();
            //mat->UpdateDisplacement(displace0);
            boll = true;
        } else {
            if ( boll == true ) {
                displace0 = displace;
                facn = FS;
                FSmin = FS;
                FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
                fac = FS;

            } else {
                displace0 = displace;
                facn = FS;
                FSmin = FS;
                FS += delta;
                fac = FS;
                //mat->UpdatePlasticStrain();
            }

        }



        //Caso tenha random field a reducao da resistencia tem que ser aplicada a todos os pontos
        if ( hhat0.nrows() != 0 ) {
            for ( int irow = 0; irow < hhat0.nrows(); irow++ ) {
                hhatcopy[irow][0] = hhat0[irow][0] / FS;
                hhatcopy[irow][1] = atan ( tan ( hhat0[irow][1] ) / FS );
            }
            meshint->SetHhat ( hhatcopy );
            c = c0 / FS;
            phi = atan ( tan ( phi0 ) / FS );
            matconsts[2] = c;
            matconsts[3] = phi;
            mat->SetMatConstants ( matconsts );
        } else {
            c = c0 / FS;
            phi = atan ( tan ( phi0 ) / FS );
            matconsts[2] = c;
            matconsts[3] = phi;
            mat->SetMatConstants ( matconsts );
        }


        counterout++;

        std::cout << " iter = " << counterout << " | FS= " << FS << " | FSmax " << FSmax << " | FSmin = " << FSmin << std::endl;



    }  while ( ( FSmax - FSmin ) / FS > tol ); //while ((FSmax - FSmin) / FS > tol && ( norm <1.e20));
    std::cout << "FOS = " << FS << std::endl;

    if ( true ) {
        MatDoub solpost23;
        solpost23.CopyFromVector ( solpost2 );
        //string names = "/home/diogo/projects/results/mathematicas/loadvsdisplacementSRM";
        //string exts = ".dat";
        //names += exts;
        string names="uxf.dat";
        std::ofstream file8 ( names );
        OutPutFile ( solpost23, file8 );

        std::vector<std::vector<double>> epsppost;
        mat->PostProcessIntegrationPointVar ( meshint->GetAllCoords(), meshint->GetMeshNodes(), meshint->GetMeshTopology(), mat->GetSolution(), epsppost );
        string name3 = "/home/diogo/projects/results/mathematicas/sqrtJ2SRM-Cho-Determx";
        string ext3 = ".dat";
        name3 += ext3;
        std::ofstream file3 ( name3 );
        OutPutPost ( epsppost, file3 );

        string filename;
        std::vector<std::vector<double>> solx, soly;
        mat->PostProcess ( meshint->GetAllCoords(), meshint->GetMeshNodes(), meshint->GetMeshTopology(), mat->GetSolution(), solx, soly );
        filename = "/home/diogo/projects/results/mathematicas/solySRM.dat";
        std::ofstream file2 ( filename );
        OutPutPost ( soly, file2 );
        filename = "/home/diogo/projects/results/mathematicas/solxSRM.dat";
        std::ofstream file22 ( filename );
        OutPutPost ( solx, file22 );
    }
    return solpost;
}





void myTreads2 ( int a, int b, slopeproject* slopeobj2, string traedN )
{
    string namefolder3 = "D:/slope-results/THREADS/GI-cho-field-Lx20-Ly4" + traedN;
    slopeobj2->MonteCarloGIM ( a, b, false, namefolder3 );

}

typedef Eigen::SparseMatrix<double> SpMat;
void slopeproject::SolveEigen2 ( MatDoub A, MatDoub b, MatDoub& x )
{
    x.assign ( A.nrows(), 1, 0. );
    SpMat Aa ( A.nrows(), A.nrows() );
    for ( int i = 0; i < A.nrows(); i++ ) {
        for ( int j = 0; j < A.ncols(); j++ ) {
            if ( fabs ( A[i][j] ) >1.e-12 ) {
                Aa.coeffRef ( i, j ) = A[i][j];
            }
        }
    }
    Aa.makeCompressed();
    VectorXd bb ( A.nrows() ), xx ( A.nrows() );
    for ( int i = 0; i < A.nrows(); i++ ) {
        bb ( i ) = b[i][0];
    }

    Eigen::SimplicialLLT<SpMat, Eigen::Lower> solver;
    // Eigen::SimplicialLLT<SpMat> solver;
    //  Eigen::SparseLU<SpMat> solver;
    // Eigen::ConjugateGradient<SpMat> solver;
    //Eigen::BiCGSTAB<SpMat> solver;
    //Eigen::SimplicialLDLT<SpMat,Eigen::Upper> solver;
    //  solver.analyzePattern(Aa);
    solver.compute ( Aa );
    if ( solver.info() !=Eigen::Success ) {
        std::cout << " Decomposition Failed " << std::endl;
        getchar();

    }
    xx = solver.solve ( bb );      // use the factorization to solve for the given right hand side
    for ( int i = 0; i < A.nrows(); i++ ) {
        x[i][0] = xx ( i );
    }
}

//#define EIGEN_USE_MKL_ALL

std::vector<std::vector<double>>   slopeproject::IterativeProcess ( int ndesi, Doub dlamb0,Doub alphatol, int niter )
{

    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    std::vector<std::vector<int>> idsvector;

    findbcids ( fmesh, idsvector );

    idsbottom = idsvector[0];
    idsleft = idsvector[1];
    idsright = idsvector[2];
    iddisplace = idsvector[3];

    Int sz = 2 * fmesh->GetMeshNodes().nrows();
    MatDoub KG ( sz, sz, 0. ), FBODY ( sz, 1, 0. ), FINT ( sz, 1, 0. ), ptsweigths;
    Doub rtemp = 10.;

    //material* mat = fmesh->fmaterial;

    fmesh->fmaterial->ResetPlasticStrain();
    fmesh->fmaterial->ResetDisplacement();
    fmesh->fmaterial->ResetMat();

    MatDoub displace, displace0;
    displace.assign ( sz, 1, 0. );
    //Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
    Doub l = 0, l0 = 0, lamb = 1., lambn =0.1, dlamb = 0.1, lamb3, diff = 100, diff2 = 100;
    Int counterout = 0, maxcountout = 30;
    std::vector<double> solcount ( 7, 0. ), uvf ( 2, 0. );
    std::vector<std::vector<double>> solpost, solpost2;
    solpost.push_back ( solcount );


    MatDoub dws ( sz, 1, 0. ), dwb ( sz, 1, 0. ), dww ( sz, 1, 0. ), dw ( sz, 1, 0. ), R;
    //mesh->fmaterial->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
    fmesh->Assemble ( KG, FINT, FBODY );
    R = FBODY;
    R *= lamb;
    R -= FINT;
    
    InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );
    SolveEigen ( KG, R, dws );
    SolveEigen ( KG, FBODY, dwb );
    MatDoub dwbt, mult;
    dwb.Transpose ( dwbt );
    dwbt.Mult ( dwb, mult );
    l0 = sqrt ( dlamb0 * dlamb0 * mult[0][0] );
    l = l0;
    //l=0.05;
    dlamb = computelamda ( dwb, dws, dw, l );
    lamb = 0;
	bool check = false;
    Doub scalar = 1.5;
    Doub rtol = 0.5;
    //cout << " \n SYSTEM SIZE  = " << KG.nrows() << std::endl;
    Doub rnorm = 10., rnormn=10.,lambn0=0.;
    Int counter = 0, maxcount = 30;
    do {
        Doub err1 = 10., err2 = 10., tol = 1.e-5;
        if ( counter<maxcount ) {
            lambn0 = lamb;
        }
        diff = 10;
        displace0 = displace;
        rnorm = 10.;
        Doub meantime=0.;
        counter = 0;
        do {

            chrono::steady_clock sc,sc2;
            // start timer
            auto start2 = sc2.now();
            fmesh->Assemble ( KG, FINT, FBODY );
            auto end2 = sc.now();
            auto time_span_assemble = static_cast<chrono::duration<double>> ( end2 - start2 );
            R = FBODY;
            R *= lamb;
            R -= FINT;

            InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );
            //SolveNR3(KG, R, dws);
            //SolveNR3(KG, FBODY, dwb);

            if ( false ) {
                string matrixout = "/home/diogocecilio/projects/dcproj/matrix.dat";
                std::ofstream matrixoutfile ( matrixout );
                PrintMathematicaFormat ( KG,matrixoutfile );
                string vecout = "/home/diogocecilio/projects/dcproj/vector.dat";
                std::ofstream vecoutfile ( vecout );
                PrintMathematicaFormat ( FBODY,vecoutfile );


                double startsolve = std::clock();
                SolveEigen ( KG, FBODY, dwb );
                double endsolve = ( std::clock() - startsolve ) / ( double ) CLOCKS_PER_SEC;
                cout << "\n time solving \n " << endsolve <<endl;

                startsolve = std::clock();
                SolveEigen2 ( KG, FBODY, dws );
                endsolve = ( std::clock() - startsolve ) / ( double ) CLOCKS_PER_SEC;

                cout << "\n time solving \n " << endsolve <<endl;
                // cout << "\n dws-=dwb \n "  <<endl;
                // dws-=dwb;
                // dws.Print();
                DebugStop();
                //return;

            }
            auto start = sc.now();
            bool llt=true;
            SolveEigenSparse ( llt, KG, R, dws );
            SolveEigenSparse ( llt, KG, FBODY, dwb );
            auto end = sc.now();
            auto time_span_solver = static_cast<chrono::duration<double>> ( end - start );
            //SolveEigen(KG, R, dws);
            //SolveEigen(KG, FBODY, dwb);
            Doub lastlambd = dlamb;
            dlamb = computelamda ( dwb, dws, dw, l );
            lamb += dlamb;
            if ( isnan ( dlamb ) == 1 || lamb<0 ) {
                std::cout << "NAN" << endl;
                std::cout << dlamb << endl;
                //fmesh->fmaterial->SetTangentMatrixType(false);
                break;
            }

            dww = dwb;
            dww *= dlamb;
            dww += dws;
            dw += dww;

            displace += dww;
            fmesh->fmaterial->UpdateDisplacement ( displace );

            rnorm = R.NRmatrixNorm();
            Doub normdw = dww.NRmatrixNorm();
            Doub unorm = displace.NRmatrixNorm();
            FBODY *= lamb;
            err1 = rnorm / FBODY.NRmatrixNorm();
            err2 = normdw / unorm;


            //cout << "Operation took: " << time_span.count() << " seconds !!!";

            //Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

            std::cout << " Time assemblimg :  " << time_span_assemble.count() <<  " Time inverting system :  " << time_span_solver.count() << " Iteration number = " << counter << " |  |R| = " << rnorm <<" | lamb  = " << lamb <<std::endl;
            //std::cout << " | time =  <<" << duration1 << std::endl;
            counter++;
            //meantime+=time_span;
            rtemp = rnorm;
            //if (counter == 1 || counter==0)rnorm = 10;


        } while ( counter < maxcount && rnorm > rtol );
        meantime/=counter;
        // Doub outtime = (std::clock() - startout) / (double)CLOCKS_PER_SEC;
        // diff2= fabs(lambn0- lamb);
        std::cout << " exter iter = " << counterout << "  | newton iters = " << counter  << " |  |R| = " << rnorm << " |  lamb  = " << lamb << std::endl;

        if ( rnorm > rtol ) {

            //std::cout << "Convergence failed. \n";

            counterout++;
            dws.assign ( sz, 1, 0. ), dwb.assign ( sz, 1, 0. ), dww.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), dw.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), FBODY.assign ( sz, 1, 0. ), FINT.assign ( sz, 1, 0. );

            fmesh->fmaterial->UpdateDisplacement ( displace0 );
            displace = displace0;
            //ndesi--;//74 s 0.497378 1.35437
            if ( check==true ) {
                dlamb0 *= 0.1/scalar;
                scalar *= 2.;
                //check = false;
                //if (rnorm> rnormn)check = false;
                // rnormn = rnorm;
            } else {
                dlamb0 *= 1/scalar;
                scalar *= 1.5;
                //check = true;
                //if (rnorm > rnormn)check = true;
                //rnormn = rnorm;
            }
            //MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
            dws.assign ( sz, 1, 0. ), dwb.assign ( sz, 1, 0. ), dww.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), dw.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), FBODY.assign ( sz, 1, 0. ), FINT.assign ( sz, 1, 0. );

            //mat->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
            fmesh->Assemble ( KG, FINT, FBODY );
            R = FBODY;
            R *= lamb;
            R -= FINT;
            InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );
            SolveEigen ( KG, R, dws );
            SolveEigen ( KG, FBODY, dwb );
            dwb.Transpose ( dwbt );
            dwbt.Mult ( dwb, mult );
            l0 = sqrt ( dlamb0 * dlamb0 * mult[0][0] );
            l = l0;
            dlamb = computelamda ( dwb, dws, dw, l );
            lamb = 0;
        } else {
            scalar = 1.;
            l *= Doub ( ndesi ) / Doub ( counter );
            if ( l > 10. ) {
                l = 10.;
            }
            diff2 = fabs ( lambn0 - lamb );
            fmesh->fmaterial->UpdatePlasticStrain();

            std::vector<double> sol ( 2 );
            sol[0]= fabs ( displace[2 * iddisplace[0] + 1][0] );
            sol[1] = lamb;
            //solpost.push_back(sol);
            std::vector<string> scalar_names;
            std::vector<string> vector_names;

            //  TPZStack<std::string> scalar_names,vector_names, tensor_names;
            vector_names.push_back ( "Displacement" );
            vector_names.push_back ( "Strain" );
            //vector_names.push_back("phi");
            //vector_names.push_back("Stress");
            Int dim=2;
            string slopestr="slope-IterativeProcessGIMArc";
            VTKGraphMesh vtkobj ( fmesh,dim,scalar_names,vector_names,slopestr );
            vtkobj.DrawSolution ( counterout );

        }

        counterout++;
        dws.assign ( sz, 1, 0. ), dwb.assign ( sz, 1, 0. ), dww.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), dw.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), FBODY.assign ( sz, 1, 0. ), FINT.assign ( sz, 1, 0. );

        solcount[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
        solcount[1] = lamb;
        solcount[2] = err1;
        solcount[3] = rnorm;
        solcount[4] = diff;
        solcount[5] = diff2;
        solcount[6] = counterout;

        uvf[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
        uvf[1] = lamb;

        solpost2.push_back ( uvf );
        solpost.push_back ( solcount );

        fmesh->fmaterial->UpdateDisplacement ( displace0 );
        std::vector<string> scalar_names;
        std::vector<string> vector_names;
        //  TPZStack<std::string> scalar_names,vector_names, tensor_names;
        vector_names.push_back ( "Displacement" );
        vector_names.push_back ( "Strain" );
        //vector_names.push_back("SqrtJ2(EPSP)");
        //vector_names.push_back("Stress");
        Int dim=2;
        string slopestr="/home/diogo/projects/dcproj/data/slope-saida-fina";
        VTKGraphMesh vtkobj ( fmesh,dim,scalar_names,vector_names,slopestr );
        vtkobj.DrawSolution ( counterout );

    } while ( counterout <= niter && fabs ( diff2 ) > alphatol ); // while (counterout <= maxcountout && fabs(diff2) > 0.05);


    if ( false ) {

        //string names = "fxu";
        //auto sss = std::to_string(cccccccc);
        //names += sss;
        //MatDoub solpost23;
        //solpost23.CopyFromVector(solpost2);
        //string exts = ".txt";
        //names += exts;
        //std::ofstream file8(names);
        //OutPutFile(solpost23, file8);

        std::cout << " GIM = " << std::endl;
        MatDoub solpost23;
        solpost23.CopyFromVector ( solpost2 );
        //string names = "/home/diogocecilio/Dropbox/slope-reliability/results/mathematicas-new/loadvsdisplacement2222";
        //string exts = ".dat";
        //names += exts;
        string names = "uxf.dat";
        std::ofstream file8 ( names );
        OutPutFile ( solpost23, file8 );


        std::vector<std::vector<double>> epsppost;
        fmesh->fmaterial->PostProcessIntegrationPointVar ( fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), epsppost );
        //string name3 = "/home/diogocecilio/Dropbox/slope-reliability/results/mathematicas-new/sqrtJ2GIM-deter2222";
        //string ext3 = ".dat";
        //name3 += ext3;
        string name3 = "/home/diogocecilio/projects/dcproj/data/sqrtJ2GIM.dat";
        std::ofstream file3 ( name3 );
        OutPutPost ( epsppost, file3 );

        string filename;
        std::vector<std::vector<double>> solx, soly;
        fmesh->fmaterial->PostProcess ( fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), solx, soly );
        filename = "/home/diogocecilio/projects/dcproj/data/soly.dat";
        std::ofstream file2 ( filename );
        OutPutPost ( soly, file2 );
        filename = "/home/diogocecilio/projects/dcproj/data/solx.dat";
        std::ofstream file22 ( filename );
        OutPutPost ( solx, file22 );
    }



//	std::cout << "  Iteration number = " << counterout <<  " |  |R| = " << rnorm << " | lamb  = " << lamb << std::endl;
    return solpost;

}


std::vector<std::vector<double>> slopeproject::IterativeProcessNew ( Int  ndesi, Doub dlamb0, Doub maxlfac, Int imc )
{
    std::vector<double> solcount ( 7, 0. ), uvf ( 2, 0. );
    std::vector<std::vector<double>> solpost, solpost2;
    solpost.push_back ( solcount );


    auto s = std::to_string ( imc );
    string filename = "/home/diogo/projects/dcproj/output/saida-gim-mc/saida-gim-mc";
    filename += s;
    filename+=".nb";
    bool print = false;
    std::ofstream file8 ( filename );
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    std::vector<std::vector<int>> idsvector;

	idsvector=fidsvector;
    //findbcids45 ( fmesh, idsvector );
   // findbcidsfatsn ( fmesh, idsvector );
    //findbcids ( fmesh, idsvector );
    //findbcidsfat ( fmesh, idsvector );
    //findbcidsfat5022 ( fmesh, idsvector );
    idsbottom = idsvector[0];
    idsleft = idsvector[1];
    idsright = idsvector[2];
    iddisplace = idsvector[3];

    Int sz = 2 * fmesh->GetMeshNodes().nrows();
    MatDoub KG ( sz, sz, 0. ), FBODY ( sz, 1, 0. ), FINT ( sz, 1, 0. ), ptsweigths;




    fmesh->fmaterial->ResetPlasticStrain();
    fmesh->fmaterial->ResetDisplacement();
    fmesh->fmaterial->ResetMat();




    MatDoub dws ( sz, 1, 0. ), dwb ( sz, 1, 0. ), dww ( sz, 1, 0. ), dw ( sz, 1, 0. ), R ( sz,1,0. ),displace ( sz,1,0. ),displace0 ( sz,1,0. );

    Doub rtol = 1.e-6;//tol newton
    Doub tol = 0.01;//tol arc length
    Doub rnorm = 10.;
    Doub lambn0 = 0.;
    //lambda  = load factor
    Doub lamb=0;
    Doub dlamb=0.;
    Doub diff2 =10.;
    //prescribed incremental solution length
    Doub l=0.;
    Int counter = 0;
    Int maxcount = 20;
    Int maxniter =20;
    Int counterout=0;

    fmesh->Assemble ( KG, FINT, FBODY );
    InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );
    Int type=0;
    SolveEigenSparse ( type,KG, FBODY, dwb );
    //initial length solution
    NRmatrix<Doub> dwbt,inner;
    dwb.Transpose ( dwbt );
    dwbt.Mult ( dwb,inner );
    Doub l0=sqrt ( dlamb0*dlamb0*inner[0][0] );
    Doub maxl =maxlfac*l0;
    l=l0;
    dlamb = dlamb0;
    chrono::steady_clock sc;
    auto start = sc.now();
//	fmesh->fmaterial->SetTangentMatrixType ( false );
    //cout << "| l0 = " << l0 << " | dlamb = " << dlamb << " | maxl = " << maxl <<endl;
    do {
        displace0 = displace;
        dws.assign ( sz, 1, 0. ), dwb.assign ( sz, 1, 0. ), dww.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), dw.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), FBODY.assign ( sz, 1, 0. ), FINT.assign ( sz, 1, 0. ),KG.assign ( sz,sz,0. );
        counter = 0;
        Doub err1=10.,err2=10.;
        rnorm=10.;
       // fmesh->Assemble ( KG, FINT, FBODY );
        do {
			

			fmesh->Assemble ( KG, FINT, FBODY );
     
           // 
            R = FBODY;
            R *= lamb;
            R -= FINT;

            InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );
            //NRmatrix<Doub> KGcp=KG;
            SolveEigenSparse ( type,KG, R, dws );
            //KG=KGcp;
            SolveEigenSparse ( type,KG, FBODY, dwb );
            //SolveEigen ( KG, R, dws );
            //SolveEigen ( KG, FBODY, dwb );

            if ( (counter==0 && counterout!=0)/*|| check==true*/ ) {
                dlamb = computelamda0 ( dwb,  dw, l );
                //dlamb = computelamda ( dwb, dws, dw, l );
                //std::cout << "l = "<< l << " | dlamb = "<< dlamb << " | lamb = "<< lamb  <<endl;
            } else {
                dlamb = computelamda ( dwb, dws, dw, l );
            }


            if ( isnan ( dlamb ) == 1 ) {
                std::cout << "isnan(dlamb)" << endl;
                std::cout << dlamb << endl;
                break;
            }
            lamb += dlamb;
            dww = dwb;
            dww *= dlamb;
            dww += dws;
            dw += dww;
            //dw delta maiusculo u
            //dww = delta minusculo u
            //dwb delta minusculo barra
            //dws delta minusculo star



            displace += dww;
            fmesh->fmaterial->UpdateDisplacement ( displace );

            rnorm = R.NRmatrixNorm();
            Doub normdw = dww.NRmatrixNorm();
            Doub unorm = displace.NRmatrixNorm();
            FBODY *= lamb;
            err1 = rnorm / FBODY.NRmatrixNorm();
            err2 = normdw / unorm;

          //  std::cout << " Iteration number = " << counter << " |  |R| = " << rnorm <<" | lamb  = " << lamb << " | dlamb = "<< dlamb  << "l = "<< l<< std::endl;

            counter++;
            //dlamb = computelamda0(dwb, dws, dw, l,lamb0);
            if ( rnorm>1000000 || lamb<0 ) {
                //fmesh->fmaterial->SetTangentMatrixType ( false );
                break;
            }

        } while ( counter < maxcount && err1 > rtol );


        //std::cout << " Iteration number = " << counter << " |  |R| = " << rnorm <<" | lamb  = " << lamb << " | dlamb = "<< dlamb  << "l = "<< l<< std::endl;
        if ( err1 > rtol || counter >= maxcount ) {
            std::cout << "  Reducing l... "  << " Iteration number = " << counter << " |  |R| = " << rnorm <<" | lamb  = " << lamb << " | dlamb = "<< dlamb  << "l = "<< l<< std::endl ;
            displace.assign(sz,1,0.);
            displace = displace0;
            fmesh->fmaterial->UpdateDisplacement ( displace0 );
            l*=0.5;
            //lamb=0;
 
           // counterout++;
            dws.assign ( sz, 1, 0. ), dwb.assign ( sz, 1, 0. ), dww.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), dw.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), FBODY.assign ( sz, 1, 0. ), FINT.assign ( sz, 1, 0. );

        } else {
 
            //fmesh->fmaterial->SetTangentMatrixType ( true );
            diff2 = fabs ( lambn0 - lamb );
            fmesh->fmaterial->UpdatePlasticStrain();
			//displace0 = displace;
			//fmesh->fmaterial->UpdateDisplacement ( displace0 );
            //fmesh->fmaterial->UpdateDisplacement(displace);
            l *= Doub ( ndesi ) / Doub ( counter );
            if ( l > maxl ) {
                l = maxl;
            }
            if ( l<0.1 ) {
                //l=0.1;
            }
            //l*=1.1;
            lambn0 = lamb;

            if ( print==true ) {
                PostVtk ( counter );
            }
            vector<Doub> uvf ( 2 );
            solcount[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
            solcount[1] = lamb;
            solcount[2] = err1;
            solcount[3] = rnorm;
            solcount[4] = 0;
            solcount[5] = diff2;
            solcount[6] = counterout;

            uvf[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
            uvf[1] = lamb;

            solpost2.push_back ( uvf );
            solpost.push_back ( solcount );
        }

        std::cout << " exter iter = " << counterout << "  | newton iters = " << counter  << " |  |R| = " << rnorm << " |  lamb  = " << lamb << "| l = "<< l  << "| //fabs(lambn0 - lamb) "<<diff2 <<std::endl;

        counterout++;
        fmesh->fmaterial->UpdateDisplacement ( displace0 );
    } while ( counterout <= maxniter && fabs ( diff2 ) > tol ); // while (counterout <= maxcountout && fabs(diff2) > 0.05);

    std::cout << " exter iter = " << counterout << "  | newton iters = " << counter  << " |  |R| = " << rnorm << " |  lamb  = " << lamb << "| l = "<< l  << "| fabs(lambn0 - lamb) "<<diff2 <<std::endl;
    auto end = sc.now();
    auto time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "| Total time taken " << time_span.count();
    std::cout << " Sucess exiting...  " <<endl;
    if ( print==false ) {
        OutPutPost ( solpost2,file8 );
    }
    return solpost;
}

void slopeproject::PostVtk ( Int step )
{

    std::vector<string> scalar_names;
    std::vector<string> vector_names;

    //  TPZStack<std::string> scalar_names,vector_names, tensor_names;
    vector_names.push_back ( "Displacement" );
    vector_names.push_back ( "Strain" );

    //vector_names.push_back("phi");
    //vector_names.push_back("Stress");
    Int dim=2;
    string slopestr="/home/diogo/projects/dcproj/output/saida-gim-mc/saida-gim-mc";
    VTKGraphMesh vtkobj ( fmesh,dim,scalar_names,vector_names,slopestr );
    vtkobj.DrawSolution ( step );
}

void  slopeproject::IterativeProcess2()
{

    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    std::vector<std::vector<int>> idsvector;
    std::vector<std::vector<double>> solpost2;


    findbcids ( fmesh, idsvector );

    idsbottom = idsvector[0];
    idsleft = idsvector[1];
    idsright = idsvector[2];
    iddisplace = idsvector[3];

    MatDoub KG,Fint,Fbody,FG;

    fmesh->Assemble ( KG, Fint, Fbody );

    int sz = KG.nrows();

    NRmatrix<Doub> u;
    Doub finalload =20.;
    u.assign ( sz,1,0. );
    // vector<double> fac = {0.1,1.,2.,3.,3.6,3.8,3.90,3.95,3.99,4.0};
    vector<double> fac = {0.1,0.4,0.6,0.8,1.};
    Int steps = 25;
    Int counterout = 1;
    MatDoub solpost ( 1000, 2, 0. );
    for ( Int iload = 0; iload < fac.size(); iload++ ) {
        std::cout << "load step = " << iload << " | fac = " <<fac[iload] << std::endl;
        Int counter = 0, maxcount = 30;
        Doub err1 = 10., err2 = 10., tol = 10.e-6;
        MatDoub dw ( sz, 1, 0. ), res ( sz, 1, 0. ), FINT,FBODY, R;
        fmesh->fmaterial->SetTangentMatrixType ( true );
        while ( counter <  maxcount && err1 > tol ) {


            fmesh->Assemble ( KG, FINT, FBODY );

            if ( false ) {
                string matrixout = "/home/diogo/projects/dcproj/matrix.dat";
                std::ofstream matrixoutfile ( matrixout );
                PrintMathematicaFormat ( KG,matrixoutfile );
            }
            // FBODY.Print();
            R = FBODY;
            R *= fac[iload];
            R -= FINT;

            InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );
            //R.Print();
            //bool llt=true;
            Int lu=1;
            SolveEigenSparse ( lu, KG, R, dw );
            //SolveEigen( KG, R, dw);
            u += dw;

            fmesh->fmaterial->UpdateDisplacement ( u );

            Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;

            rnorm = R.NRmatrixNorm();

            normdw = dw.NRmatrixNorm();

            normfg = FBODY.NRmatrixNorm();

            unorm = u.NRmatrixNorm();

            err1 = rnorm / ( fac[iload] *normfg );

            err2 = normdw / unorm;

            std::cout << " Iteration number = " << counter << " |  |R|/|FE| = " << err1 << " | deltau/u " << err2 << std::endl;
            counter++;
            fmesh->fmaterial->SetTangentMatrixType ( true );
        }
        fmesh->fmaterial->UpdatePlasticStrain();
        counterout++;
        solpost[iload][0] = fabs ( u[0][0] );
        solpost[iload][1] = fabs ( fac[iload] * finalload );

        vector<double> uvf ( 2 );
        uvf[0] = fabs ( u[2 * iddisplace[0] + 1][0] );
        uvf[1] = fac[iload];
        solpost2.push_back ( uvf );


        std::vector<double> sol ( 2 );
        sol[0]= fabs ( u[2 * iddisplace[0] + 1][0] );
        sol[1] = fabs ( fac[iload] * finalload );
        //solpost.push_back(sol);
        std::vector<string> scalar_names;
        std::vector<string> vector_names;

        //  TPZStack<std::string> scalar_names,vector_names, tensor_names;
        vector_names.push_back ( "Displacement" );
        vector_names.push_back ( "Strain" );
        //vector_names.push_back("phi");
        //vector_names.push_back("Stress");
        Int dim=2;
        string slopestr="slope-IterativeProcessGIMArc";
        VTKGraphMesh vtkobj ( fmesh,dim,scalar_names,vector_names,slopestr );
        vtkobj.DrawSolution ( counterout );

    }

    std::ofstream file8 ( "u-Fx-mc2.nb" );
    OutPutPost ( solpost2,file8 );

}
std::vector<std::vector<double>>    slopeproject::IterativeProcessGIMBinarySearch()
{
    fmesh->fmaterial->ResetPlasticStrain();
    fmesh->fmaterial->ResetDisplacement();
    fmesh->fmaterial->ResetMat();
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    std::vector<std::vector<int>> idsvector;

    findbcids ( fmesh, idsvector );

    idsbottom = idsvector[0];
    idsleft = idsvector[1];
    idsright = idsvector[2];
    iddisplace = idsvector[3];

    MatDoub KG,Fint,Fbody,FG;

    fmesh->Assemble ( KG, Fint, Fbody );

    int sz = KG.nrows();

    NRmatrix<Doub> u,u0;
    u.assign ( sz,1,0. );
    u0.assign ( sz,1,0. );
    Int counterout = 1;
    std::vector<std::vector<double>> solpost,solpost2;
    std::vector<double> solcount ( 7, 0. ), uvf ( 2, 0. );
    solpost.push_back ( solcount );
    bool converge=false;
    Doub fac = 0.1, facmin = 0., facmax =5.,factemp=0.,facmaxtemp=0.,facmintemp=0.;
    Int count2=0;
    do {
        std::cout << "load step x= " << counterout << " | fac = " <<fac  << std::endl;
        Int counter = 0, maxcount = 30,countercheck=0;
        Doub err1t=20.,err1 = 10., err2 = 10., tol = 1.e-6;
        MatDoub dw ( sz, 1, 0. ), res ( sz, 1, 0. ), FINT ( sz, 1, 0. ),FBODY ( sz, 1, 0. ), R ( sz, 1, 0. );
        KG.assign ( sz,sz,0. );
        //fmesh->fmaterial->SetTangentMatrixType(true);
        while ( counter <  maxcount && err1 > tol ) {




            fmesh->Assemble ( KG, FINT, FBODY );




            if ( false ) {
                string matrixout = "/home/diogocecilio/projects/dcproj/matrix.dat";
                std::ofstream matrixoutfile ( matrixout );
                PrintMathematicaFormat ( KG,matrixoutfile );
            }

            R = FBODY;
            R *= fac;
            R -= FINT;

            InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial );
            //R.Print();
            //SolveEigenSparse( KG, R, dw);
            //SolveEigenSparse( KG, R, dw);
            SolveEigen ( KG, R, dw );

            u += dw;

            fmesh->fmaterial->UpdateDisplacement ( u );

            Doub  normdw = 0., normfg = 0.,rnorm=0.,unorm=0.;

            rnorm = R.NRmatrixNorm();

            normdw = dw.NRmatrixNorm();

            normfg = FBODY.NRmatrixNorm();

            unorm = u.NRmatrixNorm();

            err1t=err1;
            err1 = rnorm / ( fac *normfg );

            // cout << "rnorm" << rnorm << endl;
            // cout << "normfg" << normfg << endl;
            //  cout << "normF|I" << FINT.NRmatrixNorm() << endl;
            err2 = normdw / unorm;

            //std::cout << " Iteration number = " << counter << " | rnorm "<<rnorm << " | unorm = "<<unorm <<" | dwnorm = "<<normdw << " |  |R|/|FE| = " << err1 << " | deltau/u " << err2 << "| factor = "<<fac<< std::endl;
            counter++;
            //if(err1>err1t)countercheck++;
            //if(countercheck>6)break;
            if ( rnorm>1000 ) {
                break;
            }
            if ( rnorm>100 ) {
                fmesh->fmaterial->SetTangentMatrixType ( true );
            } else  {
                fmesh->fmaterial->SetTangentMatrixType ( true ); //set the consistent tangent
            }

        }

        factemp = fac;
        facmaxtemp=facmax;
        facmintemp=facmin;
        if ( err1 > tol ||counter>=maxcount ) {
            cout <<" | facmax = "<<facmax <<" |facmin  " << facmin << " | fac = " << fac << " | (max-min)/fac="<< ( facmax - facmin ) / fac <<endl;
            //cout << " No Convergence  " << " | counter = " << counter <<endl;
            facmax = fac;
            fac = ( facmin + facmax ) / 2.;
            u=u0;
            fmesh->fmaterial->UpdateDisplacement ( u0 );
            //cout <<" | facmax = "<<facmax <<" |facmin  " << facmin << " | fac = " << fac <<endl;

            converge=false;

        } else {
            std::cout << " Iteration number = " << counter << " |  |R|/|FE| = " << err1 << " | deltau/u " << err2 << "| factor = "<<fac<< std::endl;
            cout <<"---CONV "<<" | facmax = "<<facmax <<" |facmin  " << facmin << " | fac = " << fac << " | (max-min)/fac="<< ( facmax - facmin ) / fac <<endl;
            if ( true ) {

                /*auto s2 = std::to_string(counterout);
                string filename = "/home/diogocecilio/Dropbox/slope-reliability/results/mesh-287-v2";
                std::vector<std::vector<double>> epsppost;
                fmesh->fmaterial->PostProcessIntegrationPointVar(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(),fmesh->fmaterial->GetSolution(), epsppost);
                string name3 = "/plasticsqrtj2";
                string ext3 = ".txt";
                filename += name3;
                filename += s2;
                filename += ext3;
                std::ofstream file3(filename);
                OutPutPost2(epsppost, file3);*/


                std::vector<double> sol ( 2 );
                sol[0]= fabs ( u[2 * iddisplace[0] + 1][0] );
                sol[1] = fac;
                solpost.push_back ( sol );
                std::vector<string> scalar_names;
                std::vector<string> vector_names;

                //  TPZStack<std::string> scalar_names,vector_names, tensor_names;
                vector_names.push_back ( "Displacement" );
                vector_names.push_back ( "Strain" );
                //vector_names.push_back("phi");
                //vector_names.push_back("Stress");
                Int dim=2;
                string slopestr="slope-IterativeProcessGIMBinarySearch";
                VTKGraphMesh vtkobj ( fmesh,dim,scalar_names,vector_names,slopestr );
                vtkobj.DrawSolution ( counterout );

            }
            solcount[0] = fabs ( u[2 * iddisplace[0] + 1][0] );
            solcount[1] = fac;
            solcount[2] = err1;
            solcount[3] = err2;
            solcount[4] = 0;
            solcount[5] = 0;
            solcount[6] = counterout;

            uvf[0] = fabs ( u[2 * iddisplace[0] + 1][0] );
            uvf[1] = fac;
            solpost2.push_back ( uvf );

            solpost.push_back ( solcount );

            u0=u;
            facmin=fac;
            fac = 1. / ( ( 1. / facmin + 1. / facmax ) / 2. );
            if ( fac>facmax ) {
                facmax=fac+0.5;
            }
            //if(facmax<=facmin)facmax=facmax+0.1;
            //fac +=0.5;
            fmesh->fmaterial->UpdatePlasticStrain();
            converge=true;

        }

        counterout++;







    }  while ( ( counterout<50&& ( facmax  - facmin ) / fac  > 0.001 ) || converge==false );

    std::ofstream file8 ( "u-Fx.nb" );
    OutPutPost ( solpost2,file8 );
    return solpost;
}




std::vector<std::vector<double>>   slopeproject::IterativeProcessArcLengthSRM ( int ndesi, Doub dlamb0, Doub alphatol, int niter )
{
    mesh* mesh = fmesh;
    std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
    std::vector<std::vector<int>> idsvector;

    findbcids ( mesh, idsvector );

    idsbottom = idsvector[0];
    idsleft = idsvector[1];
    idsright = idsvector[2];
    iddisplace = idsvector[3];

    Int sz = 2 * mesh->GetMeshNodes().nrows();
    MatDoub KG ( sz, sz, 0. ), FBODY ( sz, 1, 0. ), FINT ( sz, 1, 0. ), ptsweigths;
    Doub rtemp = 10.;

    material* mat = mesh->fmaterial;

    mat->ResetPlasticStrain();
    mat->ResetDisplacement();
    mat->ResetMat();

    MatDoub displace, displace0;
    displace.assign ( sz, 1, 0. );
    //Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
    Doub l = 0, l0 = 0, lamb = 1., lambn = 0, dlamb = 0., lamb3, diff = 100, diff2 = 100;
    Int counterout = 0, maxcountout = 20;
    std::vector<double> solcount ( 7, 0. ), uvf ( 2, 0. );
    std::vector<std::vector<double>> solpost, solpost2;
    solpost.push_back ( solcount );



    NRvector<Doub> matconsts ( 4, 0. );
    mat->GetMatConstants ( matconsts );
    Doub  young = matconsts[0];
    Doub nu = matconsts[1];
    Doub c = matconsts[2];
    Doub phi = matconsts[3];

    Doub c0 = c;
    Doub phi0 = phi;



    MatDoub dws ( sz, 1, 0. ), dwb ( sz, 1, 0. ), dww ( sz, 1, 0. ), dw ( sz, 1, 0. ), R;
    //mesh->fmaterial->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
    mesh->Assemble ( KG, FINT, FBODY );

    R = FINT;
    R -= FBODY;
    bool check = false;
    InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, mat );
    SolveEigen ( KG, R, dws );
    SolveEigen ( KG, FBODY, dwb );
    MatDoub dwbt, mult;
    dwb.Transpose ( dwbt );
    dwbt.Mult ( dwb, mult );
    l0 = sqrt ( dlamb0 * dlamb0 * mult[0][0] );
    l = l0;
    dlamb = computelamda ( dwb, dws, dw, l );
    lamb = 0;
    Doub scalar = 1.5;
    Doub rtol = 0.5;
    cout << " \n SYSTEM SIZE  = " << KG.nrows() << std::endl;
    cout << " \n dlamb  = " << dlamb << "\t l = "<< l << std::endl;
    Doub rnorm = 10., rnormn = 10.;
    do {
        std::clock_t start;
        double duration;
        start = std::clock();
        std::cout << "load step = " << counterout << " | l = " << l << " | diff2 = " << diff2 << std::endl;
        Int counter = 0, maxcount = 20;
        Doub err1 = 10., err2 = 10., tol = 1.e-5;
        Doub  lambn0 = lamb;
        diff = 10;
        displace0 = displace;
        rnorm = 10.;
        do {
            //c = c0 / lamb;
            //phi = atan(tan(phi0) / lamb);
            //matconsts[2] = c;
            //matconsts[3] = phi;
            //mat->SetMatConstants(matconsts);

            mesh->Assemble ( KG, FINT, FBODY );
            R = FINT;
            //R *= lamb;
            R -= FBODY;

            InserBC ( KG, R, FBODY, idsbottom, idsright, idsleft, mat );
            SolveEigen ( KG, R, dws );
            SolveEigen ( KG, FBODY, dwb );
            dlamb = computelamda ( dwb, dws, dw, l );
            lamb += dlamb;
            dww = dwb;
            dww *= dlamb;
            dww += dws;
            dw += dww;

            displace += dww;
            mat->UpdateDisplacement ( displace );

            rnorm = R.NRmatrixNorm();
            Doub normdw = dww.NRmatrixNorm();
            Doub unorm = displace.NRmatrixNorm();
            err2 = normdw / unorm;

            Doub duration1 = ( std::clock() - start ) / ( double ) CLOCKS_PER_SEC;
            c = c0 / lamb;
            phi = atan ( tan ( phi0 ) / lamb );
            matconsts[2] = c;
            matconsts[3] = phi;
            mat->SetMatConstants ( matconsts );
            std::cout << " Iteration number = " << counter << " | phi0/phi = " << tan ( phi0 ) / tan ( phi ) << " | c0/c = " << c0 / c << " | c = " << c << " | phi = " << phi << " | lamb  = " << lamb << " | time =  <<" << duration1 << std::endl;
            counter++;

            rtemp = rnorm;
            if ( counter == 1 ) {
                rnorm = 10;
            }


        } while ( counter < maxcount && rnorm > rtol );


        c = c0 / lamb;
        phi = atan ( tan ( phi0 ) / lamb );
        matconsts[2] = c;
        matconsts[3] = phi;
        mat->SetMatConstants ( matconsts );
        mat->UpdatePlasticStrain();
        counterout++;

        std::cout << " iter = " << counterout << " | FS= " << lamb << std::endl;

        dws.assign ( sz, 1, 0. ), dwb.assign ( sz, 1, 0. ), dww.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), dw.assign ( sz, 1, 0. ), R.assign ( sz, 1, 0. ), FBODY.assign ( sz, 1, 0. ), FINT.assign ( sz, 1, 0. );

        solcount[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
        solcount[1] = lamb;
        solcount[2] = err1;
        solcount[3] = rnorm;
        solcount[4] = diff;
        solcount[5] = diff2;
        solcount[6] = counterout;

        uvf[0] = fabs ( displace[2 * iddisplace[0] + 1][0] );
        uvf[1] = lamb;

        solpost2.push_back ( uvf );
        solpost.push_back ( solcount );

        //std::cout << " $$$$$ Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm <<  " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
    } while ( counterout <= niter && fabs ( diff2 ) > alphatol ); // while (counterout <= maxcountout && fabs(diff2) > 0.05);


    if ( true ) {

        //string names = "fxu";
        //auto sss = std::to_string(cccccccc);
        //names += sss;
        //MatDoub solpost23;
        //solpost23.CopyFromVector(solpost2);
        //string exts = ".txt";
        //names += exts;
        //std::ofstream file8(names);
        //OutPutFile(solpost23, file8);


        MatDoub solpost23;
        solpost23.CopyFromVector ( solpost2 );
        string names = "loadvsdisplacementnew";
        string exts = ".txt";
        names += exts;
        std::ofstream file8 ( names );
        OutPutFile ( solpost23, file8 );


        std::vector<std::vector<double>> epsppost;
        mat->PostProcessIntegrationPointVar ( mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), mat->GetSolution(), epsppost );
        string name3 = "epsppostnew";
        string ext3 = ".txt";
        name3 += ext3;
        std::ofstream file3 ( name3 );
        OutPutPost ( epsppost, file3 );

        string filename;
        std::vector<std::vector<double>> solx, soly;
        mat->PostProcess ( mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), mat->GetSolution(), solx, soly );
        filename = "soly.txt";
        std::ofstream file2 ( filename );
        OutPutPost ( soly, file2 );
        filename = "solx.txt";
        std::ofstream file22 ( filename );
        OutPutPost ( solx, file22 );
    }

    return solpost;

}

Doub  slopeproject::computelamda0 ( MatDoub& dwb,  MatDoub& dw, Doub& l )
{
    Doub tempsign;
    MatDoub dwt,solsig;
    dw.Transpose ( dwt );
    dwt.Mult ( dwb,solsig );
    Doub scal = solsig[0][0];
    Doub signum=0;
    //page 111, eq. 4.123 - Souza Neto //verificar sinal
    if ( scal>0 ) {
        signum=-1;
    } else {
        signum=1;
    }

    MatDoub dwbt, aparam;
    dwb.Transpose ( dwbt );
    dwbt.Mult ( dwb, aparam );

    return signum*l/sqrt ( aparam[0][0] ) ;


}




Doub  slopeproject::computelamda ( MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l )
{



    Int sz = dwb.nrows();
    Doub aa = 0.;
    for ( Int i = 0; i < sz; i++ ) {
        aa += dwb[i][0] * dwb[i][0];
    }
    Doub bb = 0.;
    MatDoub dwcopy = dw;
    dwcopy += dws;
    for ( Int i = 0; i < sz; i++ ) {
        bb += dwb[i][0] * dwcopy[i][0];
    }
    bb *= 2;
    Doub cc = 0.;
    for ( Int i = 0; i < sz; i++ ) {
        cc += dwcopy[i][0] * dwcopy[i][0];
    }

    cc -= l * l;
    Doub delta = bb * bb - 4. * aa * cc;
    Doub dlamb2;
    Doub dlamb1;

    //cout << "delta" <<delta << " aa = "<< aa << " bb = "<< bb << " cc = "<< cc << endl;

    if ( fabs ( aa ) >1.e-12 && delta>0 ) {
        dlamb2 = ( -bb + sqrt ( delta ) ) / ( 2. * aa ); //maior
        dlamb1= ( -bb - sqrt ( delta ) ) / ( 2. * aa ); //menor
        //return dlamb1;
        //cout << "dlamb1" <<dlamb1 << " dlamb2 = "<< dlamb2 << endl;
    } else if ( bb != 0 ) {
        //cout << "-cc/bb" <<-cc/bb << endl;
        return -cc/bb;
    } else {
        //cout << "(-bb ) / (2. * aa)" <<(-bb ) / (2. * aa)<< endl;
        return ( -bb ) / ( 2. * aa );
    }


    //page 111, eq. 4.118 - Souza Neto
    MatDoub temp1,temp1t,sol1,temp2,temp2t,sol2;
    temp1=dwb;
    temp1*=dlamb1;
    temp1+=dws;
    temp1+=dw;
    temp1.Transpose ( temp1t );
    temp1t.Mult ( dw,sol1 );

    temp2=dwb;
    temp2*=dlamb2;
    temp2+=dws;
    temp2+=dw;
    temp2.Transpose ( temp2t );
    temp2t.Mult ( dw,sol2 );



    if ( sol1[0][0]>sol2[0][0] ) {
        //cout << "return 1 " << endl;
        return dlamb1;
    } else {
        //cout << "return 2 " << endl;
        return dlamb2;
    }

    // }



    //return 0.1;
    /* // Doub rtsafe(T &funcd, const Doub x1, const Doub x2, const Doub xacc,bool &flag)
      quadraticfunction func;
      Doub x1=-10.;
      Doub x2=0.;
      Doub tol=1.e-6;
      func.Set(aa,bb,cc);

      Doub dlamb2 = rtsafe(func,x1,x2,tol);
      cout << " a = "<<aa<< " b = "<< bb <<  " c = "<<cc<< endl;
      cout << "dlamb = "<<dlamb<< endl;
      cout << "dlamb2 = "<<dlamb2<< endl;
      return dlamb;*/


    /*Int sz = dwb.nrows();
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
    Doub dlamb = (-bb + sqrt(delta)) / (2. * aa);
    return dlamb;*/

}

/*
void slopeproject::InserBC(SparseMatrix<double> & KG, VectorXd& R, VectorXd& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material* mat)
{

	//FBODY *= 1. / lamb;
	Int dir, val;
	dir = 1;
	val = 0;
	mat->DirichletBC(KG, R, idsbottom, dir, val);

	dir = 0;
	val = 0;
	mat->DirichletBC(KG, R, idsbottom, dir, val);


	dir = 0;
	val = 0;
	mat->DirichletBC(KG, R, idsright, dir, val);
	dir = 0;
	val = 0;
	mat->DirichletBC(KG, R, idsleft, dir, val);


	dir = 1;
	val = 0;
	mat->DirichletBC(KG, FBODY, idsbottom, dir, val);
	dir = 0;
	val = 0;
	mat->DirichletBC(KG, FBODY, idsbottom, dir, val);

	dir = 0;
	val = 0;
	mat->DirichletBC(KG, FBODY, idsright, dir, val);
	dir = 0;
	val = 0;
	mat->DirichletBC(KG, FBODY, idsleft, dir, val);
}*/


void slopeproject::InserBC ( MatDoub& KG, MatDoub& R, MatDoub& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material* mat )
{
    //FBODY *= 1. / lamb;
    Int dir, val;
    dir = 1;
    val = 0;
    mat->DirichletBC ( KG, R, idsbottom, dir, val );

    dir = 0;
    val = 0;
    mat->DirichletBC ( KG, R, idsbottom, dir, val );


    dir = 0;
    val = 0;
    mat->DirichletBC ( KG, R, idsright, dir, val );
    dir = 0;
    val = 0;
    mat->DirichletBC ( KG, R, idsleft, dir, val );


    dir = 1;
    val = 0;
    mat->DirichletBC ( KG, FBODY, idsbottom, dir, val );
    dir = 0;
    val = 0;
    mat->DirichletBC ( KG, FBODY, idsbottom, dir, val );

    dir = 0;
    val = 0;
    mat->DirichletBC ( KG, FBODY, idsright, dir, val );
    dir = 0;
    val = 0;
    mat->DirichletBC ( KG, FBODY, idsleft, dir, val );
}


void slopeproject::SolveEigen3 ( MatDoub A, MatDoub b, MatDoub& x )
{


    x.assign ( A.nrows(), 1, 0. );
    MatrixXd AA ( A.nrows(), A.nrows() );
    VectorXd bbb ( A.nrows() );
    for ( int i = 0; i < A.nrows(); i++ ) {
        for ( int j = 0; j < A.ncols(); j++ ) {
            AA ( i, j ) = A[i][j];
        }
    }
    for ( int i = 0; i < A.nrows(); i++ ) {
        bbb ( i ) = b[i][0];
    }
    //VectorXd xxx = AA.llt().solve ( bbb ); //mais r??pido
    //VectorXd xxx = AA.fullPivHouseholderQr().solve(bbb);
    //VectorXd xxx = AA.fullPivLu().solve(bbb);
    //VectorXd xxx = AA.ldlt().solve(bbb);
    VectorXd xxx = AA.lu().solve ( bbb );
    for ( int i = 0; i < A.nrows(); i++ ) {
        x[i][0] = xxx ( i );
    }
}

void slopeproject::SolveEigen ( MatDoub A, MatDoub b, MatDoub& x )
{

    x.assign ( A.nrows(), 1, 0. );
    MatrixXd AA ( A.nrows(), A.nrows() );
    VectorXd bbb ( A.nrows() );
    for ( int i = 0; i < A.nrows(); i++ ) {
        for ( int j = 0; j < A.ncols(); j++ ) {
            AA ( i, j ) = A[i][j];
        }
        bbb ( i ) = b[i][0];
    }
    LLT<MatrixXd> llt;
    //LU<MatrixXd> lu;
    llt.compute ( AA );
    VectorXd xxx=llt.solve ( bbb );
    for ( int i = 0; i < A.nrows(); i++ ) {
        x[i][0] = xxx ( i );
    }
}

#include<Eigen/IterativeLinearSolvers>
void slopeproject::SolveEigenSparse ( int type, MatDoub A, MatDoub b, MatDoub& x )
{

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int sz=A.nrows();

    tripletList.reserve ( sz*200 );
    // tripletList.reserve(80000);

    x.assign ( sz, 1, 0. );
    SparseMatrix<double> AA ( sz, sz );
    VectorXd bbb ( sz );
    for ( int i = 0; i < sz; i++ ) {
        for ( int j = 0; j < sz; j++ ) {
            if ( fabs ( A[i][j] ) >1.e-12) {
                tripletList.push_back ( T ( i,j,A[i][j] ) );
            }
        }
        bbb ( i ) = b[i][0];
    }

    AA.setFromTriplets ( tripletList.begin(), tripletList.end() );

    AA.makeCompressed();

    VectorXd xx;
    if ( type==0 ) {
        SimplicialLLT< SparseMatrix<double> > solver;
		solver.analyzePattern ( AA );
        solver.factorize ( AA );
		xx = solver.solve ( bbb );
       // xx = solver.compute ( AA ).solve ( bbb );
    } else if ( type==1 ) {
        SparseLU< SparseMatrix<double> > solver2;
        solver2.analyzePattern ( AA );
        solver2.factorize ( AA );
        xx = solver2.solve ( bbb );
    } else if ( type==2 ) {
        ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
        cg.compute ( AA );
        xx = cg.solve ( bbb );
        std::cout << "#iterations:     " << cg.iterations() << std::endl;
        std::cout << "estimated error: " << cg.error()      << std::endl;
    }
    else if ( type==3 ) {
        SimplicialLDLT<SparseMatrix<double> > solver;
		solver.analyzePattern ( AA );
        solver.factorize ( AA );
		xx = solver.solve ( bbb );
    }
    for ( int i=0; i<sz; i++ ) {
        x[i][0]=xx ( i );
    }
}


// void beam3dtools::SolveEigenSparse ( MatDoub A, MatDoub b, MatDoub& x )
// {
// 
//     typedef Eigen::Triplet<double> T;
//     std::vector<T> tripletList;
//     int sz=A.nrows();
// 
//     //30 = approx number of nonzero term in K matrix
//     tripletList.reserve ( sz*30 );
//     // tripletList.reserve(80000);
// 
//     x.assign ( sz, 1, 0. );
//     SparseMatrix<double> AA ( sz, sz );
//     VectorXd bbb ( sz );
//     for ( int i = 0; i < sz; i++ ) {
//         for ( int j = 0; j < sz; j++ ) {
//             if ( fabs ( A[i][j] ) >1.e-12 ) {
//                 tripletList.push_back ( T ( i,j,A[i][j] ) );
//             }
//         }
//         bbb ( i ) = b[i][0];
//     }
// 
//     AA.setFromTriplets ( tripletList.begin(), tripletList.end() );
// 
//     AA.makeCompressed();
//     SimplicialLLT< SparseMatrix<double> > solver;
//     VectorXd xx = solver.compute ( AA ).solve ( bbb );
//     for ( int i=0; i<sz; i++ ) x[i][0]=xx ( i );
// }

void  slopeproject::ReadMesh ( std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord )
{
    std::vector<std::vector<Int>> topol;
    string line, temp;

    ifstream myfile ( filenameel );
    //
    if ( myfile.is_open() ) {
        while ( getline ( myfile, line ) ) {
            std::vector<string> tokens;
            istringstream iss ( line );
            while ( iss >> temp ) {
                tokens.push_back ( temp );
            }
            std::vector<Int> input_int = vecstr_to_vecint ( tokens );
            for ( int k = 0; k < input_int.size(); k++ ) {
                input_int[k] = input_int[k] - 1;
            }
            topol.push_back ( input_int );
        }
        myfile.close();
    } else {
        std::cout << "Unable to open file";
    }

    meshtopology.CopyFromVector ( topol );
    //meshtopology.Print();

    std::vector<std::vector<Doub>> coords;
    string line2, temp2;
    ifstream myfile2 ( filenamecoord );
    if ( myfile2.is_open() ) {
        while ( getline ( myfile2, line2 ) ) {
            std::vector<string> tokens;
            istringstream iss ( line2 );
            while ( iss >> temp2 ) {
                tokens.push_back ( temp2 );
            }
            std::vector<Doub> input_doub = vecstr_to_vecdoub ( tokens );

            //std::vector<Doub> input_doub2(input_doub.size() - 1);
            //for (int k = 1;k < input_doub.size();k++)
            //{
            //	input_doub2[k] = input_doub[k];
            //}

            coords.push_back ( input_doub );
        }
        myfile2.close();
    } else {
        std::cout << "Unable to open file";
    }

    meshcoords.CopyFromVector ( coords );
    meshcoords.Print();


    std::vector<Doub> temp33 ( 3 );
    for ( Int i = 0; i < meshtopology.nrows(); i++ ) {
        std::vector< std::vector<Doub> > temp22;
        for ( Int j = 0; j < meshtopology.ncols(); j++ ) {
            Int top = meshtopology[i][j];
            temp33[0] = meshcoords[top][0];
            temp33[1] = meshcoords[top][1];
            temp33[2] = meshcoords[top][2];
            temp22.push_back ( temp33 );
        }
        allcoords.push_back ( temp22 );
    }



}
std::vector<Int>  slopeproject::vecstr_to_vecint ( std::vector<string> vs )
{
    std::vector<Int> ret;
    for ( std::vector<string>::iterator it = vs.begin() + 1; it != vs.end(); ++it ) {
        istringstream iss ( *it );
        Int temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}

std::vector<Doub>  slopeproject::vecstr_to_vecdoub ( std::vector<string> vs )
{
    std::vector<Doub> ret;
    for ( std::vector<string>::iterator it = vs.begin() ; it != vs.end() ; ++it ) {
        istringstream iss ( *it );
        Doub temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}

/*std::vector<Doub>  slopeproject::vecstr_to_vecdoub(std::vector<string> vs)
{
	std::vector<Doub> ret;
	for (std::vector<string>::iterator it = vs.begin() + 1; it != vs.end() - 1; ++it)
	{
		istringstream iss(*it);
		Doub temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}*/

std::vector<Doub>  slopeproject::vecstr_to_vecdoub2 ( std::vector<string> vs )
{
    std::vector<Doub> ret;
    for ( std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it ) {
        istringstream iss ( *it );
        Doub temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}

template <class T>
std::vector<T>  slopeproject::vecstr_to_vec ( std::vector<string> vs )
{
    std::vector<T> ret;
    for ( std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it ) {
        istringstream iss ( *it );
        T temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}

void  slopeproject::ReadMatDoub ( MatDoub& matdoub, std::string  file )
{
    std::vector<std::vector<Doub>> coords;
    string line2, temp2;
    ifstream myfile2 ( file );
    if ( myfile2.is_open() ) {
        while ( getline ( myfile2, line2 ) ) {
            std::vector<string> tokens;
            istringstream iss ( line2 );
            while ( iss >> temp2 ) {
                tokens.push_back ( temp2 );
            }
            std::vector<Doub> input_doub = vecstr_to_vecdoub2 ( tokens );
            coords.push_back ( input_doub );
        }
        myfile2.close();
    } else {
        std::cout << "Unable to open file";
    }
    //for (int i = 0;i < coords.size();i++)
    //{
    //	for (int j = 0;j < coords[0].size();j++)
    //	{
    //		cout << coords[i][j] << endl;
    //	}
    //	cout << endl;
    //}

    matdoub.CopyFromVector ( coords );
}
void slopeproject::GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub& elcoords )
{

    elcoords.assign ( allcoords[el].size(), 2, 0. );
    Int sz = allcoords[el].size();
    for ( Int j = 0; j < sz; j++ ) {
        Doub x = allcoords[el][j][0];
        Doub y = allcoords[el][j][1];
        elcoords[j][0] = x;
        elcoords[j][1] = y;
    }
}

void slopeproject::FindIdsInPath ( const MatDoub& path, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& idpath )
{
    MatDoub elcoords;
    int nels = allcoords.size();
    GetElCoords ( allcoords, 0, elcoords );
    Int nnodes = elcoords.nrows();
    for ( Int iel = 0; iel < nels; iel++ ) {
        GetElCoords ( allcoords, iel, elcoords );
        for ( Int inode = 0; inode < nnodes; inode++ ) {
            Doub x = elcoords[inode][0];
            Doub y = elcoords[inode][1];

            for ( Int ipath = 0; ipath < path.nrows(); ipath++ ) {
                Doub copathx = path[ipath][0];
                Doub copathy = path[ipath][1];

                if ( fabs ( x - copathx ) < 10.e-4 && fabs ( y - copathy ) < 10.e-4 ) {
                    idpath.push_back ( meshtopology[iel][inode] );
                    ipath = path.nrows();
                }
            }
        }

    }

    sort ( idpath.begin(), idpath.end() );
    idpath.erase ( unique ( idpath.begin(), idpath.end() ), idpath.end() );

}

void slopeproject::Line ( VecDoub a, VecDoub b, Int ndivs, MatDoub& path )
{
    Doub x0 = a[0];
    Doub xf = b[0];

    Doub y0 = a[1];
    Doub yf = b[1];

    Doub dx = ( xf - x0 ) / ndivs;
    Doub dy = ( yf - y0 ) / ndivs;

    path.assign ( ndivs, 2, 0. );

    for ( Int idiv = 0; idiv < ndivs; idiv++ ) {
        path[idiv][0] = x0 + idiv * dx;
        path[idiv][1] = y0 + idiv * dy;
    }


}

MatDoub slopeproject::AssembleHhationho ( Int i )
{

    NRmatrix<MatDoub> randomfield = frandomfield;
    Int nrandomvars = randomfield.nrows();
    Int nrowss = randomfield[0][0].nrows();
    VecDoub mean ( nrandomvars, 0. ), var ( nrandomvars, 0. );
    MatDoub hhatinho ( randomfield[0][0].nrows(), randomfield.nrows(), 0. ), posttest ( randomfield[0][0].nrows(), 1, 0. );
    MatDoub coesmat ( randomfield[0][0].nrows(), 1, 0. ), phimat ( randomfield[0][0].nrows(), 1, 0. );

    for ( Int ivar = 0; ivar < nrandomvars; ivar++ ) {
        for ( Int isample = 0; isample < nrowss; isample++ ) {
            hhatinho[isample][ivar] = randomfield[ivar][0][isample][i];
            //hhatinho[i][0] = 0.;
            mean[ivar] += hhatinho[isample][ivar];
            mean[ivar] /= ( isample + 1 );
            var[ivar] += ( mean[ivar] - hhatinho[isample][ivar] ) * ( mean[ivar] - hhatinho[isample][ivar] );
            posttest[isample][0] = hhatinho[isample][ivar];
        }
    }

    return hhatinho;
}

void slopeproject::MonteCarloSRM ( int iter,int iter2, bool print, string writenamefolder )
{

    MatDoub ptsweigths;

    Int order;
    Int elnodes =fmesh->GetMeshTopology().ncols();
    if ( elnodes==3 ) {
        order=1;
        shapetri shape = shapetri ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }
    if ( elnodes  ==4 ) {
        order=1;
        shapequad shape = shapequad ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }
    if ( elnodes==6 ) {
        order=2;
        shapetri shape = shapetri ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }
    if ( elnodes==8 ) {
        order=2;
        shapequad shape = shapequad ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }
    Int npts = ptsweigths.nrows();
    Int nglobalpts = fmesh->GetMeshTopology().nrows() * npts;
    Int sz = 2 * fmesh->GetMeshNodes().nrows();

    //material* mat = fmesh->fmaterial;


    NRvector<Doub> matconsts ( 4, 0. );
    fmesh->fmaterial->GetMatConstants ( matconsts );
    Doub  young = matconsts[0];
    Doub nu = matconsts[1];
    Doub c = matconsts[2];
    Doub phi = matconsts[3];
    MatDoub bodyforce = fmesh->fmaterial->GetBodyForce();

    string namefolder = writenamefolder;
    char* cstr = new char[namefolder.length() + 1];
    strcpy ( cstr, namefolder.c_str() );

#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
    int check = mkdir ( cstr, 777 );
#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
    int check = mkdir ( cstr );
#endif


    std::vector<double> solvec;
    int samples = frandomfield[0][0].ncols();
    MatDoub solpost ( samples, 2, 0. ), solpost2 ( samples, 1, 0. );
    int fail = 0;


//	std::vector<int> copyv = { 8, 27, 34, 49, 63, 80, 93, 151, 203, 229, 231, 301, 339, 346, 380,
//400, 522, 530, 542, 609, 655, 669, 682, 868, 872, 916, 939, 941, 990,
//1053, 1072, 1083, 1115, 1146, 1309, 1379, 1409, 1446, 1458, 1567,
//1568, 1638, 1655, 1787, 1829, 1927, 1949, 1955, 1959, 1967, 1983,
//2008, 2077 };

//	std::vector<int> copyv = { 22, 32, 42, 70, 72, 162, 175, 207, 210, 330, 379, 421, 440 };


    //std::cout << " \n fail size = " <<copyv.size() << endl;
    int idfail = 0;
    for ( int i = iter; i < iter2; i++ ) {



        //idfail = *find(copyv.begin(), copyv.end(), i);
        //if (idfail!=0 &&fabs(idfail)<2000)
        if ( false ) {
            std::cout << " \n fail = "  << endl;
            std::cout << " \n ID fail = " << idfail <<endl;
            idfail = 0;
            //continue;
        } else {



            //double min = *min_element(solvec.begin(), solvec.end());
            std::clock_t start;
            double duration;
            start = std::clock();

            MatDoub hhatinho = AssembleHhationho ( i );
            fmesh->fmaterial->ResetMat();
            fmesh->SetHhat ( hhatinho );
            fmesh->fmaterial->SetMemory ( nglobalpts, sz );
            fmesh->fmaterial->SetMatConstants ( matconsts );
            fmesh->fmaterial->UpdateBodyForce ( bodyforce );

            std::vector<std::vector<double>>  sol = IterativeProcessShearRed ( 0.5, 2.,0.01 );

            if ( print==true ) {
                PostVtk ( i );
            }

            MatDoub solpost23;
            solpost23.CopyFromVector ( sol );


            Int last = solpost23.nrows() - 1;
            Doub data = solpost23[last][1];
            solvec.push_back ( data );

            string  filename = namefolder;
            string datafile = "/information";
            string ext = ".txt";
            filename += datafile;
            auto s = std::to_string ( i );
            filename += s;
            filename += ext;
            std::ofstream fileinfo ( filename );
            fileinfo << "Monte Carlo Sample = " << i << std::endl;
            fileinfo << "Safety Factor = " << data << std::endl;
            fileinfo << "r/normr = " << solpost23[last][2] << std::endl;
            fileinfo << "u/normu = " << solpost23[last][3] << std::endl;
            fileinfo << "diff= " << solpost23[last][4] << std::endl;
            fileinfo << "diff2 = " << solpost23[last][5] << std::endl;
            fileinfo << "counterout = " << solpost23[last][6] << std::endl;


            if ( data <= 1. ) {
                fail++;
            }
            std::cout << " mc it = " << i << " | Current safety fator = " << data << endl;
            solpost2[i][0] = data;
            //delete mesh2;
            //delete mat;

            duration = ( std::clock() - start ) / ( double ) CLOCKS_PER_SEC;
            std::cout << "\n time in monte carlo iteration = " << duration << '\n';
            fileinfo << "CPU time = " << duration << std::endl;

        }
    }
    //delete mesh2;
    //delete objKLGalerkinRF;
    // delete mat;
}



void slopeproject::MonteCarloGIM ( int iter, int iter2, bool print, string writenamefolder )
{


    MatDoub ptsweigths;

    Int order;
    Int elnodes =fmesh->GetMeshTopology().ncols();
    if ( elnodes==3 ) {
        order=1;
        shapetri shape = shapetri ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }
    if ( elnodes  ==4 ) {
        order=1;
        shapequad shape = shapequad ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }
    if ( elnodes==6 ) {
        order=2;
        shapetri shape = shapetri ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }
    if ( elnodes==8 ) {
        order=2;
        shapequad shape = shapequad ( order, 1 );
        shape.pointsandweigths ( ptsweigths );
    }

    Int npts = ptsweigths.nrows();
    Int nglobalpts = fmesh->GetMeshTopology().nrows() * npts;
    Int sz = 2 * fmesh->GetMeshNodes().nrows();
    //material* materialdp = finemesh->fmaterial;

    NRvector<Doub> matconsts ( 4, 0. );
    fmesh->fmaterial->GetMatConstants ( matconsts );
    Doub  young = matconsts[0];
    Doub nu = matconsts[1];
    Doub c = matconsts[2];
    Doub phi = matconsts[3];
    MatDoub bodyforce = fmesh->fmaterial->GetBodyForce();
    chrono::steady_clock sc;
    auto start = sc.now();
    std::cout << "\n starting stochastic simulation " << endl;

    KLGalerkinRF* objKLGalerkinRF = fklgalerking;


    string namefolder = writenamefolder;

    char* cstr = new char[namefolder.length() + 1];
    strcpy ( cstr, namefolder.c_str() );


#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
    int	check = mkdir ( cstr, 777 );
#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
    int check = mkdir ( cstr );
#endif

    string datafile = namefolder;
    datafile += "/DATA.txt";
    std::ofstream file ( datafile );
    file << " Young = " << young << " | nu = " << nu << endl;
    file << " c = " << c << " | phi = " << phi << endl;
    file << " bodyforce = " << bodyforce[1][0] << endl;
    file << " samples = " << fklgalerking->GetSamples() << " | expansion order = " << fklgalerking->GetExpansionorder() << " | func type = " << fklgalerking->Getftype() << endl;
    file << "Lx = " << fklgalerking->GetLx() << " | Ly = " << fklgalerking->GetLy() << " variance = " << fklgalerking->GetSig() << endl;


    Int  fail = 0;
    std::cout << "\n starting Monte Carlo " << endl;


    NRmatrix<MatDoub> randomfield = frandomfield;
    int samples = randomfield[0][0].ncols();
	if(samples==0)
 	{
		cout << "Zero samples in MC "<<endl;
		randomfield[0][0].Print();
		DebugStop(); 
	}
    MatDoub solpost ( samples, 2, 0. ), solpost2 ( samples, 1, 0. );
    std::vector<double> solvec;
    std::cout << "\n samples " << samples << endl;



    for ( Int imc = iter; imc < iter2; imc++ ) {





        MatDoub hhatinho = AssembleHhationho ( imc );

        NRvector<Doub> matconsts ( 4, 0. );
        matconsts[0] = young;
        matconsts[1] = nu;
        matconsts[2] = c;
        matconsts[3] = phi;

        fmesh->fmaterial->SetMatConstants ( matconsts );
        fmesh->fmaterial->SetMemory ( nglobalpts, sz );
        fmesh->fmaterial->UpdateBodyForce ( bodyforce );
        fmesh->fmaterial->SetRandomField ( hhatinho );
        fmesh->SetHhat ( hhatinho );

		//c=30 phi=10
        //int desirediter = 7;
       //Doub dlamb0 =0.465;
        //Doub maxlfac=1.;
		//c=30 phi=20
		int desirediter = 8;
        Doub dlamb0 =0.6;
		Doub maxlfac=0.9;
		
        //10, 0.5, 0.01,20
        std::vector<std::vector<double>>  sol = IterativeProcessNew( desirediter, dlamb0,maxlfac,imc );
        // std::vector<std::vector<double>>  sol = IterativeProcessGIMBinarySearch();
        if ( print==true ) {
            PostVtk ( imc );
        }

        MatDoub solpost23;
        solpost23.CopyFromVector ( sol );
        Int last = solpost23.nrows() - 1;
        Doub data = solpost23[last][1];

        if ( data < 0.2 ) {
            continue;
        }

        solvec.push_back ( data );
        string  filename = namefolder;
        datafile = "/information";
        string ext = ".txt";
        filename += datafile;
        auto s = std::to_string ( imc );
        filename += s;
        filename += ext;
        std::ofstream fileinfo ( filename );
        fileinfo << "Monte Carlo Sample = " << imc << std::endl;
        fileinfo << "Safety Factor = " << data << std::endl;
        fileinfo << "r/normr = " << solpost23[last][2] << std::endl;
        fileinfo << "u/normu = " << solpost23[last][3] << std::endl;
        fileinfo << "diff= " << solpost23[last][4] << std::endl;
        fileinfo << "diff2 = " << solpost23[last][5] << std::endl;
        fileinfo << "counterout = " << solpost23[last][6] << std::endl;


        if ( data <= 1. ) {
            fail++;
        }

        std::cout << " ******** MC  realization  = " << imc << " | Current safety fator = " << data << endl;
        solpost2[imc][0] = data;

    }

    file << "failue probability = " << Doub ( fail ) / Doub ( samples ) << endl;
    std::cout << "failue probability = " << Doub ( fail ) / Doub ( samples );
    auto end = sc.now();
    auto time_span = static_cast<chrono::duration<double>> ( end - start );
    std::cout << "\n Monte Carlo simualtion time  " << time_span.count()  << '\n';
    file << "\n Monte Carlo simualtion time  " << time_span.count()  << '\n';

}

void slopeproject::PrintMCS ( string namefolder,int imc,bool print )
{
    MatDoub hhatinho = AssembleHhationho ( imc );


    if ( print ) {



        string filename = namefolder;
        std::vector<std::vector<double>> hhatx;
        string name = "/Coesao";
        string ext = ".txt";
        filename += name;
        auto s = std::to_string ( imc );
        filename += s;
        filename += ext;
        fmesh->fmaterial->PostProcess ( fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), 0, hhatinho, hhatx );
        std::ofstream file ( filename );
        OutPutPost2 ( hhatx, file );


        filename = namefolder;
        std::vector<std::vector<double>> hhatx2;
        string namesss = "/Phi";
        string extsss = ".txt";
        filename += namesss;
        auto sss = std::to_string ( imc );
        filename += sss;
        filename += ext;
        fmesh->fmaterial->PostProcess ( fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), 1, hhatinho, hhatx2 );
        std::ofstream filesss ( filename );
        OutPutPost2 ( hhatx2, filesss );



        filename = namefolder;
        std::vector<std::vector<double>> solx, soly;
        fmesh->fmaterial->PostProcess ( fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), solx, soly );
        string name2 = "/soly";
        string ext2 = ".txt";
        filename += name2;
        auto s2 = std::to_string ( imc );
        filename += s2;
        filename += ext2;
        std::ofstream file2 ( filename );
        OutPutPost2 ( soly, file2 );


        filename = namefolder;
        name2 = "/solx";
        ext2 = ".txt";
        filename += name2;
        s2 = std::to_string ( imc );
        filename += s2;
        filename += ext2;
        std::ofstream file22 ( filename );
        OutPutPost2 ( solx, file22 );


        filename = namefolder;
        string name3 = "/hhatinho";
        string ext3 = ".txt";
        filename += name3;
        auto s3 = std::to_string ( imc );
        filename += s3;
        filename += ext3;
        std::ofstream file222 ( filename );
        OutPutPost ( hhatinho, file222 );


        filename = namefolder;
        std::vector<std::vector<double>> epsppost;
        fmesh->fmaterial->PostProcessIntegrationPointVar ( fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), epsppost );
        string name4 = "/plasticsqrtj2";
        string ext4 = ".txt";
        filename += name4;
        auto s4 = std::to_string ( imc );
        filename += s4;
        filename += ext4;
        std::ofstream file4 ( filename );
        OutPutPost2 ( epsppost, file4 );

    }

}



void slopeproject::OutPutFile1var ( MatDoub& postdata, std::ofstream& file )
{

    file.clear();

    for ( Int i = 0; i < postdata.nrows(); i++ ) {
        file << postdata[i][0] << endl;
    }

    file.close();
}

void slopeproject::OutPutFile ( MatDoub& postdata, std::ofstream& file )
{

    file.clear();

    for ( Int i = 0; i < postdata.nrows(); i++ ) {
        file << postdata[i][0] << " " << postdata[i][1] << endl;
    }

    file.close();
}
/*TITLE = "Example: Simple XY Plot"
VARIABLES = "X", "Y", "Z", "Distance", "Temperature"
ZONE T="Only Zone", I=20, J=2, K=1, F=POINT
0.0 0.0 0.0 0.0 850.9
0.0 0.0 0.0 0.022 790.1
0.0 0.0 0.0 0.051 554.0
0.0 0.0 0.0 0.069 540.0
0.0 0.0 0.0 0.075 545.5
0.0 0.0 0.0 0.083 549.4
0.0 0.0 0.0 0.106 590.2
0.0 0.0 0.0 0.1234 535.1
0.0 0.0 0.0 0.15 601.0
0.0 0.0 0.0 0.201 664.5
0.0 0.0 0.0 0.25 635.9
0.0 0.0 0.0 0.35 599.9
0.0 0.0 0.0 0.4454 600.0
0.0 0.0 0.0 0.49 513.0
0.0 0.0 0.0 0.5656 442.0
0.0 0.0 0.0 0.6423 333.5
0.0 0.0 0.0 0.6604 300.4
0.0 0.0 0.0 0.7 245.7
0.0 0.0 0.0 0.89 221.4
0.0 0.0 0.0 1.3 200.0*/
/*void slopeproject::OutPutPost(std::vector<std::vector<double>>& postdata, std::ofstream& file)
{
	file.clear();

    file << "TITLE = " <<   "Example: Simple XY Plot" << endl;
    file << "VARIABLES = " <<   "X ,"<< "Y ," << "Z ," <<" sqrt(J2)"<< endl;
    file << "ZONE T= " <<  "Only Zone ," << " I =" << postdata.size() << " J = 1,"<< " K = 2, "<< "F=POINT "<<  endl;
//ZONE T="Only Zone", I=20, J=2, K=1, F=POINT

	for (Int i = 0; i < postdata.size(); i++)
	{
		for (Int j = 0; j < postdata[0].size(); j++)
		{
			file << postdata[i][j] << " ";
		}
		file << 0.;
		file << endl;
	}
	file.close();
}*/

/*void slopeproject::OutPutPost2(std::vector<std::vector<double>>& postdata, std::ofstream& file)
{
	file.clear();
    file << "x coord, "<< "y coord, "<< "z coord, "<< "scalar "<< endl;
	for (Int i = 0; i < postdata.size(); i++)
	{
		for (Int j = 0; j < postdata[0].size(); j++)
		{
            if(j==2){
                file << 0. << ", " << postdata[i][j];
            }else{
                file << postdata[i][j] << ", ";
            }
		}
		file << endl;
	}
	file.close();
}*/

void slopeproject::OutPutPost2 ( std::vector<std::vector<double>>& postdata, std::ofstream& file )
{
    file.clear();
    for ( Int i = 0; i < postdata.size(); i++ ) {
        for ( Int j = 0; j < postdata[0].size(); j++ ) {
            file << postdata[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

void slopeproject::OutPutPost ( std::vector<std::vector<double>>& postdata, std::ofstream& file )
{
    file.clear();
    file << "data = { { 0, 0 } ";
    for ( Int i = 0; i < postdata.size(); i++ ) {
        file << "{";
        for ( Int j = 0; j < postdata[0].size(); j++ ) {
            if ( j!=postdata[0].size()-1 ) {
                file << postdata[i][j] << "  ,";
            } else {
                file << postdata[i][j] ;
            }
        }
        if ( i!=postdata.size()-1 ) {
            file << "},";
        } else {
            file << "}";
        }

    }
    file << "};";
    file << endl;
    file <<"ListLinePlot[data,PlotRange->All,PlotMarkers->Automatic]";
    file.close();
}


void slopeproject::OutPutPost ( MatDoub& postdata, std::ofstream& file )
{
    file.clear();
    for ( Int i = 0; i < postdata.nrows(); i++ ) {
        for ( Int j = 0; j < postdata.ncols(); j++ ) {
            file << postdata[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

void  slopeproject::OutPutPost ( MatInt& postdata, std::ofstream& file )
{
    file.clear();
    for ( Int i = 0; i < postdata.nrows(); i++ ) {
        for ( Int j = 0; j < postdata.ncols(); j++ ) {
            file << postdata[i][j] << " " ;
        }
        file << endl;
    }
    file.close();
}

void slopeproject::PrintMathematicaFormat ( MatDoub postdata, std::ofstream& file )
{
    file.clear();
    for ( Int i = 0; i < postdata.nrows(); i++ ) {
        for ( Int j = 0; j < postdata.ncols(); j++ ) {
            file << postdata[i][j] << " " ;
        }
        file << endl;
    }
    file.close();
}
