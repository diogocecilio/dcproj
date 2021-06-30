


#include "slopeproject.h"

#ifdef _WIN32    /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
#include <direct.h>
#endif

slopeproject::slopeproject()
{

}

slopeproject::slopeproject(mesh* inmesh, KLGalerkinRF* inklgalerking)
{
	fmesh = inmesh;
	fklgalerking = inklgalerking;
}
slopeproject::slopeproject(mesh* inmesh, KLGalerkinRF* inklgalerking, NRmatrix<MatDoub> randomfield)
{
	fmesh = inmesh;
	fklgalerking = inklgalerking;
	frandomfield = randomfield;
}
slopeproject::slopeproject(mesh inmesh, KLGalerkinRF inklgalerking, NRmatrix<MatDoub> randomfield)
{
	fmesh = &inmesh;
	fklgalerking = &inklgalerking;
	frandomfield = randomfield;
}
slopeproject::~slopeproject()
{
	delete fmesh;
	delete fklgalerking;
}

slopeproject::slopeproject(slopeproject& copy)
{

}

void slopeproject::CreateRandomField(string namefolder)
{

	mesh* localmesh = fmesh;
	KLGalerkinRF* objKLGalerkinRF = fklgalerking;

	int check;

	char* cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());
#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
	check = mkdir(cstr,777);
#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
	check = mkdir(cstr);
#endif

	string datafile = namefolder;
	datafile += "/datarandom.txt";
	std::ofstream file(datafile);
	file << " samples = " << objKLGalerkinRF->GetSamples() <<
		" | expansion order = " << objKLGalerkinRF->GetExpansionorder() << " | func type = " << objKLGalerkinRF->Getftype() << endl;
	file << "Lx = " << objKLGalerkinRF->GetLx() << " | Ly = " << objKLGalerkinRF->GetLy() << " variance = " << objKLGalerkinRF->GetSig() << endl;
	VecComplex  val; MatDoub  vec, HHAT;
	NRmatrix<MatDoub> randomfield;
	NRmatrix<Doub> cosionfield, frictionfield;

	std::vector<std::vector<double>> errpost;
	objKLGalerkinRF->SolveGenEigValProblem(val, vec, randomfield, errpost);
	cosionfield = randomfield[0][0];
	frictionfield = randomfield[1][0];
	datafile = namefolder;
	datafile += "/coesionfield.txt";
	std::ofstream coesfile(datafile);
	OutPutPost(cosionfield, coesfile);

	datafile = namefolder;
	datafile += "/frictionfield.txt";
	std::ofstream frictionfile(datafile);
	OutPutPost(frictionfield, frictionfile);

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

void slopeproject::findbcids(mesh* gmesh, std::vector<std::vector<int>>& idsvector)
{
	Int ndivs = 100000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);
	a[0] = 0.; a[1] = 0.;
	b[0] = 50.; b[1] = 0;
	Line(a, b, ndivs, pathbottom);
	
	std::vector<std::vector<std::vector<double> > > ggmesh = gmesh->GetAllCoords();
	
	NRmatrix<int> meshtop = gmesh->GetMeshTopology();
	
	FindIdsInPath(pathbottom, ggmesh,meshtop, idsbottom);

	idsvector.push_back(idsbottom);

	a[0] = 0.; a[1] = 0.;
	b[0] = 0.; b[1] = 20.;
	Line(a, b, ndivs, pathleft);
	FindIdsInPath(pathleft, ggmesh,meshtop, idsleft);

	idsvector.push_back(idsleft);
	a[0] = 50.; a[1] = 0.;
	b[0] = 50.; b[1] = 10;
	Line(a, b, ndivs, pathright);
	FindIdsInPath(pathright, ggmesh,meshtop, idsright);

	idsvector.push_back(idsright);

	a[0] = 19.99; a[1] = 19.99;
	b[0] = 20.; b[1] = 20.;
	Line(a, b, ndivs, pathdisplace);
	FindIdsInPath(pathdisplace,ggmesh,meshtop, iddisplace);

	idsvector.push_back(iddisplace);

}

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
		do
		{
            std::clock_t start,startout;
            double duration=0.;
            start = std::clock();
			FINT.assign(sz, 1, 0.);
			meshint->Assemble(KG, FINT, FBODY);
			R = FBODY;
			//R *= fac;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);

			SolveEigen(KG, R, sol);

			displace += sol;
			Doub u = fabs(displace[2 * iddisplace[0] + 1][0]);
			mat->UpdateDisplacement(displace);
			norm = R.NRmatrixNorm();
            
            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
            
			std::cout << " R norm = " << norm << " | phi0/phi = " << tan(phi0) / tan(phi) << " | c0/c = " << c0 / c << " | c = " << c << "| TIME = "<<duration <<" | phi = " << phi << std::endl;
			counter++;
			postcounter++;

            if (isnan(norm) == 1) {
				std::cout << "NAN" << endl;
				return solpost;
			}

		} while (norm > 0.01 && counter <15 && norm < 5000.);
        
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
		if (norm>= 1.) {
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
	std::cout << "FOS = " << FS << std::endl;
	
	if (false)
	{
        MatDoub solpost23;
		solpost23.CopyFromVector(solpost2);
		string names = "/home/diogo/projects/results/mathematicas/loadvsdisplacementSRM";
		string exts = ".dat";
		names += exts;
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);

		std::vector<std::vector<double>> epsppost;
		mat->PostProcessIntegrationPointVar(meshint->GetAllCoords(), meshint->GetMeshNodes(), meshint->GetMeshTopology(), mat->GetSolution(), epsppost);
		string name3 = "/home/diogo/projects/results/mathematicas/sqrtJ2SRM-Cho-Determ";
		string ext3 = ".dat";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		mat->PostProcess(meshint->GetAllCoords(), meshint->GetMeshNodes(), meshint->GetMeshTopology(), mat->GetSolution(), solx, soly);
		filename = "/home/diogo/projects/results/mathematicas/solySRM.dat";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "/home/diogo/projects/results/mathematicas/solxSRM.dat";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}
	return solpost;
}

void myTreads2(int a, int b, slopeproject* slopeobj2, string traedN)
{
	string namefolder3 = "D:/slope-results/THREADS/GI-cho-field-Lx20-Ly4" + traedN;
	slopeobj2->MonteCarloGIM(a, b, false, namefolder3);
	
}


std::vector<std::vector<double>>   slopeproject::IterativeProcess( int ndesi, Doub dlamb0,Doub alphatol, int niter)
{
	
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	std::vector<std::vector<int>> idsvector;

	findbcids(fmesh, idsvector);

	idsbottom = idsvector[0];
	idsleft = idsvector[1];
	idsright = idsvector[2];
	iddisplace = idsvector[3];

	Int sz = 2 * fmesh->GetMeshNodes().nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	//material* mat = fmesh->fmaterial;

	fmesh->fmaterial->ResetPlasticStrain();
	fmesh->fmaterial->ResetDisplacement();
	fmesh->fmaterial->ResetMat();

	MatDoub displace, displace0;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l = 0, l0 = 0, lamb = 1., lambn = 0, dlamb = 0., lamb3, diff = 100, diff2 = 100;
	Int counterout = 0, maxcountout = 30;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);


	MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
	//mesh->fmaterial->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
	fmesh->Assemble(KG, FINT, FBODY);
	R = FBODY;
	R *= lamb;
	R -= FINT;
	bool check = false;
	InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial);
	SolveEigen(KG, R, dws);
	SolveEigen(KG, FBODY, dwb);
	MatDoub dwbt, mult;
	dwb.Transpose(dwbt);
	dwbt.Mult(dwb, mult);
	l0 = sqrt(dlamb0 * dlamb0 * mult[0][0]);
	l = l0;
	dlamb = computelamda(dwb, dws, dw, l);
	lamb = 0;
	Doub scalar = 1.5;
	Doub rtol = 0.5;
	//cout << " \n SYSTEM SIZE  = " << KG.nrows() << std::endl;
	Doub rnorm = 10., rnormn=10.,lambn0=0.;
    Int counter = 0, maxcount = 20;
	do
	{
		std::clock_t start,startout;
		double duration;
        startout = std::clock();
		Doub err1 = 10., err2 = 10., tol = 1.e-5;
        if(counter<maxcount)lambn0 = lamb;
		diff = 10;
		displace0 = displace;
		rnorm = 10.;
        Doub meantime=0.;
        counter = 0;
		do
		{
            start = std::clock();
			//mat->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
			fmesh->Assemble(KG, FINT, FBODY);
			R = FBODY;
			R *= lamb;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial);
			//SolveNR3(KG, R, dws);
			//SolveNR3(KG, FBODY, dwb);
			SolveEigen(KG, R, dws);
			SolveEigen(KG, FBODY, dwb);

			//SolveEigen2(KG, R, dws);
			//SolveEigen2(KG, FBODY, dwb);
			
			dlamb = computelamda(dwb, dws, dw, l);
			if (isnan(dlamb) == 1) {
				std::cout << "NAN" << endl;
				std::cout << dlamb << endl;
				break;
			}
			lamb += dlamb;
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;

			displace += dww;
			fmesh->fmaterial->UpdateDisplacement(displace);

			rnorm = R.NRmatrixNorm();
			Doub normdw = dww.NRmatrixNorm();
			Doub unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;

			Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

			std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lamb  = " << lamb << " | TIME =  <<" << duration1 << std::endl;
            //std::cout << " | time =  <<" << duration1 << std::endl;
			counter++;
            meantime+=duration1;
			rtemp = rnorm;
			if (counter == 1)rnorm = 10;


		} while (counter < maxcount && rnorm > rtol);
        meantime/=counter;
        Doub outtime = (std::clock() - startout) / (double)CLOCKS_PER_SEC;
       // diff2= fabs(lambn0- lamb);
        std::cout << " exter iter = " << counterout << "  | newton iters = " << counter  << " |  |R| = " << rnorm << " |  lamb  = " <<lamb << " | dlamb  = " <<diff2 << "| mean time in "<< meantime <<" | total time (s) = " << outtime << std::endl;
        
		if (rnorm > 0.5)
		{

			//std::cout << "Convergence failed. \n";

			counterout++;
			dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

			fmesh->fmaterial->UpdateDisplacement(displace0);
			displace = displace0;
			//ndesi--;//74 s 0.497378 1.35437
				if (check==true)
				{
					dlamb0 *= 0.1/scalar;
					scalar *= 2.;
					//check = false;
					//if (rnorm> rnormn)check = false;
					// rnormn = rnorm;
				}
				else
				{
					dlamb0 *= 1/scalar;
					scalar *= 1.5;
					//check = true;
					//if (rnorm > rnormn)check = true;
					//rnormn = rnorm;
				}
			//MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
			dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

			//mat->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
			fmesh->Assemble(KG, FINT, FBODY);
			R = FBODY;
			R *= lamb;
			R -= FINT;
			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, fmesh->fmaterial);
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
			fmesh->fmaterial->UpdatePlasticStrain();
		}

		counterout++;
		dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = rnorm;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		uvf[1] = lamb;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);

		//std::cout << " $$$$$ Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm <<  " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	} while (counterout <= niter && fabs(diff2) > alphatol );// while (counterout <= maxcountout && fabs(diff2) > 0.05);


	if (true)
	{

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
		solpost23.CopyFromVector(solpost2);
		string names = "/home/diogo/projects/results/mathematicas/loadvsdisplacement";
		string exts = ".dat";
		names += exts;
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);


		std::vector<std::vector<double>> epsppost;
		fmesh->fmaterial->PostProcessIntegrationPointVar(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), epsppost);
		string name3 = "/home/diogo/projects/results/mathematicas/sqrtJ2GIM-deter";
		string ext3 = ".dat";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		fmesh->fmaterial->PostProcess(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), solx, soly);
		filename = "/home/diogo/projects/results/mathematicas/soly.dat";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "/home/diogo/projects/results/mathematicas/solx.dat";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

//	std::cout << "  Iteration number = " << counterout <<  " |  |R| = " << rnorm << " | lamb  = " << lamb << std::endl;
	return solpost;

}


std::vector<std::vector<double>>   slopeproject::IterativeProcessArcLengthSRM(int ndesi, Doub dlamb0, Doub alphatol, int niter)
{
	mesh* mesh = fmesh;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	std::vector<std::vector<int>> idsvector;

	findbcids(mesh, idsvector);

	idsbottom = idsvector[0];
	idsleft = idsvector[1];
	idsright = idsvector[2];
	iddisplace = idsvector[3];

	Int sz = 2 * mesh->GetMeshNodes().nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	material* mat = mesh->fmaterial;

	mat->ResetPlasticStrain();
	mat->ResetDisplacement();
	mat->ResetMat();

	MatDoub displace, displace0;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l = 0, l0 = 0, lamb = 1., lambn = 0, dlamb = 0., lamb3, diff = 100, diff2 = 100;
	Int counterout = 0, maxcountout = 20;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);



	NRvector<Doub> matconsts(4, 0.);
	mat->GetMatConstants(matconsts);
	Doub  young = matconsts[0];
	Doub nu = matconsts[1];
	Doub c = matconsts[2];
	Doub phi = matconsts[3];

	Doub c0 = c;
	Doub phi0 = phi;



	MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
	//mesh->fmaterial->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
	mesh->Assemble(KG, FINT, FBODY);
	
	R = FINT;
	R -= FBODY;
	bool check = false;
	InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);
	SolveEigen(KG, R, dws);
	SolveEigen(KG, FBODY, dwb);
	MatDoub dwbt, mult;
	dwb.Transpose(dwbt);
	dwbt.Mult(dwb, mult);
	l0 = sqrt(dlamb0 * dlamb0 * mult[0][0]);
	l = l0;
	dlamb = computelamda(dwb, dws, dw, l);
	lamb = 0;
	Doub scalar = 1.5;
	Doub rtol = 0.5;
	cout << " \n SYSTEM SIZE  = " << KG.nrows() << std::endl;
	cout << " \n dlamb  = " << dlamb << "\t l = "<< l << std::endl;
	Doub rnorm = 10., rnormn = 10.;
	do
	{
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
		do
		{
			//c = c0 / lamb;
			//phi = atan(tan(phi0) / lamb);
			//matconsts[2] = c;
			//matconsts[3] = phi;
			//mat->SetMatConstants(matconsts);

			mesh->Assemble(KG, FINT, FBODY);
			R = FINT;
			//R *= lamb;
			R -= FBODY;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);
			SolveEigen(KG, R, dws);
			SolveEigen(KG, FBODY, dwb);
			dlamb = computelamda(dwb, dws, dw, l);
			lamb += dlamb;
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;

			displace += dww;
			mat->UpdateDisplacement(displace);

			rnorm = R.NRmatrixNorm();
			Doub normdw = dww.NRmatrixNorm();
			Doub unorm = displace.NRmatrixNorm();
			err2 = normdw / unorm;

			Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			c = c0 / lamb;
			phi = atan(tan(phi0) / lamb);
			matconsts[2] = c;
			matconsts[3] = phi;
			mat->SetMatConstants(matconsts);
			std::cout << " Iteration number = " << counter << " | phi0/phi = " << tan(phi0) / tan(phi) << " | c0/c = " << c0 / c << " | c = " << c << " | phi = " << phi << " | lamb  = " << lamb << " | time =  <<" << duration1 << std::endl;
			counter++;

			rtemp = rnorm;
			if (counter == 1)rnorm = 10;


		} while (counter < maxcount && rnorm > rtol);


		c = c0 / lamb;
		phi = atan(tan(phi0) / lamb);
		matconsts[2] = c;
		matconsts[3] = phi;
		mat->SetMatConstants(matconsts);
		mat->UpdatePlasticStrain();
		counterout++;

		std::cout << " iter = " << counterout << " | FS= " << lamb << std::endl;

		dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = rnorm;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		uvf[1] = lamb;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);

		//std::cout << " $$$$$ Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm <<  " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	} while (counterout <= niter && fabs(diff2) > alphatol);// while (counterout <= maxcountout && fabs(diff2) > 0.05);


	if (true)
	{

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
		solpost23.CopyFromVector(solpost2);
		string names = "loadvsdisplacementnew";
		string exts = ".txt";
		names += exts;
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);


		std::vector<std::vector<double>> epsppost;
		mat->PostProcessIntegrationPointVar(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), mat->GetSolution(), epsppost);
		string name3 = "epsppostnew";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		mat->PostProcess(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), mat->GetSolution(), solx, soly);
		filename = "soly.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solx.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

	return solpost;

}





Doub  slopeproject::computelamda(MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l)
{
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
	Doub dlamb = (-bb + sqrt(delta)) / (2. * aa);
	return dlamb;


}

void slopeproject::OutPutFile(MatDoub& postdata, std::ofstream& file)
{

	file.clear();

	for (Int i = 0; i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << " " << postdata[i][1] << endl;
	}

	file.close();
}

void slopeproject::OutPutPost(std::vector<std::vector<double>>& postdata, std::ofstream& file)
{
	file.clear();
	for (Int i = 0; i < postdata.size(); i++)
	{
		for (Int j = 0; j < postdata[0].size(); j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}

void slopeproject::InserBC(MatDoub& KG, MatDoub& R, MatDoub& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material* mat)
{
	//FBODY *= 1. / lamb;
	Int dir, val;
	dir = 1;
	val = 0;
	mat->DirichletBC(KG, R, idsbottom, dir, val);

	//dir = 0;
	//val = 0;
	//mat->DirichletBC(KG, R, idsbottom, dir, val);


	dir = 0;
	val = 0;
	mat->DirichletBC(KG, R, idsright, dir, val);
	dir = 0;
	val = 0;
	mat->DirichletBC(KG, R, idsleft, dir, val);


	dir = 1;
	val = 0;
	mat->DirichletBC(KG, FBODY, idsbottom, dir, val);
	//dir = 0;
	//val = 0;
	//mat->DirichletBC(KG, FBODY, idsbottom, dir, val);

	dir = 0;
	val = 0;
	mat->DirichletBC(KG, FBODY, idsright, dir, val);
	dir = 0;
	val = 0;
	mat->DirichletBC(KG, FBODY, idsleft, dir, val);
}


void slopeproject::SolveEigen(MatDoub A, MatDoub b, MatDoub& x)
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
	VectorXd xxx = AA.llt().solve(bbb);
	//VectorXd xxx = AA.fullPivHouseholderQr().solve(bbb);
	//VectorXd xxx = AA.fullPivLu().solve(bbb);
	//VectorXd xxx = AA.ldlt().solve(bbb);
	//VectorXd xxx = AA.lu().solve(bbb);
	for (int i = 0; i < A.nrows(); i++)x[i][0] = xxx(i);
}

void slopeproject::OutPutPost(MatDoub& postdata, std::ofstream& file)
{
	file.clear();
	for (Int i = 0; i < postdata.nrows(); i++)
	{
		for (Int j = 0; j < postdata.ncols(); j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}

void  slopeproject::OutPutPost(MatInt& postdata, std::ofstream& file)
{
	file.clear();
	for (Int i = 0; i < postdata.nrows(); i++)
	{
		for (Int j = 0; j < postdata.ncols(); j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}

void  slopeproject::ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord)
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
			std::vector<Int> input_int = vecstr_to_vecint(tokens);
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
			std::vector<Doub> input_doub = vecstr_to_vecdoub(tokens);

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


	std::vector<Doub> temp33(2);
	for (Int i = 0; i < meshtopology.nrows(); i++)
	{
		std::vector< std::vector<Doub> > temp22;
		for (Int j = 0; j < meshtopology.ncols(); j++)
		{
			Int top = meshtopology[i][j];
			temp33[0] = meshcoords[top][0];
			temp33[1] = meshcoords[top][1];
			temp22.push_back(temp33);
		}
		allcoords.push_back(temp22);
	}



}
std::vector<Int>  slopeproject::vecstr_to_vecint(std::vector<string> vs)
{
	std::vector<Int> ret;
	for (std::vector<string>::iterator it = vs.begin() + 1; it != vs.end(); ++it)
	{
		istringstream iss(*it);
		Int temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

std::vector<Doub>  slopeproject::vecstr_to_vecdoub(std::vector<string> vs)
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
}

std::vector<Doub>  slopeproject::vecstr_to_vecdoub2(std::vector<string> vs)
{
	std::vector<Doub> ret;
	for (std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it)
	{
		istringstream iss(*it);
		Doub temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

template <class T>
std::vector<T>  slopeproject::vecstr_to_vec(std::vector<string> vs)
{
	std::vector<T> ret;
	for (std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it)
	{
		istringstream iss(*it);
		T temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

void  slopeproject::ReadMatDoub(MatDoub& matdoub, std::string  file)
{
	std::vector<std::vector<Doub>> coords;
	string line2, temp2;
	ifstream myfile2(file);
	if (myfile2.is_open())
	{
		while (getline(myfile2, line2))
		{
			std::vector<string> tokens;
			istringstream iss(line2);
			while (iss >> temp2)
				tokens.push_back(temp2);
			std::vector<Doub> input_doub = vecstr_to_vecdoub2(tokens);
			coords.push_back(input_doub);
		}
		myfile2.close();
	}
	else std::cout << "Unable to open file";
	//for (int i = 0;i < coords.size();i++)
	//{
	//	for (int j = 0;j < coords[0].size();j++)
	//	{
	//		cout << coords[i][j] << endl;
	//	}
	//	cout << endl;
	//}

	matdoub.CopyFromVector(coords);
}
void slopeproject::GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub& elcoords)
{

	elcoords.assign(allcoords[el].size(), 2, 0.);
	Int sz = allcoords[el].size();
	for (Int j = 0; j < sz; j++)
	{
		Doub x = allcoords[el][j][0];
		Doub y = allcoords[el][j][1];
		elcoords[j][0] = x;
		elcoords[j][1] = y;
	}
}

void slopeproject::FindIdsInPath(const MatDoub& path, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& idpath)
{
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

			for (Int ipath = 0; ipath < path.nrows(); ipath++)
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

void slopeproject::Line(VecDoub a, VecDoub b, Int ndivs, MatDoub& path)
{
	Doub x0 = a[0];
	Doub xf = b[0];

	Doub y0 = a[1];
	Doub yf = b[1];

	Doub dx = (xf - x0) / ndivs;
	Doub dy = (yf - y0) / ndivs;

	path.assign(ndivs, 2, 0.);

	for (Int idiv = 0; idiv < ndivs; idiv++)
	{
		path[idiv][0] = x0 + idiv * dx;
		path[idiv][1] = y0 + idiv * dy;
	}


}

MatDoub slopeproject::AssembleHhationho(Int i)
{

NRmatrix<MatDoub> randomfield = frandomfield;
Int nrandomvars = randomfield.nrows();
Int nrowss = randomfield[0][0].nrows();
VecDoub mean(nrandomvars, 0.), var(nrandomvars, 0.);
MatDoub hhatinho(randomfield[0][0].nrows(), randomfield.nrows(), 0.), posttest(randomfield[0][0].nrows(), 1, 0.);
MatDoub coesmat(randomfield[0][0].nrows(), 1, 0.), phimat(randomfield[0][0].nrows(), 1, 0.);

for (Int ivar = 0; ivar < nrandomvars; ivar++)
{
	for (Int isample = 0; isample < nrowss; isample++)
	{
		hhatinho[isample][ivar] = randomfield[ivar][0][isample][i];
		//hhatinho[i][0] = 0.;
		mean[ivar] += hhatinho[isample][ivar];
		mean[ivar] /= (isample + 1);
		var[ivar] += (mean[ivar] - hhatinho[isample][ivar]) * (mean[ivar] - hhatinho[isample][ivar]);
		posttest[isample][0] = hhatinho[isample][ivar];
	}
}
return hhatinho;
}

void slopeproject::MonteCarloSRM(int iter,int iter2, bool print, string writenamefolder) {

	//mesh* mesh2 = fmesh;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = fmesh->GetMeshTopology().nrows() * npts;
	Int sz = 2 * fmesh->GetMeshNodes().nrows();

	//material* mat = fmesh->fmaterial;


	NRvector<Doub> matconsts(4, 0.);
	fmesh->fmaterial->GetMatConstants(matconsts);
	Doub  young = matconsts[0];
	Doub nu = matconsts[1];
	Doub c = matconsts[2];
	Doub phi = matconsts[3];
	MatDoub bodyforce = fmesh->fmaterial->GetBodyForce();

	string namefolder = writenamefolder;
	char* cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());

#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
	int check = mkdir(cstr, 777);
#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
	int check = mkdir(cstr);
#endif


	std::vector<double> solvec;
	int samples = frandomfield[0][0].ncols();
	MatDoub solpost(samples, 2, 0.), solpost2(samples, 1, 0.);
	int fail = 0;


//	std::vector<int> copyv = { 8, 27, 34, 49, 63, 80, 93, 151, 203, 229, 231, 301, 339, 346, 380, 
//400, 522, 530, 542, 609, 655, 669, 682, 868, 872, 916, 939, 941, 990, 
//1053, 1072, 1083, 1115, 1146, 1309, 1379, 1409, 1446, 1458, 1567, 
//1568, 1638, 1655, 1787, 1829, 1927, 1949, 1955, 1959, 1967, 1983, 
//2008, 2077 };

//	std::vector<int> copyv = { 22, 32, 42, 70, 72, 162, 175, 207, 210, 330, 379, 421, 440 };


	//std::cout << " \n fail size = " <<copyv.size() << endl;
	int idfail = 0;
	for (int i = iter; i < iter2; i++)
	{
			
		
		
			//idfail = *find(copyv.begin(), copyv.end(), i);
			//if (idfail!=0 &&fabs(idfail)<2000)
			if (false)
			{
				std::cout << " \n fail = "  << endl;
				std::cout << " \n ID fail = " << idfail <<endl;
				idfail = 0;
				//continue;
			}
			else
			{



				//double min = *min_element(solvec.begin(), solvec.end());
				std::clock_t start;
				double duration;
				start = std::clock();

				MatDoub hhatinho = AssembleHhationho(i);
				fmesh->fmaterial->ResetMat();
				fmesh->SetHhat(hhatinho);
				fmesh->fmaterial->SetMemory(nglobalpts, sz);
				fmesh->fmaterial->SetMatConstants(matconsts);
				fmesh->fmaterial->UpdateBodyForce(bodyforce);
				Doub tol = 0.01;
				std::vector<std::vector<double>>  sol = IterativeProcessShearRed(0.1,2.,tol);

				MatDoub solpost23;
				solpost23.CopyFromVector(sol);


				Int last = solpost23.nrows() - 1;
				Doub data = solpost23[last][1];
				solvec.push_back(data);

				string  filename = namefolder;
				string datafile = "/information";
				string ext = ".txt";
				filename += datafile;
				auto s = std::to_string(i);
				filename += s;
				filename += ext;
				std::ofstream fileinfo(filename);
				fileinfo << "Monte Carlo Sample = " << i << std::endl;
				fileinfo << "Safety Factor = " << data << std::endl;
				fileinfo << "r/normr = " << solpost23[last][2] << std::endl;
				fileinfo << "u/normu = " << solpost23[last][3] << std::endl;
				fileinfo << "diff= " << solpost23[last][4] << std::endl;
				fileinfo << "diff2 = " << solpost23[last][5] << std::endl;
				fileinfo << "counterout = " << solpost23[last][6] << std::endl;
if (print) {
				filename = namefolder;
				std::vector<std::vector<double>> solx, soly;
				fmesh->fmaterial->PostProcess(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), solx, soly);
				string name2 = "/soly";
				string ext2 = ".txt";
				filename += name2;
				auto s2 = std::to_string(i);
				filename += s2;
				filename += ext2;
				std::ofstream file2(filename);
				OutPutPost(soly, file2);


				filename = namefolder;
				name2 = "/solx";
				ext2 = ".txt";
				filename += name2;
				filename += s2;
				filename += ext2;
				std::ofstream file22(filename);
				OutPutPost(solx, file22);

				filename = namefolder;
				std::vector<std::vector<double>> epsppost;
				fmesh->fmaterial->PostProcessIntegrationPointVar(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(),fmesh->fmaterial->GetSolution(), epsppost);
				string name3 = "/plasticsqrtj2";
				string ext3 = ".txt";
				filename += name3;
				filename += s2;
				filename += ext3;
				std::ofstream file3(filename);
				OutPutPost(epsppost, file3);

				string filename2 = namefolder;
				filename2 += "/montecarlosafetyfactor.txt";
				std::ofstream file23(filename2);
				OutPutFile1var(solpost2, file23);

				
					filename = namefolder;
					std::vector<std::vector<double>> hhatx;
					string name = "/Coesao";
					ext = ".txt";
					filename += name;
					s = std::to_string(i);
					filename += s;
					filename += ext;
					fmesh->fmaterial->PostProcess(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), 0, hhatinho, hhatx);
					std::ofstream file(filename);
					OutPutPost(hhatx, file);


					filename = namefolder;
					std::vector<std::vector<double>> hhatx2;
					string namesss = "/Phi";
					string extsss = ".txt";
					filename += namesss;
					auto sss = std::to_string(i);
					filename += sss;
					filename += ext;
					fmesh->fmaterial->PostProcess(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), 1, hhatinho, hhatx2);
					std::ofstream filesss(filename);
					OutPutPost(hhatx2, filesss);
				}


				if (data <= 1.)
				{
					fail++;
				}
				std::cout << " mc it = " << i << " | Current safety fator = " << data << endl;
				solpost2[i][0] = data;
				//delete mesh2;
				//delete mat;

				duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
				std::cout << "\n time in monte carlo iteration = " << duration << '\n';
				fileinfo << "CPU time = " << duration << std::endl;

			}
	}
	//delete mesh2;
    //delete objKLGalerkinRF;
   // delete mat;
}



void slopeproject::MonteCarloGIM(int iter, int iter2, bool print, string writenamefolder)
{

	//mesh* finemesh = fmesh;

	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = fmesh->GetMeshTopology().nrows() * npts;
	Int sz = 2 * fmesh->GetMeshNodes().nrows();


	//material* materialdp = finemesh->fmaterial;

	NRvector<Doub> matconsts(4, 0.);
	fmesh->fmaterial->GetMatConstants(matconsts);
	Doub  young = matconsts[0];
	Doub nu = matconsts[1];
	Doub c = matconsts[2];
	Doub phi = matconsts[3];
	MatDoub bodyforce = fmesh->fmaterial->GetBodyForce();
	std::clock_t start;
	double duration;
	start = std::clock();
	std::cout << "\n starting stochastic simulation " << endl;

	KLGalerkinRF* objKLGalerkinRF = fklgalerking;

	
	string namefolder = writenamefolder;

	char* cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());


#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
    int	check = mkdir(cstr, 777);
#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
	int check = mkdir(cstr);
#endif

	string datafile = namefolder;
	datafile += "/DATA.txt";
	std::ofstream file(datafile);
	file << " Young = " << young << " | nu = " << nu << endl;
	file << " c = " << c << " | phi = " << phi << endl;
	file << " bodyforce = " << bodyforce[1][0] << endl;
	file << " samples = " << fklgalerking->GetSamples() << " | expansion order = " << fklgalerking->GetExpansionorder() << " | func type = " << fklgalerking->Getftype() << endl;
	file << "Lx = " << fklgalerking->GetLx() << " | Ly = " << fklgalerking->GetLy() << " variance = " << fklgalerking->GetSig() << endl;

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n  simualtion time  " << duration << '\n';

	Int  fail = 0;
	std::cout << "\n starting Monte Carlo " << endl;
	start = std::clock();
	Int postprintfreq = 50;
	Doub sum = 0.;
	NRmatrix<MatDoub> randomfield = frandomfield;
	int samples = randomfield[0][0].ncols();
	MatDoub solpost(samples, 2, 0.), solpost2(samples, 1, 0.);
	Doub soldatamin = 10.;
	Doub soldatamax = -10;
	std::vector<double> solvec;
	std::cout << "\n samples " << samples << endl;
	for (Int imc = iter; imc < iter2; imc++)
	{
		std::cout << "\n MC realization " << imc << endl;
		MatDoub hhatinho = AssembleHhationho(imc);

		std::clock_t start1;
		double duration1;
		start1 = std::clock();

		NRvector<Doub> matconsts(4, 0.);
		matconsts[0] = young;
		matconsts[1] = nu;
		matconsts[2] = c;
		matconsts[3] = phi;

		fmesh->fmaterial->SetMatConstants(matconsts);
		//material->fYC.setup(young, nu, c, phi);
		fmesh->fmaterial->SetMemory(nglobalpts, sz);
		fmesh->fmaterial->UpdateBodyForce(bodyforce);
		fmesh->fmaterial->SetRandomField(hhatinho);

		//cout << "all cc" << finemesh->GetAllCoords()[0].size() << endl;

		fmesh->SetHhat(hhatinho);
		//std::vector<std::vector<double>>  sol = IterativeProcessSlope(finemesh, hhatinho, material);//x = desloc y = loadfactor
		//std::vector<std::vector<double>>  sol = IterativeProcess(finemesh, hhatinho, materialdp,10,1.);//x = desloc y = loadfactor
		int maxiter = 30;
		Doub deltatol = 0.02;
		int desirediter = 10;
		Doub lamb0 = 0.2;
		//int ndesi, Doub dlamb0, Doub alphatol, int niter
		std::vector<std::vector<double>>  sol = IterativeProcess(desirediter, lamb0, deltatol, maxiter);//x = desloc y = loadfactor
		//std::vector<std::vector<double>>  sol = IterativeProcessShearRed(allcoordsfine, meshcoordsfine, meshtopologyfine, hhatinho, material);//x = desloc y = loadfactor



		duration1 = (std::clock() - start1) / (double)CLOCKS_PER_SEC;
		std::cout << "IterativeProcess time " << duration1 << std::endl;



		start1 = std::clock();
		MatDoub solpost23;
		solpost23.CopyFromVector(sol);


		Int last = solpost23.nrows() - 1;
		Doub data = solpost23[last][1];

		if (data < 0.2)
		{
			continue;
		}


		solvec.push_back(data);

		string  filename = namefolder;
		datafile = "/information";
		string ext = ".txt";
		filename += datafile;
		auto s = std::to_string(imc);
		filename += s;
		filename += ext;
		std::ofstream fileinfo(filename);
		fileinfo << "Monte Carlo Sample = " << imc << std::endl;
		fileinfo << "Safety Factor = " << data << std::endl;
		fileinfo << "r/normr = " << solpost23[last][2] << std::endl;
		fileinfo << "u/normu = " << solpost23[last][3] << std::endl;
		fileinfo << "diff= " << solpost23[last][4] << std::endl;
		fileinfo << "diff2 = " << solpost23[last][5] << std::endl;
		fileinfo << "counterout = " << solpost23[last][6] << std::endl;
        fileinfo << "time = " << duration1 << std::endl;

		if (false) {


			if (print) {
				filename = namefolder;
				std::vector<std::vector<double>> hhatx;
				string name = "/Coesao";
				ext = ".txt";
				filename += name;
				s = std::to_string(imc);
				filename += s;
				filename += ext;
				fmesh->fmaterial->PostProcess(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), 0, hhatinho, hhatx);
				std::ofstream file(filename);
				OutPutPost(hhatx, file);


				filename = namefolder;
				std::vector<std::vector<double>> hhatx2;
				string namesss = "/Phi";
				string extsss = ".txt";
				filename += namesss;
				auto sss = std::to_string(imc);
				filename += sss;
				filename += ext;
				fmesh->fmaterial->PostProcess(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), 1, hhatinho, hhatx2);
				std::ofstream filesss(filename);
				OutPutPost(hhatx2, filesss);
			}
			//	filename = namefolder;
			//	string names = "/FxU";
			//	string exts = ".txt";
				//filename += names;
				//auto ss = std::to_string(imc);
				//filename += ss;
				////filename += exts;
				//std::ofstream file8(filename);
				//OutPutFile(uvf, file8);

			filename = namefolder;
			std::vector<std::vector<double>> solx, soly;
			fmesh->fmaterial->PostProcess(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), solx, soly);
			string name2 = "/soly";
			string ext2 = ".txt";
			filename += name2;
			auto s2 = std::to_string(imc);
			filename += s2;
			filename += ext2;
			std::ofstream file2(filename);
			OutPutPost(soly, file2);


			filename = namefolder;
			name2 = "/solx";
			ext2 = ".txt";
			filename += name2;
			s2 = std::to_string(imc);
			filename += s2;
			filename += ext2;
			std::ofstream file22(filename);
			OutPutPost(solx, file22);


			filename = namefolder;
			string name3 = "/hhatinho";
			string ext3 = ".txt";
			filename += name3;
			auto s3 = std::to_string(imc);
			filename += s3;
			filename += ext3;
			std::ofstream file222(filename);
			OutPutPost(hhatinho, file222);


			filename = namefolder;
			std::vector<std::vector<double>> epsppost;
			fmesh->fmaterial->PostProcessIntegrationPointVar(fmesh->GetAllCoords(), fmesh->GetMeshNodes(), fmesh->GetMeshTopology(), fmesh->fmaterial->GetSolution(), epsppost);
			string name4 = "/plasticsqrtj2";
			string ext4 = ".txt";
			filename += name4;
			auto s4 = std::to_string(imc);
			filename += s4;
			filename += ext4;
			std::ofstream file4(filename);
			OutPutPost(epsppost, file4);
		}

		string filename2 = namefolder;
		filename2 += "/montecarlosafetyfactor.txt";
		std::ofstream file23(filename2);
		OutPutFile1var(solpost2, file23);

		duration1 = (std::clock() - start1) / (double)CLOCKS_PER_SEC;

		std::cout << " Postprocess time " << duration1 << std::endl;

		if (data <= 1.)
		{
			fail++;
		}
		std::cout << " mc it = " << imc << " | Current safety fator = " << data << endl;
		solpost2[imc][0] = data;
		//delete mat;
		//delete mesh2;
	}
	string filename = namefolder;
	filename += "/montecarlosafetyfactor.txt";
	std::ofstream file23(filename);
	OutPutFile1var(solpost2, file23);
	file << "failue probability = " << Doub(fail) / Doub(samples) << endl;
	std::cout << "failue probability = " << Doub(fail) / Doub(samples);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n Monte Carlo simualtion time  " << duration << '\n';
	file << "\n Monte Carlo simualtion time  " << duration << '\n';

	//delete finemesh;
   // delete objKLGalerkinRF;
   // delete materialdp;
}





void slopeproject::OutPutFile1var(MatDoub& postdata, std::ofstream& file)
{

	file.clear();

	for (Int i = 0; i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << endl;
	}

	file.close();
}
