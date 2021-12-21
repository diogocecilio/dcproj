#include "KLGalerkinRF.h"
#include "shapequad.h"
#include <ctime>
#include <chrono>
#include <random>

KLGalerkinRF::KLGalerkinRF( Int order, Doub Lx, Doub Ly,  Int type, Int samples, Int expansionorder) {
//	fyoung = young;
//	fnu = nu;
//	fbodyforce = bodyforce;
//	fplanestress = planestress;
//	fthickness = thickness;
	fOrder = order;
	fLx = Lx;
	fLy = Ly;
	//fsig = sig;
	ftype = type;
	fsamples = samples;
	fexpansionorder = expansionorder;
	//cout << "sz = " << fmesh.fallcoords[0].size() << endl;
}


KLGalerkinRF::~KLGalerkinRF()
{
    /*	delete fyoung;
	delete fnu;
	delete fbodyforce;
	delete fplanestress;
	delete fthickness;
	delete fOrder;
	delete fmesh;
	delete fHHAT;
	delete fhhatvel;*/
    
}

void KLGalerkinRF::ContributeB(MatDoub &BE, Doub xi, Doub eta, Doub w, MatDoub elcoords)
{
	MatDoub psis, GradPsi, Jac, GradPhi, N, NT, psist, xycoords;

	int type = 1;
	shapequad objshapes(fOrder, type);

	objshapes.shapes(psis, GradPsi, xi, eta);
	psis.Transpose(psist);
	psist.Mult(elcoords, xycoords);
	GradPsi.Mult(elcoords, Jac);
	Int nnodes = psis.nrows();
	Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
	if (DetJ <= 0)
	{
		std::cout << "\n DetJ < 0 " << std::endl;
		return;
	}

	psis.Mult(psist, BE);
	BE *= DetJ*w;

}
void KLGalerkinRF::CacStiffB(MatDoub &BE, const MatDoub  &elcoords)
{
	MatDoub ptsweigths, BEt;
	Doub xi, eta, w;
	Int nnodes = elcoords.nrows();
	BE.assign(nnodes, nnodes, 0.);

	int type = 1;
	shapequad objshapes(fOrder, type);

	objshapes.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();

	for (Int ipt = 0;ipt < npts;ipt++)
	{
		xi = ptsweigths[ipt][0];
		eta = ptsweigths[ipt][1];
		w = ptsweigths[ipt][2];
		ContributeB(BEt, xi, eta, w, elcoords);
		BE += BEt;
	}


}
void KLGalerkinRF::AssembleB(MatDoub &B)
{
	std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
	MatDoub meshnodes = fmesh.GetMeshNodes();
	MatInt meshtopology = fmesh.GetMeshTopology();

	MatDoub BE, elcoords, eltopology;
	GetElCoords(allcoords, 0, elcoords);
	Int rows = elcoords.nrows();
	Int sz = meshnodes.nrows();
	Int cols = rows;
	B.assign(sz, sz, 0.);
	Int nels = allcoords.size();
	Int fu = 0;
	for (Int iel = 0;iel < nels;iel++)
	{
		GetElCoords(allcoords, iel, elcoords);
		CacStiffB(BE, elcoords);
		for (Int irow = 0;irow < rows;irow++)
		{
			Int rowglob = meshtopology[iel][irow];
			for (Int icol = 0;icol < cols;icol++)
			{
				Int colglob = meshtopology[iel][icol];
				B[rowglob][colglob] += BE[irow][icol];
			}
		}

	}
}

void KLGalerkinRF::ContributeC(MatDoub &CE, MatDoub psis1, MatDoub GradPsi1, MatDoub elcoords1, Doub w1, MatDoub psis2, MatDoub GradPsi2, MatDoub elcoords2, Doub w2)
{
	MatDoub psis1t, psis2t, xycoords1, xycoords2, Jac1, Jac2, temp;
	psis1.Transpose(psis1t);
	psis2.Transpose(psis2t);
	psis1t.Mult(elcoords1, xycoords1);
	GradPsi1.Mult(elcoords1, Jac1);
	psis2t.Mult(elcoords2, xycoords2);
	GradPsi2.Mult(elcoords2, Jac2);
	Doub DetJ1 = -Jac1[0][1] * Jac1[1][0] + Jac1[0][0] * Jac1[1][1];
	Doub DetJ2 = -Jac2[0][1] * Jac2[1][0] + Jac2[0][0] * Jac2[1][1];

	if (DetJ1 <= 0 || DetJ2 <= 0)
	{
		std::cout << "\n DetJ < 0 " << endl;
		return;
	}
	Doub CXX = AutocorrelationFunc(xycoords1, xycoords2);
	psis1.Mult(psis2t, CE);
	CE *= CXX*DetJ1*DetJ2*w1*w2;
}
void KLGalerkinRF::CacStiffC(MatDoub &CE, const MatDoub  &elcoords1, const MatDoub  &elcoords2)
{
	int type = 1;
	shapequad objshapes(fOrder, type);
	MatDoub intrule, psis1, psis2, GradPsi1, GradPsi2, CEt;
	Int nnodes = elcoords1.nrows();
	objshapes.pointsandweigths(intrule);
	Doub xi1, eta1, w1;
	Doub xi2, eta2, w2;
	CE.assign(nnodes, nnodes, 0.);

	Int npts = intrule.nrows();
	for (Int ipt = 0;ipt < npts;ipt++)
	{
		xi1 = intrule[ipt][0];
		eta1 = intrule[ipt][1];
		w1 = intrule[ipt][2];
		objshapes.shapes(psis1, GradPsi1, xi1, eta1);
		for (Int jpt = 0;jpt < npts;jpt++)
		{
			xi2 = intrule[jpt][0];
			eta2 = intrule[jpt][1];
			w2 = intrule[jpt][2];
			objshapes.shapes(psis2, GradPsi2, xi2, eta2);
			ContributeC(CEt, psis1, GradPsi1, elcoords1, w1, psis2, GradPsi2, elcoords2, w2);
			CE += CEt;
		}
	}
}
void KLGalerkinRF::AssembleC(MatDoub &C)
{
	std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
	MatDoub meshnodes = fmesh.GetMeshNodes();
	MatInt meshtopology = fmesh.GetMeshTopology();

	MatDoub elcoords, elcoords1, elcoords2, CE;
	Int rowglob, colglob;
	Int nels = allcoords.size();
	GetElCoords(allcoords, 0, elcoords);
	Int rows = elcoords.nrows();
	Int sz = meshnodes.nrows();
	Int cols = rows;
	C.assign(sz, sz, 0.);
	for (Int iel = 0;iel < nels;iel++)
	{
		for (Int jel = 0;jel < nels;jel++)
		{
			GetElCoords(allcoords, iel, elcoords1);
			GetElCoords(allcoords, jel, elcoords2);
			CacStiffC(CE, elcoords1, elcoords2);
			for (Int irow = 0;irow < rows;irow++)
			{
				rowglob = meshtopology[iel][irow];
				for (Int icol = 0;icol < cols;icol++)
				{
					colglob = meshtopology[jel][icol];
					C[rowglob][colglob] += CE[irow][icol];
				}
			}
		}
	}

}

//Doub KLGalerkinRF::AutocorrelationFunc(MatDoub  x1, MatDoub  x2)
//{
//	Doub val = 0, xx1, xx2, yy1, yy2, dist;
//
//	xx1 = x1[0][0];
//	yy1 = x1[0][1];
//
//	xx2 = x2[0][0];
//	yy2 = x2[0][1];
//
//	//dist = sqrt((xx2 - xx1)*(xx2 - xx1)*(yy2 - yy1)*(yy2 - yy1));
//	dist = sqrt(pow(xx1 - xx2, 2) + pow(yy1 - yy2, 2));
//
//	switch (ftype)
//	{
//	case 1:
//		val = fsig*fsig* exp(-fabs(xx1 - xx2) *fabs(xx1 - xx2) / (fLx*fLx) - fabs(yy1 - yy2) *fabs(yy1 - yy2) / (fLy*fLy));
//		break;
//	case 2:
//		val = fsig*fsig* exp(-((dist / fLx)*(dist / fLx)));
//	case 3:
//		val = fsig * fsig * exp(-fabs(xx1 - xx2)   / (fLx ) - fabs(yy1 - yy2)/ (fLy ));
//		break;
//	}
//
//	return val;
//}


Doub KLGalerkinRF::AutocorrelationFunc(MatDoub  x1, MatDoub  x2)
{
	Doub val = 0, xx1, xx2, yy1, yy2, dist;

	xx1 = x1[0][0];
	yy1 = x1[0][1];

	xx2 = x2[0][0];
	yy2 = x2[0][1];

	//dist = sqrt((xx2 - xx1)*(xx2 - xx1)*(yy2 - yy1)*(yy2 - yy1));
	dist = sqrt(pow(xx1 - xx2, 2) + pow(yy1 - yy2, 2));

	switch (ftype)
	{
	case 1:
		val =  exp(-fabs(xx1 - xx2) * fabs(xx1 - xx2) / (fLx * fLx) - fabs(yy1 - yy2) * fabs(yy1 - yy2) / (fLy * fLy));
		break;
	case 2:
		val = exp(-((dist / fLx) * (dist / fLx)));
	case 3:
		val = exp(-fabs(xx1 - xx2) / (fLx)-fabs(yy1 - yy2) / (fLy));
		break;
	}

	return val;
}

void KLGalerkinRF::SolveGenEigValProblem(VecComplex & val, MatDoub & vec, NRmatrix<MatDoub>  & HHAT, std::vector<std::vector<double>> &errpost)
{
	std::clock_t start;
	double duration;
	std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
	MatDoub meshnodes = fmesh.GetMeshNodes();
	MatInt meshtopology = fmesh.GetMeshTopology();

	//cout << "sz = " << allcoords[0].size() << endl;
	/* Your algorithm here */


	std::cout << "\n Assembling the correlation matrix (C) and the deformation matrix (B)" << endl;
	start = std::clock();
	MatDoub C, B, vect, internaleigenvectors;
	AssembleC(C);
	AssembleB(B);

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "\n time assembling  =  " << duration << '\n';

	std::cout << "\n Solving the generalized eigenvalue prolem" << endl;
	start = std::clock();
	MatDoub invB, ibvBC, PHIt, PHI;
	Cholesky* chol = new Cholesky(B);
	chol->inverse(invB);
	delete chol;
	invB.Mult(C, ibvBC);

	//Jacobi* Jaco = new Jacobi(ibvBC);
	//VecDoub val1 =Jaco->d;

	Unsymmeig* Hessenberg = new Unsymmeig(ibvBC);

	VecComplex internaleigenvalues;

	internaleigenvalues = Hessenberg->wri;


	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "\n time solving the eigenvalue problem  =  " << duration << '\n';

	//Verifica se a ordem de expansao é maior que o numero de autovalores disponiveis
	Int M = fexpansionorder;
	if (internaleigenvalues.size() < M)
	{
		M = internaleigenvalues.size();
	}

	//autovettores nas colunas
	internaleigenvectors = Hessenberg->zz;

	Int degreesfredom = internaleigenvectors.nrows();

	if (meshnodes.nrows() != degreesfredom)
	{
		throw("meshnodes.nrows() != nternaleigenvectors.nrows()");
	}

	std::cout << "\n Integrating the eigenfunctions over the domain to normilize them..." << endl;
	start = std::clock();
	std::vector<double> vecint;
	val.assign(M, 0.);
	vec.assign(degreesfredom, M, 0.);
	MatDoub vectocomputeintegral(degreesfredom, 1);
	for (Int iexp = 0;iexp < M;iexp++)
	{
		val[iexp] = internaleigenvalues[iexp];
		for (Int idegree = 0;idegree < degreesfredom;idegree++)
		{
			vec[idegree][iexp] = internaleigenvectors[idegree][iexp];
			vectocomputeintegral[idegree][0] = vec[idegree][iexp];
		}
		Doub integral2 = PerfomIntegralOfListconst(vectocomputeintegral);
		vecint.push_back(integral2);

	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "\n time integrating  =  " << duration << '\n';

	cout << " \n EIGENVALUES Hessenberg  AFTER CHOP" << endl;
	val.Print();

	std::ofstream fileval("/home/diogo/Dropbox/slope-reliability/results/mathematicas/eigenvalues-Lx40-Ly2.dat");

    fileval.clear();
    for(int  i =0;i<val.size();i++)
    {
        fileval  << i  << " "<< val[i].real() << std::endl;
    }
    fileval.close();


	for (Int j = 0;j < M;j++)
	{
		for (Int i = 0;i < vec.nrows();i++)
		{
			vec[i][j] *= 1. / vecint[j];
		}
	}

	MatDoub error;
	ComputeVarianceError(val, vec, error);
	//GenerateGaussinRandomField(val, vec, HHAT, errpost);
	GenerateNonGaussinRandomField(val, vec, HHAT, errpost);

}


void KLGalerkinRF::GenerateNonGaussinRandomField(VecComplex& val, MatDoub& vec, NRmatrix<MatDoub>& HHAT, std::vector<std::vector<double>>& errpost)
{
	//	std::vector<std::vector<double>> errpost;
//	PostProcess( error, errpost);
//	std::ofstream file("error.txt");
//	OutPutPost(errpost, file);
	Int M = fexpansionorder;
	MatDoub  PHIt, PHI, vect;
	vec.Transpose(vect);

	MatDoub I(M, M, 0.);
	for (Int i = 0; i < M; i++)I[i][i] = sqrt(fabs(val[i].real()));

	I.Mult(vect, PHI);

	PHI.Transpose(PHIt);


	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::default_random_engine generator(seed);



	std::normal_distribution<Doub> distribution(0., 1.);

	Int nsamples = fsamples;
	MatDoub THETA(M, nsamples, 0.), THETA2(M, nsamples, 0.);
	Doub correlation = -0.5;
	for (int n = 0; n < nsamples; n++)
	{
		for (Int iexp = 0; iexp < M; iexp++)
		{
            std::random_device rd{};
            std::mt19937 generator{ rd() };
			Doub xic = distribution(generator);
			Doub xiphi = distribution(generator);
			THETA[iexp][n] = xic;
			THETA2[iexp][n] = xic * correlation + xiphi * sqrt(1 - correlation * correlation);
		}
	}

	MatDoub hhatphi, hhatcoes;

	PHIt.Mult(THETA, hhatcoes);
	PHIt.Mult(THETA2, hhatphi);

	HHAT.resize(2, 1);
	//437 x 5000
	//em cada coluna da hhatcoes tem um random field
	cout << " n PRECISA MUDAR AQUIx!!" << endl;
	//Distribuição log-normal
	Doub mean = 23.;
	Doub sdev = 0.3 * mean;
	Doub xi = sqrt(log(1 + pow((sdev / mean) ,2) ));
	Doub lambda = log(mean) - xi * xi / 2.;
	for (int i = 0; i < hhatcoes.nrows();i++)for (int j = 0; j < hhatcoes.ncols(); j++)hhatcoes[i][j] = exp(lambda + xi * hhatcoes[i][j]);

    mean = 0.000000001 * M_PI/180.;
	sdev = 0. * mean;
	xi = sqrt(log(1 + pow((sdev / mean), 2)));
	lambda = log(mean) - xi * xi / 2.;
	for (int i = 0; i < hhatphi.nrows(); i++)for (int j = 0; j < hhatphi.ncols(); j++)hhatphi[i][j] = exp(lambda + xi * hhatphi[i][j]);

	//Distribuição normal
	//Doub meanphi = 20.*M_PI/180.;
	//Doub sdevphi = 0.3 * meanphi;
	//for (int i = 0; i < hhatphi.nrows(); i++)for (int j = 0; j < hhatphi.ncols(); j++)hhatphi[i][j] = meanphi+ sdevphi * hhatphi[i][j];

	HHAT[0][0] = hhatcoes;
	HHAT[1][0] = hhatphi;



	std::cout << "\n Exiting  generalized eigenvalue prolem" << endl;
}


void KLGalerkinRF::GenerateGaussinRandomField(VecComplex& val, MatDoub& vec, NRmatrix<MatDoub>& HHAT, std::vector<std::vector<double>>& errpost)
{
	//	std::vector<std::vector<double>> errpost;
//	PostProcess( error, errpost);
//	std::ofstream file("error.txt");
//	OutPutPost(errpost, file);
	Int M = fexpansionorder;
	MatDoub  PHIt, PHI, vect;
	vec.Transpose(vect);

	MatDoub I(M, M, 0.);
	for (Int i = 0; i < M; i++)I[i][i] = sqrt(fabs(val[i].real()));

	I.Mult(vect, PHI);

	PHI.Transpose(PHIt);


	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::default_random_engine generator(seed);

	std::random_device rd{};
	std::mt19937 generator{ rd() };

	std::normal_distribution<Doub> distribution(0., 1.);

	Int nsamples = fsamples;
	MatDoub THETA(M, nsamples, 0.), THETA2(M, nsamples, 0.);
	Doub correlation = -0.5;
	for (int n = 0; n < nsamples; n++)
	{
		for (Int iexp = 0; iexp < M; iexp++)
		{
			Doub xic = distribution(generator);
			Doub xiphi = distribution(generator);
			THETA[iexp][n] = xic;
			THETA2[iexp][n] = xic * correlation + xiphi * sqrt(1 - correlation * correlation);
		}
	}

	MatDoub hhatphi, hhatcoes;

	PHIt.Mult(THETA, hhatcoes);
	PHIt.Mult(THETA2, hhatphi);

	HHAT.resize(2, 1);
	//437 x 5000
	//em cada coluna da hhatcoes tem um random field

	Doub mean = 18.5633;
	Doub sdev = 0.3 * mean;
	for (int i = 0; i < hhatcoes.nrows(); i++)for (int j = 0; j < hhatcoes.ncols(); j++)hhatcoes[i][j] = mean + sdev * hhatcoes[i][j];
	Doub meanphi = 20. * M_PI / 180.;
	Doub sdevphi =0.3* meanphi;
	for (int i = 0; i < hhatphi.nrows(); i++)for (int j = 0; j < hhatphi.ncols(); j++)hhatphi[i][j] = meanphi + sdevphi * hhatphi[i][j];

	HHAT[0][0] = hhatcoes;
	HHAT[1][0] = hhatphi;



	std::cout << "\n Exiting  generalized eigenvalue prolem" << endl;
}


//void KLGalerkinRF::SolveGenEigValProblem(VecComplex& val, MatDoub& vec, NRmatrix<MatDoub>& HHAT, std::vector<std::vector<double>>& errpost)
//{
//	std::clock_t start;
//	double duration;
//	std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
//	MatDoub meshnodes = fmesh.GetMeshNodes();
//	MatInt meshtopology = fmesh.GetMeshTopology();
//
//	//cout << "sz = " << allcoords[0].size() << endl;
//	/* Your algorithm here */
//
//
//	std::cout << "\n Assembling the correlation matrix (C) and the deformation matrix (B)" << endl;
//	start = std::clock();
//	MatDoub C, B, vect, internaleigenvectors;
//	AssembleC(C);
//	AssembleB(B);
//
//	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//
//	std::cout << "\n time assembling  =  " << duration << '\n';
//
//	std::cout << "\n Solving the generalized eigenvalue prolem" << endl;
//	start = std::clock();
//	MatDoub invB, ibvBC, PHIt, PHI;
//	Cholesky* chol = new Cholesky(B);
//	chol->inverse(invB);
//	delete chol;
//	invB.Mult(C, ibvBC);
//
//	//Jacobi* Jaco = new Jacobi(ibvBC);
//	//VecDoub val1 =Jaco->d;
//
//	Unsymmeig* Hessenberg = new Unsymmeig(ibvBC);
//
//	VecComplex internaleigenvalues;
//
//	internaleigenvalues = Hessenberg->wri;
//
//
//	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//
//	std::cout << "\n time solving the eigenvalue problem  =  " << duration << '\n';
//
//	//Verifica se a ordem de expansao é maior que o numero de autovalores disponiveis
//	Int M = fexpansionorder;
//	if (internaleigenvalues.size() < M)
//	{
//		M = internaleigenvalues.size();
//	}
//
//
//	//autovettores nas colunas
//	internaleigenvectors = Hessenberg->zz;
//	Int degreesfredom = internaleigenvectors.nrows();
//	MatDoub normalizedvecs(internaleigenvectors);
//	vec.assign(degreesfredom, M, 0.);
//	val.assign(M, 0.);
//	for (Int icol = 0; icol < M; icol++)
//	{
//		Doub sqrsum = 0.;
//		for (Int irow = 0; irow < internaleigenvectors.nrows(); irow++)
//		{
//			sqrsum += internaleigenvectors[irow][icol] * internaleigenvectors[irow][icol];
//		}
//		Doub norm = sqrt(sqrsum);
//		for (Int irow = 0; irow < vec.nrows(); irow++)vec[irow][icol] = internaleigenvectors[irow][icol] / norm;
//		val[icol] = internaleigenvalues[icol];
//	}
//
//	//std::cout << "normilized vecs " << std::endl;
//	Doub norm = 0.;
//	for (Int irow = 0; irow < vec.nrows(); irow++) norm += vec[irow][0] * vec[irow][0];
//	cout << "norm vec = " << sqrt(norm) << std::endl;
//	//normalizedvecs.Print();
//
//
//
//	if (meshnodes.nrows() != degreesfredom)
//	{
//		throw("meshnodes.nrows() != nternaleigenvectors.nrows()");
//	}
//
//	std::cout << "\n Integrating the eigenfunctions over the domain to normilize them..." << endl;
//	start = std::clock();
//	std::vector<double> vecint;
//
//
//	/*MatDoub vectocomputeintegral(degreesfredom, 1);
//	for (Int iexp = 0;iexp < M;iexp++)
//	{
//		val[iexp] = internaleigenvalues[iexp];
//		for (Int idegree = 0;idegree < degreesfredom;idegree++)
//		{
//			vec[idegree][iexp] = internaleigenvectors[idegree][iexp];
//			vectocomputeintegral[idegree][0] = vec[idegree][iexp];
//		}
//		Doub integral2 = PerfomIntegralOfListconst(vectocomputeintegral);
//		vecint.push_back(integral2);
//
//	}*/
//
//
//
//	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//
//	std::cout << "\n time integrating  =  " << duration << '\n';
//
//	cout << " \n EIGENVALUES Hessenberg  AFTER CHOP" << endl;
//	val.Print();
//
//	//for (Int j = 0;j < M;j++)
//	//{
//	//	for (Int i = 0;i < vec.nrows();i++)
//	//	{
//	//		vec[i][j] *= 1. / vecint[j];
//	//	}
//	//}
//
//	//norm = 0.;
//	//for (Int irow = 0; irow < vec.nrows(); irow++) norm += vec[irow][0]* vec[irow][0];
//	//cout << "norm vec 2 = " << sqrt(norm) << std::endl;
//
//
//	MatDoub error;
//	ComputeVarianceError(val, vec, error);
//
//	//	std::vector<std::vector<double>> errpost;
//	PostProcess(error, errpost);
//	//	std::ofstream file("error.txt");
//	//	OutPutPost(errpost, file);
//
//	vec.Transpose(vect);
//
//	MatDoub I(M, M, 0.);
//	for (Int i = 0; i < M; i++)I[i][i] = sqrt(fabs(val[i].real()));
//
//	I.Mult(vect, PHI);
//
//	PHI.Transpose(PHIt);
//
//
//	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//	//std::default_random_engine generator(seed);
//
//	std::random_device rd{};
//	std::mt19937 generator{ rd() };
//
//	std::normal_distribution<Doub> distribution(0., 1.);
//
//	Int nsamples = fsamples;
//	MatDoub THETA(M, nsamples, 0.), THETA2(M, nsamples, 0.);
//	Doub correlation = -0.5;
//	for (int n = 0; n < nsamples; n++)
//	{
//		for (Int iexp = 0; iexp < M; iexp++)
//		{
//			Doub xic = distribution(generator);
//			Doub xiphi = distribution(generator);
//			THETA[iexp][n] = xic;
//			THETA2[iexp][n] = xic * correlation + xiphi * sqrt(1 - correlation * correlation);
//		}
//	}
//
//	MatDoub hhatphi, hhatcoes;
//
//	PHIt.Mult(THETA, hhatcoes);
//	PHIt.Mult(THETA2, hhatphi);
//
//	HHAT.resize(2, 1);
//	HHAT[0][0] = hhatcoes;
//	HHAT[1][0] = hhatphi;
//
//
//
//	std::cout << "\n Exiting  generalized eigenvalue prolem" << endl;
//
//}


Doub KLGalerkinRF::PerfomIntegralOfListconst2( const MatDoub &Vec)
{
	std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
	MatDoub meshnodes = fmesh.GetMeshNodes();
	MatInt meshtopology = fmesh.GetMeshTopology();

	MatDoub xycoords1, xycoords2, sol1, sol2;
	Doub delta = 0.01, sum = 0., solmean, dx, dy, solu, xi, eta;
	Int nels = allcoords.size(), iel;
	for (iel = 0;iel < nels;iel++)
	{
		for (xi = -1.; xi < 1. - delta; xi += delta)
		{
			for (eta = -1.;eta < 1. - delta;eta += delta)
			{
				SolPt( iel, Vec, xi, eta, xycoords1, sol1);
				SolPt(iel, Vec, xi + delta, eta + delta, xycoords2, sol2);
				solmean = (sol1[0][0] * sol1[0][0] + sol2[0][0] * sol2[0][0])*0.5;
				dx = xycoords1[0][0] - xycoords2[0][0];
				dy = xycoords1[0][1] - xycoords2[0][1];
				sum += (dx*dy*solmean);
			}

		}
	}
	solu = sqrt(sum);
	cout << "\n Inte = " << solu << endl;

	return solu;
}

//Compute the nodal error 
void KLGalerkinRF::ComputeVarianceError(VecComplex &val, MatDoub &vec, MatDoub &error)
{
	//err2 = (sig ^ 2 - (Sum[val[[i]] vec[[i]] ^ 2, { i, 1, Length[val] }]));
	int ndof = vec.nrows();
	int ivec = vec.ncols();
	MatDoub sig(ndof, 1, 0.);
	error.assign(ndof, 1,0.);
	for (Int i = 0;i < ivec; i++)
	{
		for (Int j = 0;j < ndof;j++)
		{
			error[j][0] += fabs(val[i].real())*pow(vec[j][i], 2);
		}
	}
//	error *= 1 / (fsig*fsig);
	//sig -= error;
	//error = sig;
	Doub err = 0.;
	for (Int i = 0; i < ivec; i++)
	{
		err += fabs(val[i].real()) ;
	}
	//std::cout << " err / (fsig * fsig) = " <<err / (fsig * fsig) << std::endl;
	std::cout << "mean error do 50x20 mesh = " <<1. - 1./(50.*20.-(30.+20.)*10./2.)*  err << std::endl;

}


Doub KLGalerkinRF::PerfomIntegralOfListconst( const MatDoub &Vec)
{

	std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
	MatDoub meshnodes = fmesh.GetMeshNodes();
	MatInt meshtopology = fmesh.GetMeshTopology();

	MatDoub xycoords1, xycoords2, sol1, sol2, intrule, psis, elcoords, GradPsi, Jac;
	Doub  sum = 0., solu, xi, eta, w;
	int type = 1;
	shapequad objshapes(fOrder, type);


	objshapes.pointsandweigths(intrule);
	Int nels = allcoords.size(), iel;
	for (iel = 0;iel < nels;iel++)
	{
		Int npts = intrule.nrows();
		for (Int ipt = 0;ipt < npts;ipt++)
		{
			xi = intrule[ipt][0];
			eta = intrule[ipt][1];
			w = intrule[ipt][2];
			objshapes.shapes(psis, GradPsi, xi, eta);
			GetElCoords(allcoords, iel, elcoords);
			GradPsi.Mult(elcoords, Jac);
			Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
			for (Int inode = 0;inode < psis.nrows();inode++)
			{
				sum += psis[inode][0] * pow(Vec[meshtopology[iel][inode]][0], 2)*DetJ*w;
			}
		}
	}
	solu = sqrt(sum);
	cout << "\n Inte = " << solu << endl;

	return solu;
}
