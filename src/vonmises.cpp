#include "vonmises.h"
#include "eigen_sym.h"
#include <math.h>
#include <cmath>
#include "cholesky.h"
#include "ludcmp.h"

vonmises::vonmises(Doub young, Doub nu, Doub sigy)
{
	fyoung = young;
	fsigy = sigy;
	fnu = nu;
	fG = young / ( 2.* (1. + nu) );
	fK = young / ( 3.* (1. - 2.* nu) );
}

vonmises::vonmises()
{
}


vonmises::~vonmises()
{
}

void vonmises::closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, MatDoub & Dep, Doub & projgamma)
{
	NRtensor<Doub>  epse = epst - epsp;
	MatDoub C = GetElasticMatrix();
	//C.Print();
	//epse.Print();

	MatDoub tempepsemat, stresstrial, Dept;
	epse.FromTensorToNRmatrix(tempepsemat);
	//std::cout << "tempepsemat (elastic strain MatDoub) = " << std::endl;
	//tempepsemat.Print();
	C.Mult(tempepsemat, stresstrial);
	//stresstrial.Print();
	NRtensor<Doub>  stresstrialtensor;
	epse.FromNRmatrixToTensor(stresstrial, stresstrialtensor);
	Doub I1, J2;
	J2= stresstrialtensor.J2();
	I1 = stresstrialtensor.I1();
	Doub yieldcr = yield(J2);
	if (yieldcr < 0){

		projstress = stresstrialtensor;
		projstrain = epse;
		Dept = C;
		Dep.assign(3, 3, 0.);
		Dep[0][0] = Dept[0][0];Dep[0][1] = Dept[0][1];Dep[0][2] = Dept[0][5];
		Dep[1][0] = Dept[1][0];Dep[1][1] = Dept[1][1];Dep[1][2] = Dept[1][5];
		Dep[2][0] = Dept[5][0];Dep[2][1] = Dept[5][1];Dep[2][2] = Dept[5][5];
		projgamma = 0.;
	}
	else {
		Doub xisol = I1 / sqrt(3.);
		Doub rhotrial = sqrt(2. * J2);
		//std::cout << "rhotrial = " << rhotrial << std::endl;
		MatDoub pt, vec,vect;
		stresstrialtensor.EigenSystem(pt, vec);
		//std::cout << "autovalores  \n";
		//pt.Print();
		//std::cout << "autovetores  \n";
		//vec.Print();
		Doub sig1=pt[0][0], sig2=pt[0][1], sig3=pt[0][2];
		Doub betasol = atan((sqrt(3.)*(-sig2 + sig3)) / (-2.* sig1 + sig2 + sig3));
		//std::cout << "betasol  = "<< betasol <<std::endl;

		MatDoub sig = HW(F1HWCylVonMises(xisol, rhotrial, betasol));

		//std::cout << "sig  = "  << std::endl;
		//sig.Print();

		MatDoub fulltensorproj = stressrecosntruction(sig, vec);
		//std::cout << "fullprojtensor  = "  << std::endl;
		//fulltensorproj.Print();
		MatDoub projVoigtMat, invCe = GetInverseElasticMatrix(),epsemat;
		fulltensorproj.FromFullToVoigt(projVoigtMat);
		//std::cout << "projVoigtMat  = " << std::endl;
		//projVoigtMat.Print();

		NRtensor<Doub>  voigtProjTensor,nvec, epsptensor;
		nvec.FromNRmatrixToTensor(projVoigtMat, voigtProjTensor);
		projstress = voigtProjTensor;

		voigtProjTensor.ComputeS(nvec);

		//std::cout << "ComputeS  = " << std::endl;
		//nvec.Print();

		nvec *= sqrt(3.) / (2.* sqrt(voigtProjTensor.J2()));

		//std::cout << "nvec  = " << std::endl;
		//nvec.Print();

		invCe.Mult(projVoigtMat, epsemat);
		MatDoub epspmat = tempepsemat;
		epspmat -= epsemat;
		nvec.FromNRmatrixToTensor(epspmat, epsptensor);

		nvec.FromNRmatrixToTensor(epsemat, projstrain);
		Doub gamma = epsptensor.Norm() / nvec.Norm();
		projgamma = gamma;
		//std::cout << "gamma  = " << gamma  <<std::endl;

		MatDoub dnvecdsig = dadsig(voigtProjTensor);

		//std::cout << "dnvecdsig  = " << std::endl;

		//dnvecdsig.Print();

		MatDoub Q(6, 6, 0.),invQ,R;for (Int i = 0;i < 6;i++)Q[i][i] = 1.;
		MatDoub sec;
		C.Mult(dnvecdsig, sec);
		sec *= gamma;
		Q += sec;
		//std::cout << "Q  = " << std::endl;
		//Q.Print();

		LUdcmp *lu = new LUdcmp(Q);
		lu->inverse(invQ);
		//std::cout << "invQ  = " << std::endl;
		//invQ.Print();

		invQ.Mult(C, R);

		//std::cout << "R  = " << std::endl;
		//R.Print();

		Dept = R;
		MatDoub asol, tempprod,tempprodT,temp2;
		nvec.FromTensorToNRmatrix(asol);
		R.Mult(asol, tempprod);
		Doub sum = 0.;
		for (Int i = 0;i < asol.nrows();i++)sum += asol[i][0] * tempprod[i][0];

		//std::cout << "sum = " << sum  <<std::endl;
		tempprod.Transpose(tempprodT);
		tempprod.Mult(tempprodT, temp2);
		//std::cout << "Outer[a.R,a.R] = " << std::endl;
		//temp2.Print();

		temp2*= 1./sum;
		Dept -= temp2;

		Dep.assign(3, 3, 0.);
		Dep[0][0] = Dept[0][0];Dep[0][1] = Dept[0][1];Dep[0][2] = Dept[0][5];
		Dep[1][0] = Dept[1][0];Dep[1][1] = Dept[1][1];Dep[1][2] = Dept[1][5];
		Dep[2][0] = Dept[5][0];Dep[2][1] = Dept[5][1];Dep[2][2] = Dept[5][5];

		/*Dep2D = { { Dep[[1, 1]], Dep[[1, 2]], Dep[[1, 6]] },{ Dep[[2, 1]],
			Dep[[2, 2]], Dep[[2, 6]] },{ Dep[[6, 1]], Dep[[6, 2]],
			Dep[[6, 6]] } };
		sigproj2D = { sigprojvoigth[[1]], sigprojvoigth[[2]],
			sigprojvoigth[[6]] };
*/
		//std::cout << "DEP = " << std::endl;
		//Dep.Print();

	}

}
MatDoub vonmises::dadsig(NRtensor<Doub>  sigprojvoigt)
{
	Doub I1 = sigprojvoigt.I1();
	Doub J2 = sigprojvoigt.J2();

	MatDoub  first = P();
	first *= sqrt(3.) / (2 * sqrt(J2));
	MatDoub second,Smat,SmatT;
	NRtensor<Doub>  S;
	sigprojvoigt.ComputeS(S);
	S.FromTensorToNRmatrix(Smat);
	Smat.Transpose(SmatT);
	Smat.Mult(SmatT, second);
	second *= sqrt(3.) / (4.*pow(J2, 3. / 2.) );
	first -= second;
	return first;

}

MatDoub vonmises::stressrecosntruction(MatDoub val, MatDoub vec)
{
	MatDoub sol(3,3,0.);
	for (Int i = 0;i < 3;i++)
	{
		MatDoub colvec(3, 1, 0.),colvect,temp;
		for (Int j = 0;j < 3;j++)
		{
			colvec[j][0] = vec[j][i];
		}

		colvec.Transpose(colvect);
		colvec.Mult(colvect, temp);
		temp *= val[i][0];
		sol += temp;
	}
	//sol.Print();
	return sol;
}

Doub vonmises::yield(Doub J2)
{
	return sqrt(3 * J2) - fsigy;
}

MatDoub vonmises::GetElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G=fG, K=fK;
	C[0][0] = (4 * G) / 3 + K;    C[0][1] = -((2 * G) / 3) + K;C[0][2] = -((2 * G) / 3) + K;C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -((2 * G) / 3) + K; C[1][1] = (4 * G) / 3 + K;C[1][2] = -((2 * G) / 3) + K;C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -((2 * G) / 3) + K; C[2][1] = -((2 * G) / 3) + K;C[2][2] = (4 * G) / 3 + K;C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                  C[3][1] = 0;                 C[3][2] = 0;                 C[3][3] = G; C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                  C[4][1] = 0;                 C[4][2] = 0;                 C[4][3] = 0.;C[4][4] = G; C[4][5] = 0.;
	C[5][0] = 0;                  C[5][1] = 0;                 C[5][2] = 0;                 C[5][3] = 0.;C[5][4] = 0.;C[5][5] = G;
	return C;
}

MatDoub vonmises::GetInverseElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G = fG, K = fK;
	C[0][0] = (G + 3 * K) / (9.*G*K);    C[0][1] = -1 / (6.*G) + 1 / (9.*K);C[0][2] = -1 / (6.*G) + 1 / (9.*K);C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -1 / (6.*G) + 1 / (9.*K);  C[1][1] = (G + 3 * K) / (9.*G*K);  C[1][2] = -1 / (6.*G) + 1 / (9.*K);C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -1 / (6.*G) + 1 / (9.*K);  C[2][1] = -1 / (6.*G) + 1 / (9.*K);C[2][2] = (G + 3 * K) / (9.*G*K);  C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                         C[3][1] = 0;                       C[3][2] = 0;                       C[3][3] =1./G;C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                         C[4][1] = 0;                       C[4][2] = 0;                       C[4][3] = 0.;C[4][4] = 1./G;C[4][5] = 0.;
	C[5][0] = 0;                         C[5][1] = 0;                       C[5][2] = 0;                       C[5][3] = 0.;C[5][4] = 0.;C[5][5] = 1./G;
	return C;
}

MatDoub vonmises::P()
{
	MatDoub P(6, 6);
	P[0][0] = 2./3.;    P[0][1] =  -1./3.; P[0][2] = -1./ 3.; P[0][3] = 0.;P[0][4] = 0.;P[0][5] = 0.;
	P[1][0] = -1./3.;   P[1][1] =  2./3.;  P[1][2] = -1./ 3.; P[1][3] = 0.;P[1][4] = 0.;P[1][5] = 0.;
	P[2][0] = -1./3.;   P[2][1] = -1./3.;  P[2][2] =  2./ 3.; P[2][3] = 0.;P[2][4] = 0.;P[2][5] = 0.;
	P[3][0] = 0. ;      P[3][1] =  0.;     P[3][2] = 0. ;     P[3][3] = 2.;P[3][4] = 0.;P[3][5] = 0.;
	P[4][0] = 0.;       P[4][1] = 0.;      P[4][2] = 0.;      P[4][3] = 0.;P[4][4] = 2.;P[4][5] = 0.;
	P[5][0] = 0.;       P[5][1] = 0.;      P[5][2] = 0.;      P[5][3] = 0.;P[5][4] = 0.;P[5][5] = 2.;

	return P;
}

MatDoub vonmises::F1HWCylVonMises(Doub xi, Doub rho, Doub beta)
{
	MatDoub solproj(3,1);
	solproj[0][0] = xi;
	solproj[1][0] = sqrt(2. / 3.)* fsigy;
	solproj[2][0] = beta;
	return solproj;
}

MatDoub vonmises::HW(MatDoub sig)
{
	MatDoub sol(3,1);
	Doub xi, rho, beta;
	xi = sig[0][0];
	rho = sig[1][0];
	beta = sig[2][0];
	sol[0][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta);
	sol[1][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta - 2. *  M_PI/3.);
	sol[2][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta + 2. *  M_PI/3.);
	return sol;
}
