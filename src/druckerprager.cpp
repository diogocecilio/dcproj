#include "druckerprager.h"
#include <math.h>
#include "ludcmp.h"
#include "eigen_sym.h"
#include <math.h>
#include <cmath>
#include "cholesky.h"
class mins_ndim;
//extern bool globalfail;
druckerprager::druckerprager()
{
}

druckerprager::druckerprager(Doub young, Doub nu, Doub coesion, Doub frictionangle)
{
	setup(young, nu,coesion,frictionangle);
}


druckerprager::~druckerprager()
{

}

void druckerprager::closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep, Doub & projgamma)
{
	NRtensor<Doub>  epse = epst - epsp;
	NRmatrix<Doub>  C = GetElasticMatrix();
	//C.Print();
	//epse.Print();

	NRmatrix<Doub>  tempepsemat, stresstrial, Dept;
	epse.FromTensorToNRmatrix(tempepsemat);
	//std::cout << "tempepsemat (elastic strain MatDoub) = " << std::endl;
	//tempepsemat.Print();
	C.Mult(tempepsemat, stresstrial);
	//stresstrial.Print();
	NRtensor<Doub>  stresstrialtensor;
	epse.FromNRmatrixToTensor(stresstrial, stresstrialtensor);
	Doub I1, J2;
	J2 = stresstrialtensor.J2();
	I1 = stresstrialtensor.I1();
	Doub xi = I1 / sqrt(3.);
	Doub rho = sqrt(2. * J2);
	Doub yieldcr = yield(xi, rho);
	if (yieldcr < 0) {

		projstress = stresstrialtensor;
		projstrain = epse;
		Dept = C;
        Dep=Dept;
		projgamma = 0.;
	}
	else {
		Doub xitrial = I1 / sqrt(3.);
		NRmatrix<Doub>  pt, vec, vect;
		Doub rhotrial = sqrt(2. * J2);
		if (fabs(J2)<1.e-3 )
		{
			//cout << " fail true " << endl;
			//cout << " J2 " << J2 << endl;
			projstress = stresstrialtensor;
			projstrain = epse;
			Dept = C;
            Dep=Dept;
			projgamma = 0.;
			
		}
		else {

			stresstrialtensor.EigenSystem(pt, vec);
			Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
			Doub betasol = atan((sqrt(3.)*(-sig2 + sig3)) / (-2.* sig1 + sig2 + sig3));
			Doub xitrial = I1 / sqrt(3.);
			Doub xisol = FindMinimum(pt, xitrial, fflag);
			NRmatrix<Doub>  sig = HW(F1HWCylDruckerPragerSmoothPSMATCH(xisol, rhotrial, betasol));

			NRmatrix<Doub>  nvec;
			NRmatrix<Doub>  fulltensorproj = stressrecosntruction(sig, vec);
			NRmatrix<Doub>  Q = ComputeQ(fulltensorproj, tempepsemat, projstress, projstrain, projgamma, nvec);
			Doub checkdet = 1.e-8;
			//Doub detQ = fabs(Det(Q));
			Doub detQ =fabs( Q.DetEig());

			//if (detQ < checkdet )
			if (detQ < checkdet || stresstrialtensor.isfailed() || fflag ==false)
			{
				//cout << " fail true " << endl;
				//cout << " J2 " << J2 << endl;
				projstress = stresstrialtensor;
				projstrain = epse;
				Dept = C;
                Dep=Dept;
				projgamma = 0.;
				return;
			}
			else {
				//	cout << "detQ = " << detQ << endl;
				//	Q.Print();
				NRmatrix<Doub>  invQ, R;
				//Cholesky *chol = new Cholesky(Q);
				//chol->inverse(invQ);
				//if (chol->fail)
				//{
				//LUdcmp *lu = new LUdcmp(Q);
				//lu->inverse(invQ);
				Q.ComputeInverse(invQ);
				//}
				//if (lu->fail)
				//if (chol->fail)
				if(false)
				{
					cout << " LU fail true " << endl;
					projstress = stresstrialtensor;
					projstrain = epse;
					Dept = C;
                    Dep=Dept;
					projgamma = 0.;
				//	delete chol;
				//	delete lu;
					return;
				}
				else {

					//delete lu;

					//delete chol;

					invQ.Mult(C, R);
					Dept = R;
					NRmatrix<Doub>  asol, tempprod, tempprodT, temp2;
					//cout << "nvec = " << endl;
					//nvec.Print();
					R.Mult(nvec, tempprod);
					Doub sum = 0.;
					for (Int i = 0;i < nvec.nrows();i++)sum += nvec[i][0] * tempprod[i][0];
					//std::cout << "sum = " << sum  <<std::endl;
					tempprod.Transpose(tempprodT);
					tempprod.Mult(tempprodT, temp2);

					temp2 *= 1. / sum;
					Dept -= temp2;
                    Dep=Dept;
					fflag = true;
					//Dep.Print();
				}
			}
		}
	}
}

Doub druckerprager::FindMinimum(NRmatrix<Doub>  pt, Doub xitrial , bool flag)
{
	Doub xisol;
	Funcd func;
	func.setpt(pt, fyoung, fnu, fcoesion, fphi);
	Doub a = xitrial*10.e7;
	xisol = rtsafe(func, -a, a, 10.e-5,flag);
	return xisol;
}

NRmatrix<Doub>  druckerprager::ComputeQ(NRmatrix<Doub>  fulltensorproj, NRmatrix<Doub>  tempepsemat, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, Doub & projgamma, NRmatrix<Doub>  &nvec)
{
	NRmatrix<Doub>  C = GetElasticMatrix();

	NRmatrix<Doub>  projVoigtMat, invCe = GetInverseElasticMatrix(), epsemat;
	fulltensorproj.FromFullToVoigt(projVoigtMat);

	NRtensor<Doub>  voigtProjTensor, epsptensor;
	voigtProjTensor.FromNRmatrixToTensor(projVoigtMat, voigtProjTensor);
	projstress = voigtProjTensor;

	nvec = avec(voigtProjTensor);

	//cout << "nvec = " << endl;
	//nvec.Print();

	invCe.Mult(projVoigtMat, epsemat);
	NRmatrix<Doub>  epspmat = tempepsemat;
	epspmat -= epsemat;
	voigtProjTensor.FromNRmatrixToTensor(epspmat, epsptensor);

	voigtProjTensor.FromNRmatrixToTensor(epsemat, projstrain);
	Doub gamma = epsptensor.Norm() / nvec.NRmatrixNorm();
	projgamma = gamma;
	//std::cout << "gamma  = " << gamma  <<std::endl;

	MatDoub dnvecdsig = dadsig(voigtProjTensor);


	NRmatrix<Doub>  Q(6, 6, 0.);for (Int i = 0;i < 6;i++)Q[i][i] = 1.;
	NRmatrix<Doub>  sec;
	C.Mult(dnvecdsig, sec);
	sec *= gamma;
	Q += sec;

	return Q;
}

Doub druckerprager::yield(Doub xi,Doub rho)
{
	//return 1 + pow(rho, 2) / (2.*pow(fb, 2)) - pow((fapex - xi)/ sqrt(3), 2) / pow(fa, 2);
	if (xi > fapex) {
		xi = fapex;
	}
	return-(pow(((fapex - xi / sqrt(3.)) / fa), 2) - pow(((rho / sqrt(2.)) / fb), 2) - 1);
}


NRmatrix<Doub>  druckerprager::GetElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G = fG, K = fK;
	C[0][0] = (4 * G) / 3 + K;    C[0][1] = -((2 * G) / 3) + K;C[0][2] = -((2 * G) / 3) + K;C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -((2 * G) / 3) + K; C[1][1] = (4 * G) / 3 + K;C[1][2] = -((2 * G) / 3) + K;C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -((2 * G) / 3) + K; C[2][1] = -((2 * G) / 3) + K;C[2][2] = (4 * G) / 3 + K;C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                  C[3][1] = 0;                 C[3][2] = 0;                 C[3][3] = G; C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                  C[4][1] = 0;                 C[4][2] = 0;                 C[4][3] = 0.;C[4][4] = G; C[4][5] = 0.;
	C[5][0] = 0;                  C[5][1] = 0;                 C[5][2] = 0;                 C[5][3] = 0.;C[5][4] = 0.;C[5][5] = G;
	return C;
}
NRmatrix<Doub>  druckerprager::GetInverseElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G = fG, K = fK;
	C[0][0] = (G + 3 * K) / (9.*G*K);    C[0][1] = -1 / (6.*G) + 1 / (9.*K);C[0][2] = -1 / (6.*G) + 1 / (9.*K);C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -1 / (6.*G) + 1 / (9.*K);  C[1][1] = (G + 3 * K) / (9.*G*K);  C[1][2] = -1 / (6.*G) + 1 / (9.*K);C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -1 / (6.*G) + 1 / (9.*K);  C[2][1] = -1 / (6.*G) + 1 / (9.*K);C[2][2] = (G + 3 * K) / (9.*G*K);  C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                         C[3][1] = 0;                       C[3][2] = 0;                       C[3][3] = 1. / G;C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                         C[4][1] = 0;                       C[4][2] = 0;                       C[4][3] = 0.;C[4][4] = 1. / G;C[4][5] = 0.;
	C[5][0] = 0;                         C[5][1] = 0;                       C[5][2] = 0;                       C[5][3] = 0.;C[5][4] = 0.;C[5][5] = 1. / G;
	return C;
}
NRmatrix<Doub>  druckerprager::F1HWCylDruckerPragerSmoothPSMATCH(Doub xi, Doub rho, Doub beta)
{
	NRmatrix<Doub>  solproj(3, 1);
	Doub rhoint = sqrt(0.6666666666666666)*sqrt(-((pow(fb, 2)*(3 * pow(fa, 2) - 3 * pow(fapex, 2) + 2 * sqrt(3)*fapex*xi - pow(xi, 2))) / pow(fa, 2)));
	solproj[0][0] = xi;
	solproj[1][0] = rhoint;
	solproj[2][0] = beta;
	return solproj;
}
NRmatrix<Doub>  druckerprager::HW(NRmatrix<Doub>  sig)
{
	MatDoub sol(3, 1);
	Doub xi, rho, beta;
	xi = sig[0][0];
	rho = sig[1][0];
	beta = sig[2][0];
	sol[0][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta);
	sol[1][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta - 2. *  M_PI / 3.);
	sol[2][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta + 2. *  M_PI / 3.);
	return sol;
}
NRmatrix<Doub>  druckerprager::stressrecosntruction(NRmatrix<Doub>  val, NRmatrix<Doub>  vec)
{
	NRmatrix<Doub>  sol(3, 3, 0.);
	for (Int i = 0;i < 3;i++)
	{
		MatDoub colvec(3, 1, 0.), colvect, temp;
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
NRmatrix<Doub>  druckerprager::dadsig(NRtensor<Doub>  sigprojvoigt)
{
	//(P / b ^ 2) - (2 Outer[Times, Ii, Ii] / (9 a ^ 2))
	NRmatrix<Doub>   first = P();
	first *= 1/ (fb*fb);
	NRmatrix<Doub>  second, I(6,1,0.),It(1,6,0.);
	I[0][0] = 1.;I[1][0] = 1.;I[2][0] = 1.;
	I.Transpose(It);
	I.Mult(It, second);
	second *= 2. / (9 * fa*fa);
	first -= second;
	return first;
}

NRmatrix<Doub>  druckerprager::avec(NRtensor<Doub>  sigprojvoigt)
{
	//nvec=(6 apex b ^ 2 Ii - 2 b ^ 2 ComputeI1[sigma] Ii + 9 a ^ 2 ComputeS[sigma]) / (9 a ^ 2 b ^ 2) //. subst2
	NRtensor<Doub>  S;
	Doub I1 = sigprojvoigt.I1();
	NRmatrix<Doub>  first,second, third,I(6, 1, 0.), It(1, 6, 0.);
	I[0][0] = 1.;I[1][0] = 1.;I[2][0] = 1.;
	first = I, second = I;
	sigprojvoigt.ComputeS(S);
	S.FromTensorToNRmatrix(third);

	first *=6.*fapex*(fb*fb);
	second *= 2. * fb*fb*I1;
	third *= 9. * fa*fa;
	first -= second;
	first += third;
	first *= 1. / (9.*fa*fa*fb*fb);
	return first;
}

NRmatrix<Doub>  druckerprager::P()
{
	NRmatrix<Doub>  P(6, 6);
	P[0][0] = 2. / 3.;    P[0][1] = -1. / 3.; P[0][2] = -1. / 3.; P[0][3] = 0.;P[0][4] = 0.;P[0][5] = 0.;
	P[1][0] = -1. / 3.;   P[1][1] = 2. / 3.;  P[1][2] = -1. / 3.; P[1][3] = 0.;P[1][4] = 0.;P[1][5] = 0.;
	P[2][0] = -1. / 3.;   P[2][1] = -1. / 3.;  P[2][2] = 2. / 3.; P[2][3] = 0.;P[2][4] = 0.;P[2][5] = 0.;
	P[3][0] = 0.;      P[3][1] = 0.;     P[3][2] = 0.;     P[3][3] = 2.;P[3][4] = 0.;P[3][5] = 0.;
	P[4][0] = 0.;       P[4][1] = 0.;      P[4][2] = 0.;      P[4][3] = 0.;P[4][4] = 2.;P[4][5] = 0.;
	P[5][0] = 0.;       P[5][1] = 0.;      P[5][2] = 0.;      P[5][3] = 0.;P[5][4] = 0.;P[5][5] = 2.;

	return P;
}

void druckerprager::updateatributes(NRvector<MatDoub> mult)
{
	Doub newcoesion = 0.;
	Doub newphi = 0.;
	//Doub multcoes = 0., multphi = 0.;
	fyoung0 = fyoung;
	fnu0 = fnu;
	fcoesion0 = fcoesion;
	fphi0 = fphi;
	//Doub newyoung = fyoung + mult*fyoung;
	//Doub newcoesion = fcoesion + mult[0][0][0]*fcoesion;
	//Doub newphi = fphi + mult[1][0][0] * fphi;

	newcoesion = mult[0][0][0];
	newphi = mult[1][0][0];
	setup(fyoung, fnu, newcoesion, newphi);
}
