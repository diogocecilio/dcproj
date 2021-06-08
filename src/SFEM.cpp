#include <iostream>
#include <Eigen/Core>
#include <direct.h>
#include <iomanip>
#include <cmath>
#include "shapequad.h"
#include "gridmesh.h"
#include "elastmat2D.h"
#include "KLGalerkinRF.h"
#include "vonmises.h"
#include "druckerprager.h"
#include "elastoplastic2D.h"
#include "ludcmp.h"
#include <math.h>
#include <cmath>
#include "mesh.h"
#include "postprocess.h"
#include <algorithm>
#include <chrono>  // for high_resolution_clock
#include<Eigen/SparseCholesky>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <math.h>
#include "elastoplasticbase.h"
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

using namespace std;
NRvector<MatDoub> FindEl(const std::vector<std::vector< std::vector<Doub > > > &allcoords, const MatInt meshtopology, const VecDoub Vec,  MatDoub hhat);
void OutPutFile(MatDoub & postdata, std::ofstream &file);
void OutPutFile1var(MatDoub & postdata, std::ofstream &file);
void OutPutFile4var(MatDoub & postdata, std::ofstream &file);
void OutPutPost(std::vector<std::vector<double>> & postdata, std::ofstream &file);
void OutPutPost(MatDoub & postdata, std::ofstream &file);
void OutPutPost(MatInt & postdata, std::ofstream &file);
void ReadMesh(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, string filenameel, string filenamecoord);
std::vector<Int> vecstr_to_vecint(std::vector<string> vs);
std::vector<Doub> vecstr_to_vecdoub(std::vector<string> vs);
std::vector<std::vector<int>>  LineTopology(std::vector<int> ids, Int order);
void ToMatInt(std::vector<std::vector<int>> in, MatInt & out);
template <class T>
std::vector<T> vecstr_to_vec(std::vector<string> vs);
void MonteCarloTalude2();
void MonteCarlo();
//std::vector<std::vector<double>>   IterativeProcess(mesh * gmesh, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material);
std::vector<std::vector<double>>   IterativeProcess(mesh* gmesh, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material, int ndesi, Doub dlamb0);
void ReadMatDoub(MatDoub & matdoub, std::string  file);
void SolveEigen(MatDoub A, MatDoub b, MatDoub& x);
void SolveNR3(MatDoub A, MatDoub b, MatDoub& x);
void SolveEigen2(MatDoub A, MatDoub b, MatDoub& x);
std::vector<std::vector<double>>   IterativeSlopeStabilityEig(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material);
std::vector<std::vector<int>>  LineTopology(std::vector<int> ids, Int order);
void ToMatInt(std::vector<std::vector<int>> in, MatInt& out);
bool flag = false;
void solvedeterm(string nodestr, string elsstr, MatDoub hhatinho);
void IterativeCho(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material);
void solvedetermtest(string nodestr, string elsstr, MatDoub hhatinho, Doub c, Doub phi, Doub gamma);
void InserBC(MatDoub& KG, MatDoub& R, MatDoub& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material* mat);
Doub computelamda(MatDoub &dwb, MatDoub &dws, MatDoub &dw, Doub &l);
std::vector<std::vector<double>>    IterativeProcessShearRed(mesh* gmesh, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material);
void   IterativeProcessShearRed2(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material);
std::vector<std::vector<double>>     IterativeProcessShearRed3(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material);

void SolveShearRed(string nodestr, string elsstr, MatDoub hhatinho,Doub c, Doub phi, Doub gamma);

std::vector<std::vector<double>>  IterativeProcessShearRed(mesh* mesh);

void findbcids(mesh* gmesh, std::vector<std::vector<int>> & idsvector);

std::vector<std::vector<double>>   IterativeProcess(mesh* gmesh, Int ndesi, Doub dlamb0);

void CreateRandomField();

void ExecuteMonteCarlo();

#include "slopeproject.h"

int main2()
{


	//CreateRandomField();
	Eigen::initParallel();
	Int n = Eigen::nbThreads();
	///string nodestr = "D:/DClib/meshes/nos-297.txt";
	///string elsstr = "D:/DClib/meshes/els-297.txt";
	string nodestr = "D:/DClib/meshes/nos-414.txt";
	string elsstr = "D:/DClib/meshes/els-414.txt";

	//string nodestr = "D:/DClib/meshes/nos-40vs80-413.txt";
	//string elsstr = "D:/DClib/meshes/els-40vs80-413.txt";

	//string nodestr = "D:/DClib/meshes/nos-513-2.txt";
	//string elsstr = "D:/DClib/meshes/els-513-2.txt";


//	string nodestr = "D:/DClib/meshes/nos-gian-155.txt";
//	string elsstr = "D:/DClib/meshes/els-gian-155.txt";

	//string nodestr = "D:/DClib/meshes/nos-870.txt";
	//string elsstr = "D:/DClib/meshes/els-870.txt";


	//string nodestr = "D:/DClib/meshes/nos-298-2.txt";
	//string elsstr = "D:/DClib/meshes/els-298-2.txt";


//string nodestr = "D:/DClib/meshes/nos-513-2.txt";
//	string elsstr = "D:/DClib/meshes/els-513-2.txt";


	//string nodestr = "D:/DClib/meshes/nos-237.txt";
//	string elsstr = "D:/DClib/meshes/els-237.txt";
//	string nodestr = "D:/DClib/meshes/nos-228.txt";
//	string elsstr = "D:/DClib/meshes/els-228.txt";
 	   //	string nodestr = "D:/DClib/meshes/nos-184-c.txt";
//	string elsstr = "D:/DClib/meshes/els-184-c.txt";
//string nodestr = "D:/DClib/meshes/nos-132-c3.txt";
//string elsstr = "D:/DClib/meshes/els-132-c3.txt";
//	string nodestr = "D:/DClib/meshes/nos-132-c2.txt";
//	string elsstr = "D:/DClib/meshes/els-132-c2.txt";
//string nodestr = "D:/DClib/meshes/nos-132-c.txt";
//string elsstr = "D:/DClib/meshes/els-132-c.txt";
// 
//solvedeterm( nodestr,  elsstr);
//Doub c = 8.66286;//FS = 0.7
//Doub c = 12.3755; //FS =1.0
//Doub c = 16.161;//FS = 1.3
//Doub c = 19.8008;//FS = 1.6
 
//Doub c = 16.25;//	1.31307
//Doub phi = 20 * M_PI / 180.;
//string nodestr = "D:/DClib/meshes/nos-150.txt";
//string elsstr = "D:/DClib/meshes/els-150.txt";

	MatDoub hhatinho;
	//string file = "D:/compile-results/hhatinho722.txt";
	//string file = "D:/compile-results/hhatinho100.txt";
	//string file = "D:/compile-results/hhatinho569.txt";
	//ReadMatDoub(hhatinho, file);
	//hhatinho.Print();

	//string nodestr = "D:/DClib/meshes/liu-nodes.txt";
	//string elsstr = "D:/DClib/meshes/liu-els.txt";
	//Doub c = 40,phi=32*M_PI/180.,gamma =-17.;//2

	// SolveShearRed(nodestr, elsstr, hhatinho, c, phi, gamma);//1.4199
	//solvedetermtest(nodestr, elsstr, hhatinho, c, phi, gamma);

	//MonteCarlo();//5.83
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);




	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	//Doub c = 5.62746, phi = 30 * M_PI / 180., gamma = -20.;//2
	//Doub c = 8.44119, phi = 30 * M_PI / 180., gamma = -20.;//2
	//Doub c = 11.2549, phi = 30 * M_PI / 180., gamma = -20.;//2

	//Doub c = 24.7511, phi = 20 * M_PI / 180., gamma = -20.;//2
	Doub c =18.5633, phi = 20 * M_PI / 180., gamma = -20.;//1.5
	//Doub c = 12.3755, phi = 20 * M_PI / 180., gamma = -20.;//1

	//Doub c = 34.07, phi = 0.000001 * M_PI / 180., gamma = -20.;//1
	//Doub c = 51.1073, phi = 0.000001 * M_PI / 180., gamma = -20.;//1.5
	//Doub c = 68.14, phi = 0.000001 * M_PI / 180., gamma = -20.;//2


	//Doub c = 22.5095, phi = 30. * M_PI / 180., gamma = -20.;//2
	//Doub c = 16.8822, phi = 30 * M_PI / 180., gamma = -20.;//1.5
	//Doub c = 11.2548, phi = 30 * M_PI / 180., gamma = -20.;//1

	//Doub c = 49.5022, phi = 20 * M_PI / 180., gamma = -20.;//2
	//Doub c =37.1266, phi = 20 * M_PI / 180., gamma = -20.;//1.5
	//Doub c = 24.7511, phi = 20 * M_PI / 180., gamma = -20.;//1

	//Doub c = 136.286, phi = 0.000001 * M_PI / 180., gamma = -20.;//2
	//Doub c =102.215, phi = 0.000001 * M_PI / 180., gamma = -20.;//1.5
	//Doub c = 68.1431, phi = 0.00001 * M_PI / 180., gamma = -20.;//1

	

	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Int planestress = 0;


//	Doub c = 24.7511, phi = 20 * M_PI / 180., gamma = -20.;//2
////Doub c = 44.7511, phi = 20 * M_PI / 180., gamma = -20.;//2
//	Doub thickness = 1.;
//	Doub young = 20000.;
//	Doub nu = 0.49;
//	Int planestress = 0;


	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = gamma;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows() * npts;
	Int sz = 2 * meshcoords.nrows();

	mesh* gmesh = new mesh(allcoords, meshcoords, meshtopology);
	elastoplastic2D< druckerprager >* mat = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	mesh* mesh2 = new mesh(mat,allcoords, meshcoords, meshtopology, hhatinho);
	mat->fYC.setup(young, nu, c, phi);
	mat->SetMemory(nglobalpts, sz);
	mat->UpdateBodyForce(bodyforce);


	Doub Lx = 20.;//(*Correlation length in x direction*)
	Doub Ly = 2.;//(*Correlation length in y direction*)

	Int nsamples = 5000, expansionorder = 60;
	Doub sig = 0.3;
	Int type = 1;
	KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF(*mesh2, order, Lx, Ly, sig, type, nsamples, expansionorder);




	int ndesirediters =8,  niter=50;
	Doub dlamb0=0.2, alphatol=0.0001;
	slopeproject slopeobj = slopeproject(mesh2, objKLGalerkinRF);
	Doub tol = 0.0001;
	std::vector<std::vector<double>> soll;
	//soll = slopeobj.IterativeProcessShearRed(0.5, 1., tol);
	mat->fYC.setup(young, nu, c, phi);
	mat->SetMemory(nglobalpts, sz);
	mat->UpdateBodyForce(bodyforce);
	//soll = slopeobj.IterativeProcess(ndesirediters,  dlamb0,  alphatol,  niter);

	//soll = slopeobj.IterativeProcessArcLengthSRM(ndesirediters, dlamb0, alphatol, niter);
	

	//std::vector<std::vector<double>> soll = slopeobj.IterativeProcessShearRed(0.5,1., tol);

	//string randomfieldfolder = "D:/DClib/randomfielddataLx20-Ly2-finemesh";
	//slopeobj.CreateRandomField(randomfieldfolder);


	MatDoub coesionrandomfield, frictionrandomfield;
	string filerf = "D:/DClib/randomfielddataLx20-Ly2-finemesh/coesionfield.txt";
	ReadMatDoub(coesionrandomfield, filerf);
	string filerff = "D:/DClib/randomfielddataLx20-Ly2-finemesh/frictionfield.txt";
	ReadMatDoub(frictionrandomfield, filerff);

	NRmatrix<MatDoub> randomfield(2,1);
	randomfield[0][0] = coesionrandomfield;
	randomfield[1][0] = frictionrandomfield;

	slopeproject slopeobj2 = slopeproject(mesh2, objKLGalerkinRF, randomfield);
	//std::vector<std::vector<double>> soltl = slopeobj2.IterativeProcess(10, 0.2,0.05,10);
	//std::vector<std::vector<double>> soltl = slopeobj2.IterativeProcessShearRed(0.5,1.);



	//string namefolder2 = "D:/DClib/SRM-Lx1000-Ly1000";

	//string namefolder3 = "D:/DClib/GIM-Lx1000-Ly1000";

	string namefolder2 = "D:/DClib/SRM-Lx20-Ly2-fine2";

	string namefolder3 = "D:/DClib/GIM-Lx20-Ly2-fine2";

	bool print = true;
	slopeobj2.MonteCarloGIM(405, 406, print, namefolder3);

	slopeobj2.MonteCarloSRM(405, 406, print, namefolder2);


	// slopeobj2.MonteCarloGIM(956, nsamples, namefolder3);

	//slopeobj2.MonteCarloSRM(956, nsamples,namefolder2);



	//int ndesirediters =8;
	//Doub dlamb0 = 0.2;
	//std::vector<std::vector<double>>  solu = IterativeProcess(gmesh, hhatinho, mat, ndesirediters, dlamb0);
	

	std::clock_t start;
	double duration;
	start = std::clock();

	std::vector<std::vector<double>>  solu = IterativeProcess(mesh2, ndesirediters, dlamb0);

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n time in deterministc GIM simulation = " << duration << '\n';

	//MonteCarlo();


	string namefolder = "D:/DClib/SRM";
	char* cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());
	int check = mkdir(cstr);


	MatDoub FSGI;
	string filemc = "D:/DClib/GI/montecarlosafetyfactor.txt";
	ReadMatDoub(FSGI, filemc);



	std::vector<double> solvec;
	int samples =FSGI.nrows();
	MatDoub solpost(samples, 2, 0.), solpost2(samples, 1, 0.);
	int fail = 0;
	for (int i = 1340; i < samples; i++)
	{

		if (FSGI[i][0] < 0.2)
		{
			std::cout << " \n GI FS too SMALL = " << FSGI[i][0] << std::endl;
		}
		else {

		std::cout << " \n GI FS = " << FSGI[i][0] << std::endl;
		std::clock_t start;
		double duration;
		start = std::clock();

		string file = "D:/DClib/GI/hhatinho";
		string ext3 = ".txt";
		auto s3 = std::to_string(i);
		file += s3;
		file += ext3;
		ReadMatDoub(hhatinho, file);

		mat->ResetMat();
		mat->fYC.setup(young, nu, c, phi);
		mat->SetMemory(nglobalpts, sz);
		mat->UpdateBodyForce(bodyforce);

		mesh2->SetHhat(hhatinho);
		std::vector<std::vector<double>>  sol = IterativeProcessShearRed(mesh2);

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
		
		filename = namefolder;
		std::vector<std::vector<double>> solx, soly;
		mat->PostProcess(allcoords, meshcoords, meshtopology, mat->GetSolution(), solx, soly);
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
		s2 = std::to_string(i);
		filename += s2;
		filename += ext2;
		std::ofstream file22(filename);
		OutPutPost(solx, file22);

		filename = namefolder;
		std::vector<std::vector<double>> epsppost;
		mat->PostProcessIntegrationPointVar(allcoords, meshcoords, meshtopology, mat->GetSolution(), epsppost);
		string name3 = "/plasticsqrtj2";
		 ext3 = ".txt";
		filename += name3;
		filename += s3;
		filename += ext3;
		std::ofstream file3(filename);
		OutPutPost(epsppost, file3);

		string filename2 = namefolder;
		filename2 += "/montecarlosafetyfactor.txt";
		std::ofstream file23(filename2);
		OutPutFile1var(solpost2, file23);


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
	//std::vector<std::vector<double>>  solu = IterativeProcessShearRed(mesh2);//x = desloc y = loadfactor
	


	//MonteCarlo();


	return 0;

}

void ExecuteMonteCarlo()
{

}

void CreateRandomField()
{

	string nodestr = "D:/DClib/meshes/nos-132-c3.txt";
	string elsstr = "D:/DClib/meshes/els-132-c3.txt";

	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);


	material *dummymat;

	mesh* localmesh = new mesh(dummymat, allcoords, meshcoords, meshtopology);

	Doub Lx = 20.;//(*Correlation length in x direction*)
	Doub Ly = 2.0;//(*Correlation length in y direction*)

	Int samples = 5000, expansionorder = 40;
	Doub sig = 0.3;
	Int type = 1,order=2;
	KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF(*localmesh, order, Lx, Ly, sig, type, samples, expansionorder);

	int check;
	string namefolder = "D:/DClib/randomfielddata";

	char* cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());

	check = mkdir(cstr);

	string datafile = namefolder;
	datafile += "/datarandom.txt";
	std::ofstream file(datafile);
	file << " samples = " << samples << " | expansion order = " << expansionorder << " | func type = " << type << endl;
	file << "Lx = " << Lx << " | Ly = " << Ly << " variance = " << sig << endl;
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

int cccccccc = 0;


void findbcids(mesh* gmesh, std::vector<std::vector<int>>& idsvector)
{
	//Int ndivs = 100000;
	//MatDoub pathbottom, pathleft, pathright, pathdisplace;
	//std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	//VecDoub a(2), b(2);


	//a[0] = 0.; a[1] = 0.;
	//b[0] = 36.; b[1] = 0;
	//gridmesh::Line(a, b, ndivs, pathbottom);
	//gridmesh::FindIdsInPath(pathbottom, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsbottom);

	//idsvector.push_back(idsbottom);


	//a[0] = 0.; a[1] = 0.;
	//b[0] = 0.; b[1] = 10.;
	//gridmesh::Line(a, b, ndivs, pathleft);
	//gridmesh::FindIdsInPath(pathleft, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsleft);


	//idsvector.push_back(idsleft);

	//a[0] = 36.; a[1] = 0.;
	//b[0] = 36.; b[1] = 18.;
	//gridmesh::Line(a, b, ndivs, pathright);
	//gridmesh::FindIdsInPath(pathright, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsright);

	//idsvector.push_back(idsright);

	//a[0] = 23.99; a[1] = 17.99;
	//b[0] = 24.; b[1] = 18.;
	//gridmesh::Line(a, b, ndivs, pathdisplace);
	//gridmesh::FindIdsInPath(pathdisplace, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), iddisplace);

	//idsvector.push_back(iddisplace);



	
Int ndivs = 100000;
MatDoub pathbottom, pathleft, pathright, pathdisplace;
std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
VecDoub a(2), b(2);
a[0] = 0.; a[1] = 0.;
b[0] = 50.; b[1] = 0;
gridmesh::Line(a, b, ndivs, pathbottom);
gridmesh::FindIdsInPath(pathbottom, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsbottom);

idsvector.push_back(idsbottom);

a[0] = 0.; a[1] = 0.;
b[0] = 0.; b[1] = 20.;
gridmesh::Line(a, b, ndivs, pathleft);
gridmesh::FindIdsInPath(pathleft, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsleft);

idsvector.push_back(idsleft);
a[0] = 50.; a[1] = 0.;
b[0] = 50.; b[1] = 10;
gridmesh::Line(a, b, ndivs, pathright);
gridmesh::FindIdsInPath(pathright, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsright);

idsvector.push_back(idsright);

a[0] = 19.99; a[1] = 19.99;
b[0] = 20.; b[1] = 20.;
gridmesh::Line(a, b, ndivs, pathdisplace);
gridmesh::FindIdsInPath(pathdisplace, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), iddisplace);

idsvector.push_back(iddisplace);

}


std::vector<std::vector<double>>  IterativeProcessShearRed(mesh *gmesh, elastoplastic2D< druckerprager >* material) {

	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	std::vector<std::vector<int>> idsvector;

	findbcids(gmesh, idsvector);

	idsbottom = idsvector[0];
	idsleft = idsvector[1];
	idsright = idsvector[2];
	iddisplace = idsvector[3];


	Int sz = 2 * gmesh->GetMeshNodes().nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	MatDoub displace, displace0,R,sol;
	displace.assign(sz, 1, 0.);
	displace0.assign(sz, 1, 0.);

	Doub l = 0, l0 = 0, lamb = 1., lambn = 0, dlamb = 0., lamb3, diff = 100, diff2 = 100;
	Int counterout = 0, maxcountout = 80;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);

	Int niter =500;
	Doub fac =0.5;
	Doub delta = 0.1;
	Int postcounter = 0;
	Doub c = material->fYC.GetCoes();
	Doub phi = material->fYC.GetPhi();
	Doub  young = material->fYC.GetYoung();
	Doub nu = material->fYC.GetNu();

	Doub c0 =c;
	Doub phi0 =phi;
	c = c0 / fac;
	phi = atan(tan(phi0 / fac));
	material->fYC.reset();
	material->fYC.setup(young, nu, c, phi);
	Doub res = 10, facn = 0,FS= fac,FSmin=0,FSmax=100000;
	Doub norm = 100000.;
	bool boll = false;
	do
	{
		Int counter = 0;
		//displace.assign(sz, 1, 0.);
		//material->UpdateDisplacement(displace);
		sol.assign(sz,1,0.);
		norm = 1000.;
		do
		{
			FINT.assign(sz, 1, 0.);
			material->Assemble(gmesh->GetAllCoords(), gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), KG, FINT, FBODY);
			R = FBODY;
			//R *= fac;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);

			SolveEigen2(KG, R, sol);

			displace += sol;
			Doub u = fabs(displace[2 * iddisplace[0] + 1][0]);
			material->UpdateDisplacement(displace);
			norm = R.NRmatrixNorm();
			std::cout << " R norm = " << norm <<" | phi0/phi = " << tan(phi0) / tan(phi) << " | c0/c = " << c0 / c << " | c = "<< c <<  " | phi = "<< phi<< std::endl;
			counter++;
			postcounter++;
			

		} while (norm > 0.01 && counter <30 && norm < 3000. );
	
		res = (fac - facn)/fac;
		if (R.NRmatrixNorm() > 1) {
			displace = displace0;
			//R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);
			//facn = FS;
			FSmax =  FS;
			FS =  (FSmin + FSmax) / 2.;
			//fac = FS;
			//material->ResetMat();
			material->UpdateDisplacement(displace0);
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
			}

		}



		c = c0 / FS;
		phi = atan(tan(phi0 ) / FS);
		material->fYC.reset();
		material->fYC.setup(young, nu, c, phi);
		
		counterout++;

		std::cout << " iter = " << counterout << " | FS= " << FS  << " | FSmax " << FSmax << " | FSmin = " << FSmin << std::endl;


		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = FS;
		solcount[2] = 0;
		solcount[3] = 0;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = 0;
		uvf[1] = FS;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);

	} while ((FSmax-FSmin) / FS > 0.001);
	std::cout << "FOS = " <<FS << std::endl;
	material->UpdatePlasticStrain();
	if (false)
	{

		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(gmesh->GetAllCoords(), gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), material->fdisplace, epsppost);
		string name3 = "epsppostnewHUM";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		material->PostProcess(gmesh->GetAllCoords(), gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), material->fdisplace, solx, soly);
		filename = "solyHUM.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solxHUM.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}
	return solpost;
}


std::vector<std::vector<double>>  IterativeProcessShearRed(mesh* mesh) {

	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	std::vector<std::vector<int>> idsvector;

	findbcids(mesh, idsvector);

	material *mat = mesh->fmaterial;

	idsbottom = idsvector[0];
	idsleft = idsvector[1];
	idsright = idsvector[2];
	iddisplace = idsvector[3];

	mat->ResetMat();
	Int sz = 2 * mesh->GetMeshNodes().nrows();
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
	Doub fac = 0.5;
	Doub delta = 1.;
	Int postcounter = 0;

	NRvector<Doub> matconsts(4, 0.);
	mat->GetMatConstants(matconsts);
	Doub  young = matconsts[0];
	Doub nu = matconsts[1];
	Doub c = matconsts[2];
	Doub phi = matconsts[3];

	Doub c0 = c;
	Doub phi0 = phi;
	c = c0 / fac;
	phi = atan(tan(phi0 / fac));
	matconsts[2]=c;
	matconsts[3]=phi;
	mat->SetMatConstants(matconsts);


	Doub res = 10, facn = 0, FS = fac, FSmin = 0, FSmax = 100000;
	Doub norm = 100000.;
	bool boll = false;
	do
	{
		Int counter = 0;
		//displace.assign(sz, 1, 0.);
		//material->UpdateDisplacement(displace);
		sol.assign(sz, 1, 0.);
		norm = 1000.;
		do
		{
			FINT.assign(sz, 1, 0.);
			mesh->Assemble(KG, FINT, FBODY);
			R = FBODY;
			//R *= fac;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);

			SolveEigen(KG, R, sol);

			displace += sol;
			Doub u = fabs(displace[2 * iddisplace[0] + 1][0]);
			mat->UpdateDisplacement(displace);
			norm = R.NRmatrixNorm();
			std::cout << " R norm = " << norm << " | phi0/phi = " << tan(phi0) / tan(phi) << " | c0/c = " << c0 / c << " | c = " << c << " | phi = " << phi << std::endl;
			counter++;
			postcounter++;


		} while (norm > 0.01 && counter < 10 && norm < 3000.);

		res = (fac - facn) / fac;
		if (R.NRmatrixNorm() > 1) {
			displace = displace0;
			//R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);
			//facn = FS;
			FSmax = FS;
			FS = (FSmin + FSmax) / 2.;
			//fac = FS;
			//material->ResetMat();
			mat->UpdateDisplacement(displace0);
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
			}

		}



		c = c0 / FS;
		phi = atan(tan(phi0) / FS);
		//YC->reset();
		//YC->setup(young, nu, c, phi);
		matconsts[2] = c;
		matconsts[3] = phi;
		mat->SetMatConstants(matconsts);


		counterout++;

		std::cout << " iter = " << counterout << " | FS= " << FS << " | FSmax " << FSmax << " | FSmin = " << FSmin << std::endl;


		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = FS;
		solcount[2] = 0;
		solcount[3] = 0;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = 0;
		uvf[1] = FS;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);

	} while ((FSmax - FSmin) / FS > 0.05);
	std::cout << "FOS = " << FS << std::endl;
	mat->UpdatePlasticStrain();
	if (false)
	{

		std::vector<std::vector<double>> epsppost;
		mat->PostProcessIntegrationPointVar(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), mat->GetSolution(), epsppost);
		string name3 = "epsppostnewHUM";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		mat->PostProcess(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), mat->GetSolution(), solx, soly);
		filename = "solyHUM.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solxHUM.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}
	return solpost;
}



std::vector<std::vector<double>>   IterativeProcess(mesh* mesh, int ndesi, Doub dlamb0)
{


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


	MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
	//mesh->fmaterial->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
	mesh->Assemble( KG, FINT, FBODY);
	R = FBODY;
	R *= lamb;
	R -= FINT;

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
	cout << " \n SYSTEM SIZE  = " << KG.nrows() << std::endl;
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
		Doub rnorm = 10.;
		do
		{

			//mat->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
			mesh->Assemble(KG, FINT, FBODY);
			R = FBODY;
			R *= lamb;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);
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
			mat->UpdateDisplacement(displace);

			rnorm = R.NRmatrixNorm();
			Doub normdw = dww.NRmatrixNorm();
			Doub unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;

			Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

			std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lamb  = " << lamb << " | time =  <<" << duration1 << std::endl;
			counter++;

			rtemp = rnorm;
			if (counter == 1)rnorm = 10;


		} while (counter < maxcount && rnorm > 0.5);


		if (rnorm > 10)
		{
			
			std::cout << "Convergence failed. \n";

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


			return solpost;
			mat->UpdateDisplacement(displace0);
			displace = displace0;
			//ndesi--;//74 s 0.497378 1.35437
			if (isnan(dlamb) == 1) {
				dlamb0 *= 0.8;
			}
			else {
				dlamb0 *= 0.8;
			}


			//MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
			dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

			//mat->Assemble(mesh->GetAllCoords(), mesh->GetMeshNodes(), mesh->GetMeshTopology(), KG, FINT, FBODY);
			mesh->Assemble(KG, FINT, FBODY);
			R = FBODY;
			R *= lamb;
			R -= FINT;
			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, mat);
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

			l *= Doub(ndesi) / Doub(counter);
			if (l > 2.)l = 2.;
			diff2 = fabs(lambn0 - lamb);
			mat->UpdatePlasticStrain();
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

		//std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	}while (counterout <= 4);// while (counterout <= maxcountout && fabs(diff2) > 0.05);

	cccccccc++;
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


std::vector<std::vector<double>>   IterativeProcess( mesh*gmesh, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material,int ndesi,Doub dlamb0)
{


	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	std::vector<std::vector<int>> idsvector;

	findbcids(gmesh, idsvector);

	idsbottom = idsvector[0];
	idsleft = idsvector[1];
	idsright = idsvector[2];
	iddisplace = idsvector[3];

	Int sz = 2 * gmesh->GetMeshNodes().nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	material->ResetPlasticStrain();
	material->ResetDisplacement();
	material->ResetMat();

	MatDoub displace,displace0;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l = 0, l0 = 0, lamb = 1., lambn = 0,dlamb=0., lamb3, diff = 100, diff2 = 100;
	Int counterout = 0,maxcountout =20;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);
	//int ndesi =10;//65 s u = 0.312275 lamb=1.35255
	//Doub dlamb0 = 1.1;
	//int ndesi = 5;//52 s u = 0.0969881 lamb = 1.34115
	//Doub dlamb0 = 1.1;

	//int ndesi = 5;//42 s u =  0.121154 lamb =1.34468
	//Doub dlamb0 = 1.5;

	//int ndesi = 3;//35 s u = u = 0.127228 lamb =1.34519
	//Doub dlamb0 = 2.;

	//int ndesi = 10;//60 s 0.21547 1.35037
	//Doub dlamb0 = 0.5;

	//int ndesi = 7;//63 s 0.0977655 1.34152
	//Doub dlamb0 = 0.5;
	// 
	
//	int ndesi = 8;//74 s 0.497378 1.35437  //FUNCIONA PARA MALHA 228 para cima e FS 1.3
//	Doub dlamb0 =	0.8;


	//	while (counterout < 15 && fabs(diff)>0.1)

	MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.),dw(sz, 1, 0.), R;
	material->Assemble(gmesh->GetAllCoords(),gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), KG, FINT, FBODY);
	R = FBODY;
	R *= lamb;
	R -= FINT;
	InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);
	SolveEigen(KG, R, dws);
	SolveEigen(KG, FBODY, dwb);
	MatDoub dwbt, mult;
	dwb.Transpose(dwbt);
	dwbt.Mult(dwb, mult);
	l0 = sqrt(dlamb0 * dlamb0 * mult[0][0]);
	l = l0;
	dlamb = computelamda(dwb, dws, dw, l);
	lamb = 0;
	do
	{
		std::clock_t start;
		double duration;
		start = std::clock();
		std::cout << "load step = " << counterout <<  " | l = " << l << " | diff2 = " << diff2 << std::endl;
		Int counter = 0, maxcount = 30;
		Doub err1 = 10., err2 = 10., tol = 1.e-5;
		Doub  lambn0 = lamb;
		diff = 10;
		displace0 = displace;
		Doub rnorm = 10.;
		do
		{
			
			material->Assemble(gmesh->GetAllCoords(), gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), KG, FINT, FBODY);

			R = FBODY;
			R *= lamb;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);

			SolveEigen(KG, R, dws);
			SolveEigen(KG, FBODY, dwb);



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
			material->UpdateDisplacement(displace);

			rnorm = R.NRmatrixNorm();
			Doub normdw = dww.NRmatrixNorm();
			Doub unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;

			Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

			std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm <<" | lamb  = " << lamb << " | time =  <<" << duration1 << std::endl;
			counter++;

			rtemp = rnorm;
			if (counter == 1)rnorm = 10;


		} while (counter < maxcount && rnorm > 0.1);


		if (rnorm > 1)
		{
			std::cout << "Convergence failed. \n";
			material->UpdateDisplacement(displace0);
			displace = displace0;
			//ndesi--;//74 s 0.497378 1.35437
			if (isnan(dlamb) == 1) {
				dlamb0 *= 0.5;
			}else{
				dlamb0 *= 0.9;
			}
			

			//MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
			dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

			material->Assemble(gmesh->GetAllCoords(), gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), KG, FINT, FBODY);
			R = FBODY;
			R *= lamb;
			R -= FINT;
			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);
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

			l *= Doub(ndesi) / Doub(counter);
			if (l > 2.)l = 2.;
			diff2 = fabs(lambn0 - lamb);
			material->UpdatePlasticStrain();
		}
		
		counterout++;
		dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.),dw.assign(sz,1,0.),R.assign(sz, 1, 0.),FBODY.assign(sz,1,0.),FINT.assign(sz,1,0.);

		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = err2;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		uvf[1] = lamb;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);
		
		//std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	} while (counterout <= maxcountout && fabs(diff2) > 0.001);
	string names = "fxUsn";
	auto sss = std::to_string(cccccccc);
	names += sss;
	MatDoub solpost23;
	solpost23.CopyFromVector(solpost2);
	string exts = ".txt";
	names += exts;
	std::ofstream file8(names);
	OutPutFile(solpost23, file8);
	cccccccc++;

	if (true)
	{


		MatDoub solpost23;
		solpost23.CopyFromVector(solpost2);
		string names = "loadvsdisplacementnew";
		string exts = ".txt";
		names += exts;
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);


		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(gmesh->GetAllCoords(), gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), material->fdisplace, epsppost);
		string name3 = "epsppostnew";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		material->PostProcess(gmesh->GetAllCoords(), gmesh->GetMeshNodes(), gmesh->GetMeshTopology(), material->fdisplace, solx, soly);
		filename = "soly.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solx.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

	return solpost;

}
Doub computelamda(MatDoub &dwb, MatDoub &dws, MatDoub &dw, Doub &l )
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
void InserBC(MatDoub &KG, MatDoub & R, MatDoub & FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material* mat)
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
}

void MonteCarlo()
{

//	string nodestrfine = "D:/DClib/src/nos-228.txt";
//	string elsstrfine = "D:/DClib/src/els-228.txt";

	//string nodestrfine = "D:/DClib/meshes/nos-297.txt";
	//string elsstrfine = "D:/DClib/meshes/els-297.txt";
//	string nodestrfine = "D:/DClib/meshes/nos-237.txt";
//	string elsstrfine = "D:/DClib/meshes/els-237.txt";

	string nodestrfine = "D:/DClib/meshes/nos-132-c3.txt";
	string elsstrfine = "D:/DClib/meshes/els-132-c3.txt";

//	string nodestrfine = "D:/DClib/meshes/nos-150.txt";
//	string elsstrfine = "D:/DClib/meshes/els-150.txt";

	MatDoub  meshcoordsfine, elcoordsfine;
	MatInt meshtopologyfine;
	std::vector<std::vector<std::vector<Doub>>> allcoordsfine, test;
	ReadMesh(allcoordsfine, meshcoordsfine, meshtopologyfine, elsstrfine, nodestrfine);


	mesh* finemesh = new  mesh(allcoordsfine, meshcoordsfine, meshtopologyfine);


	std::clock_t start;
	double duration;
	start = std::clock();

	std::cout << "\n starting stochastic simulation " << endl;


	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	//Doub c = 8.66286;//FS = 0.7
	//Doub c = 12.3755; //FS =1.0
	//Doub c = 16.161;//FS = 1.3
	//Doub c = 50;//FS = 1.3
	//Doub c = 19.8008;//FS = 1.6
   // Doub phi = 20 * M_PI / 180.;
	//Doub c = 16.0882, phi = 20 * M_PI / 180., gamma = -20.;//2

	Doub c = 24.7511, phi = 20 * M_PI / 180., gamma = -20.;//2




	//Doub young =100000.;
	//Doub nu = 0.3;
	//Doub c = 23;
	//Doub phi = 0.01 * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = -20.;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = finemesh->GetMeshTopology().nrows() * npts;
	Int sz = 2 * finemesh->GetMeshNodes().nrows();

	Doub Lx = 25.;//(*Correlation length in x direction*)
	Doub Ly = 2.5;//(*Correlation length in y direction*)

	Int samples = 2300, expansionorder = 40;
	Doub sig = 0.3;
	Int type = 1;
	KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF(*finemesh, order, Lx, Ly, sig, type, samples, expansionorder);

	//system("mkdir -p D:\DClib\results");

	//boost::filesystem::create_directories("D:\DClib\results");
	int check;
	//char* dirname = "D:\DClib\results";
	//auto s = std::to_string( rand() % 30 + 1985);
	//auto s = std::to_string(4);
	string namefolder = "D:/DClib/GIM-2300";
	//namefolder += s;

	char* cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());

	check = mkdir(cstr);

	string datafile = namefolder;
	datafile += "/DATA.txt";
	std::ofstream file(datafile);
	file << " Young = " << young << " | nu = " << nu << endl;
	file << " c = " << c << " | phi = " << phi << endl;
	file << " bodyforce = " << bodyforce[1][0] << endl;
	file << " Mesh  = " << nodestrfine << endl;
	file << " samples = " << samples << " | expansion order = " << expansionorder << " | func type = " << type << endl;
	file << "Lx = " << Lx << " | Ly = " << Ly << " variance = " << sig << endl;
	VecComplex  val; MatDoub  vec, HHAT;
	NRmatrix<MatDoub> randomfield;
	//objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, randomfield);


	std::vector<std::vector<double>> errpost;
	objKLGalerkinRF->SolveGenEigValProblem(val, vec, randomfield, errpost);
	datafile = namefolder;
	datafile += "/error2.txt";
	std::ofstream fileerro(datafile);
	OutPutPost(errpost, fileerro);

	delete objKLGalerkinRF;


	datafile = namefolder;
	datafile += "/vec.txt";
	std::ofstream filevec(datafile);
	OutPutPost(vec, filevec);

	datafile = namefolder;
	datafile += "/val.txt";
	std::ofstream fileval(datafile);
	for (Int j = 0; j < val.size(); j++)
	{
		fileval << val[j].real() << endl;
	}

	datafile = namefolder;
	datafile += "/meshcoords.txt";
	std::ofstream filemesh1(datafile);
	OutPutPost(finemesh->GetMeshNodes(), filemesh1);

	datafile = namefolder;
	datafile += "/meshtopology.txt";
	std::ofstream filemesh2(datafile);

	OutPutPost(finemesh->GetMeshTopology(), filemesh2);


	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n  simualtion time  " << duration << '\n';

	Int  fail = 0;
	std::cout << "\n starting Monte Carlo " << endl;
	start = std::clock();
	Int postprintfreq = 50;
	Doub sum = 0.;
	MatDoub solpost(samples, 2, 0.), solpost2(samples, 1, 0.);
	Doub soldatamin = 10.;
	Doub soldatamax = -10;
	std::vector<double> solvec;


	delete finemesh;
	//postprocess post;

	elastoplastic2D< druckerprager >* materialdp = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order);


	mesh* mesh2 = new mesh(materialdp, allcoordsfine, meshcoordsfine, meshtopologyfine);
	//material* mat = mesh2->fmaterial;

	for (Int imc = 0; imc < samples; imc++)
	{



		Int nrandomvars = randomfield.nrows();
		Int nrowss = randomfield[0][0].nrows();
		VecDoub mean(nrandomvars, 0.), var(nrandomvars, 0.);
		MatDoub hhatinho(randomfield[0][0].nrows(), randomfield.nrows(), 0.), posttest(randomfield[0][0].nrows(), 1, 0.);
		MatDoub coesmat(randomfield[0][0].nrows(), 1, 0.), phimat(randomfield[0][0].nrows(), 1, 0.);
		for (Int ivar = 0; ivar < nrandomvars; ivar++)
		{
			for (Int isample = 0; isample < nrowss; isample++)
			{
				hhatinho[isample][ivar] = randomfield[ivar][0][isample][imc];
				//hhatinho[i][0] = 0.;
				mean[ivar] += hhatinho[isample][ivar];
				mean[ivar] /= (isample + 1);
				var[ivar] += (mean[ivar] - hhatinho[isample][ivar]) * (mean[ivar] - hhatinho[isample][ivar]);
				posttest[isample][0] = hhatinho[isample][ivar];
			}
		}


		
		std::clock_t start1;
		double duration1;
		start1 = std::clock();



		NRvector<Doub> matconsts(4, 0.);
		matconsts[0] = young;
		matconsts[1] = nu;
		matconsts[2] = c;
		matconsts[3] = phi;

		materialdp->SetMatConstants(matconsts);
		//material->fYC.setup(young, nu, c, phi);
		materialdp->SetMemory(nglobalpts, sz);
		materialdp->UpdateBodyForce(bodyforce);
		materialdp->SetRandomField(hhatinho);

		//cout << "all cc" << finemesh->GetAllCoords()[0].size() << endl;

		mesh2->SetHhat(hhatinho);
		start1 = std::clock();
		//std::vector<std::vector<double>>  sol = IterativeProcessSlope(finemesh, hhatinho, material);//x = desloc y = loadfactor
		//std::vector<std::vector<double>>  sol = IterativeProcess(finemesh, hhatinho, materialdp,10,1.);//x = desloc y = loadfactor
		std::vector<std::vector<double>>  sol = IterativeProcess(mesh2,10,0.2);//x = desloc y = loadfactor
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

		if (true) {



			//filename = namefolder;
			//std::vector<std::vector<double>> hhatx;
			//string name = "/Coesao";
			//ext = ".txt";
			//filename += name;
			//s = std::to_string(imc);
			//filename += s;
			//filename += ext;
			//materialdp->PostProcess(allcoordsfine, meshcoordsfine, meshtopologyfine, 0, hhatinho, hhatx);
			//std::ofstream file(filename);
			//OutPutPost(hhatx, file);


			//filename = namefolder;
			//std::vector<std::vector<double>> hhatx2;
			//string namesss = "/Phi";
			//string extsss = ".txt";
			//filename += namesss;
			//auto sss = std::to_string(imc);
			//filename += sss;
			//filename += ext;
			//materialdp->PostProcess(allcoordsfine, meshcoordsfine, meshtopologyfine, 1, hhatinho, hhatx2);
			//std::ofstream filesss(filename);
			//OutPutPost(hhatx2, filesss);

			//	filename = namefolder;
			//	string names = "/FxU";
			//	string exts = ".txt";
				//filename += names;
				//auto ss = std::to_string(imc);
				//filename += ss;
				////filename += exts;
				//std::ofstream file8(filename);
				//OutPutFile(uvf, file8);

			//filename = namefolder;
			//std::vector<std::vector<double>> solx, soly;
			//materialdp->PostProcess(allcoordsfine, meshcoordsfine, meshtopologyfine, materialdp->GetSolution(), solx, soly);
			//string name2 = "/soly";
			//string ext2 = ".txt";
			//filename += name2;
			//auto s2 = std::to_string(imc);
			//filename += s2;
			//filename += ext2;
			//std::ofstream file2(filename);
			//OutPutPost(soly, file2);


			//filename = namefolder;
			//name2 = "/solx";
			//ext2 = ".txt";
			//filename += name2;
			//s2 = std::to_string(imc);
			//filename += s2;
			//filename += ext2;
			//std::ofstream file22(filename);
			//OutPutPost(solx, file22);


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
			materialdp->PostProcessIntegrationPointVar(allcoordsfine, meshcoordsfine, meshtopologyfine, materialdp->GetSolution(), epsppost);
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

}



void OutPutFile(MatDoub & postdata, std::ofstream &file)
{

	file.clear();

	for (Int i = 0;i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << " " << postdata[i][1] << endl;
	}

	file.close();
}

void OutPutFile1var(MatDoub & postdata, std::ofstream &file)
{

	file.clear();

	for (Int i = 0;i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << endl;
	}

	file.close();
}

void OutPutFile4var(MatDoub & postdata, std::ofstream &file)
{

	file.clear();

	for (Int i = 0;i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << " " << postdata[i][1] << " " << postdata[i][2] << " " << postdata[i][3] << endl;
	}

	file.close();
}


void OutPutPost(std::vector<std::vector<double>> & postdata, std::ofstream &file)
{
	file.clear();
	for (Int i = 0;i < postdata.size(); i++)
	{
		for (Int j = 0;j < postdata[0].size();j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}

void OutPutPost(MatDoub & postdata, std::ofstream &file)
{
	file.clear();
	for (Int i = 0;i < postdata.nrows(); i++)
	{
		for (Int j = 0;j < postdata.ncols();j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}
void OutPutPost(MatInt & postdata, std::ofstream &file)
{
	file.clear();
	for (Int i = 0;i < postdata.nrows(); i++)
	{
		for (Int j = 0;j < postdata.ncols();j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}

void ReadMesh(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, string filenameel, string filenamecoord)
{
	std::vector<std::vector<Int>> topol;
	string line,temp;

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
			for (int k = 0;k < input_int.size();k++)
			{
				input_int[k] = input_int[k]-1;
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
	for (Int i = 0; i < meshtopology.nrows();i++)
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
std::vector<Int> vecstr_to_vecint(std::vector<string> vs)
{
	std::vector<Int> ret;
	for (std::vector<string>::iterator it = vs.begin()+1;it != vs.end();++it)
	{
		istringstream iss(*it);
		Int temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

std::vector<Doub> vecstr_to_vecdoub(std::vector<string> vs)
{
	std::vector<Doub> ret;
	for (std::vector<string>::iterator it = vs.begin()+1;it != vs.end()-1;++it)
	{
		istringstream iss(*it);
		Doub temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

std::vector<Doub> vecstr_to_vecdoub2(std::vector<string> vs)
{
	std::vector<Doub> ret;
	for (std::vector<string>::iterator it = vs.begin() ;it != vs.end() ;++it)
	{
		istringstream iss(*it);
		Doub temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

template <class T>
std::vector<T> vecstr_to_vec(std::vector<string> vs)
{
	std::vector<T> ret;
	for (std::vector<string>::iterator it = vs.begin();it != vs.end();++it)
	{
		istringstream iss(*it);
		T temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

void ReadMatDoub(MatDoub & matdoub, std::string  file)
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




//// REMEMBER to update "lu.hpp" header includes from boost-CVS
//
//
///* Matrix inversion routine.
//Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
//bool InvertMatrix(const matrix<double>& input, matrix<double>& inverse) {
//	using namespace boost::numeric::ublas;
//	typedef permutation_matrix<std::size_t> pmatrix;
//	// create a working copy of the input
//	matrix<double> A(input);
//	// create a permutation matrix for the LU-factorization
//	pmatrix pm(A.size1());
//
//	// perform LU-factorization
//	int res = lu_factorize(A, pm);
//	if (res != 0) return false;
//
//	// create identity matrix of "inverse"
//	inverse.assign(identity_matrix<double>(A.size1()));
//
//	// backsubstitute to get the inverse
//	lu_substitute(A, pm, inverse);
//
//	return true;
//}
//void repruducecriticalcase()
//{
//	MatDoub hhatinho;
//	ReadMatDoub(hhatinho, "D:\\DClib\\results-75els-011996\\hhatinho1.145088.txt");
//	hhatinho.Print();
//
//	string nodestr = "nos-75.txt";
//	string elsstr = "els-75.txt";
//	MatDoub  meshcoords, elcoords;
//	MatInt meshtopology;
//	std::vector<std::vector<std::vector<Doub>>> allcoords;
//	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);
//
//	std::ofstream filemesh1("meshcoords.txt");
//	OutPutPost(meshcoords, filemesh1);
//
//	std::ofstream filemesh2("meshtopology.txt");
//
//	OutPutPost(meshtopology, filemesh2);
//
//	std::clock_t start;
//	double duration;
//	start = std::clock();
//
//	cout << "\n starting stochastic simulation " << endl;
//
//
//	Doub thickness = 1.;
//	Doub young = 20000.;
//	Doub nu = 0.49;
//	Doub c = 16.25;
//	Doub phi = 20 * M_PI / 180.;
//	Int planestress = 0;
//	MatDoub bodyforce(2, 1, 0.), newbodyforce;
//	bodyforce[1][0] = -20.;
//	MatDoub ptsweigths;
//	int order = 2;
//	shapequad shape = shapequad(order, 1);
//	shape.pointsandweigths(ptsweigths);
//	Int npts = ptsweigths.nrows();
//	Int nglobalpts = meshtopology.nrows()* npts;
//	Int sz = 2 * meshcoords.nrows();
//
//	elastoplastic2D< druckerprager > *  material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
//	material->fYC.setup(young, nu, c, phi);
//	material->SetMemory(nglobalpts, sz);
//	material->UpdateBodyForce(bodyforce);
//
//	std::vector<std::vector<double>>  sol = IterativeSlopeStabilityNew(allcoords, meshcoords, meshtopology, hhatinho, material);//x = desloc y = loadfactor
//
//	int check;
//	auto s = std::to_string(rand() % 30 + 1985);
//	string namefolder = "D:/DClib/criticalcases";
//	namefolder += s;
//
//	char *cstr = new char[namefolder.length() + 1];
//	strcpy(cstr, namefolder.c_str());
//
//	check = mkdir(cstr);
//
//	string datafile = namefolder;
//
//	string filename = namefolder;
//	std::vector<std::vector<double>> hhatx;
//	string name = "/Coesao";
//	string ext = ".txt";
//	filename += name;
//	filename += ext;
////	material->PostProcess(0, allcoords, meshtopology, hhatinho, hhatx);
//	std::ofstream file(filename);
//	OutPutPost(hhatx, file);
//
//
//	filename = namefolder;
//	std::vector<std::vector<double>> hhatx2;
//	string namesss = "/Phi";
//	string extsss = ".txt";
//	filename += namesss;
//	filename += ext;
////	material->PostProcess(1, allcoords, meshtopology, hhatinho, hhatx2);
//	std::ofstream filesss(filename);
//	OutPutPost(hhatx2, filesss);
//
//
//
//	filename = namefolder;
//	std::vector<std::vector<double>> solx, soly;
////	material->PostProcess(allcoords, meshtopology, material->fdisplace, solx, soly);
//	string name2 = "/soly";
//	string ext2 = ".txt";
//	filename += name2;
//	filename += ext2;
//	std::ofstream file2(filename);
//	OutPutPost(soly, file2);
//
//
//	filename = namefolder;
//	name2 = "/solx";
//	ext2 = ".txt";
//	filename += name2;
//	filename += ext2;
//	std::ofstream file22(filename);
//	OutPutPost(solx, file22);
//
//
//	filename = namefolder;
//	name2 = "/hhatinho";
//	ext2 = ".txt";
//	filename += name2;
//	filename += ext2;
//	std::ofstream file222(filename);
//	OutPutPost(hhatinho, file222);
//
//
//	filename = namefolder;
//	std::vector<std::vector<double>> epsppost;
////	material->PostProcessIntegrationPointVar(allcoords, meshtopology, material->fdisplace, epsppost);
//	string name3 = "/plasticsqrtj2";
//	string ext3 = ".txt";
//	filename += name3;
//	filename += ext3;
//	std::ofstream file3(filename);
//	OutPutPost(epsppost, file3);
//
//}

//MatDoub boostsolve(MatDoub KG, MatDoub FG)
//{
//	int sz = KG.nrows();
//	matrix<double> A;
//	boost::numeric::ublas::vector<double> b(sz);
//	KG.Toboost(A);
//	for (int i = 0;i < sz;i++)b(i) = FG[i][0];
//	//compressed_matrix<double, column_major, 0> A(sz,sz,sz*sz);
//	//A = m;
//	permutation_matrix<size_t> pm(A.size1());
//	lu_factorize(A, pm);
//	lu_substitute(A, pm, b);
//	MatDoub x(sz,1, 0.);
//	for (int i = 0;i < sz;i++)x[i][0] = b(i);
//	return x;
//}



NRvector<MatDoub> FindEl(const std::vector<std::vector< std::vector<Doub > > > &allcoords, const MatInt meshtopology, const VecDoub Vec, MatDoub hhat)
{
	Int el =0;
	Doub xi = 0.;
	Doub eta = 0.;
	MatDoub psis, GradPsi, elcoords, psist, solel,xycoords,sol;
	shapequad shape = shapequad(2, 1);
	shape.shapes(psis, GradPsi, xi, eta);
	Doub l0 = 100;
	Doub l = 100;
	std::vector<int> possible;
	int counter = 0;
	bool  breaktrue = false;
	while (counter<20)
	{
	//	if (breaktrue == true)
		//{
		//	possible.resize(0);
	//	}
	
	for (int iel = 0;iel < allcoords.size();iel++)
	{
		gridmesh::GetElCoords(allcoords, iel, elcoords);
		Int nodes = psis.nrows();
		Int nstatevars = 1;
		solel.assign(nodes, 1, 0);
		psis.Transpose(psist);
		psist.Mult(elcoords, xycoords);
		Doub dx  = xycoords[0][0]-Vec[0];
		Doub dy = xycoords[0][1]-Vec[1];
		Doub dist = sqrt(dx*dx + dy*dy);

		if (dist<l0)
		{
			possible.push_back(iel);
			//breaktrue = true;
		}
		else {
			//breaktrue = false;
			//break;
		}
	}
	cout << "l0 = " << l0 << endl;
	l0 *= 0.5;
	counter++;
	}

	el = possible.size() -1;
	cout << possible[el] << std::endl;
	gridmesh::GetElCoords(allcoords, possible[el], elcoords);
	elcoords.Print();
	Int rows = elcoords.nrows();
	NRvector<MatDoub> fhhatvel;
	fhhatvel.resize(hhat.ncols());
	for (Int ivar = 0;ivar < hhat.ncols();ivar++) {
		fhhatvel[ivar].assign(rows, 1, 0.);
		for (Int inode = 0;inode < rows;inode++)
		{
			fhhatvel[ivar][inode][0] = hhat[meshtopology[possible[el]][inode]][ivar];
		}
	}


	Doub tol = 10e-3;
	for (Doub xi = -1.;xi <= 1.;xi += 0.01)
	{
		for (Doub eta = -1.;eta <= 1.;eta += 0.01)
		{
			shape.shapes(psis, GradPsi, xi, eta);
			psis.Transpose(psist);
			//psis.Print();
			//psist.Print();
			
			psist.Mult(elcoords, xycoords);
			Doub dx = xycoords[0][0] - Vec[0];
			Doub dy = xycoords[0][1] - Vec[1];
			Doub dist = sqrt(dx*dx + dy*dy);

			NRvector<MatDoub> hhat2(fhhatvel.size());

			if (dist<tol)
			{
				for (Int ivar = 0;ivar < fhhatvel.size();ivar++)
				{
					psist.Mult(fhhatvel[ivar], hhat2[ivar]);
				}
				cout << "XY SOL = " << endl;
				xycoords.Print();
				return hhat2;
			}
		}

	}

}
void SolveEigen2(MatDoub A, MatDoub b, MatDoub& x)
{
	x.assign(A.nrows(), 1, 0.);
	SpMat Aa(A.nrows(), A.nrows());
	for (int i = 0; i < A.nrows(); i++)
	{
		for (int j = 0; j < A.ncols(); j++)
		{
			if (fabs(A[i][j]) > 1.e-6)
			{
				Aa.insert(i, j) = A[i][j];
			}
		}
	}
	Aa.makeCompressed();
	VectorXd bb(A.nrows()), xx(A.nrows());
	for (int i = 0; i < A.nrows(); i++)bb(i) = b[i][0];

	Eigen::SimplicialCholesky<SpMat> solver;
	//	solver.setMode(SimplicialCholeskyLDLT);//31.54
	solver.setMode(SimplicialCholeskyLLT);//32.01 malha 221
	solver.compute(Aa);
	xx = solver.solve(bb);         // use the factorization to solve for the given right hand side

	for (int i = 0; i < A.nrows(); i++)x[i][0] = xx(i);
}

void SolveNR3(MatDoub A, MatDoub b, MatDoub& x)
{
	Int sz = A.nrows();
	VecDoub bv(sz,0.), xv(sz,0.);
	for (Int i = 0; i < sz; i++)
	{
		bv[i] = b[i][0];
	}
	Cholesky chol(A);
	chol.solve(bv, xv);
	x.assign(sz, 1, 0.);
	for (Int i = 0; i < sz; i++)
	{
		x[i][0] = xv[i];
	}
}

void SolveEigen(MatDoub A, MatDoub b, MatDoub& x)
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

std::vector<std::vector<int>>  LineTopology(std::vector<int> ids, Int order)
{
	Int k = 0;

	std::vector<std::vector<int>> vg;
	for (int j = 0; j < ids.size() / order; j++)
	{
		std::vector<int> v;
		for (int i = 0; i < order + 1; i++)
		{
			v.push_back(ids[i + k]);
		}
		vg.push_back(v);
		k += order;
	}
	return vg;
}

void ToMatInt(std::vector<std::vector<int>> in, MatInt& out)
{
	Int  rows = in.size();
	Int cols = in[0].size();
	out.assign(rows, cols, 0.);
	for (Int i = 0; i < rows; i++)for (Int j = 0; j < cols; j++)out[i][j] = in[i][j];
}

void IterativeCho(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material)
{
	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);
	a[0] = 0.; a[1] = 0.;
	b[0] = 30.; b[1] = 0;
	gridmesh::Line(a, b, ndivs, pathbottom);
	gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	//cout << "IDS BOTTOM " << endl;
	//for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

	a[0] = 0.; a[1] = 0.;
	b[0] = 0.; b[1] = 15.;
	gridmesh::Line(a, b, ndivs, pathleft);
	gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	//cout << "IDS idsleft " << endl;
	//for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
	a[0] = 30.; a[1] = 0.;
	b[0] = 30.; b[1] = 5.;
	gridmesh::Line(a, b, ndivs, pathright);
	gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	//cout << "IDS idsright " << endl;
	//for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
	a[0] = 9.99; a[1] = 14.99;
	b[0] = 10.; b[1] = 15.;
	gridmesh::Line(a, b, ndivs, pathdisplace);
	gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	cout << "IDS iddisplace " << endl;
	for (Int i = 0; i < iddisplace.size(); i++)cout << iddisplace[i] << endl;




	Int sz = 2 * meshnodes.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	material->ResetPlasticStrain();
	material->ResetDisplacement();
	material->ResetMat();

	MatDoub displace;
	displace.assign(sz, 1, 0.);
	Doub l = 0, l0 = 0, lamb = 1., lambn = 0, lamb3, diff = 100, diff2 = 100;
	Int counterout = 0;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);
	int ndesi = 10;
	while (counterout < ndesi && fabs(diff2)>0.001)
	{
		std::clock_t start;
		double duration;
		start = std::clock();
		std::cout << "load step = " << counterout << std::endl;
		std::cout << "l= " << l << std::endl;
		Int counter = 0, maxcount = 10;
		Doub err1 = 10., err2 = 10., tol = 1.e-3;
		MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), R;
		Doub dlamb = 1., lambn0 = lamb;
		MatDoub dw(sz, 1, 0.);
		std::cout << "diff = " << diff << std::endl;
		diff = 10;
		while (counter <  maxcount && err1 > tol)
		{
			lambn = lamb;
			material->Assemble(allcoords, meshnodes, meshtopology, KG, FINT, FBODY);
			R = FBODY;
			R *= lamb;
			R -= FINT;

			Int dir, val;
			dir = 1;
			val = 0;
			material->DirichletBC(KG, R, idsbottom, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsbottom, dir, val);

			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsleft, dir, val);


			dir = 1;
			val = 0;
			material->DirichletBC(KG, FBODY, idsbottom, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsbottom, dir, val);

			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsleft, dir, val);

			SolveEigen(KG, R, dws);
			SolveEigen(KG, FBODY, dwb);

			Doub aa = 0.;
			for (int i = 0; i < sz; i++)aa += dwb[i][0] * dwb[i][0];
			Doub bb = 0.;
			MatDoub dwcopy = dw;
			dwcopy += dws;
			for (int i = 0; i < sz; i++)bb += dwb[i][0] * dwcopy[i][0];
			bb *= 2;
			Doub cc = 0.;
			for (int i = 0; i < sz; i++)cc += dwcopy[i][0] * dwcopy[i][0];

			if (counter == 0 && counterout == 0) {
				MatDoub dwbt, mult;
				dwb.Transpose(dwbt);
				dwbt.Mult(dwb, mult);
				mult.Print();
				l0 = sqrt(mult[0][0]);
				l = l0;
			}
			cc -= l * l;
			Doub delta = bb * bb - 4. * aa * cc;
			dlamb = (-bb + sqrt(delta)) / (2. * aa);
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;
			lamb += dlamb;
			displace += dww;
			material->UpdateDisplacement(displace);
			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;
			rnorm = R.NRmatrixNorm();
			normdw = dww.NRmatrixNorm();
			unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;

			Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

			std::cout << " Iteration number = " << counter << " |  |du|/|u| = " << err2 << " |  |R| = " << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << " | time =  <<" << duration1 << std::endl;
			counter++;
			diff = fabs(lamb) - fabs(lambn);
			rtemp = rnorm;
			if (counter == 1)err1 = 10;
		}

		if (rtemp > 10.) {
			l = l0 * Doub(ndesi) / Doub(counter);
			l *= 0.1;
		}
		else {
			l = l0 * Doub(ndesi) / Doub(counter);
		}
		//l = 0.2344;
		diff2 = fabs(lambn0 - lamb);
		cout << " diff2 = " << diff2 << endl;
		material->UpdatePlasticStrain();


		solcount[0] = diff;
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = err2;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		uvf[1] = lamb;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);
		counterout++;
		//std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	}


	if (false)
	{


		MatDoub solpost23;
		solpost23.CopyFromVector(solpost2);
		string names = "loadvsdisplacementnew";
		string exts = ".txt";
		names += exts;
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);


		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(allcoords, meshnodes, meshtopology, material->fdisplace, epsppost);
		string name3 = "epsppostnew";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		material->PostProcess(allcoords, meshnodes, meshtopology, material->fdisplace, solx, soly);
		filename = "soly.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solx.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}


}

void solvedetermtest(string nodestr, string elsstr,  MatDoub hhatinho, Doub c, Doub phi, Doub gamma)
{

	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);


	mesh* gmesh = new mesh(allcoords, meshcoords, meshtopology);


	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(gmesh->GetMeshNodes(), filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(gmesh->GetMeshTopology(), filemesh2);

	std::clock_t start;
	double duration;
	start = std::clock();
	//Doub c = 12.3755;//1.0
	//Doub c = 15.4694; // 1.25
	//Doub c = 18.5633; // 1.5
//	Doub c = 21.6572;//1.75
	//Doub c = 24.7511;//2
	//Doub c = 27.845;//2,0
	//Doub c = 6;//2.5
	//Doub c = 9.28166;//0.75
	//Doub c = 12.3755;
	//Doub c = 6.18777;//0.5
	//Doub c = 16.25;
	Doub thickness = 1.;
	Doub young = 50000.;
	Doub nu = 0.32;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = gamma;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows() * npts;
	Int sz = 2 * meshcoords.nrows();
	mesh* localmesh = new  mesh(allcoords, meshcoords, meshtopology);
	elastoplastic2D< druckerprager >* material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	material->fYC.setup(young, nu, c, phi);
	material->SetMemory(nglobalpts, sz);
	material->UpdateBodyForce(bodyforce);

	//IterativeCho(allcoords, meshcoords, meshtopology, hhatinho, material);//x = desloc y = loadfactor
	IterativeProcess(gmesh, hhatinho, material,10,1.);//x = desloc y = loadfactor


	if (false)
	{

		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(allcoords, meshcoords, meshtopology, material->fdisplace, epsppost);
		string name3 = "epsppostnew2HUM";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		material->PostProcess(allcoords, meshcoords, meshtopology, material->fdisplace, solx, soly);
		filename = "soly2HUM.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solx2HUM.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n Tempo Eigen  " << duration << '\n';
}

void solvedeterm(string nodestr, string elsstr)
{

	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);




	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	std::clock_t start;
	double duration;
	start = std::clock();




	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 10.;
	//Doub young = 100000.;
	//Doub nu = 0.3;
	//Doub c = 23.;
	Doub phi = 30. * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = -20.;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows() * npts;
	Int sz = 2 * meshcoords.nrows();
	MatDoub hhatinho;
	mesh* localmesh = new  mesh(allcoords, meshcoords, meshtopology);
	elastoplastic2D< druckerprager >* material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	material->fYC.setup(young, nu, c, phi);
	material->SetMemory(nglobalpts, sz);
	material->UpdateBodyForce(bodyforce);

	IterativeCho(allcoords, meshcoords, meshtopology, hhatinho, material);//x = desloc y = loadfactor
	//IterativeProcess(allcoords, meshcoords, meshtopology, hhatinho, material);//x = desloc y = loadfactor

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n Tempo Eigen  " << duration << '\n';
}


std::vector<std::vector<double>>  IterativeSlopeStabilityEig(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material)
{

	//std::vector<std::vector< std::vector<Doub > > > allcoords = inmesh->GetAllCoords();
	//MatDoub meshcoords = inmesh ->GetMeshNodes();
	//MatInt meshtopology = inmesh ->GetMeshTopology();


	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);
	//a[0] = 0.;a[1] = 0.;
	//b[0] = 75.;b[1] = 0;
	//gridmesh::Line(a, b, ndivs, pathbottom);
	//gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	////cout << "IDS BOTTOM " << endl;
	////for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

	//a[0] = 0.;a[1] = 0.;
	//b[0] = 0.;b[1] = 40.;
	//gridmesh::Line(a, b, ndivs, pathleft);
	//gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	////cout << "IDS idsleft " << endl;
	////for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
	//a[0] = 75.;a[1] = 0.;
	//b[0] = 75.;b[1] = 30;
	//gridmesh::Line(a, b, ndivs, pathright);
	//gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	////cout << "IDS idsright " << endl;
	////for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
	//a[0] = 34.99;a[1] = 39.99;
	//b[0] = 35.;b[1] = 40.;
	//gridmesh::Line(a, b, ndivs, pathdisplace);
	//gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	////cout << "IDS iddisplace " << endl;
	////for (Int i = 0;i < iddisplace.size();i++)cout << iddisplace[i] << endl;
	//Int sz = 2 * meshcoords.nrows();
	//MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;




	//a[0] = 0.;a[1] = 0.;
	//b[0] = 30.;b[1] = 0;
	//gridmesh::Line(a, b, ndivs, pathbottom);
	//gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	////cout << "IDS BOTTOM " << endl;
	////for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

	//a[0] = 0.;a[1] = 0.;
	//b[0] = 0.;b[1] = 15.;
	//gridmesh::Line(a, b, ndivs, pathleft);
	//gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	////cout << "IDS idsleft " << endl;
	////for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
	//a[0] = 30.;a[1] = 0.;
	//b[0] = 30.;b[1] = 15;
	//gridmesh::Line(a, b, ndivs, pathright);
	//gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	////cout << "IDS idsright " << endl;
	////for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
	//a[0] = 9.99;a[1] = 14.99;
	//b[0] = 10.;b[1] =15.;
	//gridmesh::Line(a, b, ndivs, pathdisplace);
	//gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	//cout << "IDS iddisplace " << endl;
	//for (Int i = 0;i < iddisplace.size();i++)cout << iddisplace[i] << endl;
	//Int sz = 2 * meshcoords.nrows();
	//MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;



	a[0] = 0.; a[1] = 0.;
	b[0] = 50.; b[1] = 0;
	gridmesh::Line(a, b, ndivs, pathbottom);
	gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	//cout << "IDS BOTTOM " << endl;
	//for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

	a[0] = 0.; a[1] = 0.;
	b[0] = 0.; b[1] = 20.;
	gridmesh::Line(a, b, ndivs, pathleft);
	gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	//cout << "IDS idsleft " << endl;
	//for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
	a[0] = 50.; a[1] = 0.;
	b[0] = 50.; b[1] = 10;
	gridmesh::Line(a, b, ndivs, pathright);
	gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	//cout << "IDS idsright " << endl;
	//for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
	a[0] = 19.99; a[1] = 19.99;
	b[0] = 20.; b[1] = 20.;
	gridmesh::Line(a, b, ndivs, pathdisplace);
	gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	cout << "IDS iddisplace " << endl;
	for (Int i = 0; i < iddisplace.size(); i++)cout << iddisplace[i] << endl;
	Int sz = 2 * meshnodes.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	material->ResetPlasticStrain();
	material->ResetDisplacement();
	material->ResetMat();

	MatDoub displace;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l = 0, l0 = 0, lamb = 1., lambn = 0, lamb3, diff = 100, diff2 = 100;
	Int counterout = 0;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);
	int ndesi = 10;
	//	while (counterout < 15 && fabs(diff)>0.1)
	while (counterout < ndesi && fabs(diff2)>0.001)
	{
		std::clock_t start;
		double duration;
		start = std::clock();
		std::cout << "load step = " << counterout << std::endl;
		std::cout << "l= " << l << std::endl;
		Int counter = 0, maxcount = 10;
		Doub err1 = 10., err2 = 10., tol = 1.e-3;
		MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), R;

		//lamb = 1.;
		Doub dlamb = 1., lambn0 = lamb;
		MatDoub dw(sz, 1, 0.);
		//while (counter <  maxcount && err1 > tol && fabs(dlamb)>1.e-3)
		std::cout << "diff = " << diff << std::endl;
		diff = 10;
		//while (counter <  maxcount && err1 > tol && fabs(diff)>0.001)

		while (counter <  maxcount && err1 > tol)
		{
			lambn = lamb;
			//newbodyforce = bodyforce;
			//newbodyforce *= lamb;
			//material->UpdateBodyForce(newbodyforce);

			material->Assemble(allcoords, meshnodes, meshtopology, KG, FINT, FBODY);
			//material->Assemble( *inmesh,KG, FINT, FBODY);



			R = FBODY;
			R *= lamb;
			R -= FINT;

			//FBODY *= 1. / lamb;
			Int dir, val;
			dir = 1;
			val = 0;
			material->DirichletBC(KG, R, idsbottom, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsbottom, dir, val);

			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsleft, dir, val);


			dir = 1;
			val = 0;
			material->DirichletBC(KG, FBODY, idsbottom, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsbottom, dir, val);

			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsleft, dir, val);

			SolveEigen(KG, R, dws);
			SolveEigen(KG, FBODY, dwb);

			Doub aa = 0.;
			for (int i = 0; i < sz; i++)aa += dwb[i][0] * dwb[i][0];
			Doub bb = 0.;
			MatDoub dwcopy = dw;
			dwcopy += dws;
			for (int i = 0; i < sz; i++)bb += dwb[i][0] * dwcopy[i][0];
			bb *= 2;
			Doub cc = 0.;
			for (int i = 0; i < sz; i++)cc += dwcopy[i][0] * dwcopy[i][0];

			if (counter == 0 && counterout == 0) {
				MatDoub dwbt, mult;
				dwb.Transpose(dwbt);
				dwbt.Mult(dwb, mult);
				mult.Print();
				l0 = sqrt(mult[0][0]);
				l = l0;
			}
			cc -= l * l;
			Doub delta = bb * bb - 4. * aa * cc;
			dlamb = (-bb + sqrt(delta)) / (2. * aa);
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;
			lamb += dlamb;
			displace += dww;
			material->UpdateDisplacement(displace);
			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;
			rnorm = R.NRmatrixNorm();
			normdw = dww.NRmatrixNorm();
			unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;

			Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

			std::cout << " Iteration number = " << counter << " |  |du|/|u| = " << err2 << " |  |R| = " << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << " | time =  <<" << duration1 << std::endl;
			counter++;
			diff = fabs(lamb) - fabs(lambn);
			rtemp = rnorm;
			if (counter == 1)err1 = 10;
		}

		if (rtemp > 10.) {
			l = l0 * Doub(ndesi) / Doub(counter);
			l *= 0.1;
		}
		else {
			l = l0 * Doub(ndesi) / Doub(counter);
		}
		//l = 0.2344;
		diff2 = fabs(lambn0 - lamb);
		cout << " diff2 = " << diff2 << endl;
		material->UpdatePlasticStrain();
		//	solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		//	solcount[1] = lamb;

		solcount[0] = diff;
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = err2;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		uvf[1] = lamb;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);
		counterout++;
		//std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	}

	//MatDoub out(2, 1, 0.);
	//out[0][0] = solpost[counterout - 1][0];
	//out[1][0] = solpost[counterout - 1][1];

	if (false)
	{


		MatDoub solpost23;
		solpost23.CopyFromVector(solpost2);
		string names = "loadvsdisplacementnew";
		string exts = ".txt";
		names += exts;
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);


		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(allcoords, meshnodes, meshtopology, material->fdisplace, epsppost);
		string name3 = "epsppostnew";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		material->PostProcess(allcoords, meshnodes, meshtopology, material->fdisplace, solx, soly);
		filename = "soly.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solx.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

	return solpost;

}



void MonteCarloTalude2()
{



	//string nodestr = "nodes-mais-fina.txt";
	//string elsstr = "elements-mais-fina.txt";
	//string nodestr = "nodes-fina.txt";
	//string elsstr = "elements-fina.txt";
	//string nodestrfine = "nodes-extrema-fina.txt";
	//string elsstrfine = "elements-extrema-fina.txt";
//	string nodestrfine = "concentrated-nos.txt";
//	string elsstrfine = "concentrated-els.txt";
	//string nodestr = "nodes-good.txt";
	//string elsstr = "els-good.txt";

	//string nodestr = "nos-110.txt";
	//string elsstr = "els-110.txt";

	//string nodestr = "nos-74.txt";
	//string elsstr = "els-74.txt";

	//string nodestr = "nos-75.txt";
	//string elsstr = "els-75.txt";

	//string nodestr = "sn-nos-91.txt";
	//string elsstr = "sn-els-91.txt";

	//string nodestrfine = "sn-nos-91.txt";
//	string elsstrfine = "sn-els-91.txt";

	//string nodestr = "sn-nos-182.txt";
	//string elsstr = "sn-els-182.txt";

	//string nodestr = "sn-nos-200.txt";//1.33
	//string elsstr = "sn-els-200.txt";

	//string nodestr = "nos-sz.txt";
	//string elsstr = "els-sz.txt";

	//string nodestr = "nos-sz-v2.txt";
	//string elsstr = "els-sz-v2.txt";
	//string nodestr = "nos-101.txt";
	//string elsstr = "els-101.txt";

	//string nodestr = "nos-156.txt";
	//string elsstr = "els-156.txt";

	//string nodestrfine = "nos-igual-221.txt";
	//string elsstrfine = "els-igual-221.txt";

	//string nodestrfine = "sn-nos-121-b.txt";//1.337 //16.6 s
	//string elsstrfine = "sn-els-121-b.txt";

	//string nodestr = "nos-sz.txt";
	//string elsstr = "els-sz.txt";

	//string nodestrfine = "nos-sz.txt";
	//string elsstrfine = "els-sz.txt";

	//MatDoub  meshcoords, elcoords;
	//MatInt meshtopology;
	//std::vector<std::vector<std::vector<Doub>>> allcoords;
	//ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);
	//mesh  coarsemesh =   mesh(allcoords, meshcoords, meshtopology);

	//string nodestrfine = "nos-132-c3.txt";
	//string elsstrfine = "els-132-c3.txt";


	//string nodestrfine = "nos-184-c.txt";
	//string elsstrfine = "els-184-c.txt";

	string nodestrfine = "D:/DClib/src/nos-228.txt";
	string elsstrfine = "D:/DClib/src/els-228.txt";


	MatDoub  meshcoordsfine, elcoordsfine;
	MatInt meshtopologyfine;
	std::vector<std::vector<std::vector<Doub>>> allcoordsfine, test;
	ReadMesh(allcoordsfine, meshcoordsfine, meshtopologyfine, elsstrfine, nodestrfine);

	mesh* finemesh = new  mesh(allcoordsfine, meshcoordsfine, meshtopologyfine);
	//meshcoords = finemesh.GetMeshNodes();
	//meshtopology = finemesh.GetMeshTopology();
	//test = finemesh.GetAllCoords();
//	cout << "sz = " << test[0].size() << endl;
	//mesh * finemesh = new  mesh(allcoords, meshcoords, meshtopology);

	std::clock_t start;
	double duration;
	start = std::clock();

	cout << "\n starting stochastic simulation " << endl;


	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 16.25;
	Doub phi = 20 * M_PI / 180.;

	//Doub young =100000.;
	//Doub nu = 0.3;
	//Doub c = 23;
	//Doub phi = 0.01 * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = -20.;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = finemesh->GetMeshTopology().nrows() * npts;
	Int sz = 2 * finemesh->GetMeshNodes().nrows();

	Doub Lx = 40;//(*Correlation length in x direction*)
	Doub Ly = 4;//(*Correlation length in y direction*)

	Int samples = 1000, expansionorder = 100;
	Doub sig = 0.3;
	Int type = 1;
	KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF(*finemesh,  order, Lx, Ly, sig, type, samples, expansionorder);

	//system("mkdir -p D:\DClib\results");

	//boost::filesystem::create_directories("D:\DClib\results");
	int check;
	//char* dirname = "D:\DClib\results";
	//auto s = std::to_string( rand() % 30 + 1985);
	auto s = std::to_string(4);
	string namefolder = "D:/posprocess-SlopeProject-40-4";

	namefolder += s;

	char* cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());

	check = mkdir(cstr);

	string datafile = namefolder;
	datafile += "/DATA.txt";
	std::ofstream file(datafile);
	file << " Young = " << young << " | nu = " << nu << endl;
	file << " c = " << c << " | phi = " << phi << endl;
	file << " bodyforce = " << bodyforce[1][0] << endl;
	file << " Mesh  = " << nodestrfine << endl;
	file << " samples = " << samples << " | expansion order = " << expansionorder << " | func type = " << type << endl;
	file << "Lx = " << Lx << " | Ly = " << Ly << " variance = " << sig << endl;
	VecComplex  val; MatDoub  vec, HHAT;
	NRmatrix<MatDoub> randomfield;
	//objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, randomfield);


	std::vector<std::vector<double>> errpost;
	objKLGalerkinRF->SolveGenEigValProblem(val, vec, randomfield, errpost);
	datafile = namefolder;
	datafile += "/error2.txt";
	std::ofstream fileerro(datafile);
	OutPutPost(errpost, fileerro);

	delete objKLGalerkinRF;


	datafile = namefolder;
	datafile += "/vec.txt";
	std::ofstream filevec(datafile);
	OutPutPost(vec, filevec);

	datafile = namefolder;
	datafile += "/val.txt";
	std::ofstream fileval(datafile);
	for (Int j = 0; j < val.size(); j++)
	{
		fileval << val[j].real() << endl;
	}

	datafile = namefolder;
	datafile += "/meshcoords.txt";
	std::ofstream filemesh1(datafile);
	OutPutPost(finemesh->GetMeshNodes(), filemesh1);

	datafile = namefolder;
	datafile += "/meshtopology.txt";
	std::ofstream filemesh2(datafile);

	OutPutPost(finemesh->GetMeshTopology(), filemesh2);


	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n  simualtion time  " << duration << '\n';

	Int  fail = 0;
	cout << "\n starting Monte Carlo " << endl;
	start = std::clock();
	Int postprintfreq = 50;
	Doub sum = 0.;
	MatDoub solpost(samples, 2, 0.), solpost2(samples, 1, 0.);
	Doub soldatamin = 10.;
	Doub soldatamax = -10;
	std::vector<double> solvec;


	delete finemesh;
	//postprocess post;

	elastoplastic2D< druckerprager >* material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order);
	for (Int imc = 0; imc < samples; imc++)
	{


		Int nrandomvars = randomfield.nrows();
		Int nrowss = randomfield[0][0].nrows();
		VecDoub mean(nrandomvars, 0.), var(nrandomvars, 0.);
		MatDoub hhatinho(randomfield[0][0].nrows(), randomfield.nrows(), 0.), posttest(randomfield[0][0].nrows(), 1, 0.);
		MatDoub coesmat(randomfield[0][0].nrows(), 1, 0.), phimat(randomfield[0][0].nrows(), 1, 0.);
		for (Int ivar = 0; ivar < nrandomvars; ivar++)
		{
			for (Int isample = 0; isample < nrowss; isample++)
			{
				hhatinho[isample][ivar] = randomfield[ivar][0][isample][imc];
				//hhatinho[i][0] = 0.;
				mean[ivar] += hhatinho[isample][ivar];
				mean[ivar] /= (isample + 1);
				var[ivar] += (mean[ivar] - hhatinho[isample][ivar]) * (mean[ivar] - hhatinho[isample][ivar]);
				posttest[isample][0] = hhatinho[isample][ivar];
			}
		}
		//VecDoub coordsvec(2, 0.);
		///coordsvec[0] = 32.5;
		//coordsvec[1] = 39.5;
		//FindEl(allcoords, meshtopology, coordsvec, hhatinho);
		//finemesh->FindSolution(coordsvec, posttest);

		//std::vector<std::vector<double>> finefield, coarsefield;
		//post.PostProcess(finemesh, hhatinho, finefield);
		//std::ofstream file("finefield.txt");
		//OutPutPost(finefield, file);

		//posttest = finemesh.TransferSolution(coarsemesh, hhatinho);

		//post.PostProcess(coarsemesh, posttest, coarsefield);
		//std::ofstream file222222("coarsefield.txt");
		//OutPutPost(coarsefield, file222222);



		for (int i = 0; i < randomfield[0][0].nrows(); i++)coesmat[i][0] = hhatinho[i][0];
		for (int i = 0; i < randomfield[0][0].nrows(); i++)phimat[i][0] = hhatinho[i][1];

		std::clock_t start1;
		double duration1;
		start1 = std::clock();


		material->fYC.setup(young, nu, c, phi);
		material->SetMemory(nglobalpts, sz);
		material->UpdateBodyForce(bodyforce);
		material->SetRandomField(hhatinho);

		//cout << "all cc" << finemesh->GetAllCoords()[0].size() << endl;

		start1 = std::clock();
		//std::vector<std::vector<double>>  sol = IterativeProcessSlope(finemesh, hhatinho, material);//x = desloc y = loadfactor
		std::vector<std::vector<double>>  sol = IterativeSlopeStabilityEig(allcoordsfine, meshcoordsfine, meshtopologyfine, hhatinho, material);//x = desloc y = loadfactor
		duration1 = (std::clock() - start1) / (double)CLOCKS_PER_SEC;
		std::cout << "IterativeProcessSlope time " << duration1 << std::endl;



		start1 = std::clock();

		MatDoub solpost23;
		solpost23.CopyFromVector(sol);

		Int last = solpost23.nrows() - 1;
		Doub data = solpost23[last][1];
		solvec.push_back(data);

		double min = *min_element(solvec.begin(), solvec.end());
		double max = *max_element(solvec.begin(), solvec.end());


		//solpost.Print();

		//if (soldatamin>min  || soldatamax<max)
	//	{

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

		//	solcount[0] = diff;
		//	solcount[1] = lamb;
		//	solcount[2] = err1;
		//	solcount[3] = err2;
		//	solcount[4] = diff;
		///	solcount[5] = diff2;
		//	solcount[6] = counterout;

		if (false) {


			filename = namefolder;
			std::vector<std::vector<double>> hhatx;
			string name = "/Coesao";
			ext = ".txt";
			filename += name;
			s = std::to_string(imc);
			filename += s;
			filename += ext;
			material->PostProcess(allcoordsfine, meshcoordsfine, meshtopologyfine, 0, hhatinho, hhatx);
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
			material->PostProcess(allcoordsfine, meshcoordsfine, meshtopologyfine, 1, hhatinho, hhatx2);
			std::ofstream filesss(filename);
			OutPutPost(hhatx2, filesss);

			filename = namefolder;
			string names = "/FxU";
			string exts = ".txt";
			filename += names;
			auto ss = std::to_string(imc);
			filename += ss;
			filename += exts;
			std::ofstream file8(filename);
			OutPutFile(solpost23, file8);

			filename = namefolder;
			std::vector<std::vector<double>> solx, soly;
			material->PostProcess(allcoordsfine, meshcoordsfine, meshtopologyfine, material->fdisplace, solx, soly);
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


			//filename = namefolder;
			//std::vector<std::vector<double>> hc, hphi;
			//material->PostProcessIntegrationPointVar(allcoords, meshtopology, coesmat, hc);
			//material->PostProcessIntegrationPointVar(allcoords, meshtopology, phimat, hphi);
			//name2 = "/CoesionIntPoint";
			//ext2 = ".txt";
			//filename += name2;
			//s2 = std::to_string(data);
			//filename += s2;
			//filename += ext2;
			//std::ofstream filehc(filename);
			//OutPutPost(hc, filehc);

			//filename = namefolder;
			//name2 = "/PhiIntPoint";
			//ext2 = ".txt";
			//filename += name2;
			//s2 = std::to_string(data);
			//filename += s2;
			//filename += ext2;
			//std::ofstream filephi(filename);
			//OutPutPost(hphi, filephi);


			filename = namefolder;
			name2 = "/hhatinho";
			ext2 = ".txt";
			filename += name2;
			s2 = std::to_string(imc);
			filename += s2;
			filename += ext2;
			std::ofstream file222(filename);
			OutPutPost(hhatinho, file222);
		}

		filename = namefolder;
		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(allcoordsfine, meshcoordsfine, meshtopologyfine, material->fdisplace, epsppost);
		string name3 = "/plasticsqrtj2";
		string ext3 = ".txt";
		filename += name3;
		auto s3 = std::to_string(imc);
		filename += s3;
		filename += ext3;
		std::ofstream file3(filename);
		OutPutPost(epsppost, file3);


		if (false) {


			std::cout << "min = " << min << std::endl;
			std::cout << "max = " << max << std::endl;
			std::cout << "soldatamax = " << soldatamax << std::endl;
			std::cout << "soldatamin = " << soldatamin << std::endl;

			if (soldatamin > min)
			{
				soldatamin = min;
			}
			else {
				soldatamax = max;
			}
		}

		//	}

		string filename2 = namefolder;
		filename2 += "/montecarlosafetyfactor.txt";
		std::ofstream file23(filename2);
		OutPutFile1var(solpost2, file23);

		if (data <= 1.)
		{
			fail++;
		}
		cout << " mc it = " << imc << " | Current safety fator = " << data << " time on plastic fem integration = " << duration1 << endl;
		solpost2[imc][0] = data;
		//duration1 = (std::clock() - start1) / (double)CLOCKS_PER_SEC;
		//std::cout << " post time time " << duration1 << std::endl;
		//
	}
	string filename = namefolder;
	filename += "/montecarlosafetyfactor.txt";
	std::ofstream file23(filename);
	OutPutFile1var(solpost2, file23);
	file << "failue probability = " << Doub(fail) / Doub(samples) << endl;
	cout << "failue probability = " << Doub(fail) / Doub(samples);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n Monte Carlo simualtion time  " << duration << '\n';
	file << "\n Monte Carlo simualtion time  " << duration << '\n';

	delete material;

}

void IterativeProcessShearRed2(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material) {
	Doub c = material->fYC.GetCoes();
	Doub phi = material->fYC.GetPhi();
	Doub  young = material->fYC.GetYoung();
	Doub nu = material->fYC.GetNu();
	//Doub  young = 14000, nu = 0.3, phi = 35 * M_PI / 180.;
	//Doub c = 12.3755;
	//Doub c = 19.8008;//FS = 1.6
	//Doub c = 16.25;//	1.31307
	//Doub c = 6.18777;//0.5
	//Doub c = 9.28166;//0.75
	// Doub c = 12.3755; // 1.0
	//Doub c = 15.4694; // 1.25
	//Doub c = 20.; // 1.5
	//Doub c = 21.6572;//1.75
	//Doub c = 24.7511;//2,0
	//Doub c = 27.875;//2,5
	//Doub c = 30.9389;//2,5
	//Doub c = 10;

	//Doub c = 2.81369;//0.5
	//Doub c = 4.22054;//0.75
	//Doub c = 5.62739;//1
	//Doub c = 7.03423;//1.25
	//Doub c = 8.44108;//1.5
	//Doub c = 9.84793;//1.75
	//Doub c = 11.2548;//2
	//Doub c = 14.06685;//2.5


	//Doub c = 17.03;//0.5
	//Doub c = 25.5522;//0.75
	//Doub c = 34.0695;//1
	//Doub c = 42.5869;//1.25
	//Doub c = 51.1043;//1.5
	//Doub c = 59.6217;//1.75
	// Doub c = 68.1391;//2
	//Doub c = 85.1739;//2.5
	MatDoub bodyforce(2, 1, 0.);
	bodyforce[1][0] = -20.;
	material->UpdateBodyForce(bodyforce);

	material->fYC.setup(young, nu, c, phi);
	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);
	a[0] = 0.; a[1] = 0.;
	b[0] = 50.; b[1] = 0;
	gridmesh::Line(a, b, ndivs, pathbottom);
	gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	a[0] = 0.; a[1] = 0.;
	b[0] = 0.; b[1] = 20.;
	gridmesh::Line(a, b, ndivs, pathleft);
	gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	a[0] = 50.; a[1] = 0.;
	b[0] = 50.; b[1] = 10;
	gridmesh::Line(a, b, ndivs, pathright);
	gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	a[0] = 19.99; a[1] = 19.99;
	b[0] = 20.; b[1] = 20.;
	gridmesh::Line(a, b, ndivs, pathdisplace);
	gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	Int sz = 2 * meshnodes.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	MatDoub displace, displace0, R, sol;
	displace.assign(sz, 1, 0.);

	Int niter = 30;
	//Doub fac =  1.46;
	//Doub fac = 0.79;
	//Doub fac = 1.258;//1.6
	Doub fac = 0.5;
	Doub delta = 0.1;
	Int postcounter = 0;



	Doub FSmin = 0.;
	Doub FSmax = 1000000000.;
	Doub FS = 0.01;
	Doub c0 = c;
	Doub phi0 = phi;

	c = c0 / FS;
	phi = atan(tan(phi0 / FS));
	material->fYC.reset();
	material->fYC.setup(young, nu, c, phi);

	displace0 = displace;
	while ((FSmax - FSmin) / FS > 1.e-2)
	{
		Int counter = 0;
		R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.), displace.assign(sz, 1, 0.);
		sol.assign(sz, 1, 0.);
		std::cout << " FS = " << FS << " FSmin " << FSmin << " FSmax " << FSmax << std::endl;

		do
		{
			Doub rnorm = R.NRmatrixNorm();
			material->Assemble(allcoords, meshnodes, meshtopology, KG, FINT, FBODY);
			R = FBODY;
			//R *= fac;
			R -= FINT;

			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);

			SolveEigen2(KG, R, sol);

			displace += sol;
			Doub u = fabs(displace[2 * iddisplace[0] + 1][0]);
			material->UpdateDisplacement(displace);

			std::cout << " R norm = " << R.NRmatrixNorm() << " u " << u << " | phi0/phi = " << phi0 / phi << " | c0/c = " << c0 / c << " | c = " << c << " | phi = " << phi << std::endl;
			counter++;
			postcounter++;

			if ((R.NRmatrixNorm() > rnorm && counter >= 20) || R.NRmatrixNorm() > 10.e3)break;
			//if ( R.NRmatrixNorm() > 10.e3) {
			//	FS *= 0.9;
			//	c = c0 / FS;
			//	phi = atan(tan(phi0) / FS);
			//	material->fYC.reset();
			//	material->fYC.setup(young, nu, c, phi);
			//}
		} while (R.NRmatrixNorm() > 0.1 && counter <= 20);

		if (R.NRmatrixNorm() > 1)
		{
			FSmax = 1.0*FS;
			FS =1.* (FSmin + FSmax) / 2.;
			//material->UpdateDisplacement(displace0);
			material->ResetMat();
		}
		else {
			FSmin =FS;
			FS = 1. / ((1. / FSmin + 1. / FSmax) / 2.);
			//FS *= 1.01;
			material->fYC.reset();
			material->fYC.setup(young, nu, c, phi);
			material->UpdatePlasticStrain();
			displace0 = displace;
			//material->ResetMat();
		}

		c = c0 / FS;
		phi = atan(tan(phi0 ) / FS);
		material->fYC.reset();
		material->fYC.setup(young, nu, c, phi);


	}
	material->UpdatePlasticStrain();
	if (false)
	{

		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(allcoords, meshnodes, meshtopology, material->fdisplace, epsppost);
		string name3 = "epsppostnewHUM";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		material->PostProcess(allcoords, meshnodes, meshtopology, material->fdisplace, solx, soly);
		filename = "solyHUM.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solxHUM.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

}


std::vector<std::vector<double>>   IterativeProcessShearRed3(std::vector<std::vector< std::vector<Doub > > > allcoords, MatDoub meshnodes, MatInt meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager >* material)
{

	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);

	a[0] = 0.; a[1] = 0.;
	b[0] = 50.; b[1] = 0;
	gridmesh::Line(a, b, ndivs, pathbottom);
	gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	//cout << "IDS BOTTOM " << endl;
	//for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

	a[0] = 0.; a[1] = 0.;
	b[0] = 0.; b[1] = 20.;
	gridmesh::Line(a, b, ndivs, pathleft);
	gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	//cout << "IDS idsleft " << endl;
	//for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
	a[0] = 50.; a[1] = 0.;
	b[0] = 50.; b[1] = 10;
	gridmesh::Line(a, b, ndivs, pathright);
	gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	//cout << "IDS idsright " << endl;
	//for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
	a[0] = 19.99; a[1] = 19.99;
	b[0] = 20.; b[1] = 20.;
	gridmesh::Line(a, b, ndivs, pathdisplace);
	gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	cout << "IDS iddisplace " << endl;
	for (Int i = 0; i < iddisplace.size(); i++)cout << iddisplace[i] << endl;
	Int sz = 2 * meshnodes.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
	Doub rtemp = 10.;

	material->ResetPlasticStrain();
	material->ResetDisplacement();
	material->ResetMat();

	MatDoub displace, displace0;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l = 0, l0 = 0, lamb = 1., lambn = 0, dlamb = 0., lamb3, diff = 100, diff2 = 100;
	Int counterout = 0, maxcountout = 80;
	std::vector<double> solcount(7, 0.), uvf(2, 0.);
	std::vector<std::vector<double>> solpost, solpost2;
	solpost.push_back(solcount);
	//int ndesi =10;//65 s u = 0.312275 lamb=1.35255
	//Doub dlamb0 = 1.1;
	//int ndesi = 5;//52 s u = 0.0969881 lamb = 1.34115
	//Doub dlamb0 = 1.1;

	//int ndesi = 5;//42 s u =  0.121154 lamb =1.34468
	//Doub dlamb0 = 1.5;

	//int ndesi = 3;//35 s u = u = 0.127228 lamb =1.34519
	//Doub dlamb0 = 2.;

	//int ndesi = 10;//60 s 0.21547 1.35037
	//Doub dlamb0 = 0.5;

	//int ndesi = 7;//63 s 0.0977655 1.34152
	//Doub dlamb0 = 0.5;
	// 

//	int ndesi = 8;//74 s 0.497378 1.35437  //FUNCIONA PARA MALHA 228 para cima e FS 1.3
//	Doub dlamb0 =	0.8;




	Doub c = material->fYC.GetCoes();
	Doub phi = material->fYC.GetPhi();
	Doub  young = material->fYC.GetYoung();
	Doub nu = material->fYC.GetNu();







	int ndesi = 10;//74 s 0.497378 1.35437  //FUNCIONA PARA MALHA 228 para cima e FS 1.3
	Doub dlamb0 = 0.2;
	//	while (counterout < 15 && fabs(diff)>0.1)

	MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
	material->Assemble(allcoords, meshnodes, meshtopology, KG, FINT, FBODY);
	R = FBODY;
	R -= FINT;
	InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);

	Doub c0 = c;
	Doub phi0 = phi;
	c = c0 / lamb;
	phi = atan(tan(phi0 / lamb));
	material->fYC.reset();
	material->fYC.setup(young, nu, c, phi);

	SolveEigen(KG, R, dws);




	SolveEigen(KG, FBODY, dwb);
	MatDoub dwbt, mult;
	dwb.Transpose(dwbt);
	dwbt.Mult(dwb, mult);
	l0 = sqrt(dlamb0 * dlamb0 * mult[0][0]);
	l = l0;
	dlamb = computelamda(dwb, dws, dw, l);
	lamb = 0.001;
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
		Doub rnorm = 10.;
		do
		{

			c = c0 / lamb;
			phi = atan(tan(phi0 / lamb));
			material->fYC.reset();
			material->fYC.setup(young, nu, c, phi);
			material->Assemble(allcoords, meshnodes, meshtopology, KG, FINT, FBODY);
			std::cout << "1 KG.Norm()" << KG.NRmatrixNorm() << endl;
			R = FBODY;
			R -= FINT;
			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);
			SolveEigen(KG, R, dws);

			material->fYC.reset();
			material->fYC.setup(young, nu, c0, phi0);
			material->Assemble(allcoords, meshnodes, meshtopology, KG, FINT, FBODY);
			std::cout << "KG.Norm()" << KG.NRmatrixNorm() << endl;
			R = FBODY;
			R -= FINT;
			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);

			SolveEigen(KG, R, dwb);

			std::cout << "R.Norm()" << R.NRmatrixNorm() << endl;
			std::cout << "dws.Norm()" << dws.NRmatrixNorm() << endl;
			std::cout << "dwb.Norm()" << dwb.NRmatrixNorm() << endl;
			std::cout << "dw.Norm()" << dw.NRmatrixNorm() << endl;
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
			material->UpdateDisplacement(displace);

			rnorm = R.NRmatrixNorm();
			Doub normdw = dww.NRmatrixNorm();
			Doub unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;

			Doub duration1 = (std::clock() - start) / (double)CLOCKS_PER_SEC;

			std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lamb  = " << lamb << " | time =  <<" << duration1 << std::endl;
			counter++;

			rtemp = rnorm;
			if (counter == 1)rnorm = 10;


		} while (counter < maxcount && rnorm > 0.1);


		if (rnorm > 1)
		{
			std::cout << "Convergence failed. \n";
			material->UpdateDisplacement(displace0);
			displace = displace0;
			//ndesi--;//74 s 0.497378 1.35437
			if (isnan(dlamb) == 1) {
				dlamb0 *= 0.5;
			}
			else {
				dlamb0 *= 0.9;
			}


			//MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), dw(sz, 1, 0.), R;
			dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

			material->Assemble(allcoords, meshnodes, meshtopology, KG, FINT, FBODY);
			R = FBODY;
			R *= lamb;
			R -= FINT;
			InserBC(KG, R, FBODY, idsbottom, idsright, idsleft, material);
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

			l *= Doub(ndesi) / Doub(counter);
			//if (l > 10)l = 10;
			diff2 = fabs(lambn0 - lamb);
			material->UpdatePlasticStrain();
		}

		counterout++;
		dws.assign(sz, 1, 0.), dwb.assign(sz, 1, 0.), dww.assign(sz, 1, 0.), R.assign(sz, 1, 0.), dw.assign(sz, 1, 0.), R.assign(sz, 1, 0.), FBODY.assign(sz, 1, 0.), FINT.assign(sz, 1, 0.);

		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = err2;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;

		uvf[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		uvf[1] = lamb;

		solpost2.push_back(uvf);
		solpost.push_back(solcount);

		//std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	} while (counterout <= maxcountout && fabs(diff2) > 0.001);
	string names = "fxUu";
	auto sss = std::to_string(cccccccc);
	names += sss;
	MatDoub solpost23;
	solpost23.CopyFromVector(solpost2);
	string exts = ".txt";
	names += exts;
	std::ofstream file8(names);
	OutPutFile(solpost23, file8);
	cccccccc++;

	if (false)
	{


		MatDoub solpost23;
		solpost23.CopyFromVector(solpost2);
		string names = "loadvsdisplacementnew";
		string exts = ".txt";
		names += exts;
		std::ofstream file8(names);
		OutPutFile(solpost23, file8);


		std::vector<std::vector<double>> epsppost;
		material->PostProcessIntegrationPointVar(allcoords, meshnodes, meshtopology, material->fdisplace, epsppost);
		string name3 = "epsppostnew";
		string ext3 = ".txt";
		name3 += ext3;
		std::ofstream file3(name3);
		OutPutPost(epsppost, file3);

		string filename;
		std::vector<std::vector<double>> solx, soly;
		material->PostProcess(allcoords, meshnodes, meshtopology, material->fdisplace, solx, soly);
		filename = "soly.txt";
		std::ofstream file2(filename);
		OutPutPost(soly, file2);
		filename = "solx.txt";
		std::ofstream file22(filename);
		OutPutPost(solx, file22);
	}

	return solpost;

}


//Int ndivs = 100000;
//MatDoub pathbottom, pathleft, pathright, pathdisplace;
//std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
//VecDoub a(2), b(2);
//a[0] = 0.; a[1] = 0.;
//b[0] = 50.; b[1] = 0;
//gridmesh::Line(a, b, ndivs, pathbottom);
//gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
//a[0] = 0.; a[1] = 0.;
//b[0] = 0.; b[1] = 20.;
//gridmesh::Line(a, b, ndivs, pathleft);
//gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
//a[0] = 50.; a[1] = 0.;
//b[0] = 50.; b[1] = 10;
//gridmesh::Line(a, b, ndivs, pathright);
//gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
//a[0] = 19.99; a[1] = 19.99;
//b[0] = 20.; b[1] = 20.;
//gridmesh::Line(a, b, ndivs, pathdisplace);
//gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
//Int sz = 2 * meshnodes.nrows();
//MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
//Doub rtemp = 10.;


//Int ndivs = 100000;
//MatDoub pathbottom, pathleft, pathright, pathdisplace;
//std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
//VecDoub a(2), b(2);
//
//a[0] = 0.; a[1] = 0.;
//b[0] = 36.; b[1] = 0;
//gridmesh::Line(a, b, ndivs, pathbottom);
//gridmesh::FindIdsInPath(pathbottom, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsbottom);
////cout << "IDS BOTTOM " << endl;
////for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;
//
//a[0] = 0.; a[1] = 0.;
//b[0] = 0.; b[1] = 10.;
//gridmesh::Line(a, b, ndivs, pathleft);
//gridmesh::FindIdsInPath(pathleft, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsleft);
////cout << "IDS idsleft " << endl;
////for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
//a[0] = 36.; a[1] = 0.;
//b[0] = 36.; b[1] = 18.;
//gridmesh::Line(a, b, ndivs, pathright);
//gridmesh::FindIdsInPath(pathright, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), idsright);
////cout << "IDS idsright " << endl;
////for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
//a[0] = 23.99; a[1] = 17.99;
//b[0] = 24.; b[1] = 18.;
//gridmesh::Line(a, b, ndivs, pathdisplace);
//gridmesh::FindIdsInPath(pathdisplace, gmesh->GetAllCoords(), gmesh->GetMeshTopology(), iddisplace);
//cout << "IDS iddisplace " << endl;
//for (Int i = 0; i < iddisplace.size(); i++)cout << iddisplace[i] << endl;
//Int sz = 2 * gmesh->GetMeshNodes().nrows();
//MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
//Doub rtemp = 10.;



//Int ndivs = 1000;
//MatDoub pathbottom, pathleft, pathright, pathdisplace;
//std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
//VecDoub a(2), b(2);

//a[0] = 0.; a[1] = 0.;
//b[0] = 20.; b[1] = 0;
//gridmesh::Line(a, b, ndivs, pathbottom);
//gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
////cout << "IDS BOTTOM " << endl;
////for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

//a[0] = 0.; a[1] = 0.;
//b[0] = 0.; b[1] = 4.;
//gridmesh::Line(a, b, ndivs, pathleft);
//gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
////cout << "IDS idsleft " << endl;
////for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
//a[0] = 20.; a[1] = 0.;
//b[0] = 20.; b[1] = 10;
//gridmesh::Line(a, b, ndivs, pathright);
//gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
////cout << "IDS idsright " << endl;
////for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
//a[0] = 9.99; a[1] = 9.99;
//b[0] = 10.; b[1] = 10.;
//gridmesh::Line(a, b, ndivs, pathdisplace);
//gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
//cout << "IDS iddisplace " << endl;
//for (Int i = 0; i < iddisplace.size(); i++)cout << iddisplace[i] << endl;
//Int sz = 2 * meshnodes.nrows();
//MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
//Doub rtemp = 10.;


//Int ndivs = 1000;
//MatDoub pathbottom, pathleft, pathright, pathdisplace;
//std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
//VecDoub a(2), b(2);

//a[0] = 0.; a[1] = 0.;
//b[0] = 0.3; b[1] = 0;
//gridmesh::Line(a, b, ndivs, pathbottom);
//gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
////cout << "IDS BOTTOM " << endl;
////for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

//a[0] = 0.; a[1] = 0.;
//b[0] = 0.; b[1] = 0.06;
//gridmesh::Line(a, b, ndivs, pathleft);
//gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
////cout << "IDS idsleft " << endl;
////for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
//a[0] = 0.3; a[1] = 0.;
//b[0] = 0.3; b[1] = 0.3;
//gridmesh::Line(a, b, ndivs, pathright);
//gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
////cout << "IDS idsright " << endl;
////for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
//a[0] = 0.1199; a[1] = 0.30;
//b[0] = 0.12; b[1] = 0.30;
//gridmesh::Line(a, b, ndivs, pathdisplace);
//gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
////	cout << "IDS iddisplace " << endl;
//	//for (Int i = 0; i < iddisplace.size(); i++)cout << iddisplace[i] << endl;
//Int sz = 2 * meshnodes.nrows();
//MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;
//Doub rtemp = 10.;

