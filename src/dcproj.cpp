
// dcproj.cpp : Defines the entry point for the application.
//

#include "dcproj.h"
#include "mesh.h"
#include <thread>
//#include <boost/thread.hpp>
//using namespace boost; 
using namespace std;

void ReadMatDoub(MatDoub& matdoub, std::string  file);
void ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord);
template <class T>
void OutPutPost(NRmatrix<T>& postdata, std::ofstream& file);
template <class T>
void vecstr_to_vec(std::vector<std::string> vs, std::vector<T>& ret);
std::vector<Int> vecstr_to_vecint(std::vector<string> vs);
std::vector<Doub> vecstr_to_vecdoub(std::vector<string> vs);
void myTreads(int a, int b, slopeproject* slopeobj2, string traedN);
void leakcraw();

class SynchronizedFile {
public:
    SynchronizedFile (const string& path) : _path(path) {
        // Open file for writing...
    }

    void write (const string& dataToWrite) {
        // Write to the file in a synchronized manner (described below)...
    }

private:
    string _path;
};

class Writer {
public:
    Writer (std::shared_ptr<SynchronizedFile> sf) : _sf(sf) {}

    void someFunctionThatWritesToFile () {
        // Do some work...
        _sf->write("Some data to write...");
    }
private:
    std::shared_ptr<SynchronizedFile> _sf;
};

//<<<<<<< HEAD
void myTreads(int a, int b, slopeproject* slopeobj2,string namefolder3)
{
	slopeobj2->MonteCarloGIM(a, b, false, namefolder3);
    delete slopeobj2;
}

void myTreadsSRM(int a, int b, slopeproject* slopeobj2,string namefolder3)
{
	slopeobj2->MonteCarloSRM(a, b, false, namefolder3);
    delete slopeobj2;
}

void mainlinux(int simtype,int comeco,int fim)
{
    
    //cout<< "simtype"<< simtype;
	//string nodestr = "/home/diogo/projects/dcproj/nos-132-c3.txt";
	//string elsstr = "/home/diogo/projects/dcproj/els-132-c3.txt";
    
    //string nodestr = "/home/diogo/projects/dcproj/nos-cho.txt";
	//string elsstr = "/home/diogo/projects/dcproj/els-cho.txt";
    
    string nodestr = "/home/diogo/projects/dcproj/nos-414.txt";
	string elsstr = "/home/diogo/projects/dcproj/els-414.txt";


	MatDoub hhatinho;
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);
	std::ofstream filemesh2("meshtopology.txt");
	OutPutPost(meshtopology, filemesh2);

	//Doub c = 18.5633, phi = 20 * M_PI / 180., gamma = -20.;//1.5
	Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5

	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	//Doub young = 100000.;
	//Doub nu = 0.3;
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

	elastoplastic2D< druckerprager >* mat0 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat1 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat2 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat3 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat4 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat5 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat6 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat7 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat8 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat9 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat10 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
    elastoplastic2D< druckerprager >* mat11 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	
	mesh* mesh0 = new mesh(mat0, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh1 = new mesh(mat1, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh2 = new mesh(mat2, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh3 = new mesh(mat3, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh4 = new mesh(mat4, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh5 = new mesh(mat5, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh6 = new mesh(mat6, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh7 = new mesh(mat7, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh8 = new mesh(mat8, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh9 = new mesh(mat9, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh10 = new mesh(mat10, allcoords, meshcoords, meshtopology, hhatinho);
    mesh* mesh11 = new mesh(mat11, allcoords, meshcoords, meshtopology, hhatinho);

	mat0->fYC.setup(young, nu, c, phi);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);
    
    mat1->fYC.setup(young, nu, c, phi);
	mat1->SetMemory(nglobalpts, sz);
	mat1->UpdateBodyForce(bodyforce);

    mat2->fYC.setup(young, nu, c, phi);
	mat2->SetMemory(nglobalpts, sz);
	mat2->UpdateBodyForce(bodyforce);

    mat3->fYC.setup(young, nu, c, phi);
	mat3->SetMemory(nglobalpts, sz);
	mat3->UpdateBodyForce(bodyforce);

    mat4->fYC.setup(young, nu, c, phi);
	mat4->SetMemory(nglobalpts, sz);
	mat4->UpdateBodyForce(bodyforce);

    mat5->fYC.setup(young, nu, c, phi);
	mat5->SetMemory(nglobalpts, sz);
	mat5->UpdateBodyForce(bodyforce);

    mat6->fYC.setup(young, nu, c, phi);
	mat6->SetMemory(nglobalpts, sz);
	mat6->UpdateBodyForce(bodyforce);

    mat7->fYC.setup(young, nu, c, phi);
	mat7->SetMemory(nglobalpts, sz);
	mat7->UpdateBodyForce(bodyforce);

    mat8->fYC.setup(young, nu, c, phi);
	mat8->SetMemory(nglobalpts, sz);
	mat8->UpdateBodyForce(bodyforce);

    mat9->fYC.setup(young, nu, c, phi);
	mat9->SetMemory(nglobalpts, sz);
	mat9->UpdateBodyForce(bodyforce);

    mat10->fYC.setup(young, nu, c, phi);
	mat10->SetMemory(nglobalpts, sz);
	mat10->UpdateBodyForce(bodyforce);

        mat11->fYC.setup(young, nu, c, phi);
	mat11->SetMemory(nglobalpts, sz);
	mat11->UpdateBodyForce(bodyforce);
    
    Doub Lx = 40.;//(*Correlation length in x direction*)
	Doub Ly = 2.;//(*Correlation length in y direction*)
	Int nsamples = 50000, expansionorder = 150;
	Int type = 3;
	KLGalerkinRF* objKLGalerkinRF0 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF1 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF2 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF3 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF4 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF5 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF6 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF7 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF8 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF9 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF10 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
    KLGalerkinRF* objKLGalerkinRF11 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	objKLGalerkinRF0->SetMesh(mesh0);
    objKLGalerkinRF1->SetMesh(mesh1);
    objKLGalerkinRF2->SetMesh(mesh2);
    objKLGalerkinRF3->SetMesh(mesh3);
    objKLGalerkinRF4->SetMesh(mesh4);
    objKLGalerkinRF5->SetMesh(mesh5);
    objKLGalerkinRF6->SetMesh(mesh6);
    objKLGalerkinRF7->SetMesh(mesh7);
    objKLGalerkinRF8->SetMesh(mesh8);
    objKLGalerkinRF9->SetMesh(mesh9);
    objKLGalerkinRF10->SetMesh(mesh10);
    objKLGalerkinRF11->SetMesh(mesh11);

    if(false)
    {
        slopeproject* slopeobj0 = new slopeproject(mesh1, objKLGalerkinRF1);
        string filename = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2";
        slopeobj0->CreateRandomField(filename);

        return;
    }




    MatDoub coesionrandomfield, frictionrandomfield;
	string filerf = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2/coesionfield.txt";
	ReadMatDoub(coesionrandomfield, filerf);
	string filerff = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2/frictionfield.txt";
	ReadMatDoub(frictionrandomfield, filerff);

	NRmatrix<MatDoub> randomfield(2, 1);
	randomfield[0][0] = coesionrandomfield;
	randomfield[1][0] = frictionrandomfield;

    slopeproject* slopeobj0 = new slopeproject(mesh0, objKLGalerkinRF0,randomfield);
    slopeproject* slopeobj1 = new slopeproject(mesh1, objKLGalerkinRF1,randomfield);
    slopeproject* slopeobj2 = new slopeproject(mesh2, objKLGalerkinRF2,randomfield);
    slopeproject* slopeobj3 = new slopeproject(mesh3, objKLGalerkinRF3,randomfield);
    slopeproject* slopeobj4 = new slopeproject(mesh4, objKLGalerkinRF4,randomfield);
    slopeproject* slopeobj5 = new slopeproject(mesh5, objKLGalerkinRF5,randomfield);
    slopeproject* slopeobj6 = new slopeproject(mesh6, objKLGalerkinRF6,randomfield);
    slopeproject* slopeobj7 = new slopeproject(mesh7, objKLGalerkinRF7,randomfield);
    slopeproject* slopeobj8 = new slopeproject(mesh8, objKLGalerkinRF8,randomfield);
    slopeproject* slopeobj9 = new slopeproject(mesh9, objKLGalerkinRF9,randomfield);
    slopeproject* slopeobj10 = new slopeproject(mesh10, objKLGalerkinRF10,randomfield);
    slopeproject* slopeobj11 = new slopeproject(mesh11, objKLGalerkinRF11,randomfield);
    
    if(false)//Print a single simulation
    {
    
    MatDoub hhatinho2 = slopeobj0->AssembleHhationho(1156);//pior Lx =20, Ly =4
    //MatDoub hhatinho2 = slopeobj0->AssembleHhationho(1055);//melhor Lx =20, Ly =4
    mesh0->SetHhat(hhatinho2);
    
    string filename = "/home/diogo/projects/results/mathematicas/Coesao.dat";
    std::vector<std::vector<double>> hhatx;
    mesh0->fmaterial->PostProcess(mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 0, hhatinho2, hhatx);
    std::ofstream file(filename);
    slopeobj0->OutPutPost(hhatx, file);
    
    string filename2 = "/home/diogo/projects/results/mathematicas/Phi.dat";
    mesh0->fmaterial->PostProcess(mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 1, hhatinho2, hhatx);
    std::ofstream file2(filename);
    slopeobj0->OutPutPost(hhatx, file2);
    
    return;
    
    bool deterministicsol = true;
	if (deterministicsol == true)
	{
		int ndesirediters = 8, niter = 50;
		Doub dlamb0 = 0.2, alphatol = 0.0001;

		Doub tol = 0.001;
		std::vector<std::vector<double>> soll;
		//soll = slopeobj2->IterativeProcessShearRed(0.0001, 1., tol);
		mat0->fYC.setup(young, nu, c, phi);
		mat0->SetMemory(nglobalpts, sz);
		mat0->UpdateBodyForce(bodyforce);
        soll = slopeobj0->IterativeProcess(20, 0.1, 0.0001,10);
	}
	cout <<"\n initializing  ---->" << endl;
    return;
    }
    //int simtype,int comeco,int fim
    int GIMorSRM,begin,end;
    if ( simtype< 0)
    {
    std::cout << " \n Please input the simulation type: GIM = 0, SMR = 1: ";
    std::cin >> GIMorSRM;


    std::cout << " \n Please set the initial sample: ";
   std::cin >> begin;
    
    std::cout << " \n Please set the last sample: ";
    std::cin >> end;
    }else
    {
        std::cout << " \n simulation type: " << simtype <<endl;
        std::cout << " \n initial sample: "<< comeco <<endl;
        std::cout << " \n last sample: "<< fim <<endl;

        GIMorSRM=simtype;
        begin=comeco;
        end=fim;




    }
    if(GIMorSRM==0)
    {
        string namefolder = "/home/diogo/projects/results/gim-414els-Lx40-Ly2";
        int a=begin,b,c,d,e,f,g,h,i,j,l,m,n;
        int delta=int((end-begin)/10);
        b=a+delta;
        c=b+delta;
        d=c+delta;
        e=d+delta;
        f=e+delta;
        g=f+delta;
        h=g+delta;
        i=h+delta;
        j=i+delta;
        l=j+delta;
        m=l+delta;
        n=m+delta;
        a=250,b=1250,c=2250,d=3250,e=4250,f=5250,g=6250,h=7250,i=8250,j=9250,l=10000;

        //multthread
        std::thread thread0(myTreads,a,b, slopeobj0,namefolder);
        std::thread thread1(myTreads,b,c, slopeobj1,namefolder);
        std::thread thread2(myTreads,c,d, slopeobj2,namefolder);
        std::thread thread3(myTreads,d,e, slopeobj3,namefolder);
        std::thread thread4(myTreads,e,f, slopeobj4,namefolder);
        std::thread thread5(myTreads,f,g, slopeobj5,namefolder);
        std::thread thread6(myTreads,g,h, slopeobj6,namefolder);
        std::thread thread7(myTreads,h,i, slopeobj7,namefolder);
        std::thread thread8(myTreads,i,j, slopeobj8,namefolder);
        std::thread thread9(myTreads,j,l, slopeobj9,namefolder);
       // std::thread thread10(myTreads,l,m, slopeobj10,namefolder);
       // std::thread thread11(myTreads,m,n, slopeobj11,namefolder);

        thread0.join();
        thread1.join();
        thread2.join();
        thread3.join();
        thread4.join();
        thread5.join();
        thread6.join();
        thread7.join();
        thread8.join();
        thread9.join();
      //  thread10.join();
     //   thread11.join();
        //serial
        //slopeobj0->MonteCarloGIM(begin,end, false, namefolder);
       // slopeobj0->MonteCarloGIM(0,5, false, namefolder);



        return;
    }
    
    if(GIMorSRM==1)
    {
        string namefolder = "/home/diogo/projects/results/srm-414els-Lx40-Ly2";
        //slopeobj2->MonteCarloSRM(begin,end, false, namefolder3);
        int a=begin,b,c;
        int delta=int((end-begin)/2);
        b=a+delta;
        c=b+delta;
        cout << "a = "<< a <<endl;
        cout << "b = "<< b <<endl;
        cout << "c = "<< c <<endl;
        cout << "delta = "<< delta <<endl;
       // std::thread thread0(myTreadsSRM,a,b, slopeobj0,namefolder3);
      //  std::thread thread1(myTreadsSRM,b,c, slopeobj1,namefolder3);
      //  thread0.join();
      //  thread1.join();

        //serial
        slopeobj0->MonteCarloSRM(begin,end, false, namefolder);
        return;
    }

    
    std::cout << " \n wrong type input!";
    return;
     
}

void mainwindows()
{
	cout << "HIasd " << endl;

	string nodestr = "D:/DClibrary/meshes/nos-132-c3.txt";
	string elsstr = "D:/DClibrary/meshes/els-132-c3.txt";
    
    
    

	MatDoub hhatinho;
	cout << "HI " << endl;
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	cout << "HI " << endl;
	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);
	std::ofstream filemesh2("meshtopology.txt");
	OutPutPost(meshtopology, filemesh2);

	//Doub c = 18.5633, phi = 20 * M_PI / 180., gamma = -20.;//1.5
	Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5

	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	//Doub young = 100000.;
	//Doub nu = 0.3;
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

	elastoplastic2D< druckerprager >* mat2 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	elastoplastic2D< druckerprager >* mat3 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	elastoplastic2D< druckerprager >* mat4 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	elastoplastic2D< druckerprager >* mat5 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	elastoplastic2D< druckerprager >* mat6 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	elastoplastic2D< druckerprager >* mat7 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	elastoplastic2D< druckerprager >* mat8 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	elastoplastic2D< druckerprager >* mat9 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);

	mesh* mesh2 = new mesh(mat2, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh3 = new mesh(mat3, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh4 = new mesh(mat4, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh5 = new mesh(mat5, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh6 = new mesh(mat6, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh7 = new mesh(mat7, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh8 = new mesh(mat8, allcoords, meshcoords, meshtopology, hhatinho);
	mesh* mesh9 = new mesh(mat9, allcoords, meshcoords, meshtopology, hhatinho);

	int szdebug = mesh2->GetAllCoords().size();

	mat2->fYC.setup(young, nu, c, phi);
	mat2->SetMemory(nglobalpts, sz);
	mat2->UpdateBodyForce(bodyforce);


	mat3->fYC.setup(young, nu, c, phi);
	mat3->SetMemory(nglobalpts, sz);
	mat3->UpdateBodyForce(bodyforce);

	mat4->fYC.setup(young, nu, c, phi);
	mat4->SetMemory(nglobalpts, sz);
	mat4->UpdateBodyForce(bodyforce);

	mat5->fYC.setup(young, nu, c, phi);
	mat5->SetMemory(nglobalpts, sz);
	mat5->UpdateBodyForce(bodyforce);

	mat6->fYC.setup(young, nu, c, phi);
	mat6->SetMemory(nglobalpts, sz);
	mat6->UpdateBodyForce(bodyforce);

	mat7->fYC.setup(young, nu, c, phi);
	mat7->SetMemory(nglobalpts, sz);
	mat7->UpdateBodyForce(bodyforce);

	mat8->fYC.setup(young, nu, c, phi);
	mat8->SetMemory(nglobalpts, sz);
	mat8->UpdateBodyForce(bodyforce);

	mat9->fYC.setup(young, nu, c, phi);
	mat9->SetMemory(nglobalpts, sz);
	mat9->UpdateBodyForce(bodyforce);

	Doub Lx = 20.;//(*Correlation length in x direction*)
	Doub Ly = 2.;//(*Correlation length in y direction*)

	Int nsamples = 5000, expansionorder = 150;
	Int type = 3;
	KLGalerkinRF* objKLGalerkinRF2 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	KLGalerkinRF* objKLGalerkinRF3 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	KLGalerkinRF* objKLGalerkinRF4 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	KLGalerkinRF* objKLGalerkinRF5 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	KLGalerkinRF* objKLGalerkinRF6 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	KLGalerkinRF* objKLGalerkinRF7 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	KLGalerkinRF* objKLGalerkinRF8 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	KLGalerkinRF* objKLGalerkinRF9 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);

	objKLGalerkinRF2->SetMesh(mesh2);
	objKLGalerkinRF3->SetMesh(mesh3);
	objKLGalerkinRF4->SetMesh(mesh4);
	objKLGalerkinRF5->SetMesh(mesh5);
	objKLGalerkinRF6->SetMesh(mesh6);
	objKLGalerkinRF7->SetMesh(mesh7);
	objKLGalerkinRF8->SetMesh(mesh8);
	objKLGalerkinRF9->SetMesh(mesh9);
	slopeproject* slopeobj = new slopeproject(mesh2, objKLGalerkinRF2);
	//slopeproject* slopeobj3 = new slopeproject(mesh3, objKLGalerkinRF3);
	//	
	bool deterministicsol = false;
	if (deterministicsol == true)
	{
		int ndesirediters = 8, niter = 50;
		Doub dlamb0 = 0.2, alphatol = 0.0001;

		Doub tol = 0.001;
		std::vector<std::vector<double>> soll;
		soll = slopeobj->IterativeProcessShearRed(0.01, 2., tol);
		mat2->fYC.setup(young, nu, c, phi);
		mat2->SetMemory(nglobalpts, sz);
		mat2->UpdateBodyForce(bodyforce);

		//soll = slopeproject->IterativeProcess(10, 0.2, 0.001, 30);
	}

	MatDoub coesionrandomfield, frictionrandomfield;
	string filerf = "D:/slope-results/cho-field-Lx20-Ly4/coesionfield.txt";
	ReadMatDoub(coesionrandomfield, filerf);
	string filerff = "D:/slope-results/cho-field-Lx20-Ly4/frictionfield.txt";
	ReadMatDoub(frictionrandomfield, filerff);

	NRmatrix<MatDoub> randomfield(2, 1);
	randomfield[0][0] = coesionrandomfield;
	randomfield[1][0] = frictionrandomfield;

	slopeproject* slopeobj2 = new slopeproject(mesh2, objKLGalerkinRF2, randomfield);
	slopeproject* slopeobj3 = new slopeproject(mesh3, objKLGalerkinRF3, randomfield);
	slopeproject* slopeobj4 = new slopeproject(mesh4, objKLGalerkinRF4, randomfield);
	slopeproject* slopeobj5 = new slopeproject(mesh5, objKLGalerkinRF5, randomfield);
	slopeproject* slopeobj6 = new slopeproject(mesh6, objKLGalerkinRF6, randomfield);
	slopeproject* slopeobj7 = new slopeproject(mesh7, objKLGalerkinRF7, randomfield);
	slopeproject* slopeobj8 = new slopeproject(mesh8, objKLGalerkinRF8, randomfield);
	slopeproject* slopeobj9 = new slopeproject(mesh9, objKLGalerkinRF9, randomfield);

	string namefolder2 = "D:/slope-results/SRM-cho-field-Lx20-Ly4";

	string namefolder3 = "D:/slope-results/GI-cho-field-Lx20-Ly4";

	bool print = false;
	//slopeobj2->MonteCarloSRM(1900, 5000, print, namefolder2);

	// slopeobj2->MonteCarloGIM(1290, 5000, print, namefolder3);
    string number= "2";
	std::thread thread2(myTreads, 0,250, slopeobj2, number);
    number= "3";
	std::thread thread3(myTreads, 251, 500, slopeobj3, number);
	std::thread thread4(myTreads, 501, 750, slopeobj4, number);
	std::thread thread5(myTreads, 751,1000, slopeobj5, number);
	std::thread thread6(myTreads, 1001, 1250, slopeobj6, number);
	std::thread thread7(myTreads, 1251,1500, slopeobj7, number);
	std::thread thread8(myTreads, 1501, 1750, slopeobj8, number);
	std::thread thread9(myTreads, 1751, 2000, slopeobj9, number);

	thread2.join();
	thread3.join();
	thread4.join();
	thread5.join();
	thread6.join();
	thread7.join();
	thread8.join();
	thread9.join();
}

void createrandomfield()
{
        //cout<< "simtype"<< simtype;
	//string nodestr = "/home/diogo/projects/dcproj/nos-132-c3.txt";
	//string elsstr = "/home/diogo/projects/dcproj/els-132-c3.txt";

    //string nodestr = "/home/diogo/projects/dcproj/nos-cho.txt";
	//string elsstr = "/home/diogo/projects/dcproj/els-cho.txt";

    string nodestr = "/home/diogo/projects/dcproj/nos-414.txt";
	string elsstr = "/home/diogo/projects/dcproj/els-414.txt";

	MatDoub hhatinho;
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);
	std::ofstream filemesh2("meshtopology.txt");
	OutPutPost(meshtopology, filemesh2);

	Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5
	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
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

	elastoplastic2D< druckerprager >* mat0 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);

	mesh* mesh0 = new mesh(mat0, allcoords, meshcoords, meshtopology, hhatinho);

	mat0->fYC.setup(young, nu, c, phi);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);

    Doub Lx = 40.;//(*Correlation length in x direction*)
	Doub Ly = 2.;//(*Correlation length in y direction*)
	Int nsamples = 50000, expansionorder = 150;
	Int type = 3;

	KLGalerkinRF* objKLGalerkinRF0 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);

	objKLGalerkinRF0->SetMesh(mesh0);


    if(false)
    {
        slopeproject* slopeobj0 = new slopeproject(mesh0, objKLGalerkinRF0);
        string filename = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2";
        slopeobj0->CreateRandomField(filename);

        return;
    }
}

void mainlinuxserial(int simtype,int comeco,int fim)
{

    //cout<< "simtype"<< simtype;
	//string nodestr = "/home/diogo/projects/dcproj/nos-132-c3.txt";
	//string elsstr = "/home/diogo/projects/dcproj/els-132-c3.txt";

    //string nodestr = "/home/diogo/projects/dcproj/nos-cho.txt";
	//string elsstr = "/home/diogo/projects/dcproj/els-cho.txt";

    string nodestr = "/home/diogo/projects/dcproj/nos-414.txt";
	string elsstr = "/home/diogo/projects/dcproj/els-414.txt";

	MatDoub hhatinho;
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);
	std::ofstream filemesh2("meshtopology.txt");
	OutPutPost(meshtopology, filemesh2);

	Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5
	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
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

	elastoplastic2D< druckerprager >* mat0 = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);

	mesh* mesh0 = new mesh(mat0, allcoords, meshcoords, meshtopology, hhatinho);

	mat0->fYC.setup(young, nu, c, phi);
	mat0->SetMemory(nglobalpts, sz);
	mat0->UpdateBodyForce(bodyforce);

    Doub Lx = 40.;//(*Correlation length in x direction*)
	Doub Ly = 2.;//(*Correlation length in y direction*)
	Int nsamples = 50000, expansionorder = 150;
	Int type = 3;

	KLGalerkinRF* objKLGalerkinRF0 = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);

	objKLGalerkinRF0->SetMesh(mesh0);


    if(false)
    {
        slopeproject* slopeobj0 = new slopeproject(mesh0, objKLGalerkinRF0);
        string filename = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2";
        slopeobj0->CreateRandomField(filename);

        return;
    }

    MatDoub coesionrandomfield, frictionrandomfield;
	string filerf = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2/coesionfield.txt";
	ReadMatDoub(coesionrandomfield, filerf);
	string filerff = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2/frictionfield.txt";
	ReadMatDoub(frictionrandomfield, filerff);

	NRmatrix<MatDoub> randomfield(2, 1);
	randomfield[0][0] = coesionrandomfield;
	randomfield[1][0] = frictionrandomfield;

    slopeproject* slopeobj0 = new slopeproject(mesh0, objKLGalerkinRF0,randomfield);

    if(false)//Print a single simulation
    {
        MatDoub hhatinho2 = slopeobj0->AssembleHhationho(1156);//pior Lx =20, Ly =4
        //MatDoub hhatinho2 = slopeobj0->AssembleHhationho(1055);//melhor Lx =20, Ly =4
        mesh0->SetHhat(hhatinho2);

        string filename = "/home/diogo/projects/results/mathematicas/Coesao.dat";
        std::vector<std::vector<double>> hhatx;
        mesh0->fmaterial->PostProcess(mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 0, hhatinho2, hhatx);
        std::ofstream file(filename);
        slopeobj0->OutPutPost(hhatx, file);

        string filename2 = "/home/diogo/projects/results/mathematicas/Phi.dat";
        mesh0->fmaterial->PostProcess(mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 1, hhatinho2, hhatx);
        std::ofstream file2(filename);
        slopeobj0->OutPutPost(hhatx, file2);

        return;

        bool deterministicsol = true;
        if (deterministicsol == true)
        {
            int ndesirediters = 8, niter = 50;
            Doub dlamb0 = 0.2, alphatol = 0.0001;

            Doub tol = 0.001;
            std::vector<std::vector<double>> soll;
            //soll = slopeobj2->IterativeProcessShearRed(0.0001, 1., tol);
            mat0->fYC.setup(young, nu, c, phi);
            mat0->SetMemory(nglobalpts, sz);
            mat0->UpdateBodyForce(bodyforce);
            soll = slopeobj0->IterativeProcess(20, 0.1, 0.0001,10);
        }
        return;
    }

    int GIMorSRM,begin,end;
    if ( simtype< 0)
    {
        std::cout << " \n Please input the simulation type: GIM = 0, SMR = 1: ";
        std::cin >> GIMorSRM;
        std::cout << " \n Please set the initial sample: ";
        std::cin >> begin;
        std::cout << " \n Please set the last sample: ";
        std::cin >> end;
    }
    else
    {
        std::cout << " \n simulation type: " << simtype <<endl;
        std::cout << " \n initial sample: "<< comeco <<endl;
        std::cout << " \n last sample: "<< fim <<endl;

        GIMorSRM=simtype;
        begin=comeco;
        end=fim;
    }
    if(GIMorSRM==0)
    {
        string namefolder = "/home/diogo/projects/results/gim-414els-Lx40-Ly2";
        slopeobj0->MonteCarloGIM(begin,end, false, namefolder);
        return;
    }

    if(GIMorSRM==1)
    {
        string namefolder = "/home/diogo/projects/results/srm-414els-Lx40-Ly2";
        slopeobj0->MonteCarloSRM(begin,end, false, namefolder);
        return;
    }

    std::cout << " \n wrong type input!";
    return;

}

int main(int argc, char *argv[])
{

#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
    //leakcraw();
  //  cout <<argc<<"\n";

   // mainlinuxserial(-1,-1,-1);
    Eigen::initParallel();
    setNbThreads(5);
    if (argc > 3) {
        mainlinux(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

    }
    else
    {
        mainlinux(-1,-1,-1);
    }


#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
	mainwindows();
#endif

	cout << "Hello CMake." << endl;
	return 0;
}

template <class T>
void OutPutPost(NRmatrix<T>& postdata, std::ofstream& file)
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

//void ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord)
//{
//	std::vector<std::vector<Int>> topol;
//	string line, temp;
//
//	ifstream myfile(filenameel);
//	//
//	if (myfile.is_open())
//	{
//		while (getline(myfile, line))
//		{
//			std::vector<string> tokens;
//			istringstream iss(line);
//			while (iss >> temp)
//				tokens.push_back(temp);
//			std::vector<Int> input_int;
//			vecstr_to_vec(tokens, input_int);
//			for (int k = 0; k < input_int.size(); k++)
//			{
//				input_int[k] = input_int[k] - 1;
//			}
//			topol.push_back(input_int);
//		}
//		myfile.close();
//	}
//	else std::cout << "Unable to open file";
//
//	meshtopology.CopyFromVector(topol);
//	//meshtopology.Print();
//
//	std::vector<std::vector<Doub>> coords;
//	string line2, temp2;
//	ifstream myfile2(filenamecoord);
//	if (myfile2.is_open())
//	{
//		while (getline(myfile2, line2))
//		{
//			std::vector<string> tokens;
//			istringstream iss(line2);
//			while (iss >> temp2)
//				tokens.push_back(temp2);
//			std::vector<Doub> input_doub;
//			vecstr_to_vec(tokens, input_doub);
//
//			//std::vector<Doub> input_doub2(input_doub.size() - 1);
//			//for (int k = 1;k < input_doub.size();k++)
//			//{
//			//	input_doub2[k] = input_doub[k];
//			//}
//
//			coords.push_back(input_doub);
//		}
//		myfile2.close();
//	}
//	else std::cout << "Unable to open file";
//
//	meshcoords.CopyFromVector(coords);
//	//meshcoords.Print();
//
//
//	std::vector<Doub> temp33(2);
//	for (Int i = 0; i < meshtopology.nrows(); i++)
//	{
//		std::vector< std::vector<Doub> > temp22;
//		for (Int j = 0; j < meshtopology.ncols(); j++)
//		{
//			Int top = meshtopology[i][j];
//			temp33[0] = meshcoords[top][0];
//			temp33[1] = meshcoords[top][1];
//			temp22.push_back(temp33);
//		}
//		allcoords.push_back(temp22);
//	}
//
//
//
//}
//
//
//void ReadMatDoub(MatDoub& matdoub, std::string  file)
//{
//	std::vector<std::vector<Doub>> coords;
//	string line2, temp2;
//	ifstream myfile2(file);
//	if (myfile2.is_open())
//	{
//		while (getline(myfile2, line2))
//		{
//			std::vector<string> tokens;
//			istringstream iss(line2);
//			while (iss >> temp2)
//				tokens.push_back(temp2);
//			std::vector<Doub> input_doub;
//			vecstr_to_vec(tokens, input_doub);
//			coords.push_back(input_doub);
//		}
//		myfile2.close();
//	}
//	else std::cout << "Unable to open file";
//	//for (int i = 0;i < coords.size();i++)
//	//{
//	//	for (int j = 0;j < coords[0].size();j++)
//	//	{
//	//		cout << coords[i][j] << endl;
//	//	}
//	//	cout << endl;
//	//}
//
//	matdoub.CopyFromVector(coords);
//}
//
//template <class T> 
//void vecstr_to_vec(std::vector<std::string> vs, std::vector<T> &ret)
//{
//
//	for (std::vector<std::string>::iterator it = vs.begin(); it != vs.end(); ++it)
//	{
//		istringstream iss(*it);
//		T temp;
//		iss >> temp;
//		ret.push_back(temp);
//	}
//
//}



void ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord)
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
std::vector<Int> vecstr_to_vecint(std::vector<string> vs)
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

std::vector<Doub> vecstr_to_vecdoub(std::vector<string> vs)
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

std::vector<Doub> vecstr_to_vecdoub2(std::vector<string> vs)
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
std::vector<T> vecstr_to_vec(std::vector<string> vs)
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

void ReadMatDoub(MatDoub& matdoub, std::string  file)
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
