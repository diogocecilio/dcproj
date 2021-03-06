
// dcproj.cpp : Defines the entry point for the application.
//
#include "readgidmesh.h"
#include "dcproj.h"
#include "mesh.h"
#include "elastoplastic3D.h"
#include <thread>
#include "beam3dtools.h"
#include "mohrcoulomb.h"
#include "shapetri.h"
//#include <boost/thread.hpp>
//using namespace boost;
//#define EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL
//#define EIGEN_NO_DEBUG
//#define EIGEN_DONT_PARALLELIZE
//#include <mkl.h>
#define EIGEN_USE_MKL_ALL
//#include <lapacke.h>
#include "pressurizedhole.h"
#include "readgidmesh.h"
using namespace std;

void ReadMatDoub ( MatDoub& matdoub, std::string  file );
void ReadMesh ( std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord );
template <class T>
void OutPutPost ( NRmatrix<T>& postdata, std::ofstream& file );
template <class T>
void vecstr_to_vec ( std::vector<std::string> vs, std::vector<T>& ret );
std::vector<Int> vecstr_to_vecint ( std::vector<string> vs );
std::vector<Doub> vecstr_to_vecdoub ( std::vector<string> vs );
void myTreads ( int a, int b, slopeproject* slopeobj2, string traedN );
void leakcraw();

class SynchronizedFile
{
public:
    SynchronizedFile ( const string& path ) : _path ( path )
    {
        // Open file for writing...
    }

    void write ( const string& dataToWrite )
    {
        // Write to the file in a synchronized manner (described below)...
    }

private:
    string _path;
};

class Writer
{
public:
    Writer ( std::shared_ptr<SynchronizedFile> sf ) : _sf ( sf ) {}

    void someFunctionThatWritesToFile ()
    {
        // Do some work...
        _sf->write ( "Some data to write..." );
    }
private:
    std::shared_ptr<SynchronizedFile> _sf;
};

//<<<<<<< HEAD
void myTreads ( int a, int b, slopeproject* slopeobj2,string namefolder3 )
{
    bool print =  true;
    slopeobj2->MonteCarloGIM ( a, b, print, namefolder3 );

	
    delete slopeobj2;
}
void  PrintMathematicaFormat ( MatDoub postdata, std::ofstream& file )
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
void myTreadsSRM ( int a, int b, slopeproject* slopeobj2,string namefolder3 )
{
    slopeobj2->MonteCarloSRM ( a, b, false, namefolder3 );
    delete slopeobj2;
}

void solvequd();
int main ( int argc, char *argv[] )
{
	//string file ="/home/diogo/Desktop/mesh.msh";
	//readgidmesh read = readgidmesh(file);
	//read.ReadMesh();
	//return 0;

//     STRAT(1)xx  -9.2299616122518157E-004
//      STRAT(2)yy   3.6101399878573409E-005
//       STRAT(3)xy  -1.3917205303761122E-005
//          STRAT(4)zz   1.7273057237631743E-003  */



/*    Doub Phi=20*M_PI/180.;
    Doub Psi=20*M_PI/180.;
    Doub c=50.;
    Doub young=20000.;
    Doub nu=0.49;
    mohrcoulomb *mohr = new mohrcoulomb ( Phi,  Psi,  c, young,  nu );

    mohr->SetUp ( Phi,  Psi,  c, young,  nu );

    Doub exx=-2.9780912404287953E-003 ,exy=3.0354177217632269E-002,exz=0.,eyy=4.2960017010575877E-003 ,eyz=0.;
    Doub ezz= -1.9353008778866254E-017  */;
	
	Doub Phi=20*M_PI/180.;
    Doub Psi=20*M_PI/180.;
    Doub c=490.;
    Doub young=0.1e8;
    Doub nu=0.48;
    mohrcoulomb *mohr = new mohrcoulomb ( young,  nu, c ,Phi,Psi );

    mohr->SetUp ( young,  nu, c ,Phi,Psi );

    Doub exx=-1.2560424036224354E-003 ,exy=-3.0306513266847217E-004,exz=0.,eyy=7.5770078874204600E-004 ,eyz=0.;
    Doub ezz= 1.6664811388913960E-003  ;
	
        
    NRtensor<Doub> epst,epsp,projstress,projstrain,test;
    epst.XX() =exx;
    epst.XY() =exy;
    epst.XZ() =exz;
    epst.YY() =eyy;
    epst.YZ() =eyz;
    epst.ZZ() =ezz;

    NRmatrix<Doub>   Dep;
    Doub projgamma;
	
    mohr->closestpointproj( epst, epsp, projstress, projstrain, Dep,projgamma );

    //projstress.Print();

    Dep.Print();
	
    // cout<<"calculou!"<<endl;
    //return 0;
	//solvequd();
   slope2x1( );
	cout << "Hello CMake." << endl;
    return 0;
    // beam3dtools beam0bj = beam3dtools();
    // beam0bj.SolveElasticBeam();
    //beam0bj.IterativeProcess();

    //beam0bj.SolveElasticCube();
    //pressurizedhole hole0bj = pressurizedhole();
    //hole0bj.SolveElasticHole();
    //hole0bj.IterativeProcess();
    //  return 0;

#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
    //leakcraw();
    //  cout <<argc<<"\n";

    // mainlinuxserial(-1,-1,-1);
    Eigen::initParallel();
    setNbThreads ( 1 );
    if ( argc > 3 ) {
        //mainlinux(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
    } else {
        // mainlinux(-1,-1,-1);
        //mainlinux2(-1,-1,-1);
    }


#elif defined(_WIN32) || defined(WIN32)     /* _Win32 is usually defined by compilers targeting 32 or   64 bit Windows systems */
    mainwindows();
#endif

    cout << "Hello CMake." << endl;
    return 0;
}
 


void slope2x1( )
{


   //el padroa 9 sn 10
   //no padrao 4 sn 3
   
     MatDoub hhatinho;
    MatDoub  meshcoords, elcoords;
    MatInt meshtopology;
    std::vector<std::vector<std::vector<Doub>>> allcoords;
	//string file ="/home/diogo/projects/dcproj/data/mesh-el952-no2k.msh";//GIM 1.40 22.4s   com 7 pts
	//string file ="/home/diogo/projects/dcproj/data/meshtri-el1k-no3k.msh";//GIM 1.40 22.4s   com 7 pts
//	string file ="/home/diogo/projects/dcproj/data/meshtri-el1k-no2k.msh";//GIM 1.40 30s   com 7 pts
	//string file ="/home/diogo/projects/dcproj/data/mesh-el890-no1k.msh";//GIM 1.39 18s 20 s com 7 pts
	//string file ="/home/diogo/projects/dcproj/data/mesh-el816-no1k.msh";//GIM 1.39 18.34 com 7 pts
	//string file ="/home/diogo/projects/dcproj/data/mesh-el558-no1k.msh";//GIM 1.39 7.75 s 9.4s com 7 pts
	//string file ="/home/diogo/projects/dcproj/data/mesh-el1k-no650.msh";//GIM 1.54
	//string file ="/home/diogo/projects/dcproj/data/mesh-el558-no314.msh";//GIM 1.59 2.7s
	//string file ="/home/diogo/projects/dcproj/data/mesh-el6k-no3k.msh";
	//string file ="/home/diogo/projects/dcproj/data/mesh-el1k-no980.msh";//GIM 1.58
	//string file ="/home/diogo/projects/dcproj/data/mesh-el3k-no6k-p2.msh";//GIM 1.43
	//string file ="/home/diogo/projects/dcproj/data/mesh-el1k-no3k-quad.msh";
	//string file ="/home/diogo/projects/dcproj/data/mesh-el138-no467.msh";
	//string file ="/home/diogo/projects/dcproj/data/mesh-el551-no1k.msh";//GIM 1.44, 
	//string file ="/home/diogo/projects/dcproj/data/mesh-el1k-no3k.msh";//SRM 1.23
	//string file ="/home/diogo/projects/dcproj/data/mesh-el1k-no3k-b.msh";//GIM 1.44
	//string file ="/home/diogo/projects/dcproj/data/quad.msh";//SRM 1.23
	string file ="/home/diogo/projects/dcproj/data/mesh-el287-no922.msh";
	//string file ="/home/diogo/projects/dcproj/data/mesh-el606-no1k.msh";
	//string file ="/home/diogo/projects/dcproj/data/snmesh-3knodes-1kels.msh";
	//string file ="/home/diogo/projects/dcproj/data/snmesh-2knodes-725els.msh";
	//string file ="/home/diogo/projects/dcproj/data/snmesh-1knodes-360els.msh";
	//string file ="/home/diogo/projects/dcproj/data/mesh-el1k-no2k-renumber.msh";
	//string file ="/home/diogo/projects/dcproj/data/snmesh-tri-1knodes-900els.msh";
	//string file ="/home/diogo/projects/dcproj/data/snmesh-tri-2knodes-1kels-dir.msh";
	//string file ="/home/diogo/projects/dcproj/data/quad.msh";
	readgidmesh read = readgidmesh(file);
	read.ReadMesh();
	meshtopology = read.GetTopology();
	meshcoords = read.GetCoords();
	allcoords = read.GetAllCoords();
    //ReadMesh ( allcoords, meshcoords, meshtopology, elsstr, nodestr );


	Doub thickness = 1.;
	//malha 1x1
   // Doub c = 50., phi = 20. * M_PI / 180., gamma = -20.;
	
	
	 Doub c = 10., phi =30.* M_PI / 180., gamma = -20.;//1.5
	
    //Doub c = 23., phi =0.000000001* M_PI / 180., gamma = -20.;//1.5
    
    //Doub c = 5., phi =20.* M_PI / 180., gamma = -20.;//1.5
    
    //Doub c = 36.1664, phi =0.00001* M_PI / 180., gamma = -20.;//1.5, phi =20.* M_PI / 180., gamma = -20.;//1.5
	
	
	//para 45 graus LE preve FS =1.2 com NS=5.53
	//Doub c = 30., phi = 10. * M_PI / 180., gamma = -20.;
	
	//Doub c = 10., phi = 30. * M_PI / 180., gamma = -20.;
   
    Doub Phi=phi;
    Doub Psi=phi;


    
    Doub young = 20000.;
    Doub nu = 0.49;
    //Doub young = 100000.;
    //Doub nu = 0.3;
    Int planestress = 0;

    MatDoub bodyforce ( 2, 1, 0. ), newbodyforce;
    bodyforce[1][0] = gamma;
    MatDoub ptsweigths;
	int order;
	Int elnodes = meshtopology.ncols();
	cout << "elnodes"<<endl;
	cout << elnodes << endl;
	if(elnodes==3)
 	{
	  order=1;
	  shape* shapelocal = new shapetri(order,1);
	  shapelocal->pointsandweigths ( ptsweigths );
	}
	if(elnodes  ==4)
 	{
	  order=1;
	  shape* shapelocal = new shapequad(order,1);
	  shapelocal->pointsandweigths ( ptsweigths );
	}
	if(elnodes==6)
 	{
	  order=2;
	  shape* shapelocal = new shapetri(order,1);
	  shapelocal->pointsandweigths ( ptsweigths );
	}
	if(elnodes==8)
 	{
	  order=2;
	  shape* shapelocal = new shapequad(order,1);
	  shapelocal->pointsandweigths ( ptsweigths );
	}

	cout << "order  = "<< order << endl;
	
    Int npts = ptsweigths.nrows();
    Int nglobalpts = meshtopology.nrows() * npts;
    Int sz = 2 * meshcoords.nrows();
    int nthreads =15;
    std::vector <std::thread> threadsmat;

    Doub Lx = 20.;//(*Correlation length in x direction*)
    Doub Ly = 2.;//(*Correlation length in y direction*)
    Int nsamples = 2000, expansionorder = 300;
    Int type = 3;
    NRmatrix<MatDoub> randomfield ( 2, 1 );
    Int dim =2;

    elastoplastic2D< druckerprager >* mat = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, elnodes, hhatinho );
    mat->fYC.setup ( young, nu, c, phi );
    mat->SetMemory ( nglobalpts, sz );
    mat->UpdateBodyForce ( bodyforce );

    elastoplastic2D< mohrcoulomb >* matmohr = new elastoplastic2D< mohrcoulomb > ( thickness, bodyforce, planestress, elnodes, hhatinho );
    matmohr->fYC.SetUp (  young, nu,c,Phi,Psi );
    matmohr->SetMemory ( nglobalpts, sz );
    matmohr->UpdateBodyForce ( bodyforce );
	
// 	bodyforce.assign(2,1,0.);
// 	thickness=1.;
// 	planestress=1;
// 	young = 10000.;
// 	nu = 0.3;
 	Doub sigy=10000.;
	elastoplastic2D< vonmises >* matvon = new elastoplastic2D< vonmises > ( thickness, bodyforce, planestress, elnodes, hhatinho );
	matvon->fYC.setup( young,  nu ,sigy);
    matvon->SetMemory ( nglobalpts, sz );
    matvon->UpdateBodyForce ( bodyforce );

 //	mesh* meshs = new mesh ( dim,matmohr, allcoords, meshcoords, meshtopology, hhatinho );
     mesh* meshs = new mesh ( dim,mat , allcoords, meshcoords, meshtopology, hhatinho );

//	mesh* meshs = new mesh(dim,mat, allcoords, meshcoords, meshtopology, hhatinho);
// 	NRmatrix<Doub> KG,FINT,FBODY;
// 	meshs->Assemble ( KG, FINT, FBODY );
// 	KG.Print();
// 	return ;

    KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF ( elnodes, Lx, Ly, type, nsamples, expansionorder );
    objKLGalerkinRF->SetMesh ( meshs );

    slopeproject* slopeobj = new slopeproject ( meshs, objKLGalerkinRF,randomfield );

    //SRM nao bate com artigos para phi = 0, verificar norma do residuo. Para esta exemple convergencia estabiliza em 20. Proble instável.
    //DETERMINISTC SOL. Phi=0 aproxima tresca?
	
    //SRM  com tangente numerica converge no MC e com tangente analitica nao.
	
	if (false ) {

        //last filed created c =23 e phi=0
        string filename = "mesh-secondexample";
        filename+="/field-Lx";
        filename+=to_string ( Int ( Lx ) );
        filename+="-Ly";
        filename+=to_string ( Int ( Ly ) );
        slopeobj->CreateRandomField ( filename );
		
		return;
		
    }




//     randomfield[0][0] = coesionrandomfield;
//     randomfield[1][0] = frictionrandomfield;
// 	slopeproject* slopeobj0 = new slopeproject ( meshs, objKLGalerkinRF,randomfield );
// 	MatDoub hhatinho2 = slopeobj0->AssembleHhationho ( 1003 ); //pior Lx =20, Ly =4
// 	//MatDoub hhatinho2 = slopeobj0->AssembleHhationho(1055);//melhor Lx =20, Ly =4
// 	meshs->SetHhat ( hhatinho2 );

    if ( false ) {
        cout <<"\n  DETERMINISTC  " << endl;
		

		
		
        slopeproject* slopeobj = new slopeproject ( meshs, objKLGalerkinRF );

   
        std::vector<std::vector<double>> soll;
        mat->fYC.setup ( young, nu, c, phi );
        mat->SetMemory ( nglobalpts, sz );
        mat->UpdateBodyForce ( bodyforce );

        matmohr->fYC.SetUp (young,  nu,c,Phi,Psi );
        matmohr->SetMemory ( nglobalpts, sz );
        matmohr->UpdateBodyForce ( bodyforce );

		matvon->fYC.setup( young,  nu ,c);
		matvon->SetMemory ( nglobalpts, sz );
		matvon->UpdateBodyForce ( bodyforce );
		//24 s
		//int desirediter = 10;
       // Doub dlamb0 =0.25;
		//Doub maxlfac=2;
		//28.4 s
       // int desirediter = 10;
      //  Doub dlamb0 =0.25;
		//Doub maxlfac=1.5;
		//22.72 s
		//int desirediter = 7;
        //Doub dlamb0 =0.25;
		//Doub maxlfac=2;

		//18.48 s
		//int desirediter = 9;
      //  Doub dlamb0 =0.6;
	//	Doub maxlfac=0.9;
		
		int desirediter = 10;
        Doub dlamb0 =0.6;
		Doub maxlfac=0.9;
		
        slopeobj->IterativeProcessNew ( desirediter, dlamb0,maxlfac,0);
		//slopeobj->PostVtk(0);
       // slopeobj->IterativeProcess2( );
        //soll = slopeobj->IterativeProcessShearRed( 0.1, 2.,0.01);
		//cout <<"\n  exit  " << endl;
      //  soll = slopeobj->IterativeProcessGIMBinarySearch();
        return;

    }


    std::cout << "lendo fields = "<<std::endl;
    MatDoub coesionrandomfield, frictionrandomfield;
    string filerf = "mesh-secondexample/field-Lx20-Ly2/coesionfield.txt";
    ReadMatDoub ( coesionrandomfield, filerf );
    string filerff = "mesh-secondexample/field-Lx20-Ly2/frictionfield.txt";
    ReadMatDoub ( frictionrandomfield, filerff );

	randomfield[0][0] = coesionrandomfield;
    randomfield[1][0] = frictionrandomfield;
	//slopeproject* slopeobj0 = new slopeproject ( meshs, objKLGalerkinRF,randomfield );
	
    string namefolder = "mesh-secondexample/gim-Lx20-Ly2-dp";
    //int delta = int ( nsamples/nthreads );
	int delta = 100;
    int a=0,b=delta;


    for ( int i=0; i<nthreads; i++ ) {
//	for ( int i=0; i<1; i++ ) {
      //    elastoplastic2D< druckerprager >* mat = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
        std::cout << "criando malha "<<std::endl;

        //elastoplastic2D< mohrcoulomb >* matmohr = new elastoplastic2D< mohrcoulomb > ( thickness, bodyforce, planestress, elnodes, hhatinho );
		elastoplastic2D< druckerprager >* mat = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, elnodes, hhatinho );
        mesh* meshs = new mesh ( dim,mat, allcoords, meshcoords, meshtopology, hhatinho );
		
       // matmohr->fYC.SetUp ( young, nu,c,Phi,Psi );
      // matmohr->SetMemory ( nglobalpts, sz );
      //  matmohr->UpdateBodyForce ( bodyforce );
		
		
		mat->fYC.setup ( young, nu, c, phi );
        mat->SetMemory ( nglobalpts, sz );
       mat->UpdateBodyForce ( bodyforce );
		
        KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF ( elnodes, Lx, Ly, type, nsamples, expansionorder );
        objKLGalerkinRF->SetMesh ( meshs );
        std::cout << "criando slopeproject "<<std::endl;
        slopeproject* slopeobj = new slopeproject ( meshs, objKLGalerkinRF,randomfield );
        std::cout << "a = "<< a <<std::endl;
        std::cout << "b = "<< b <<std::endl;
        std::thread threadx ( myTreads,a,b, slopeobj,namefolder );
		//slopeobj->MonteCarloSRM ( a, b, false, namefolder );
    //   std::thread threadx(myTreadsSRM,a,b, slopeobj,namefolder);
        a=b+1;
        b+=delta;

        threadsmat.push_back ( std::move ( threadx ) );
    }

    for ( auto &threadx: threadsmat ) threadx.join();



}


void mainlinux ( int simtype,int comeco,int fim )
{




    //cout<< "simtype"<< simtype;
    //string nodestr = "/home/diogo/projects/dcproj/nos-132-c3.txt";
    //string elsstr = "/home/diogo/projects/dcproj/els-132-c3.txt";

    // string nodestr = "/home/diogo/projects/dcproj/nos-cho.txt";
    //string elsstr = "/home/diogo/projects/dcproj/els-cho.txt";

    //string nodestr = "/home/diogo/projects/dcproj/nos-912.txt";
    //string elsstr = "/home/diogo/projects/dcproj/els-912.txt";

    //   string nodestr = "/home/diogo/projects/dcproj/nodes.dat";//mathematica
    //string elsstr = "/home/diogo/projects/dcproj/topol.dat";//mathematica

    //  string nodestr = "/home/diogocecilio/projects/dcproj/nodes-unstructured.txt";
    //string elsstr = "/home/diogocecilio/projects/dcproj/els-unstructured.txt";

//   string nodestr = "/home/diogocecilio/projects/dcproj/nodes1k.txt";
//	string elsstr = "/home/diogocecilio/projects/dcproj/els1k.txt";
//

    //string nodestr = "/home/diogocecilio/projects/dcproj/nodes-2k.txt";
    //string elsstr = "/home/diogocecilio/projects/dcproj/els-2k.txt";

    //string nodestr = "/home/diogocecilio/projects/dcproj/data/nos-287.txt";
    //string elsstr = "/home/diogocecilio/projects/dcproj/data/els-287.txt";

    //string nodestr = "/home/diogocecilio/projects/dcproj/nos-381.txt";
    //string elsstr = "/home/diogocecilio/projects/dcproj/els-381.txt";
	
	   string nodestr = "/home/diogocecilio/projects/dcproj/data/triumdoido-nodes.txt";
    string elsstr = "/home/diogocecilio/projects/dcproj/data/triumdoido-els.txt";

    // string nodestr = "/home/diogocecilio/projects/dcproj/nodes-606.txt";
    //string elsstr = "/home/diogocecilio/projects/dcproj/els-606.txt";

    //string nodestr = "/home/diogo/projects/dcproj/nos-445.txt";
    //string elsstr = "/home/diogo/projects/dcproj/els-445.txt";


    int dim=2;
    MatDoub hhatinho;
    MatDoub  meshcoords, elcoords;
    MatInt meshtopology;
    std::vector<std::vector<std::vector<Doub>>> allcoords;
    ReadMesh ( allcoords, meshcoords, meshtopology, elsstr, nodestr );

    std::ofstream filemesh1 ( "meshcoords.txt" );
    OutPutPost ( meshcoords, filemesh1 );
    std::ofstream filemesh2 ( "meshtopology.txt" );
    OutPutPost ( meshtopology, filemesh2 );

    //Doub c = 18.5633, phi = 20 * M_PI / 180., gamma = -20.;//1.5
    Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5

    Doub thickness = 1.;
    Doub young = 20000.;
    Doub nu = 0.49;
    //Doub young = 100000.;
    //Doub nu = 0.3;
    Int planestress = 0;

    MatDoub bodyforce ( 2, 1, 0. ), newbodyforce;
    bodyforce[1][0] = gamma;
    MatDoub ptsweigths;
    int order = 2;
    shapequad shape = shapequad ( order, 1 );
    shape.pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();
    Int nglobalpts = meshtopology.nrows() * npts;
    Int sz = 2 * meshcoords.nrows();

    Doub Lx = 20.;//(*Correlation length in x direction*)
    Doub Ly = 2.;//(*Correlation length in y direction*)
    Int nsamples = 50000, expansionorder = 150;
    Int type = 3;

    elastoplastic2D< druckerprager >* mat0 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat1 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat2 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat3 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat4 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat5 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat6 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat7 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat8 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat9 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat10 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat11 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat12= new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat13 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat14 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat15 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat16 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat17 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat18 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );

    mesh* mesh0 = new mesh ( dim,mat0, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh1 = new mesh ( dim,mat1, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh2 = new mesh ( dim,mat2, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh3 = new mesh ( dim,mat3, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh4 = new mesh ( dim,mat4, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh5 = new mesh ( dim,mat5, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh6 = new mesh ( dim,mat6, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh7 = new mesh ( dim,mat7, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh8 = new mesh ( dim,mat8, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh9 = new mesh ( dim,mat9, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh10 = new mesh ( dim,mat10, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh11 = new mesh ( dim,mat11, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh12 = new mesh ( dim,mat12, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh13 = new mesh ( dim,mat13, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh14 = new mesh ( dim,mat14, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh15 = new mesh ( dim,mat15, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh16 = new mesh ( dim,mat16, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh17 = new mesh ( dim,mat17, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh18 = new mesh ( dim,mat18, allcoords, meshcoords, meshtopology, hhatinho );


    mat0->fYC.setup ( young, nu, c, phi );
    mat0->SetMemory ( nglobalpts, sz );
    mat0->UpdateBodyForce ( bodyforce );

    mat1->fYC.setup ( young, nu, c, phi );
    mat1->SetMemory ( nglobalpts, sz );
    mat1->UpdateBodyForce ( bodyforce );

    mat2->fYC.setup ( young, nu, c, phi );
    mat2->SetMemory ( nglobalpts, sz );
    mat2->UpdateBodyForce ( bodyforce );

    mat3->fYC.setup ( young, nu, c, phi );
    mat3->SetMemory ( nglobalpts, sz );
    mat3->UpdateBodyForce ( bodyforce );

    mat4->fYC.setup ( young, nu, c, phi );
    mat4->SetMemory ( nglobalpts, sz );
    mat4->UpdateBodyForce ( bodyforce );

    mat5->fYC.setup ( young, nu, c, phi );
    mat5->SetMemory ( nglobalpts, sz );
    mat5->UpdateBodyForce ( bodyforce );

    mat6->fYC.setup ( young, nu, c, phi );
    mat6->SetMemory ( nglobalpts, sz );
    mat6->UpdateBodyForce ( bodyforce );

    mat7->fYC.setup ( young, nu, c, phi );
    mat7->SetMemory ( nglobalpts, sz );
    mat7->UpdateBodyForce ( bodyforce );

    mat8->fYC.setup ( young, nu, c, phi );
    mat8->SetMemory ( nglobalpts, sz );
    mat8->UpdateBodyForce ( bodyforce );

    mat9->fYC.setup ( young, nu, c, phi );
    mat9->SetMemory ( nglobalpts, sz );
    mat9->UpdateBodyForce ( bodyforce );

    mat10->fYC.setup ( young, nu, c, phi );
    mat10->SetMemory ( nglobalpts, sz );
    mat10->UpdateBodyForce ( bodyforce );

    mat11->fYC.setup ( young, nu, c, phi );
    mat11->SetMemory ( nglobalpts, sz );
    mat11->UpdateBodyForce ( bodyforce );

    mat12->fYC.setup ( young, nu, c, phi );
    mat12->SetMemory ( nglobalpts, sz );
    mat12->UpdateBodyForce ( bodyforce );

    mat13->fYC.setup ( young, nu, c, phi );
    mat13->SetMemory ( nglobalpts, sz );
    mat13->UpdateBodyForce ( bodyforce );

    mat14->fYC.setup ( young, nu, c, phi );
    mat14->SetMemory ( nglobalpts, sz );
    mat14->UpdateBodyForce ( bodyforce );

    mat15->fYC.setup ( young, nu, c, phi );
    mat15->SetMemory ( nglobalpts, sz );
    mat15->UpdateBodyForce ( bodyforce );

    mat16->fYC.setup ( young, nu, c, phi );
    mat16->SetMemory ( nglobalpts, sz );
    mat16->UpdateBodyForce ( bodyforce );

    mat17->fYC.setup ( young, nu, c, phi );
    mat17->SetMemory ( nglobalpts, sz );
    mat17->UpdateBodyForce ( bodyforce );

    mat18->fYC.setup ( young, nu, c, phi );
    mat18->SetMemory ( nglobalpts, sz );
    mat18->UpdateBodyForce ( bodyforce );

    //Doub Lx = 20.;//(*Correlation length in x direction*)
    //Doub Ly = 2.;//(*Correlation length in y direction*)
    //Int nsamples = 50000, expansionorder = 150;
    //Int type = 3;
    KLGalerkinRF* objKLGalerkinRF0 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF1 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF2 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF3 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF4 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF5 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF6 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF7 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF8 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF9 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF10 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF11 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF12 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF13 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF14 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF15 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF16 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF17 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF18 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );

    objKLGalerkinRF0->SetMesh ( mesh0 );
    objKLGalerkinRF1->SetMesh ( mesh1 );
    objKLGalerkinRF2->SetMesh ( mesh2 );
    objKLGalerkinRF3->SetMesh ( mesh3 );
    objKLGalerkinRF4->SetMesh ( mesh4 );
    objKLGalerkinRF5->SetMesh ( mesh5 );
    objKLGalerkinRF6->SetMesh ( mesh6 );
    objKLGalerkinRF7->SetMesh ( mesh7 );
    objKLGalerkinRF8->SetMesh ( mesh8 );
    objKLGalerkinRF9->SetMesh ( mesh9 );
    objKLGalerkinRF10->SetMesh ( mesh10 );
    objKLGalerkinRF11->SetMesh ( mesh11 );
    objKLGalerkinRF12->SetMesh ( mesh12 );
    objKLGalerkinRF13->SetMesh ( mesh13 );
    objKLGalerkinRF14->SetMesh ( mesh14 );
    objKLGalerkinRF15->SetMesh ( mesh15 );
    objKLGalerkinRF16->SetMesh ( mesh16 );
    objKLGalerkinRF17->SetMesh ( mesh17 );
    objKLGalerkinRF18->SetMesh ( mesh18 );

    string mathematicapath = "/home/diogocecilio/projects/results/mathematica-287";
    string resultfolderpath = "/home/diogocecilio/Dropbox/slope-reliability/results/mesh-287-novostestes";

    if ( true ) {
        cout <<"\n initializing hhhh ---->" << endl;
        slopeproject* slopeobj0 = new slopeproject ( mesh0, objKLGalerkinRF0 );
        int ndesirediters = 8, niter = 50;
        Doub dlamb0 = 0.2, alphatol = 0.0001;
        Doub tol = 0.001;
        std::vector<std::vector<double>> soll;
        mat0->fYC.setup ( young, nu, c, phi );
        mat0->SetMemory ( nglobalpts, sz );
        mat0->UpdateBodyForce ( bodyforce );
        soll = slopeobj0->IterativeProcess ( 10, 0.5, 0.01,30 );
        //soll = slopeobj0->IterativeProcessShearRed(0.1,2.,tol);
        //'return;
        string filename = resultfolderpath;
        filename+="/field-Lx";
        filename+=to_string ( Int ( Lx ) );
        filename+="-Ly";
        filename+=to_string ( Int ( Ly ) );
        slopeobj0->CreateRandomField ( filename );

        return;
    }




    MatDoub coesionrandomfield, frictionrandomfield;
    string filerf=resultfolderpath;
    filerf+="/field-Lx";
    filerf+=to_string ( Int ( Lx ) );
    filerf+="-Ly";
    filerf+=to_string ( Int ( Ly ) );
    filerf += "/coesionfield.txt";
    ReadMatDoub ( coesionrandomfield, filerf );
    string filerff =resultfolderpath;
    filerff+="/field-Lx";
    filerff+=to_string ( Int ( Lx ) );
    filerff+="-Ly";
    filerff+=to_string ( Int ( Ly ) );
    filerff+= "/frictionfield.txt";
    ReadMatDoub ( frictionrandomfield, filerff );

    NRmatrix<MatDoub> randomfield ( 2, 1 );
    randomfield[0][0] = coesionrandomfield;
    randomfield[1][0] = frictionrandomfield;

    slopeproject* slopeobj0 = new slopeproject ( mesh0, objKLGalerkinRF0,randomfield );
    slopeproject* slopeobj1 = new slopeproject ( mesh1, objKLGalerkinRF1,randomfield );
    slopeproject* slopeobj2 = new slopeproject ( mesh2, objKLGalerkinRF2,randomfield );
    slopeproject* slopeobj3 = new slopeproject ( mesh3, objKLGalerkinRF3,randomfield );
    slopeproject* slopeobj4 = new slopeproject ( mesh4, objKLGalerkinRF4,randomfield );
    slopeproject* slopeobj5 = new slopeproject ( mesh5, objKLGalerkinRF5,randomfield );
    slopeproject* slopeobj6 = new slopeproject ( mesh6, objKLGalerkinRF6,randomfield );
    slopeproject* slopeobj7 = new slopeproject ( mesh7, objKLGalerkinRF7,randomfield );
    slopeproject* slopeobj8 = new slopeproject ( mesh8, objKLGalerkinRF8,randomfield );
    slopeproject* slopeobj9 = new slopeproject ( mesh9, objKLGalerkinRF9,randomfield );
    slopeproject* slopeobj10 = new slopeproject ( mesh10, objKLGalerkinRF10,randomfield );
    slopeproject* slopeobj11 = new slopeproject ( mesh11, objKLGalerkinRF11,randomfield );

    slopeproject* slopeobj12 = new slopeproject ( mesh12, objKLGalerkinRF12,randomfield );
    slopeproject* slopeobj13 = new slopeproject ( mesh13, objKLGalerkinRF13,randomfield );
    slopeproject* slopeobj14 = new slopeproject ( mesh14, objKLGalerkinRF14,randomfield );
    slopeproject* slopeobj15 = new slopeproject ( mesh15, objKLGalerkinRF15,randomfield );
    slopeproject* slopeobj16 = new slopeproject ( mesh16, objKLGalerkinRF16,randomfield );
    slopeproject* slopeobj17 = new slopeproject ( mesh17, objKLGalerkinRF17,randomfield );
    slopeproject* slopeobj18 = new slopeproject ( mesh18, objKLGalerkinRF18,randomfield );
    MatDoub hhatinho2 = slopeobj0->AssembleHhationho ( 1156 );
    if ( false ) { //Print a single simulation

        MatDoub hhatinho2 = slopeobj0->AssembleHhationho ( 1156 ); //pior Lx =20, Ly =4
        //MatDoub hhatinho2 = slopeobj0->AssembleHhationho(1055);//melhor Lx =20, Ly =4
        mesh0->SetHhat ( hhatinho2 );
        string filename = mathematicapath;
        filename += "/Coesao.dat";
        std::vector<std::vector<double>> hhatx;
        mesh0->fmaterial->PostProcess ( mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 0, hhatinho2, hhatx );
        std::ofstream file ( filename );
        slopeobj0->OutPutPost ( hhatx, file );

        string filename2 = mathematicapath;
        filename2+="/Phi.dat";
        mesh0->fmaterial->PostProcess ( mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 1, hhatinho2, hhatx );
        std::ofstream file2 ( filename );
        slopeobj0->OutPutPost ( hhatx, file2 );



        bool deterministicsol = true;
        if ( deterministicsol == true ) {
            int ndesirediters = 8, niter = 50;
            Doub dlamb0 = 0.2, alphatol = 0.0001;

            Doub tol = 0.001;
            std::vector<std::vector<double>> soll;
            //soll = slopeobj2->IterativeProcessShearRed(0.0001, 1., tol);
            mat0->fYC.setup ( young, nu, c, phi );
            mat0->SetMemory ( nglobalpts, sz );
            mat0->UpdateBodyForce ( bodyforce );
            soll = slopeobj0->IterativeProcess ( 30, 0.1, 0.01,20 );
        }

        return;
    }
    //int simtype,int comeco,int fim
    int GIMorSRM,begin,end;
    /*   if ( simtype< 0)
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




       }*/
    GIMorSRM=0;
    begin=0;
    end=50000;
    if ( GIMorSRM==0 ) {
        /*  std::cout << " \n GIM ";
          string namefolder = "/home/diogo/Dropbox/slope-reliability/results/mesh-287/gim-Lx20-Ly2";
          int a=begin,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,u,v;
          int delta=int((end-begin)/19);
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
          o=n+delta;
          p=o+delta;
          q=p+delta;
          r=q+delta;
          s=r+delta;
          t=s+delta;
          u=t+delta;
          std::cout << " \n a ";
          std::cout << " \n u ";
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
         std::thread thread10(myTreads,l,m, slopeobj10,namefolder);
          std::thread thread11(myTreads,m,n, slopeobj11,namefolder);
          std::thread thread12(myTreads,n,o, slopeobj12,namefolder);
          std::thread thread13(myTreads,o,p, slopeobj13,namefolder);
          std::thread thread14(myTreads,p,q, slopeobj14,namefolder);
          std::thread thread15(myTreads,q,r, slopeobj15,namefolder);
          std::thread thread16(myTreads,r,s, slopeobj16,namefolder);
          std::thread thread17(myTreads,s,t, slopeobj17,namefolder);
          std::thread thread18(myTreads,t,u, slopeobj18,namefolder);

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
          thread10.join();
          thread11.join();
          thread12.join();
          thread13.join();
          thread14.join();
          thread15.join();
          thread16.join();
          thread17.join();
          thread18.join();
        */
        //serial
        //slopeobj0->MonteCarloGIM(begin,end, false, namefolder);
        std::cout << " \n teste vtk ";
        auto namefolder=resultfolderpath;
        namefolder += "/gim-Lx";
        namefolder+=to_string ( Int ( Lx ) );
        namefolder+="-Ly";
        namefolder+=to_string ( Int ( Ly ) );
        slopeobj0->MonteCarloGIM ( 0,1, false, namefolder );
        return;
    }

    if ( GIMorSRM==1 ) {
        std::cout << " \n GIM ";
        string namefolder = "/home/diogo/Dropbox/slope-reliability/results/mesh-287/srm-Lx20-Ly2";
        int a=begin,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,u,v;
        int delta=int ( ( end-begin ) /19 );
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
        o=n+delta;
        p=o+delta;
        q=p+delta;
        r=q+delta;
        s=r+delta;
        t=s+delta;
        u=t+delta;
        std::cout << " \n a " << a <<std::endl;
        std::cout << " \n u "<< u <<std::endl;
        //multthread
        std::thread thread0 ( myTreadsSRM,a,b, slopeobj0,namefolder );
        std::thread thread1 ( myTreadsSRM,b,c, slopeobj1,namefolder );
        std::thread thread2 ( myTreadsSRM,c,d, slopeobj2,namefolder );
        std::thread thread3 ( myTreadsSRM,d,e, slopeobj3,namefolder );
        std::thread thread4 ( myTreadsSRM,e,f, slopeobj4,namefolder );
        std::thread thread5 ( myTreadsSRM,f,g, slopeobj5,namefolder );
        std::thread thread6 ( myTreadsSRM,g,h, slopeobj6,namefolder );
        std::thread thread7 ( myTreadsSRM,h,i, slopeobj7,namefolder );
        std::thread thread8 ( myTreadsSRM,i,j, slopeobj8,namefolder );
        std::thread thread9 ( myTreadsSRM,j,l, slopeobj9,namefolder );
        std::thread thread10 ( myTreadsSRM,l,m, slopeobj10,namefolder );
        std::thread thread11 ( myTreadsSRM,m,n, slopeobj11,namefolder );
        std::thread thread12 ( myTreadsSRM,n,o, slopeobj12,namefolder );
        std::thread thread13 ( myTreadsSRM,o,p, slopeobj13,namefolder );
        std::thread thread14 ( myTreadsSRM,p,q, slopeobj14,namefolder );
        std::thread thread15 ( myTreadsSRM,q,r, slopeobj15,namefolder );
        std::thread thread16 ( myTreadsSRM,r,s, slopeobj16,namefolder );
        std::thread thread17 ( myTreadsSRM,s,t, slopeobj17,namefolder );
        std::thread thread18 ( myTreadsSRM,t,u, slopeobj18,namefolder );

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
        thread10.join();
        thread11.join();
        thread12.join();
        thread13.join();
        thread14.join();
        thread15.join();
        thread16.join();
        thread17.join();
        thread18.join();
    }


    std::cout << " \n wrong type input!";
    return;

}

void mainwindows()
{
    cout << "HIasd " << endl;

    string nodestr = "D:/DClibrary/meshes/nos-132-c3.txt";
    string elsstr = "D:/DClibrary/meshes/els-132-c3.txt";


    int dim =2;

    MatDoub hhatinho;
    cout << "HI " << endl;
    MatDoub  meshcoords, elcoords;
    MatInt meshtopology;
    std::vector<std::vector<std::vector<Doub>>> allcoords;
    ReadMesh ( allcoords, meshcoords, meshtopology, elsstr, nodestr );

    cout << "HI " << endl;
    std::ofstream filemesh1 ( "meshcoords.txt" );
    OutPutPost ( meshcoords, filemesh1 );
    std::ofstream filemesh2 ( "meshtopology.txt" );
    OutPutPost ( meshtopology, filemesh2 );

    //Doub c = 18.5633, phi = 20 * M_PI / 180., gamma = -20.;//1.5
    Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5

    Doub thickness = 1.;
    Doub young = 20000.;
    Doub nu = 0.49;
    //Doub young = 100000.;
    //Doub nu = 0.3;
    Int planestress = 0;

    MatDoub bodyforce ( 2, 1, 0. ), newbodyforce;
    bodyforce[1][0] = gamma;
    MatDoub ptsweigths;
    int order = 2;
    shapequad shape = shapequad ( order, 1 );
    shape.pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();
    Int nglobalpts = meshtopology.nrows() * npts;
    Int sz = 2 * meshcoords.nrows();

    elastoplastic2D< druckerprager >* mat2 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat3 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat4 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat5 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat6 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat7 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat8 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );
    elastoplastic2D< druckerprager >* mat9 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );

    mesh* mesh2 = new mesh ( dim,mat2, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh3 = new mesh ( dim,mat3, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh4 = new mesh ( dim,mat4, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh5 = new mesh ( dim,mat5, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh6 = new mesh ( dim,mat6, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh7 = new mesh ( dim,mat7, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh8 = new mesh ( dim,mat8, allcoords, meshcoords, meshtopology, hhatinho );
    mesh* mesh9 = new mesh ( dim,mat9, allcoords, meshcoords, meshtopology, hhatinho );

    int szdebug = mesh2->GetAllCoords().size();

    mat2->fYC.setup ( young, nu, c, phi );
    mat2->SetMemory ( nglobalpts, sz );
    mat2->UpdateBodyForce ( bodyforce );


    mat3->fYC.setup ( young, nu, c, phi );
    mat3->SetMemory ( nglobalpts, sz );
    mat3->UpdateBodyForce ( bodyforce );

    mat4->fYC.setup ( young, nu, c, phi );
    mat4->SetMemory ( nglobalpts, sz );
    mat4->UpdateBodyForce ( bodyforce );

    mat5->fYC.setup ( young, nu, c, phi );
    mat5->SetMemory ( nglobalpts, sz );
    mat5->UpdateBodyForce ( bodyforce );

    mat6->fYC.setup ( young, nu, c, phi );
    mat6->SetMemory ( nglobalpts, sz );
    mat6->UpdateBodyForce ( bodyforce );

    mat7->fYC.setup ( young, nu, c, phi );
    mat7->SetMemory ( nglobalpts, sz );
    mat7->UpdateBodyForce ( bodyforce );

    mat8->fYC.setup ( young, nu, c, phi );
    mat8->SetMemory ( nglobalpts, sz );
    mat8->UpdateBodyForce ( bodyforce );

    mat9->fYC.setup ( young, nu, c, phi );
    mat9->SetMemory ( nglobalpts, sz );
    mat9->UpdateBodyForce ( bodyforce );

    Doub Lx = 20.;//(*Correlation length in x direction*)
    Doub Ly = 2.;//(*Correlation length in y direction*)

    Int nsamples = 5000, expansionorder = 150;
    Int type = 3;
    KLGalerkinRF* objKLGalerkinRF2 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF3 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF4 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF5 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF6 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF7 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF8 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );
    KLGalerkinRF* objKLGalerkinRF9 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );

    objKLGalerkinRF2->SetMesh ( mesh2 );
    objKLGalerkinRF3->SetMesh ( mesh3 );
    objKLGalerkinRF4->SetMesh ( mesh4 );
    objKLGalerkinRF5->SetMesh ( mesh5 );
    objKLGalerkinRF6->SetMesh ( mesh6 );
    objKLGalerkinRF7->SetMesh ( mesh7 );
    objKLGalerkinRF8->SetMesh ( mesh8 );
    objKLGalerkinRF9->SetMesh ( mesh9 );
    slopeproject* slopeobj = new slopeproject ( mesh2, objKLGalerkinRF2 );
    //slopeproject* slopeobj3 = new slopeproject(mesh3, objKLGalerkinRF3);
    //
    bool deterministicsol = false;
    if ( deterministicsol == true ) {
        int ndesirediters = 8, niter = 50;
        Doub dlamb0 = 0.2, alphatol = 0.0001;

        Doub tol = 0.001;
        std::vector<std::vector<double>> soll;
        soll = slopeobj->IterativeProcessShearRed ( 0.01, 2., tol );
        mat2->fYC.setup ( young, nu, c, phi );
        mat2->SetMemory ( nglobalpts, sz );
        mat2->UpdateBodyForce ( bodyforce );

        //soll = slopeproject->IterativeProcess(10, 0.2, 0.001, 30);
    }

    MatDoub coesionrandomfield, frictionrandomfield;
    string filerf = "D:/slope-results/cho-field-Lx20-Ly4/coesionfield.txt";
    ReadMatDoub ( coesionrandomfield, filerf );
    string filerff = "D:/slope-results/cho-field-Lx20-Ly4/frictionfield.txt";
    ReadMatDoub ( frictionrandomfield, filerff );

    NRmatrix<MatDoub> randomfield ( 2, 1 );
    randomfield[0][0] = coesionrandomfield;
    randomfield[1][0] = frictionrandomfield;

    slopeproject* slopeobj2 = new slopeproject ( mesh2, objKLGalerkinRF2, randomfield );
    slopeproject* slopeobj3 = new slopeproject ( mesh3, objKLGalerkinRF3, randomfield );
    slopeproject* slopeobj4 = new slopeproject ( mesh4, objKLGalerkinRF4, randomfield );
    slopeproject* slopeobj5 = new slopeproject ( mesh5, objKLGalerkinRF5, randomfield );
    slopeproject* slopeobj6 = new slopeproject ( mesh6, objKLGalerkinRF6, randomfield );
    slopeproject* slopeobj7 = new slopeproject ( mesh7, objKLGalerkinRF7, randomfield );
    slopeproject* slopeobj8 = new slopeproject ( mesh8, objKLGalerkinRF8, randomfield );
    slopeproject* slopeobj9 = new slopeproject ( mesh9, objKLGalerkinRF9, randomfield );

    string namefolder2 = "D:/slope-results/SRM-cho-field-Lx20-Ly4";

    string namefolder3 = "D:/slope-results/GI-cho-field-Lx20-Ly4";

    bool print = false;
    //slopeobj2->MonteCarloSRM(1900, 5000, print, namefolder2);

    // slopeobj2->MonteCarloGIM(1290, 5000, print, namefolder3);
    string number= "2";
    std::thread thread2 ( myTreads, 0,250, slopeobj2, number );
    number= "3";
    std::thread thread3 ( myTreads, 251, 500, slopeobj3, number );
    std::thread thread4 ( myTreads, 501, 750, slopeobj4, number );
    std::thread thread5 ( myTreads, 751,1000, slopeobj5, number );
    std::thread thread6 ( myTreads, 1001, 1250, slopeobj6, number );
    std::thread thread7 ( myTreads, 1251,1500, slopeobj7, number );
    std::thread thread8 ( myTreads, 1501, 1750, slopeobj8, number );
    std::thread thread9 ( myTreads, 1751, 2000, slopeobj9, number );

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
    ReadMesh ( allcoords, meshcoords, meshtopology, elsstr, nodestr );

    std::ofstream filemesh1 ( "meshcoords.txt" );
    OutPutPost ( meshcoords, filemesh1 );
    std::ofstream filemesh2 ( "meshtopology.txt" );
    OutPutPost ( meshtopology, filemesh2 );

    Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5
    Doub thickness = 1.;
    Doub young = 20000.;
    Doub nu = 0.49;
    Int planestress = 0;

    MatDoub bodyforce ( 2, 1, 0. ), newbodyforce;
    bodyforce[1][0] = gamma;
    MatDoub ptsweigths;
    int order = 2;
    shapequad shape = shapequad ( order, 1 );
    shape.pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();
    Int nglobalpts = meshtopology.nrows() * npts;
    Int sz = 2 * meshcoords.nrows();

    elastoplastic2D< druckerprager >* mat0 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );

    int dim=2;
    mesh* mesh0 = new mesh ( dim,mat0, allcoords, meshcoords, meshtopology, hhatinho );

    mat0->fYC.setup ( young, nu, c, phi );
    mat0->SetMemory ( nglobalpts, sz );
    mat0->UpdateBodyForce ( bodyforce );

    Doub Lx = 40.;//(*Correlation length in x direction*)
    Doub Ly = 2.;//(*Correlation length in y direction*)
    Int nsamples = 50000, expansionorder = 150;
    Int type = 3;

    KLGalerkinRF* objKLGalerkinRF0 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );

    objKLGalerkinRF0->SetMesh ( mesh0 );


    if ( false ) {
        slopeproject* slopeobj0 = new slopeproject ( mesh0, objKLGalerkinRF0 );
        string filename = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2";
        slopeobj0->CreateRandomField ( filename );

        return;
    }
}

void mainlinuxserial ( int simtype,int comeco,int fim )
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
    ReadMesh ( allcoords, meshcoords, meshtopology, elsstr, nodestr );

    std::ofstream filemesh1 ( "meshcoords.txt" );
    OutPutPost ( meshcoords, filemesh1 );
    std::ofstream filemesh2 ( "meshtopology.txt" );
    OutPutPost ( meshtopology, filemesh2 );

    Doub c = 10., phi = 30 * M_PI / 180., gamma = -20.;//1.5
    Doub thickness = 1.;
    Doub young = 20000.;
    Doub nu = 0.49;
    Int planestress = 0;

    MatDoub bodyforce ( 2, 1, 0. ), newbodyforce;
    bodyforce[1][0] = gamma;
    MatDoub ptsweigths;
    int order = 2;
    shapequad shape = shapequad ( order, 1 );
    shape.pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();
    Int nglobalpts = meshtopology.nrows() * npts;
    Int sz = 2 * meshcoords.nrows();

    elastoplastic2D< druckerprager >* mat0 = new elastoplastic2D< druckerprager > ( thickness, bodyforce, planestress, order, hhatinho );

    int dim=2;
    mesh* mesh0 = new mesh ( dim,mat0, allcoords, meshcoords, meshtopology, hhatinho );

    mat0->fYC.setup ( young, nu, c, phi );
    mat0->SetMemory ( nglobalpts, sz );
    mat0->UpdateBodyForce ( bodyforce );

    Doub Lx = 40.;//(*Correlation length in x direction*)
    Doub Ly = 2.;//(*Correlation length in y direction*)
    Int nsamples = 50000, expansionorder = 150;
    Int type = 3;

    KLGalerkinRF* objKLGalerkinRF0 = new KLGalerkinRF ( order, Lx, Ly, type, nsamples, expansionorder );

    objKLGalerkinRF0->SetMesh ( mesh0 );


    if ( false ) {
        slopeproject* slopeobj0 = new slopeproject ( mesh0, objKLGalerkinRF0 );
        string filename = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2";
        slopeobj0->CreateRandomField ( filename );

        return;
    }

    MatDoub coesionrandomfield, frictionrandomfield;
    string filerf = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2/coesionfield.txt";
    ReadMatDoub ( coesionrandomfield, filerf );
    string filerff = "/home/diogo/projects/results/cho-field-414els-Lx40-Ly2/frictionfield.txt";
    ReadMatDoub ( frictionrandomfield, filerff );

    NRmatrix<MatDoub> randomfield ( 2, 1 );
    randomfield[0][0] = coesionrandomfield;
    randomfield[1][0] = frictionrandomfield;

    slopeproject* slopeobj0 = new slopeproject ( mesh0, objKLGalerkinRF0,randomfield );

    if ( false ) { //Print a single simulation
        MatDoub hhatinho2 = slopeobj0->AssembleHhationho ( 1156 ); //pior Lx =20, Ly =4
        //MatDoub hhatinho2 = slopeobj0->AssembleHhationho(1055);//melhor Lx =20, Ly =4
        mesh0->SetHhat ( hhatinho2 );

        string filename = "/home/diogo/projects/results/mathematicas/Coesao.dat";
        std::vector<std::vector<double>> hhatx;
        mesh0->fmaterial->PostProcess ( mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 0, hhatinho2, hhatx );
        std::ofstream file ( filename );
        slopeobj0->OutPutPost ( hhatx, file );

        string filename2 = "/home/diogo/projects/results/mathematicas/Phi.dat";
        mesh0->fmaterial->PostProcess ( mesh0->GetAllCoords(), mesh0->GetMeshNodes(), mesh0->GetMeshTopology(), 1, hhatinho2, hhatx );
        std::ofstream file2 ( filename );
        slopeobj0->OutPutPost ( hhatx, file2 );

        return;

        bool deterministicsol = true;
        if ( deterministicsol == true ) {
            int ndesirediters = 8, niter = 50;
            Doub dlamb0 = 0.2, alphatol = 0.0001;

            Doub tol = 0.001;
            std::vector<std::vector<double>> soll;
            //soll = slopeobj2->IterativeProcessShearRed(0.0001, 1., tol);
            mat0->fYC.setup ( young, nu, c, phi );
            mat0->SetMemory ( nglobalpts, sz );
            mat0->UpdateBodyForce ( bodyforce );
            soll = slopeobj0->IterativeProcess ( 20, 0.1, 0.0001,10 );
        }
        return;
    }

    int GIMorSRM,begin,end;
    if ( simtype< 0 ) {
        std::cout << " \n Please input the simulation type: GIM = 0, SMR = 1: ";
        std::cin >> GIMorSRM;
        std::cout << " \n Please set the initial sample: ";
        std::cin >> begin;
        std::cout << " \n Please set the last sample: ";
        std::cin >> end;
    } else {
        std::cout << " \n simulation type: " << simtype <<endl;
        std::cout << " \n initial sample: "<< comeco <<endl;
        std::cout << " \n last sample: "<< fim <<endl;

        GIMorSRM=simtype;
        begin=comeco;
        end=fim;
    }
    if ( GIMorSRM==0 ) {
        string namefolder = "/home/diogo/projects/results/gim-414els-Lx40-Ly2";
        slopeobj0->MonteCarloGIM ( begin,end, false, namefolder );
        return;
    }

    if ( GIMorSRM==1 ) {
        string namefolder = "/home/diogo/projects/results/srm-414els-Lx40-Ly2";
        slopeobj0->MonteCarloSRM ( begin,end, false, namefolder );
        return;
    }

    std::cout << " \n wrong type input!";
    return;

}

#include <Eigen/Core>
using namespace Eigen;
using namespace std;


int main3()
{
    ifstream datamatrix ( "matrix.txt" );
    ifstream datavec ( "vector.txt" );
    int m=2630;
    MatrixXd A ( m,m );
    SparseMatrix<double> SP ( m,m );
    VectorXd b ( m );
    chrono::steady_clock sc;
    auto start = sc.now();
    double temp =0.;
    for ( int i = 0; i < m; i++ ) {
        for ( int j = 0; j < m; j++ ) {
            datamatrix>>  temp;
            A ( i,j ) =temp;
            if ( fabs ( A ( i,j ) ) !=0. ) {
                SP.coeffRef ( i, j ) =temp;
            }
            // SP.coeffRef(i,j)=A(i,j);
        }
        datavec >> b ( i );
    }
    auto end = sc.now();
    auto time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "Operation took: " << time_span.count() << " seconds.";
    start = sc.now();
    VectorXd x = A.llt().solve ( b );
    end = sc.now();
    // measure time span between start & end
    time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "Operation took: " << time_span.count() << " seconds.";


    start = sc.now();
    LLT<MatrixXd> llt;
    llt.compute ( A );
    VectorXd xxx=llt.solve ( b );
    end = sc.now();

    time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "Operation took: " << time_span.count() << " seconds.";

    start = sc.now();
    SimplicialLLT< SparseMatrix<double> > solver;
    x = solver.compute ( SP ).solve ( b );
    end = sc.now();
    time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "Operation took: " << time_span.count() << " seconds.";
    //cout << xxx-x << endl;

}




template <class T>
void OutPutPost ( NRmatrix<T>& postdata, std::ofstream& file )
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


void ReadMesh ( std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord )
{
    std::vector<std::vector<Int>> topol;
    string line, temp;

    ifstream myfile ( filenameel );
    //
    if ( myfile.is_open() ) {
        while ( getline ( myfile, line ) ) {
            std::vector<string> tokens;
            istringstream iss ( line );
            while ( iss >> temp )
                tokens.push_back ( temp );
            std::vector<Int> input_int = vecstr_to_vecint ( tokens );
            for ( int k = 0; k < input_int.size(); k++ ) {
                input_int[k] = input_int[k] - 1;
            }
            topol.push_back ( input_int );
        }
        myfile.close();
    } else std::cout << "Unable to open file";

    meshtopology.CopyFromVector ( topol );
    meshtopology.Print();

    std::vector<std::vector<Doub>> coords;
    string line2, temp2;
    ifstream myfile2 ( filenamecoord );
    if ( myfile2.is_open() ) {
        while ( getline ( myfile2, line2 ) ) {
            std::vector<string> tokens;
            istringstream iss ( line2 );
            while ( iss >> temp2 )
                tokens.push_back ( temp2 );
            std::vector<Doub> input_doub = vecstr_to_vecdoub ( tokens );

            //std::vector<Doub> input_doub2(input_doub.size() - 1);
            //for (int k = 1;k < input_doub.size();k++)
            //{
            //	input_doub2[k] = input_doub[k];
            //}

            coords.push_back ( input_doub );
        }
        myfile2.close();
    } else std::cout << "Unable to open file";

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

  //  cout << "leu"<<endl;
   // DebugStop();
}
std::vector<Int> vecstr_to_vecint ( std::vector<string> vs )
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
// std::vector<Int> vecstr_to_vecint ( std::vector<string> vs )
// {
//     std::vector<Int> ret;
//     for ( std::vector<string>::iterator it = vs.begin() + 2; it != vs.end(); ++it ) {
//         istringstream iss ( *it );
//         Int temp;
//         iss >> temp;
//         ret.push_back ( temp );
//     }
//     return ret;
// }

std::vector<Doub> vecstr_to_vecdoub ( std::vector<string> vs )
{
    std::vector<Doub> ret;
    for ( std::vector<string>::iterator it = vs.begin() +1; it != vs.end() ; ++it ) {
        istringstream iss ( *it );
        Doub temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    //ret.push_back(0.);
    return ret;
}

std::vector<Doub> vecstr_to_vecdoub2 ( std::vector<string> vs )
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
std::vector<T> vecstr_to_vec ( std::vector<string> vs )
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

void ReadMatDoub ( MatDoub& matdoub, std::string  file )
{

    std::vector<std::vector<Doub>> coords;
    string line2, temp2;
    ifstream myfile2 ( file );
    if ( myfile2.is_open() ) {
        while ( getline ( myfile2, line2 ) ) {
            std::vector<string> tokens;
            istringstream iss ( line2 );
            while ( iss >> temp2 )
                tokens.push_back ( temp2 );
            std::vector<Doub> input_doub = vecstr_to_vecdoub2 ( tokens );
            coords.push_back ( input_doub );
        }
        myfile2.close();
    } else std::cout << "Unable to open file";
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

void solvequd()
{
	string file ="/home/diogo/projects/dcproj/data/tri-p2.msh";//SRM 1.23
	readgidmesh read = readgidmesh(file);
	read.ReadMesh();
	MatInt meshtopology = read.GetTopology();
	MatDoub meshcoords = read.GetCoords();
	std::vector<std::vector<std::vector<Doub>>>  allcoords = read.GetAllCoords();
	//Uma coordenada em cima da linha. Qualquer coordenada
	NRvector<double> constcoorddata ( 3,0. );
    constcoorddata[0]=0.;
    constcoorddata[1]=0.;
    constcoorddata[2]=0.;
    std::vector<int> idsline;
	//Direcao a buscar. 0 significa que o algoritmo e livre para buscar naquela direcao e 1 que dizer que é fixo.
    NRvector<int>constcoord2 ( 3 );
    constcoord2[0]=0;//livre
    constcoord2[1]=1;//fixo
    constcoord2[2]=1;//fixo
	read.FindIds(constcoorddata,constcoord2,idsline);
	 cout << "idsline[i] "<< endl;
	for(Int i=0;i<idsline.size();i++)cout << idsline[i] << endl;
	//Coordenada do ponto
    constcoorddata[0]=10.;
    constcoorddata[1]=10.;
    constcoorddata[2]=0.;
    std::vector<int> idspoint;
	//Direcao a buscar. 0 significa que o algoritmo e livre para buscar naquela direcao e 1 que dizer que é fixo.
    constcoord2[0]=1;//livre
    constcoord2[1]=1;//fixo
    constcoord2[2]=1;//fixo
	read.FindIds(constcoorddata,constcoord2,idspoint);
	
	
	//Uma coordenada em cima da linha. Qualquer coordenada
    constcoorddata[0]=10.;
    constcoorddata[1]=10.;
    constcoorddata[2]=0.;
    std::vector<int> idslinetop;
	//Direcao a buscar. 0 significa que o algoritmo e livre para buscar naquela direcao e 1 que dizer que é fixo.
    constcoord2[0]=0;//livre
    constcoord2[1]=1;//fixo
    constcoord2[2]=1;//fixo
	read.FindIds(constcoorddata,constcoord2,idslinetop);

	 cout << "idslinetop[i] "<< endl;
	for(Int i=0;i<idslinetop.size();i++)cout << idslinetop[i] << endl;
	cout << "end "<< endl;
	//elastmat2D::elastmat2D (Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order )
	Doub young=10000.,  nu=0.3,  thickness=1.,bodyforce=0.;
	Int planestress=1;

	Int elnodes = meshtopology.ncols();
	elastmat2D * mat = new elastmat2D ( young,  nu,  thickness,  bodyforce,  planestress,  elnodes );
	Int dim =2;
	mesh* meshs = new mesh( dim, mat, allcoords, meshcoords, meshtopology);
	NRmatrix<Doub> KG,F,Kinv,sol;
 	meshs->AssembleLinear(KG,F);
	KG.Print();
	//Aplica forca em x no no (10,10)

	
	Int dir, val;
    dir = 1;
    val = 0;
    mat->DirichletBC ( KG, F, idsline, dir, val );
	dir = 0;
    mat->DirichletBC ( KG, F, idsline, dir, val );
	KG.Print();
//     LUdcmp* lu = new LUdcmp ( KG );
//     lu->inverse ( Kinv );
	string matrixout = "/home/diogo/projects/dcproj/matrix.dat";
	std::ofstream matrixoutfile ( matrixout );
	PrintMathematicaFormat ( KG,matrixoutfile );
 
	NRmatrix<Int> toplinetopol ;
	//toplinetopol[0][1]=2;
 	std::vector<std::vector<int>> topolline=read.LineTopology(idslinetop);
 	for(Int i=0;i<topolline.size();i++)cout << topolline[i][0] << endl;
 	read.ToMatInt(topolline,toplinetopol);
	toplinetopol.Print();
	Doub force;
	force=1.;
	mat->AssembleLoadvector ( KG,  F, meshcoords,toplinetopol,  force );
	F.Print();
	matrixout = "/home/diogo/projects/dcproj/f.dat";
	std::ofstream foutfile ( matrixout );
	PrintMathematicaFormat ( F,foutfile );
	
	//F[2*2][0]=10.;
	F.Print();
	KG.ComputeInverse(Kinv);
	Kinv.Mult(F,sol);
	NRmatrix<Doub> temp;
	KG.Mult(sol,temp);
	temp-=F;
	cout << "residual"<<endl;
	temp.Print();
	cout << "sol"<<endl;
	sol.Print();
	
	meshs->fmaterial->UpdateDisplacement(sol);
	
	std::vector<string> scalar_names;
    std::vector<string> vector_names;
    vector_names.push_back ( "Displacement" );
    vector_names.push_back ( "Strain" );
    string slopestr="/home/diogo/projects/dcproj/output/quad";
    VTKGraphMesh vtkobj ( meshs,dim,scalar_names,vector_names,slopestr );
    vtkobj.DrawSolution ( 0);
	
}
