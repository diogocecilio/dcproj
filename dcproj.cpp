// dcproj.cpp : Defines the entry point for the application.
//

#include "dcproj.h"
#include "mesh.h"


using namespace std;

void ReadMatDoub(MatDoub& matdoub, std::string  file);
void ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord);
template <class T>
void OutPutPost(NRmatrix<T>& postdata, std::ofstream& file);
template <class T>
void vecstr_to_vec(std::vector<std::string> vs, std::vector<T> &ret);
std::vector<Int> vecstr_to_vecint(std::vector<string> vs);
std::vector<Doub> vecstr_to_vecdoub(std::vector<string> vs);

void mainlinux()
{
	cout << "HIasd " << endl;

	string nodestr = "/home/joaohenrique/Documents/dcproj/nos-132-c3.txt";
	string elsstr = "/home/joaohenrique/Documents/dcproj/els-132-c3.txt";

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

	elastoplastic2D< druckerprager >* mat = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);

	mesh* mesh2 = new mesh(mat, allcoords, meshcoords, meshtopology, hhatinho);

	int szdebug = mesh2->GetAllCoords().size();

	mat->fYC.setup(young, nu, c, phi);
	mat->SetMemory(nglobalpts, sz);
	mat->UpdateBodyForce(bodyforce);

	Doub Lx = 20.;//(*Correlation length in x direction*)
	Doub Ly = 2.;//(*Correlation length in y direction*)

	Int nsamples = 5000, expansionorder = 150;
	Int type = 3;
	KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	objKLGalerkinRF->SetMesh(mesh2);
	slopeproject* slopeobj = new slopeproject(mesh2, objKLGalerkinRF);
	//	
	bool deterministicsol = false;
	if (deterministicsol == true)
	{
		int ndesirediters = 8, niter = 50;
		Doub dlamb0 = 0.2, alphatol = 0.0001;

		Doub tol = 0.001;
		std::vector<std::vector<double>> soll;
		soll = slopeobj->IterativeProcessShearRed(0.01, 2., tol);
		mat->fYC.setup(young, nu, c, phi);
		mat->SetMemory(nglobalpts, sz);
		mat->UpdateBodyForce(bodyforce);

		//soll = slopeproject->IterativeProcess(10, 0.2, 0.001, 30);
	}
	cout << "HI " << endl;
	cout << "HIsdsadasd " << endl;
	MatDoub coesionrandomfield, frictionrandomfield;
	string filerf = "/home/joaohenrique/Documents/dcproj/randomfielddataLx20-Ly2-lognormal/coesionfield.txt";
	ReadMatDoub(coesionrandomfield, filerf);
	string filerff = "/home/joaohenrique/Documents/dcproj/randomfielddataLx20-Ly2-lognormal/frictionfield.txt";
	ReadMatDoub(frictionrandomfield, filerff);
	cout << "Hdasdasdasdasdasdadsadasda " << endl;

	NRmatrix<MatDoub> randomfield(2, 1);
	randomfield[0][0] = coesionrandomfield;
	randomfield[1][0] = frictionrandomfield;

	slopeproject* slopeobj2 = new slopeproject(mesh2, objKLGalerkinRF, randomfield);

	string namefolder2 = "/home/joaohenrique/Documents/dcproj/SRM-Lx20-Ly2";

	string namefolder3 = "/home/joaohenrique/Documents/dcproj/GIM-Lx20-Ly2";

	bool print = false;
	slopeobj2->MonteCarloGIM(0, 2000, print, namefolder3);
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

	elastoplastic2D< druckerprager >* mat = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);

	mesh* mesh2 = new mesh(mat, allcoords, meshcoords, meshtopology, hhatinho);

	int szdebug = mesh2->GetAllCoords().size();

	mat->fYC.setup(young, nu, c, phi);
	mat->SetMemory(nglobalpts, sz);
	mat->UpdateBodyForce(bodyforce);

	Doub Lx = 20.;//(*Correlation length in x direction*)
	Doub Ly = 2.;//(*Correlation length in y direction*)

	Int nsamples = 5000, expansionorder = 150;
	Int type = 3;
	KLGalerkinRF* objKLGalerkinRF = new KLGalerkinRF(order, Lx, Ly, type, nsamples, expansionorder);
	objKLGalerkinRF->SetMesh(mesh2);
	slopeproject* slopeobj = new slopeproject(mesh2, objKLGalerkinRF);
	//	
	bool deterministicsol = false;
	if (deterministicsol == true)
	{
		int ndesirediters = 8, niter = 50;
		Doub dlamb0 = 0.2, alphatol = 0.0001;

		Doub tol = 0.001;
		std::vector<std::vector<double>> soll;
		soll = slopeobj->IterativeProcessShearRed(0.01, 2., tol);
		mat->fYC.setup(young, nu, c, phi);
		mat->SetMemory(nglobalpts, sz);
		mat->UpdateBodyForce(bodyforce);

		//soll = slopeproject->IterativeProcess(10, 0.2, 0.001, 30);
	}
	cout << "HI " << endl;
	cout << "HIsdsadasd " << endl;
	MatDoub coesionrandomfield, frictionrandomfield;
	string filerf = "D:/slope-results/cho-field/coesionfield.txt";
	ReadMatDoub(coesionrandomfield, filerf);
	string filerff = "D:/slope-results/cho-field/frictionfield.txt";
	ReadMatDoub(frictionrandomfield, filerff);
	cout << "Hdasdasdasdasdasdadsadasda " << endl;

	NRmatrix<MatDoub> randomfield(2, 1);
	randomfield[0][0] = coesionrandomfield;
	randomfield[1][0] = frictionrandomfield;

	slopeproject* slopeobj2 = new slopeproject(mesh2, objKLGalerkinRF, randomfield);

	string namefolder2 = "D:/slope-results/SRM-cho-field";

	string namefolder3 = "D:/slope-results/GI-cho-field";

	bool print = false;
	slopeobj2->MonteCarloGIM(1900, 5000, print, namefolder3);
}

int main()
{
#ifdef __unix__                    /* __unix__ is usually defined by compilers targeting Unix systems */
	mainlinux();
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
		for (Int j = 0; j < postdata.ncols();j++)
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
