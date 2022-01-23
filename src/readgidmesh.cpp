#include "readgidmesh.h"


readgidmesh::readgidmesh(){
	
}

readgidmesh::~readgidmesh(){
	
}

readgidmesh::readgidmesh(string file ){
	ffile=file;
}

template <class T>
std::vector<T> readgidmesh::str_vec( std::vector<std::string> &vs )
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

void readgidmesh::ReadMesh (  )
{
    std::vector<std::vector<Int>> topol;
    string line, temp;

    ifstream myfile ( ffile );

	
    std::vector<std::vector<Doub>> coords;

    if ( myfile.is_open() ) {
		getline ( myfile, line );
		std::cout << line<< std::endl;
		getline ( myfile, line );
		std::cout << line<< std::endl;
        while ( getline ( myfile, line ) ) {
			if(line=="End Coordinates")break;
            std::vector<string> tokens;
            istringstream iss ( line );
            while ( iss >> temp )
                tokens.push_back ( temp );
			std::vector<Doub> input_doub_temp= str_vec<Doub> ( tokens );
            coords.push_back ( input_doub_temp );
        }
       // myfile.close();
    } else std::cout << "Unable to open file";
	Int nodes = coords.size();
	fmeshcoords.assign(nodes,3,0.);
	for(Int inode=0;inode<nodes;inode++){
		fmeshcoords[inode][0]=coords[inode][1];
		fmeshcoords[inode][1]=coords[inode][2];
		fmeshcoords[inode][2]=coords[inode][3];
	}
    fmeshcoords.Print();
	
    if ( myfile.is_open() ) {
		getline ( myfile, line );
		std::cout << line<< std::endl;
		getline ( myfile, line );
		std::cout << line<< std::endl;
        while ( getline ( myfile, line ) ) {
			if(line=="End Elements")break;
            std::vector<string> tokens;
            istringstream iss ( line );
            while ( iss >> temp )
                tokens.push_back ( temp );
            std::vector<Int> input_int = str_vec<Int> ( tokens );
            topol.push_back ( input_int );
        }
        myfile.close();
    } else std::cout << "Unable to open file";

	Int els = topol.size();
	Int elnodes = topol[0].size()-1;
	fmeshtopology.assign(els,elnodes,0);
	for(Int iel=0;iel<els;iel++){
		for(Int elnode=0;elnode<elnodes;elnode++)
		{
			fmeshtopology[iel][elnode]=topol[iel][elnode+1]-1;
		}
	}
    fmeshtopology.Print();
	
	std::vector<Doub> temp33 ( 3 );
    for ( Int i = 0; i < fmeshtopology.nrows(); i++ ) {
        std::vector< std::vector<Doub> > temp22;
        for ( Int j = 0; j < fmeshtopology.ncols(); j++ ) {
            Int top = fmeshtopology[i][j];
            temp33[0] = fmeshcoords[top][0];
            temp33[1] = fmeshcoords[top][1];
            temp33[2] = fmeshcoords[top][2];
            temp22.push_back ( temp33 );
        }
        fallcoords.push_back ( temp22 );
    }
    elnodes=topol[0].size()-1;
	
	if(elnodes==3)
 	{
	  fOrder=1;
	}
	if(elnodes  ==4)
 	{
	  fOrder=1;
	}
	if(elnodes==6)
 	{
	  fOrder=2;
	}
	if(elnodes==8)
 	{
	  fOrder=2;
	}
}
