#include <sstream>
#include "vtkmesh.h"
//using namespace std;


VTKGraphMesh::VTKGraphMesh(mesh *cmesh, int dimension,std::vector<std::string> &scalnames, const std::vector<std::string> &vecnames,  std::string FileName)
{
    fmesh=cmesh;
    fdim=dimension;
    fFileName = FileName;
    fdata=fmesh->fmaterial->GetSolution();
    fVecNames=vecnames;
    fScalNames=scalnames;
}

void VTKGraphMesh::DrawSolution(Int step, Doub time){

	NRmatrix<Int> meshtopol=GetMeshTopologyVTK();

    material * matp = fmesh->fmaterial;
	if(fOutFile.is_open())
	{
		fOutFile.close();
	}
	int n;
	/*{
		std::stringstream sout;
		sout << fFileName.substr(0,fFileName.size()-4) << ".scal_vec." << step << ".vtk";
       // sout << fFileName.substr(0,fFileName.size()-4) << "beamxvtkfile" << step << ".vtk";
		fOutFile.open(sout.str().c_str());
	}*/

	std::string extension= ".scal_vec.vtk";
	auto name =fFileName;
    name+=std::to_string(step);
    name+=extension;
    std::stringstream sout;
    sout << name;
    fOutFile.open(sout.str().c_str());

	fOutFile << "# vtk DataFile Version 3.0" << endl;
	fOutFile << "File generated from DClib" << endl;
	fOutFile << "ASCII\n" << endl;
	fOutFile << "DATASET UNSTRUCTURED_GRID" << endl;
	//DrawNodes();
	//DrawConnectivity(ECube);
    Int els=meshtopol.nrows();
    NRmatrix<Doub> nodes = fmesh->GetMeshNodes();
    Int nnodes = nodes.nrows();
    Int elnodes = meshtopol.ncols();
    fOutFile << "POINTS " << nnodes << " float"<<endl;
    for (Int inode = 0; inode < nnodes; inode++)
	{
        for(Int idim=0;idim<3;idim++)
        {
            fOutFile << nodes[inode][idim] << " ";
        }
        fOutFile<<std::endl;
	}

    Int type;
    if(elnodes==8)
    {
        type = 23;
    }
    else if(elnodes==20)
    {
        type =25;
    }
	(fOutFile) << "CELLS " << els << " "<< els*elnodes+els <<  endl;
    for (Int iel = 0; iel < els; iel++)
	{
        fOutFile << elnodes << " ";
        for(Int ielnode=0;ielnode<elnodes;ielnode++)
        {
            fOutFile << meshtopol[iel][ielnode] << " ";
        }
        fOutFile<<std::endl;
	}

	(fOutFile) << "CELL_TYPES " << els <<  endl;
    for (Int iel = 0; iel < els; iel++)
	{
        fOutFile << type <<  endl;
	}
	//Int nsol = fdata.nrows();
    fOutFile << "POINT_DATA " << nnodes << endl;
	fOutFile << "SCALARS scalars float 1" << endl;
	fOutFile << "LOOKUP_TABLE default" << endl;
    for (Int inode = 0; inode < nnodes; inode++)
	{
        fOutFile << fdata[inode*fdim][0] <<   endl;//x
        //fOutFile << sqrt(fdata[inode][0]*fdata[inode][0]+fdata[inode+1][0]*fdata[inode+1][0]) <<   endl;
	}

  /*  for (Int iel = 0; iel < els; iel++)
	{
        //fOutFile << elnodes << " ";
        for(Int ielnode=0;ielnode<elnodes;ielnode++)
        {
            fOutFile << fdata[meshtopol[iel][ielnode]][0] << " ";
        }
        fOutFile<<std::endl;
	}*/

    Int nvecnames = fVecNames.size();
    for(Int ivar=0;ivar<1;ivar++)
    {
        (fOutFile) << "VECTORS " << fVecNames[ivar] << " float" << std::endl;
        if(fVecNames[ivar] =="Displacement")
        {
            for (Int inode = 0; inode < nnodes; inode++)
            {
                for(Int idim=0;idim<fdim;idim++)
                {
                    fOutFile << fdata[inode*fdim+idim][0] << " ";
                }
                fOutFile << 0.;
                fOutFile << std::endl;
            }
        }
        else if(fVecNames[ivar] =="Strain")
        {
            for (Int inode = 0; inode < nnodes; inode++)
            {
                for(Int idim=0;idim<fdim;idim++)
                {
                    fmesh->fmaterial->GetPlas
                    fOutFile << fdata[inode*fdim+idim][0] << " ";
                }
                fOutFile << 0.;
                fOutFile << std::endl;
            }
        }
        else if(fVecNames[ivar] =="SqrtJ2(EPSP)")
        {

        }
        else if(fVecNames[ivar] =="Stress")
        {

        }*/

    }


	fOutFile.close();
}

