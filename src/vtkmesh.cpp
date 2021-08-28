#include <sstream>
#include "vtkmesh.h"
//using namespace std;
#include "druckerprager.h"

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

     NRvector<  NRmatrix<Doub>  > sol,dsol;
    fmesh->fmaterial->ComputeSolAndDSol(fmesh,sol,dsol);

    fOutFile << "POINT_DATA " << nnodes << endl;
	fOutFile << "SCALARS phi float 1" << endl;
	fOutFile << "LOOKUP_TABLE default" << endl;

    for (Int inode = 0; inode < nnodes; inode++)
    {
        NRmatrix<Doub> eps,gradu,gradut;
        NRtensor<Doub> tensor(0.);
        gradu=dsol[inode];
       // gradu.Transpose(gradut);
       // eps=gradu;
       // eps+=gradut;
      //  eps*=1/2.;
        if(fdim==3)
        {
            //tensor.XX()=eps[0][0];tensor.XY()=eps[0][1];tensor.XZ()=eps[0][2];
            //tensor.XY()=eps[1][0];tensor.YY()=eps[1][1];tensor.YZ()=eps[1][2];
            //tensor.XZ()=eps[2][0];tensor.YZ()=eps[2][1];tensor.ZZ()=eps[2][2];
        }
        else if(fdim==2)
        {
            //	Doub ex = gradu[0][0];// dudx;
	//Doub ey = gradu[1][1];
	//Doub exy = (gradu[0][1] + gradu[1][0]);
            //tensor.XX()=gradu[0][0];tensor.XY()=(gradu[1][0]+gradu[0][1]);
            //tensor.XY()=(gradu[1][0]+gradu[0][1]);tensor.YY()=gradu[1][1];
            //tensor.YY()=gradu[1][1];
            tensor.XX()=gradu[0][0];tensor.YY()=gradu[1][1];tensor.XY()=(gradu[0][1] + gradu[1][0]);
        }
        Doub valphi =fmesh->fmaterial->ComputePhi(tensor) ;
        fOutFile << valphi <<   endl;//x
    }
        //sol[inode].assign(2,1,0.);
        //dsol[inode].assign(2,2,0.);

    NRmatrix<Doub> HHAT;
    fmesh->GetHhat(HHAT);
    if(HHAT.nrows()!=0)
    {
        fOutFile << "SCALARS cohesion float 1" << endl;
        fOutFile << "LOOKUP_TABLE default" << endl;

        for (Int inode = 0; inode < nnodes; inode++)
        {
            fOutFile << HHAT[inode][0] <<   endl;//x
        }
    }




    //sol[0].Print();
    //sol[5].Print();
    //fdata.Print();
    Int nvecnames = fVecNames.size();
    for(Int ivar=0;ivar<nvecnames;ivar++)
    {
        (fOutFile) << "VECTORS " << fVecNames[ivar] << " float" << std::endl;
        if(fVecNames[ivar] =="Displacement")
        {
            for (Int inode = 0; inode < nnodes; inode++)
            {
                //fOutFile << sol[inode][0][0] << " " << sol[inode][0][1] << " ";
                fOutFile << sol[inode][0][0] << " " << sol[inode][0][1] << " ";
                fOutFile << 0;
                fOutFile << std::endl;

            }

        }
        else if(fVecNames[ivar] =="Strain")
        {
            NRmatrix<Doub> eps,gradu,gradut;
            //NRmatrix<Doub> eps,gradu,gradut;
            for (Int inode = 0; inode < nnodes; inode++)
            {
                gradu=dsol[inode];
                //gradu.Transpose(gradut);
                //eps=gradu;
                //eps+=gradut;
                //eps*=1/2.;
                	//Doub ex = gradu[0][0];// dudx;
	//Doub ey = gradu[1][1];
	//Doub exy = (gradu[0][1] + gradu[1][0]);
                fOutFile << gradu[0][0] << " " << gradu[1][1] <<  " "<< (gradu[1][0] );
                //fOutFile << gradu[0][0] << " " << gradu[0][1]<< " " << " "<< gradu[1][0];
                fOutFile << std::endl;
            }

        }
        else if(fVecNames[ivar] =="SqrtJ2(EPSP)")
        {

        }
        else if(fVecNames[ivar] =="Stress")
        {

        }

    }


	fOutFile.close();
}
