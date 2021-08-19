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

    fOutFile << "POINT_DATA " << nnodes << endl;
	fOutFile << "SCALARS DisplacementMagnitude float 1" << endl;
	fOutFile << "LOOKUP_TABLE default" << endl;
    for (Int inode = 0; inode < nnodes; inode++)
    {
        fOutFile << sqrt(fdata[inode*fdim+0][0]*fdata[inode*fdim+0][0]+fdata[inode*fdim+1][0]*fdata[inode*fdim+1][0]) <<   endl;//x
    }


    NRmatrix<Doub> HHAT;
    fmesh->GetHhat(HHAT);



   // NRvector<NRtensor<Doub> > plasticstrain=fmesh->fmaterial->GetPlasticStrain();
 //NRmatrix<Doub> sol,dsol;

 //   fmesh->fmaterial->ComputeSolAndDSol(fmesh,sol,dsol);


     NRvector<  NRmatrix<Doub>  > sol,dsol;
    fmesh->fmaterial->ComputeSolAndDSol(fmesh,sol,dsol);
    sol[0].Print();
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

                fOutFile << sol[inode][0][0] << " " << sol[inode][0][1] << " ";

                fOutFile << 0;
                fOutFile << std::endl;
               /* NRvector<int>constcoord2(3);
                constcoord2[0]=1;//livre para procurar x
                constcoord2[1]=1;//livre para procurar y
                constcoord2[2]=1;//livre para procurar z
                std::vector<Int> ids;
                NRvector<Doub> doubcoords(3);
                doubcoords[0]=nodes[inode][0];
                doubcoords[1]=nodes[inode][1];
                doubcoords[2]=nodes[inode][2];
                fmesh->FindIds(doubcoords, constcoord2 , ids);*/


               // for(Int idim=0;idim<fdim;idim++)
                //{
                    //fOutFile << sol[inode*fdim+idim][0] << " ";
               // }
               // fOutFile << 0.;
              //  fOutFile << std::endl;
            }

          /*  for (Int inode = 0; inode < nnodes; inode++)
            {
                Int count=0;
                for (Int iel = 0; iel < els; iel++)
                {
                    for(Int ielnode=0;ielnode<elnodes;ielnode++)
                    {
                        Int id;
                        id = meshtopol[iel][ielnode];
                        if(id==inode)
                        {
                            for(Int idim=0;idim<fdim;idim++)
                            {
                                fOutFile << sol[id*fdim+idim][0] << " ";
                            }
                            fOutFile << 0.;
                            fOutFile << std::endl;
                        }
                        count++;
                    }
                }

            }*/




        }
        else if(fVecNames[ivar] =="Strain")
        {
            NRmatrix<Doub> eps,gradu,gradut;
            for (Int inode = 0; inode < nnodes; inode++)
            {
                gradu=dsol[inode];
                gradu.Transpose(gradut);
                eps=gradu;
                eps+=gradut;
                eps*=1/2.;
                fOutFile << eps[0][0] << " " << eps[1][1]<< " " << " "<< eps[1][0];
                fOutFile << std::endl;
            }

            /*
            for (Int inode = 0; inode < nnodes; inode++)
            {
                for(Int idim=0;idim<fdim;idim++)
                {
                    fOutFile << dsol[inode*fdim+idim][0] << " ";//duxdx e duxdy
                }
                fOutFile << 0.;
                fOutFile << std::endl;
            }
            */


        }
        else if(fVecNames[ivar] =="SqrtJ2(EPSP)")
        {

        }
        else if(fVecNames[ivar] =="Stress")
        {

        }

    }

   /* Int nvecnames = fVecNames.size();
    for(Int ivar=0;ivar<nvecnames;ivar++)
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
            NRvector<NRvector<NRtensor<Doub>>> strain;
            fmesh->fmaterial->PostProcessStrain(fmesh,strain);

            bool repeat;
            for (Int inode = 0; inode < nnodes; inode++)
            {
                    repeat=false;
                    NRtensor<Doub> straint;
                    Int count=0;
                    for (Int iel = 0; iel < strain.size(); iel++)
                    {

                        for(Int ielnode=0;ielnode<strain[iel].size();ielnode++)
                        {
                            Int id;
                            id= meshtopol[iel][ielnode];
                            if(id==inode && repeat==false)
                            //if(id==inode)
                            {
                                straint+=strain[iel][ielnode];
                               // std::cout << " id = " << id <<std::endl;
                                //strain[iel][ielnode].Print();
                               // fOutFile << strain[iel][ielnode].XX() << " "<<strain[iel][ielnode].YY() << " "<< strain[iel][ielnode].XY();
                                //repeat=true;
                                fOutFile << std::endl;
                                count++;
                            }
                        }
                    }
                    straint*=1./Doub(count);
                    fOutFile << straint.XX() << " "<<straint.YY() << " "<< straint.XY();
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
    */


	fOutFile.close();
}
