#include <sstream>
#include "vtkmesh.h"
//using namespace std;
#include "druckerprager.h"

VTKGraphMesh::VTKGraphMesh ( mesh *cmesh, int dimension,std::vector<std::string> &scalnames, const std::vector<std::string> &vecnames,  std::string FileName )
{
    fmesh=cmesh;
    fdim=dimension;
    fFileName = FileName;
    fdata=fmesh->fmaterial->GetSolution();
    fVecNames=vecnames;
    fScalNames=scalnames;
}

void VTKGraphMesh::DrawSolution ( Int step )
{

    NRmatrix<Int> meshtopol=GetMeshTopologyVTK();

    material * matp = fmesh->fmaterial;

    std::string extension= ".scal_vec.vtk";
    auto name =fFileName;
    name+=std::to_string ( step );
    name+=extension;
    std::stringstream sout;
    sout << name;
    fOutFile.open ( sout.str().c_str() );

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
    for ( Int inode = 0; inode < nnodes; inode++ ) {
        for ( Int idim=0; idim<3; idim++ ) {
            fOutFile << nodes[inode][idim] << " ";
        }
        fOutFile<<std::endl;
    }

    Int type;
    if ( elnodes==8 ) {
        type = 23;
    } else if ( elnodes==20 ) {
        type =25;
    } else if ( elnodes==6 ) {
        type =22;
    } else if ( elnodes==3 ) {
        type =5;
    } else if ( elnodes==4 ) {
        type=9;
    }
	
    ( fOutFile ) << "CELLS " << els << " "<< els*elnodes+els <<  endl;
    for ( Int iel = 0; iel < els; iel++ ) {
        fOutFile << elnodes << " ";
        for ( Int ielnode=0; ielnode<elnodes; ielnode++ ) {
            fOutFile << meshtopol[iel][ielnode] << " ";
        }
        fOutFile<<std::endl;
    }


    ( fOutFile ) << "CELL_TYPES " << els <<  endl;
    for ( Int iel = 0; iel < els; iel++ ) {
        fOutFile << type <<  endl;
    }

    //  NRvector<  NRmatrix<Doub>  > sol,dsol;
    //  fmesh->fmaterial->ComputeSolAndDSol(fmesh,sol,dsol);

    NRmatrix<Doub>   sol2,dsol2;
    fmesh->fmaterial->ComputeSolAndDSol ( fmesh,sol2,dsol2 );

    fOutFile << "POINT_DATA " << nnodes << endl;


//     for ( Int inode = 0; inode < nnodes; inode++ ) {
//         NRmatrix<Doub> eps,gradu,gradut;
//         NRtensor<Doub> straintensor ( 0. );
// 
//         if ( fdim==3 ) {
// 			DebugStop();
//         } else if ( fdim==2 ) {
//             Doub dudx= dsol2[inode*2][0];
//             Doub dudy= dsol2[inode*2][1];
//             Doub dwdx= dsol2[inode*2+1][0];
//             Doub dwdy= dsol2[inode*2+1][1];
// 			straintensor.XX()=dudx;straintensor.XY()=(dwdx+dudy)/2.;
// 			straintensor.YY()=dwdy;
//         }
//         NRvector<Doub> valphi =fmesh->fmaterial->ComputePhi ( straintensor ) ;
//         fOutFile << valphi[0]<<   endl;//x
//     }
    NRmatrix<Doub> HHAT;
    fmesh->GetHhat ( HHAT );
    if ( HHAT.nrows() !=0 ) {
        fOutFile << "SCALARS cohesion float 1" << endl;
        fOutFile << "LOOKUP_TABLE default" << endl;
        for ( Int inode = 0; inode < nnodes; inode++ ) {
            fOutFile << HHAT[inode][0] << std::endl;
        }
		fOutFile << "SCALARS friction float 1" << endl;
        fOutFile << "LOOKUP_TABLE default" << endl;
		for ( Int inode = 0; inode < nnodes; inode++ ) {
			fOutFile << (HHAT[inode][1]*180/M_PI) << std::endl;
        }
    }

    


    Int nvecnames = fVecNames.size();
    for ( Int ivar=0; ivar<nvecnames; ivar++ ) {
        ( fOutFile ) << "VECTORS " << fVecNames[ivar] << " float" << std::endl;
        if ( fVecNames[ivar] =="Displacement" ) {
            for ( Int inode = 0; inode < nnodes; inode++ ) {
                fOutFile << sol2[inode*2][0] << " " << sol2[inode*2+1][0] << " ";
                fOutFile << 0;
                fOutFile << std::endl;

            }

        } else if ( fVecNames[ivar] =="Strain" ) {

            for ( Int inode = 0; inode < nnodes; inode++ ) {
                Doub dudx= dsol2[inode*2][0];
                Doub dudy= dsol2[inode*2][1];
                Doub dwdx= dsol2[inode*2+1][0];
                Doub dwdy= dsol2[inode*2+1][1];
                fOutFile << dudx << " " << dwdy <<  " "<< ( dudy+dwdx ) /2.;
                fOutFile << std::endl;
            }

        } else if ( fVecNames[ivar] =="sas" ) {

        } else if ( fVecNames[ivar] =="Stress" ) {

        }

    }


    fOutFile.close();
}
