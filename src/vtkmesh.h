#include "nr3.h"
#include "mesh.h"
#include <iostream>
#include <set>

//enum eltype {Ehexa =25, Equad =8};

class VTKGraphMesh
{

public:

    /** @brief Constructor for graphical mesh using VTK format with tensor variables */
    VTKGraphMesh ( mesh *cmesh, int dimension,std::vector<std::string> &scalnames, const std::vector<std::string> &vecnames, std::string FileName );

    void DrawNodes();
    void DrawConnectivity();
    void DrawSolution ( Int step, Doub time );

    MatInt GetMeshTopologyVTK()
    {
        NRmatrix<Int> meshtopology = fmesh->GetMeshTopology();
        Int els = meshtopology.nrows();
        Int elnodes =  meshtopology.ncols();
        NRmatrix<Int> vtktopol ( els,elnodes );
        if ( elnodes==8 && fdim==2 ) { //quad 8 noded
            vtktopol=meshtopology;
        } else if ( elnodes==4 && fdim==2 ) { //quad 4 noded
            vtktopol=meshtopology;
        } else if ( elnodes==8 && fdim==3 ) { //cube 8 noded
            vtktopol=meshtopology;
        } else if ( elnodes==6 && fdim==2 ) { //cube 8 noded
            vtktopol=meshtopology;
        } else if ( elnodes==3 && fdim==2 ) { //cube 8 noded
            vtktopol=meshtopology;
        }else if ( elnodes==20 && fdim==3 ) { //cube 20 noded
            for ( Int iel = 0; iel < els; iel++ ) {

                vtktopol[iel][0]=meshtopology[iel][0];
                vtktopol[iel][1]=meshtopology[iel][1];
                vtktopol[iel][2]=meshtopology[iel][2];
                vtktopol[iel][3]=meshtopology[iel][3];

                vtktopol[iel][4]=meshtopology[iel][4];
                vtktopol[iel][5]=meshtopology[iel][5];
                vtktopol[iel][6]=meshtopology[iel][6];
                vtktopol[iel][7]=meshtopology[iel][7];

                vtktopol[iel][8]=meshtopology[iel][8];
                vtktopol[iel][9]=meshtopology[iel][9];
                vtktopol[iel][10]=meshtopology[iel][10];
                vtktopol[iel][11]=meshtopology[iel][11];

                vtktopol[iel][12]=meshtopology[iel][16];
                vtktopol[iel][13]=meshtopology[iel][17];
                vtktopol[iel][14]=meshtopology[iel][18];
                vtktopol[iel][15]=meshtopology[iel][19];

                vtktopol[iel][16]=meshtopology[iel][12];
                vtktopol[iel][17]=meshtopology[iel][13];
                vtktopol[iel][18]=meshtopology[iel][14];
                vtktopol[iel][19]=meshtopology[iel][15];

            }

        }
        return vtktopol;
    }




protected:

    std::ofstream fOutFile;
    NRmatrix<Doub> fdata;
    mesh * fmesh;
    std::string fFileName;
    Int fdim;
    std::vector<std::string> fVecNames;
    std::vector<std::string> fScalNames;


};
