
//#include "mins_ndim.h"
#include <stdio.h>
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include "mesh.h"

mesh::mesh ( int dim,std::vector<std::vector< std::vector<Doub > > > &allcoords, NRmatrix<Doub> &meshnodes, NRmatrix<Int> &meshtopology )
{
    fallcoords = allcoords;
    fmeshnodes = meshnodes;
    fmeshtopology = meshtopology;
    fdim=dim;
}

mesh::mesh ( int dim,material *mat, std::vector<std::vector< std::vector<Doub > > >& allcoords, NRmatrix<Doub>& meshnodes, NRmatrix<Int>& meshtopology )
{
    fmaterial = mat;
    fallcoords = allcoords;
    fmeshnodes = meshnodes;
    fmeshtopology = meshtopology;
    fdim=dim;

}


mesh::mesh ( int dim,material* mat, std::vector<std::vector< std::vector<Doub > > >& allcoords, NRmatrix<Doub>& meshnodes, NRmatrix<Int>& meshtopology,MatDoub  HHAT )
{
    fmaterial = mat;
    fallcoords = allcoords;
    fmeshnodes = meshnodes;
    fmeshtopology = meshtopology;
    fHHAT = HHAT;
    fdim=dim;
}

mesh::mesh()
{
}

mesh::mesh ( mesh &copy )
{
    copy = mesh ( fdim,fallcoords, fmeshnodes, fmeshtopology );
}

mesh::~mesh()
{
    delete fmaterial;
}

std::vector<std::vector< std::vector<Doub > > >  mesh::GetAllCoords()
{
    int szdebug = fallcoords.size();
    return fallcoords;
}

MatDoub mesh::GetMeshNodes()
{
    return fmeshnodes;
}

MatInt mesh::GetMeshTopology()
{
    return fmeshtopology;
}


//material mesh::GetMaterial()
//{
//	return *fmaterial;
//}


void mesh::GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords )
{
    elcoords.assign ( allcoords[el].size(), 3, 0. );
    Int sz = allcoords[el].size();
    for ( Int j = 0; j <sz; j++ ) {
        for ( int i =0; i<3; i++ ) elcoords[j][i]=allcoords[el][j][i];

    }
}

MatDoub mesh::FindSolution ( VecDoub coord, MatDoub datatosearch )
{
    //MatDoub datatosearch [ data]

    MatDoub solout,returnsol ( 1,3,0. );
    Int el = 0;
    Doub xi = 0.;
    Doub eta = 0.;
    MatDoub psis, GradPsi, elcoords, psist, solel, xycoords, sol;
    shapequad shape = shapequad ( 2, 1 );
    shape.shapes ( psis, GradPsi, xi, eta );
    Doub l0 = 200;
    Doub l = 200;
    std::vector<int> possible;
    int counter = 0;
    bool  breaktrue = false;
    while ( counter<40 ) {
        for ( int iel = 0; iel < fallcoords.size(); iel++ ) {
            GetElCoords ( fallcoords, iel, elcoords );
            Int nodes = psis.nrows();
            Int nstatevars = 1;
            solel.assign ( nodes, 1, 0 );
            psis.Transpose ( psist );
            psist.Mult ( elcoords, xycoords );
            Doub dx = xycoords[0][0] - coord[0];
            Doub dy = xycoords[0][1] - coord[1];
            Doub dist = sqrt ( dx*dx + dy*dy );
            if ( dist<l0 ) {
                possible.push_back ( iel );
            }
        }
        //cout << "l0 = " << l0 << endl;
        l0 *= 0.5;
        counter++;
    }
   // std::cout << "in co  = " << endl;
   // coord.Print();
    el = possible.size() - 1;
    std::cout << "el = " <<possible[el] << endl;
    GetElCoords ( fallcoords, possible[el], elcoords );
   // std::cout << "elcoords = " <<  endl;
   // elcoords.Print();
    MatDoub elnodesol ( elcoords.nrows(), 1, 0. );
    for ( Int inode = 0; inode < elcoords.nrows(); inode++ ) {
        elnodesol[inode][0] = datatosearch[ fmeshtopology[possible[el]][inode] ][0];
    }
   // std::cout << "elnodesol = " << endl;
   // elnodesol.Print();

    VecDoub p ( 2, -0. );
    FuncdSearch func2;
    func2.set ( coord, elcoords );
    //Linemethod<FuncdSearch> line(func2);
    //Doub min = line.linmin();
    //std::cout << "min = " << min  << endl;


    //Powell<FuncdSearch> powell(func2);


    for ( Doub xi = -1.; xi <= 1.; xi += 0.0001 ) {
        for ( Doub eta = -1.; eta <= 1.; eta += 0.0001 ) {
            shape.shapes ( psis, GradPsi, xi, eta );
            psis.Transpose ( psist );
            psist.Mult ( elcoords, xycoords );
            Doub dx = xycoords[0][0] - coord[0];
            Doub dy = xycoords[0][1] - coord[1];
            Doub dist = sqrt ( dx*dx + dy*dy );
            if ( dist < 0.1 ) {
                p[0] = xi;
                p[1] = eta;
                break;
            }

        }
    }


    //p = powell.minimize(p);

    std::cout << "solu = " << endl;
    //p.Print();

    shape.shapes ( psis, GradPsi, p[0], p[1] );
    psis.Transpose ( psist );
    psist.Mult ( elcoords, xycoords );
    psist.Mult ( elnodesol, solout );
    returnsol[0][0] = xycoords[0][0];
    returnsol[0][1] = xycoords[0][1];
    returnsol[0][2] = solout[0][0];
    returnsol.Print();
    return returnsol;




}

MatDoub mesh::TransferSolution ( mesh &out, MatDoub datatosearch )
{
    Int outsize = out.fmeshnodes.nrows();//Um valor para cada n???
    MatDoub outsol ( outsize,1,0. ); // newval na malha coarse
    Int outnels = out.fallcoords.size();
    MatDoub outelcoords,returnsol;
    VecDoub co ( 2,0. );
    //preciso fornecer para a malha fina a coordenada do no da malha grossa para que a busca na malha fina seja realizada
    //portanto preciso iterar na malha grossa out
    for ( int iel = 0; iel < outnels; iel++ ) {
        GetElCoords ( fallcoords, iel, outelcoords );
        for ( Int inode = 0; inode < outelcoords.nrows(); inode++ ) {
            co[0] = outelcoords[inode][0];
            co[1] = outelcoords[inode][1];
            returnsol = FindSolution ( co, datatosearch ); //(acha coordena x y e val na malha fina a partir da coordenada
            //fornecido pela malha grossa)
            outsol[out.fmeshtopology[iel][inode]][0] = returnsol[0][2];
        }

    }
    return outsol;
}

/*
void mesh::Assemble(MatDoub& KG, MatDoub& Fint, MatDoub& Fbody)
{
	//std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
	//MatDoub meshnodes = fmesh.GetMeshNodes();
	//MatInt meshtopology = fmesh.GetMeshTopology();

	//cout << "all cc" << allcoords[0].size() <<endl;

	MatDoub ek, efint, efbody, elcoords, eltopology;
	GetElCoords(fallcoords, 0, elcoords);
	Int nnodes = fmeshnodes.nrows();
	Int rows = elcoords.nrows();
	Int sz = 2 * nnodes;
	Int cols = rows;
	KG.assign(sz, sz, 0.);
	Fint.assign(sz, 1, 0.);
	Fbody.assign(sz, 1, 0.);
	Int nels = fallcoords.size();

	//uglob = Table[
		//Table[{displacement[[2 topol[[k, j]] - 1]],
		//	displacement[[2 topol[[k, j]]]]}, { j, 1,
			//Length[topol[[k]]] }], { k, 1, nels }];
	NRmatrix<NRvector<Doub>> uglob;
	uglob.resize(nels, nnodes);
	for (int i = 0; i < nels; i++)
	{
		for (int j = 0; j < nnodes; j++) {
			uglob[i][j].assign(2, 0.);
		}
	}

	for (Int iel = 0; iel < nels; iel++)
	{
		for (Int node = 0; node < rows; node++)
		{
			uglob[iel][node][0] = fmaterial->GetSolution()[2 * fmeshtopology[iel][node]][0];
			uglob[iel][node][1] = fmaterial->GetSolution()[2 * fmeshtopology[iel][node] + 1][0];
		}
	}

	//	uglob.Print2();

	Int fu = 0;
	for (Int iel = 0; iel < nels; iel++)
	{
		if (fHHAT.nrows() != 0)
		{
			fhhatvel.resize(fHHAT.ncols());
			for (Int ivar = 0; ivar < fHHAT.ncols(); ivar++) {
				fhhatvel[ivar].assign(rows, 1, 0.);
				for (Int inode = 0; inode < rows; inode++)
				{
					fhhatvel[ivar][inode][0] = fHHAT[fmeshtopology[iel][inode]][ivar];
				}
			}
		}
		fmaterial->SetRandomField(fHHAT);
		fmaterial->SetRandomFieldLocal(fhhatvel);
		MatDoub elementdisplace(elcoords.nrows(), 2, 0.);
		for (Int i = 0; i < elcoords.nrows(); i++)for (Int j = 0; j < 2; j++)elementdisplace[i][j] = uglob[iel][i][j];
		GetElCoords(fallcoords, iel, elcoords);
		fmaterial->CacStiff(ek, efint, efbody, elcoords, elementdisplace);
		for (Int irow = 0; irow < rows; irow++)
		{
			Int rowglob = fmeshtopology[iel][irow];
			for (Int icol = 0; icol < cols; icol++)
			{
				Int colglob = fmeshtopology[iel][icol];
				KG[2 * rowglob + fu][2 * colglob + fu] += ek[2 * irow + fu][2 * icol + fu];
				KG[2 * rowglob + fu][2 * colglob + 1 + fu] += ek[2 * irow + fu][2 * icol + 1 + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + fu] += ek[2 * irow + 1 + fu][2 * icol + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + 1 + fu] += ek[2 * irow + 1 + fu][2 * icol + 1 + fu];

			}
			Fbody[2 * rowglob + fu][0] += efbody[2 * irow + fu][0];
			Fbody[2 * rowglob + 1 + fu][0] += efbody[2 * irow + 1 + fu][0];

			Fint[2 * rowglob + fu][0] += efint[2 * irow + fu][0];
			Fint[2 * rowglob + 1 + fu][0] += efint[2 * irow + 1 + fu][0];
		}

	}
	fmaterial->ResetCounter();
}

*/
void mesh::AssembleLinear ( MatDoub& KG, MatDoub& F )
{

    int ndof_per_node=fdim;
    MatDoub ek, ef,elcoords, eltopology;
    GetElCoords ( fallcoords, 0, elcoords );
    Int nnodes = fmeshnodes.nrows();
    Int rows = elcoords.nrows();
    Int sz = ndof_per_node * nnodes;
    Int cols = rows;
    KG.assign ( sz, sz, 0. );
    F.assign ( sz, 1, 0. );
    Int nels = fallcoords.size();

    NRmatrix<NRvector<Doub>> uglob;
    uglob.resize ( nels, rows );
    for ( int i = 0; i < nels; i++ ) {
        for ( int j = 0; j < rows; j++ ) {
            uglob[i][j].assign ( ndof_per_node, 0. );
        }
    }


    Int fu = 0;
    for ( Int iel = 0; iel < nels; iel++ ) {

        GetElCoords ( fallcoords, iel, elcoords );
        fmaterial->CalcStiff ( ek, ef, elcoords);
		cout << "ekekekek" <<endl;
		ek.Print();
        for ( Int irow = 0; irow < rows; irow++ ) {
            Int rowglob = fmeshtopology[iel][irow];
            for ( Int icol = 0; icol < cols; icol++ ) {
                Int colglob = fmeshtopology[iel][icol];
                int n=ndof_per_node;
                for ( int idof=0; idof<n; idof++ ) {
                    for ( int jdof=0; jdof<n; jdof++ ) {
                        Int linepos =n * rowglob - idof+n-1;
						Int colpos=n * colglob- jdof+n-1;
                        if ( linepos==4 && colpos==5 ) {
                            cout<< " | lineek = "<< n * irow- idof+n-1 ;
                            cout<< " | colek = "<< n * icol- jdof+n-1<< endl;
							cout<< " | ek[n * irow- idof+n-1][n * icol- jdof+n-1] = ";
							cout<< ek[n * irow- idof+n-1][n * icol- jdof+n-1] << endl;
                        }
                        KG[n * rowglob - idof+n-1][n * colglob- jdof+n-1] += ek[n * irow- idof+n-1][n * icol- jdof+n-1];
                    }

                }

            }

            int n=ndof_per_node;
            for ( int idof=0; idof<n; idof++ ) {
                F[n * rowglob - idof+n-1][0] += ef[n * irow- idof+n-1][0];
               
            }
        }

    }

}

void mesh::Assemble ( MatDoub& KG, MatDoub& Fint, MatDoub& Fbody )
{

    int ndof_per_node=fdim;
    MatDoub ek, efint, efbody, elcoords, eltopology;
    GetElCoords ( fallcoords, 0, elcoords );
    Int nnodes = fmeshnodes.nrows();
    Int rows = elcoords.nrows();
    Int sz = ndof_per_node * nnodes;
    Int cols = rows;
    KG.assign ( sz, sz, 0. );
    Fint.assign ( sz, 1, 0. );
    Fbody.assign ( sz, 1, 0. );
    Int nels = fallcoords.size();

    NRmatrix<NRvector<Doub>> uglob;
    uglob.resize ( nels, rows );
    for ( int i = 0; i < nels; i++ ) {
        for ( int j = 0; j < rows; j++ ) {
            uglob[i][j].assign ( ndof_per_node, 0. );
        }
    }


    for ( Int iel = 0; iel < nels; iel++ ) {
        for ( Int node = 0; node < rows; node++ ) {
            for ( int idof=0; idof<ndof_per_node; idof++ ) {
                uglob[iel][node][idof] = fmaterial->GetSolution() [ndof_per_node * fmeshtopology[iel][node]+idof][0];
            }
        }
    }

    Int fu = 0;
    for ( Int iel = 0; iel < nels; iel++ ) {
        if ( fHHAT.nrows() != 0 ) {
            fhhatvel.resize ( fHHAT.ncols() );
            for ( Int ivar = 0; ivar < fHHAT.ncols(); ivar++ ) {
                fhhatvel[ivar].assign ( rows, 1, 0. );
                for ( Int inode = 0; inode < rows; inode++ ) {
                    fhhatvel[ivar][inode][0] = fHHAT[fmeshtopology[iel][inode]][ivar];
                }
            }
        }
        fmaterial->SetRandomField ( fHHAT );
        fmaterial->SetRandomFieldLocal ( fhhatvel );
        MatDoub elementdisplace ( elcoords.nrows(),ndof_per_node, 0. );
        for ( Int i = 0; i < elcoords.nrows(); i++ ) for ( Int j = 0; j < ndof_per_node; j++ ) elementdisplace[i][j] = uglob[iel][i][j];
        GetElCoords ( fallcoords, iel, elcoords );
		//elcoords.Print();
		//DebugStop();
        fmaterial->CalcStiff ( ek, efint, efbody, elcoords, elementdisplace );

        for ( Int irow = 0; irow < rows; irow++ ) {
            Int rowglob = fmeshtopology[iel][irow];
            for ( Int icol = 0; icol < cols; icol++ ) {
                Int colglob = fmeshtopology[iel][icol];
                int n=ndof_per_node;
                for ( int idof=0; idof<n; idof++ ) {
                    for ( int jdof=0; jdof<n; jdof++ ) {
                        KG[n * rowglob - idof+n-1][n * colglob- jdof+n-1] += ek[n * irow- idof+n-1][n * icol- jdof+n-1];
                    }

                }

            }

            int n=ndof_per_node;
            for ( int idof=0; idof<n; idof++ ) {
                Fbody[n * rowglob - idof+n-1][0] += efbody[n * irow- idof+n-1][0];
                Fint[n * rowglob - idof+n-1][0]  += efint[n * irow- idof+n-1][0];
            }
        }

    }
    fmaterial->ResetCounter();
}

void mesh::Assemble ( MatDoub& Fint, MatDoub& Fbody )
{

    int ndof_per_node=fdim;
    MatDoub ek, efint, efbody, elcoords, eltopology;
    GetElCoords ( fallcoords, 0, elcoords );
    Int nnodes = fmeshnodes.nrows();
    Int rows = elcoords.nrows();
    Int sz = ndof_per_node * nnodes;
    Int cols = rows;

    Fint.assign ( sz, 1, 0. );
    Fbody.assign ( sz, 1, 0. );
    Int nels = fallcoords.size();

    NRmatrix<NRvector<Doub>> uglob;
    uglob.resize ( nels, rows );
    for ( int i = 0; i < nels; i++ ) {
        for ( int j = 0; j < rows; j++ ) {
            uglob[i][j].assign ( ndof_per_node, 0. );
        }
    }

    for ( Int iel = 0; iel < nels; iel++ ) {
        for ( Int node = 0; node < rows; node++ ) {
            for ( int idof=0; idof<ndof_per_node; idof++ ) {
                uglob[iel][node][idof] = fmaterial->GetSolution() [ndof_per_node * fmeshtopology[iel][node]+idof][0];
            }
        }
    }

    Int fu = 0;
    for ( Int iel = 0; iel < nels; iel++ ) {
        if ( fHHAT.nrows() != 0 ) {
            fhhatvel.resize ( fHHAT.ncols() );
            for ( Int ivar = 0; ivar < fHHAT.ncols(); ivar++ ) {
                fhhatvel[ivar].assign ( rows, 1, 0. );
                for ( Int inode = 0; inode < rows; inode++ ) {
                    fhhatvel[ivar][inode][0] = fHHAT[fmeshtopology[iel][inode]][ivar];
                }
            }
        }
        fmaterial->SetRandomField ( fHHAT );
        fmaterial->SetRandomFieldLocal ( fhhatvel );
        MatDoub elementdisplace ( elcoords.nrows(),ndof_per_node, 0. );
        for ( Int i = 0; i < elcoords.nrows(); i++ ) for ( Int j = 0; j < ndof_per_node; j++ ) elementdisplace[i][j] = uglob[iel][i][j];
        GetElCoords ( fallcoords, iel, elcoords );

        fmaterial->CalcStiff ( ek, efint, efbody, elcoords, elementdisplace );
        for ( Int irow = 0; irow < rows; irow++ ) {
            Int rowglob = fmeshtopology[iel][irow];
            int n=ndof_per_node;
            for ( int idof=0; idof<n; idof++ ) {
                Fbody[n * rowglob - idof+n-1][0] += efbody[n * irow- idof+n-1][0];
                Fint[n * rowglob - idof+n-1][0]  += efint[n * irow- idof+n-1][0];
            }
        }

    }
    fmaterial->ResetCounter();
}
void mesh::ComputeSolAndDSol ( NRmatrix<Doub> &sol, NRmatrix<Doub> &dsol )
{

    int ndof_per_node=fdim;
    MatDoub ek, efint, efbody, elcoords, eltopology;
    GetElCoords ( fallcoords, 0, elcoords );
    Int nnodes = fmeshnodes.nrows();
    Int rows = elcoords.nrows();
    Int sz = ndof_per_node * nnodes;
    Int cols = rows;

    Int nels = fallcoords.size();
    sol.assign ( nnodes*fdim,1,0. );
    dsol.assign ( nnodes*fdim,2,0. );

    NRmatrix<NRvector<Doub>> uglob;
    uglob.resize ( nels, rows );
    for ( int i = 0; i < nels; i++ ) {
        for ( int j = 0; j < rows; j++ ) {
            uglob[i][j].assign ( ndof_per_node, 0. );
        }
    }

    for ( Int iel = 0; iel < nels; iel++ ) {
        for ( Int node = 0; node < rows; node++ ) {
            for ( int idof=0; idof<ndof_per_node; idof++ ) {
                uglob[iel][node][idof] = fmaterial->GetSolution() [ndof_per_node * fmeshtopology[iel][node]+idof][0];
            }
        }
    }

    for ( Int iel = 0; iel < nels; iel++ ) {
        MatDoub elementdisplace ( elcoords.nrows(),ndof_per_node, 0. );
        for ( Int i = 0; i < elcoords.nrows(); i++ ) for ( Int j = 0; j < ndof_per_node; j++ ) elementdisplace[i][j] = uglob[iel][i][j];
        for ( Int irow = 0; irow < rows; irow++ ) {
            Int rowglob = fmeshtopology[iel][irow];
            int n=ndof_per_node;
            for ( int idof=0; idof<n; idof++ ) {
                sol[n * rowglob - idof+n-1][0] += efbody[n * irow- idof+n-1][0];
                dsol[n * rowglob - idof+n-1][0]  += efint[n * irow- idof+n-1][0];
            }
        }

    }
}

void mesh::Assemble ( SparseMatrix<double>  &KG, VectorXd &Fint, VectorXd &Fbody )
{
    //std::vector<std::vector< std::vector<Doub > > > allcoords = fmesh.GetAllCoords();
    //MatDoub meshnodes = fmesh.GetMeshNodes();
    //MatInt meshtopology = fmesh.GetMeshTopology();
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    //tripletList.reserve(estimation_of_entries);
    //cout << "all cc" << allcoords[0].size() <<endl;

    MatDoub ek, efint, efbody, elcoords, eltopology;

    GetElCoords ( fallcoords, 0, elcoords );
    Int nnodes = fmeshnodes.nrows();
    Int rows = elcoords.nrows();
    Int sz = 2 * nnodes;
    Int cols = rows;
    KG.resize ( sz, sz );

    Fint.resize ( sz );
    Fbody.resize ( sz );
    Int nels = fallcoords.size();
    MatrixXd KGt ( sz,sz );
    tripletList.reserve ( sz*sz );
    //uglob = Table[
    //Table[{displacement[[2 topol[[k, j]] - 1]],
    //	displacement[[2 topol[[k, j]]]]}, { j, 1,
    //Length[topol[[k]]] }], { k, 1, nels }];
    NRmatrix<NRvector<Doub>> uglob;
    uglob.resize ( nels, rows );
    for ( int i = 0; i < nels; i++ ) {
        for ( int j = 0; j < rows; j++ ) {
            uglob[i][j].assign ( 2, 0. );
        }
    }

    for ( Int iel = 0; iel < nels; iel++ ) {
        for ( Int node = 0; node < rows; node++ ) {
            uglob[iel][node][0] = fmaterial->GetSolution() [2 * fmeshtopology[iel][node]][0];
            uglob[iel][node][1] = fmaterial->GetSolution() [2 * fmeshtopology[iel][node] + 1][0];
        }
    }

    //	uglob.Print2();

    Int fu = 0;
    for ( Int iel = 0; iel < nels; iel++ ) {
        if ( fHHAT.nrows() != 0 ) {
            fhhatvel.resize ( fHHAT.ncols() );
            for ( Int ivar = 0; ivar < fHHAT.ncols(); ivar++ ) {
                fhhatvel[ivar].assign ( rows, 1, 0. );
                for ( Int inode = 0; inode < rows; inode++ ) {
                    fhhatvel[ivar][inode][0] = fHHAT[fmeshtopology[iel][inode]][ivar];
                }
            }
        }
        fmaterial->SetRandomField ( fHHAT );
        fmaterial->SetRandomFieldLocal ( fhhatvel );
        MatDoub elementdisplace ( elcoords.nrows(), 2, 0. );
        for ( Int i = 0; i < elcoords.nrows(); i++ ) for ( Int j = 0; j < 2; j++ ) elementdisplace[i][j] = uglob[iel][i][j];
        GetElCoords ( fallcoords, iel, elcoords );
        fmaterial->CalcStiff ( ek, efint, efbody, elcoords, elementdisplace );
        for ( Int irow = 0; irow < rows; irow++ ) {
            Int rowglob = fmeshtopology[iel][irow];
            for ( Int icol = 0; icol < cols; icol++ ) {
                Int colglob = fmeshtopology[iel][icol];

                /*double val1=0,val2=0,val3=0,val4=0;
                val1 = ek[2 * irow + fu][2 * icol + fu];
                val2 = ek[2 * irow + fu][2 * icol + 1 + fu];
                val3 = ek[2 * irow + 1 + fu][2 * icol + fu];
                val4 = ek[2 * irow + 1 + fu][2 * icol + 1 + fu];
                if(fabs(val1)>1.e-12)KG.coeffRef(2 * rowglob + fu,2 * colglob + fu) +=val1;
                if(fabs(val2)>1.e-12)KG.coeffRef(2 * rowglob + fu,2 * colglob + 1 + fu)+=val2;
                if(fabs(val3)>1.e-12)KG.coeffRef(2 * rowglob + 1 + fu,2 * colglob + fu)+=val3;
                if(fabs(val4)>1.e-12)KG.coeffRef(2 * rowglob + 1 + fu,2 * colglob + 1 + fu)+=val4;*/

                //tripletList.push_back(T(i,j,v_ij));

                KGt ( 2 * rowglob + fu,2 * colglob + fu ) += ek[2 * irow + fu][2 * icol + fu];
                KGt ( 2 * rowglob + fu,2 * colglob + 1 + fu ) += ek[2 * irow + fu][2 * icol + 1 + fu];
                KGt ( 2 * rowglob + 1 + fu,2 * colglob + fu ) += ek[2 * irow + 1 + fu][2 * icol + fu];
                KGt ( 2 * rowglob + 1 + fu,2 * colglob + 1 + fu ) += ek[2 * irow + 1 + fu][2 * icol + 1 + fu];

            }
            Fbody ( 2 * rowglob + fu ) += efbody[2 * irow + fu][0];
            Fbody ( 2 * rowglob + 1 + fu ) += efbody[2 * irow + 1 + fu][0];

            Fint ( 2 * rowglob + fu ) += efint[2 * irow + fu][0];
            Fint ( 2 * rowglob + 1 + fu ) += efint[2 * irow + 1 + fu][0];
        }

    }


    for ( int i=0; i<sz; i++ ) {
        for ( int j=0; j<sz; j++ ) {
            if ( fabs ( KGt ( i,j ) ) >1.e-12 ) {
                tripletList.push_back ( T ( i,j,KGt ( i,j ) ) );

            }
        }
    }

    KG.setFromTriplets ( tripletList.begin(), tripletList.end() );
    //KG.makeCompressed();
    fmaterial->ResetCounter();
}

void  mesh::FindIds ( NRvector<double> constcoorddata,NRvector<int> constcoord, std::vector<int> & ids )
{
    std::vector<std::vector< std::vector<Doub > > > allcoords =fallcoords;
    MatInt meshtopology = fmeshtopology;
//constcoorddata vector containig info about face to search id. It must contain any coodinate locate in the in the face
    MatDoub elcoords;
    int nels = allcoords.size();
    GetElCoords ( allcoords, 0, elcoords );
    Int nnodes = elcoords.nrows();
    int sum=0;
    //constcoord.size() = 1 face
    //constcoord.size() = 2 linha
    //constcoord.size() = 3 pontos

    std::vector<int> dirs;
    for ( int iconst=0; iconst<constcoord.size(); iconst++ ) {
        sum+=constcoord[iconst];
        if ( constcoord[iconst]==1 ) dirs.push_back ( iconst );

    }
    for ( Int iel = 0; iel < nels; iel++ ) {
        GetElCoords ( allcoords, iel, elcoords );
        for ( Int inode = 0; inode < nnodes; inode++ ) {

            if ( sum==1 ) {
                if ( fabs ( elcoords[inode][dirs[0]] - constcoorddata[dirs[0]] ) <1.e-4 ) {
                    ids.push_back ( meshtopology[iel][inode] );
                }
            } else if ( sum==2 ) {
                if ( fabs ( elcoords[inode][dirs[0]] - constcoorddata[dirs[0]] ) <1.e-4 && fabs ( elcoords[inode][dirs[1]] - constcoorddata[dirs[1]] ) <1.e-4 ) {
                    ids.push_back ( meshtopology[iel][inode] );
                }
            } else if ( sum==3 ) {
                //elcoords.Print();
                if ( fabs ( elcoords[inode][dirs[0]] - constcoorddata[dirs[0]] ) <1.e-4 && fabs ( elcoords[inode][dirs[1]] - constcoorddata[dirs[1]] ) <1.e-4 && fabs ( elcoords[inode][dirs[2]] - constcoorddata[dirs[2]] ) <1.e-4 ) {
                    ids.push_back ( meshtopology[iel][inode] );
                }
            }


        }


    }

    if ( ids.size() >0 ) {
        sort ( ids.begin(), ids.end() );
        ids.erase ( unique ( ids.begin(), ids.end() ), ids.end() );
    }
}



