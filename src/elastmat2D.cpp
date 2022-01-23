#include "elastmat2D.h"
#include "gridmesh.h"

using namespace std;
elastmat2D::elastmat2D (Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int elnodes )
{
    fyoung = young;
    fnu = nu;
    fbodyforce = bodyforce;
    fplanestress = planestress;
    fthickness = thickness;
	if(elnodes==3)
 	{
	  fOrder=1;
	  fshape = new shapetri(fOrder,1);
	}
	if(elnodes  ==4)
 	{
	  fOrder=1;
	  fshape = new shapequad(fOrder,1);
	}
	if(elnodes==6)
 	{
	  fOrder=2;
	  fshape = new shapetri(fOrder,1);
	}
	if(elnodes==8)
 	{
	  fOrder=2;
	  fshape = new shapequad(fOrder,1);
	}

}

elastmat2D::elastmat2D ( Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int elnodes, MatDoub  HHAT )
{
    fyoung = young;
    fnu = nu;
    fbodyforce = bodyforce;
    fplanestress = planestress;
    fthickness = thickness;
	if(elnodes==3)
 	{
	  fOrder=1;
	  fshape = new shapetri(fOrder,1);
	}
	if(elnodes  ==4)
 	{
	  fOrder=1;
	  fshape = new shapequad(fOrder,1);
	}
	if(elnodes==6)
 	{
	  fOrder=2;
	  fshape = new shapetri(fOrder,1);
	}
	if(elnodes==8)
 	{
	  fOrder=2;
	  fshape = new shapequad(fOrder,1);
	}
    fHHAT = HHAT;
}

elastmat2D::elastmat2D()
{

}

elastmat2D::~elastmat2D()
{


}


void elastmat2D::CalcStiff ( MatDoub &ek, MatDoub &ef, const MatDoub  &elcoords )
{

    MatDoub ptsweigths, ekt, eft;
    Doub xi, eta, w;
    Int nnodes = elcoords.nrows();
    ek.assign ( nnodes * 2, nnodes * 2, 0. );
    ef.assign ( nnodes * 2, 1, 0. );
    //shapequad shape = shapequad ( fOrder, 1 );
    fshape->pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();

    for ( Int ipt = 0; ipt < npts; ipt++ ) {
        xi = ptsweigths[ipt][0];
        eta = ptsweigths[ipt][1];
        w = ptsweigths[ipt][2];
        Contribute ( ekt, eft, xi, eta, w, elcoords );
        ek += ekt;
        ef += eft;
    }

}


void elastmat2D::Contribute ( MatDoub &ek, MatDoub &ef, Doub xi, Doub eta, Doub w, MatDoub elcoords )
{
    MatDoub psis, GradPsi, elcoordst, xycoords, Jac, InvJac ( 2, 2 ), GradPhi, B, BT, N, NT, psist, C, BC, BCS, stress ( 3, 1, 0. ), bodyforce ( 2, 1 ), temp, CS, KSt;
    bodyforce[0][0] = 0;
    bodyforce[1][0] = -fbodyforce;

    int type = 1;
    //shapequad objshapes ( fOrder, type );

    fshape->shapes ( psis, GradPsi, xi, eta );
    psis.Transpose ( psist );
    psist.Mult ( elcoords, xycoords );
    MatDoub hhat;
    if ( fhhatvel.nrows() != 0 ) {
        psist.Mult ( fhhatvel, hhat );
    }

    GradPsi.Mult ( elcoords, Jac );
    Int nnodes = psis.nrows();
    Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
    if ( DetJ <= 0 ) {

        std::cout << "\n DetJ < 0 " << endl;
        std::cout << "\n xi " << xi << endl;
        std::cout << "\n eta " << eta << endl;
        GradPsi.Print();
        elcoords.Print();
        GradPsi.Mult ( elcoords, Jac );
        xycoords.Print();

        psis.Print();

        fshape->shapes ( psis, GradPsi, xi, eta );
        return;
    }
    InvJac[0][0] = Jac[1][1] / DetJ;
    InvJac[0][1] = -Jac[0][1] / DetJ;
    InvJac[1][0] = -Jac[1][0] / DetJ;
    InvJac[1][1] = Jac[0][0] / DetJ;
    InvJac.Mult ( GradPsi, GradPhi );
    assembleBandN ( B, N, psis, GradPhi );
	cout << "elcoords,psis,gradps,Jac,xycords"<<endl;
	elcoords.Print();
	 psis.Print();
	 GradPsi.Print();
	 Jac.Print();
	 xycoords.Print();
    N.Transpose ( NT );
    B.Transpose ( BT );
    assembleConstitutiveMatrix ( C);
    BT.Mult ( C, BC );
    BC.Mult ( B, ek );


    ek *= w*DetJ*fthickness;
    BT.Mult ( stress, ef );
    NT.Mult ( bodyforce, temp );
    ef -= temp;
    ef *= w*DetJ;
}



void elastmat2D::assembleBandN ( MatDoub &B, MatDoub &N, const MatDoub &psis, const MatDoub &GradPhi )
{
    B.assign ( 3, psis.nrows() * 2, 0. );
    N.assign ( 2, psis.nrows() * 2, 0. );

    Int j = 0, k = 0;
    for ( Int i = 0; i < psis.nrows() * 2; i++ ) {
        if ( i % 2 == 0 ) {
            B[0][i] = GradPhi[0][j];
            B[1][i] = 0;
            B[2][i] = GradPhi[1][j];

            N[0][i] = psis[j][0];
            N[1][i] = 0;
            j++;
        } else {
            B[0][i] = 0;
            B[1][i] = GradPhi[1][k];
            B[2][i] = GradPhi[0][k];

            N[0][i] = 0;
            N[1][i] = psis[k][0];

            k++;
        }

    }
}
void elastmat2D::assembleConstitutiveMatrix ( MatDoub &ce )
{
    //plane stress
    Doub nusqr = fnu*fnu;
    Doub young = fyoung, nu = fnu;
    if ( fplanestress ==0) {

        ce.assign ( 3, 3, 0. );
        ce[0][0] = young / ( 1 - nusqr );ce[0][1] = nu*young / ( 1 - nusqr );ce[0][2] = 0.;
        ce[1][0] = nu*young / ( 1 - nusqr );ce[1][1] = young / ( 1 - nusqr );ce[1][2] = 0.;
        ce[2][0] = 0.;ce[2][1] = 0.;ce[2][2] = young / ( 1 - nusqr ) * ( 1-nu ) /2.;
    } else {

        Doub temp = young/ ( ( 1.+nu ) * ( 1.-2.*nu ) );
        ce.assign ( 3, 3, 0. );
        ce[0][0] = 1.-nu;ce[0][1] = nu;ce[0][2] = 0.;
        ce[1][0] = nu;ce[1][1] = 1.-nu;ce[1][2] = 0.;
        ce[2][0] = 0.;ce[2][1] = 0.;ce[2][2] = ( 1.-2.*nu ) /2.;
        ce*=temp;

    }

	

}



void elastmat2D::DirichletBC ( MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val )
{
    Int nodes = ids.size();
    Int sz = KG.nrows();
    Int fu = 0, cu = 0;
    for ( Int i = 0; i < nodes; i++ ) {
        Int pso = ids[i];
        if ( dir == 0 ) {
            for ( Int j = 0; j < sz; j++ ) {
                KG[2 * pso][j] = 0;
                KG[j][2 * pso] = 0;
            }
            KG[2 * pso][2 * pso] = 1;
            FG[2 * pso][0] = val;
        } else {
            for ( Int j = 0; j < sz; j++ ) {
                KG[2 * pso + 1][j] = 0;
                KG[j][2 * pso + 1] = 0;
            }
            KG[2 * pso + 1][2 * pso + 1] = 1;
            FG[2 * pso + 1][0] = val;
        }


    }

}
void elastmat2D::ContributeLineNewan ( MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val )
{
    MatDoub psis, gradpsis, ptsws;

   // int type = 1;
    //shapequad objshapes ( fOrder, type );

    fshape->pointsandweigths1D ( ptsws );
    Int npts = ptsws.nrows();

    Doub xi, w;
    MatDoub DetJ;
    for ( Int ipt = 0; ipt < npts; ipt++ ) {
        xi = ptsws[ipt][0];
        w = ptsws[ipt][2];
        fshape->shapes1D ( psis, gradpsis, xi );
        psis.Mult ( gradpsis, DetJ );

    }

}


// void elastmat2D::ComputeLoadVector( MatDoub &KG, MatDoub & FG, std::vector<int> ids,MatDoub & force )
// {
//     MatDoub psis, gradpsis, ptsws;
// 
//    // int type = 1;
//     //shapequad objshapes ( fOrder, type );
// 
//     fshape->pointsandweigths1D ( ptsws );
//     Int npts = ptsws.nrows();
// 
//     Doub xi, w;
//     MatDoub DetJ;
//     for ( Int ipt = 0; ipt < npts; ipt++ ) {
//         xi = ptsws[ipt][0];
//         w = ptsws[ipt][2];
//         fshape->shapes1D ( psis, gradpsis, xi );
//         psis.Mult ( gradpsis, DetJ );
// 
//     }
// 
// }


void elastmat2D::AssembleLoadvector ( NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology,Doub force )
{
    int type = 1;
    //shapequad objshapes ( fOrder, type );
	shape * objshapes = fshape;
    Int sz = 2*meshnodes.nrows();
    FG.assign ( sz, 1, 0. );
    //std::cout << "sz = " << sz << std::endl;
    MatDoub psis, gradpsis, ptsws, DetJ, psist;

    objshapes->pointsandweigths1D ( ptsws );

    Int  els=linetopology.nrows(),nodes= fOrder+1;
    Doub xi, w;
    for ( Int iel = 0; iel < els; iel++ ) {
        MatDoub xy ( 1, 2, 0. ),elcoords ( nodes,2,0. ), diff ( 1,2,0. ),temp ( 1, 2, 0. ), temp2 ( 1, 2, 0. ), xycoords ( 1, 2, 0. );
        Doub x=0., y=0.;
        for ( Int inode = 0; inode < nodes; inode++ ) {
            Doub xmesh = meshnodes[linetopology[iel][inode]][0];
            Doub ymesh = meshnodes[linetopology[iel][inode]][1];
            elcoords[inode][0] = xmesh;
            elcoords[inode][1] = ymesh;
        }
        NRvector<Doub> vec(2,0.),normal(2,0.);
		vec[0]=elcoords[nodes-1][0]-elcoords[0][0];
		vec[1]=elcoords[nodes-1][1]-elcoords[0][1];
		Doub elsize = sqrt(pow(vec[0],2)+pow(vec[1],2));
		normal[0]=-vec[1]/elsize;
		normal[1]=vec[0]/elsize;
		cout << " normal " <<endl;
		cout <<  normal[0] <<endl;
		cout << normal[1] <<endl;

		
		Doub jac=elsize/2.;
		cout << " jac " <<endl;
		cout << jac <<endl;
        Int npts = ptsws.nrows();
        MatDoub integral ( fOrder+1, 1, 0. ) ;
        for ( Int inode = 0; inode < fOrder + 1; inode++ ) {
			
			
			
            for ( Int ipt = 0; ipt < npts; ipt++ ) {
                xi = ptsws[ipt][0];
                w = ptsws[ipt][1];
                objshapes->shapes1D ( psis, gradpsis, xi );
                //integral[inode][0] += psis[inode][0] * jac*w;
				FG[ linetopology[iel][inode] * 2 ][0] += psis[inode][0]*force*normal[0]* jac*w ;
            	FG[ linetopology[iel][inode] * 2 + 1 ][0] += psis[inode][0] * force*normal[1]*jac*w ;

            }

            //FG[ linetopology[iel][inode] * 2 ][0] += integral[inode][0]*force*normal[0] ;
           // FG[ linetopology[iel][inode] * 2 + 1 ][0] += integral[inode][0]*force*normal[1] ;
        }

    }

    //FG.Print();
}


void elastmat2D::ComputeSolAndDSol ( mesh * inmesh,NRmatrix<Doub>&sol,NRmatrix<Doub>&dsol )
{
    Doub xi,eta,w;
    std::vector<std::vector< std::vector<Doub > > >  allcoords = inmesh->GetAllCoords();
    NRmatrix<Int> meshtopology = inmesh->GetMeshTopology();
    NRmatrix<Doub> nodalsol = GetSolution();
    Int meshnodes = inmesh->GetMeshNodes().nrows();

    NRmatrix<Doub>  elcoords, gradpsis;
    GetElCoords ( allcoords, 0, elcoords );
    Int nels = allcoords.size();
    Int nvars = 2;
    Int sz=nodalsol.nrows();

    //shapequad objshapes ( 2, 1 );
    NRmatrix<Doub> base=fshape->GetBaseNodes();

    sol.assign ( sz,2,0. );
    dsol.assign ( sz,2,0. );


    xi=0.;
    eta=0.;
    w=0.;

    NRmatrix<Doub>  psis,psist, GradPsi,elementdisplace;
    fshape->shapes ( psis, GradPsi, xi, eta );
    psis.Transpose ( psist );

    Int elnodes = elcoords.nrows();

    vector<Int> indexvec;
    NRmatrix<NRvector<Doub>> uglob;
    uglob.resize ( nels, elnodes );

    for ( int i = 0; i < nels; i++ ) {
        for ( int j = 0; j < elnodes; j++ ) {
            uglob[i][j].assign ( 2, 0. );
        }
    }

    for ( Int iel = 0; iel < nels; iel++ ) {
        for ( Int node = 0; node < elnodes; node++ ) {
            for ( int idof=0; idof<2; idof++ ) {
                uglob[iel][node][idof] = GetSolution() [2 * meshtopology[iel][node]+idof][0];
            }
        }
    }


    for ( Int iel = 0; iel < nels; iel++ ) {
        GetElCoords ( allcoords, iel, elcoords );
        elementdisplace.assign ( elnodes,2,0. );
        for ( Int i = 0; i < elnodes; i++ ) for ( Int j = 0; j < 2; j++ ) elementdisplace[i][j] = uglob[iel][i][j];
        for ( Int inode=0; inode<elnodes; inode++ ) {
            NRmatrix<Doub>  psis,GradPsi,Jac,InvJac ( 2,2 ),GradPhi,gradu;
            xi=base[inode][0];
            eta=base[inode][1];
            fshape->shapes ( psis, GradPsi, xi, eta );


            GradPsi.Mult ( elcoords, Jac );

            Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];


            InvJac[0][0] = Jac[1][1] / DetJ;
            InvJac[0][1] = -Jac[0][1] / DetJ;
            InvJac[1][0] = -Jac[1][0] / DetJ;
            InvJac[1][1] = Jac[0][0] / DetJ;

            InvJac.Mult ( GradPsi, GradPhi );

            GradPhi.Mult ( elementdisplace, gradu );

            Int index=meshtopology[iel][inode];

            sol[2 * index ][0]=  nodalsol[2 * index ][0];
            sol[2 * index +1][0]= nodalsol[2 * index +1][0];

            dsol[index * 2][0]= gradu[0][0] ;//dudx
            dsol[index * 2 + 1][0]= gradu[1][0] ;//dwdx

            dsol[index* 2][1]= gradu[0][1] ;//dudy
            dsol[index * 2 + 1][1]= gradu[1][1] ;//dwdy
        }

    }

}



