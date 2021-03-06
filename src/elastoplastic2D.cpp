#include "elastoplastic2D.h"
#include "vonmises.h"
#include "druckerprager.h"
#include "mohrcoulomb.h"
template <class YC>
elastoplastic2D<YC>::elastoplastic2D ( Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int elnodes )
{
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

//template <class YC>
//elastoplastic2D<YC>::elastoplastic2D(Doub young, Doub nu, Doub sigy, Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int order, NRmatrix<Doub>   HHAT)
//{
//	fbodyforce = bodyforce;
//	fplanestress = planestress;
//	fthickness = thickness;
//	fOrder = order;
//	fHHAT = HHAT;

//}

template <class YC>
elastoplastic2D<YC>::elastoplastic2D ( Doub thickness, NRmatrix<Doub>  bodyforce, Int planestress, Int elnodes, NRmatrix<Doub>   HHAT )
{

    fbodyforce = bodyforce;
    fplanestress = planestress;
    fthickness = thickness;
    fHHAT = HHAT;
	
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

template <class YC>
elastoplastic2D<YC>::elastoplastic2D ( elastoplastic2D & copy )
{
}

template <class YC>
elastoplastic2D<YC>::elastoplastic2D()
{
}



template <class YC>
elastoplastic2D<YC>::~elastoplastic2D()
{


    /*  delete  fbodyforce;
    delete fplanestress;
    delete fthickness;
    delete fOrder;
    delete fHHAT;
    deletefhhatvel;
    delete fYC;

    delete  fdisplace;
    delete fepspvec;
    delete fepspsolitern;
    delete fglobalcounter;*/
}

template <class YC>
void elastoplastic2D<YC>::SetMemory ( Int nglobalpts, Int systemsize )
{
    fdisplace.assign ( systemsize,1, 0. );

    fepspvec.assign ( nglobalpts, 0. );
    fepspsolitern.assign ( nglobalpts, 0. );
    fglobalcounter = 0;

}

template <class YC>
void elastoplastic2D<YC>::ResetMemory()
{
    Int systemsize = fdisplace.nrows();
    fdisplace.assign ( systemsize, 1, 0. );

    Int nglobalpts = fepspvec.size();
    fepspvec.assign ( nglobalpts, 0. );
    fepspsolitern.assign ( nglobalpts, 0. );
    fglobalcounter = 0;
}

template <class YC>
void elastoplastic2D<YC>::SetRandomField ( NRmatrix<Doub>  hhat )
{
    fHHAT = hhat;
}
template <class YC>
void elastoplastic2D<YC>::SetRandomFieldLocal ( NRvector<NRmatrix<Doub> >  hhatvel )
{
    fhhatvel = hhatvel;
}


template <class YC>
void elastoplastic2D<YC>::UpdateDisplacement ( NRmatrix<Doub>  displace )
{
    fdisplace = displace;
}

template <class YC>
void elastoplastic2D<YC>::UpdateDisplacement ( VectorXd displace )
{
    double sz =  displace.size();
    fdisplace.resize ( sz,1 );
    for ( int i=0; i<sz; i++ ) fdisplace[i][0] = displace ( i );
}

template <class YC>
void elastoplastic2D<YC>::UpdatePlasticStrain()
{
    fepspsolitern = fepspvec;
    fglobalcounter = 0;

}

template <class YC>
void elastoplastic2D<YC>::UpdateBodyForce ( NRmatrix<Doub>  newbodyforce )
{
    fbodyforce = newbodyforce;
}

template <class YC>
void elastoplastic2D<YC>::ResetPlasticStrain()
{


    //NRvector<TensorDoub> fepspvec, fepspsolitern;
    Int rows = fepspvec.size();
    for ( Int i = 0; i < rows; i++ ) {
        for ( Int j = 0; j < 6; j++ ) {
            fepspvec[i].XX() = 0. ;
            fepspvec[i].YY() = 0. ;
            fepspvec[i].ZZ() = 0. ;
            fepspvec[i].XZ() = 0. ;
            fepspvec[i].YZ() = 0. ;
            fepspvec[i].XY() = 0. ;

            fepspsolitern[i].XX() = 0.;
            fepspsolitern[i].YY() = 0.;
            fepspsolitern[i].ZZ() = 0.;
            fepspsolitern[i].XZ() = 0.;
            fepspsolitern[i].YZ() = 0.;
            fepspsolitern[i].XY() = 0.;
        }
    }

}

//Return the last converged plastic deformation
template <class YC>
NRvector<NRtensor<Doub> >   elastoplastic2D<YC>::GetPlasticStrain()
{
    return fepspsolitern;
}

template <class YC>
void elastoplastic2D<YC>::ResetMat()
{
    fglobalcounter = 0;
    ResetDisplacement();
    ResetPlasticStrain();
    ResetMemory();

}

template <class YC>
void elastoplastic2D<YC>::ResetDisplacement()
{
    Int rows = fdisplace.nrows();
    fdisplace.assign ( rows, 1,0. );
}



template <class YC>
void elastoplastic2D<YC>::Contribute ( NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody,NRvector<Doub> ptsw, NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace )
{
    Doub xi,eta,w;
    xi=ptsw[0];
    eta=ptsw[1];
    w=ptsw[2];
    int type = 1;
    //shapequad objshapes ( fOrder, type );
	shape *objshapes = fshape;

    NRmatrix<Doub>  psis, GradPsi, elcoordst, xycoords, Jac, InvJac ( 2, 2 ), GradPhi, B, BT, N, NT, psist, C, BC, BCS, stress ( 3, 1, 0. ), temp, CS, KSt;
    objshapes->shapes ( psis, GradPsi, xi, eta );
    psis.Transpose ( psist );
    psist.Mult ( elcoords, xycoords );
    NRvector<NRmatrix<Doub> > hhat ( fhhatvel.size() );

    if ( fhhatvel.size() != 0 ) {
        for ( Int ivar = 0; ivar < fhhatvel.size(); ivar++ ) {
            psist.Mult ( fhhatvel[ivar], hhat[ivar] );
        }

    }
    GradPsi.Mult ( elcoords, Jac );
    Int nnodes = psis.nrows();
    Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
    if ( DetJ <= 0 ) {

        std::cout << "\n DetJ < 0 " << endl;
        std::cout << "\n xi " << xi << endl;
        std::cout << "\n eta " << eta << endl;

        std::cout << "\n w " << w << endl;
        Jac.Print();
        GradPsi.Print();
        elcoords.Print();
        GradPsi.Mult ( elcoords, Jac );
        xycoords.Print();

        psis.Print();

        objshapes->shapes ( psis, GradPsi, xi, eta );
        std::cout << "\n Det < 0 " <<std::endl;
        DebugStop();
    }
    InvJac[0][0] = Jac[1][1] / DetJ;
    InvJac[0][1] = -Jac[0][1] / DetJ;
    InvJac[1][0] = -Jac[1][0] / DetJ;
    InvJac[1][1] = Jac[0][0] / DetJ;
    InvJac.Mult ( GradPsi, GradPhi );
    NRmatrix<Doub>  gradu;
    GradPhi.Mult ( eldisplace, gradu );

    //{ {dudx, dudy}, { dvdx, dvdy }} = gradprevsol;
    Doub ex = gradu[0][0];// dudx;
    Doub ey = gradu[1][1];
    Doub exy = ( gradu[0][1] + gradu[1][0] );


    NRtensor<Doub>   epst ( 0. ), epsp ( 0. ), projstress ( 0. ), projstrain ( 0. ), epspeint ( 0. );

    NRmatrix<Doub>  Dep,Dept;
    Doub  projgamma=0.;
    epst.XX() = ex;
    epst.YY() = ey;
    epst.XY() = exy;
    //epsp = epspsoliternGLOBAL[globalcounter];
    epsp = fepspsolitern[fglobalcounter];

    //std::cout << "epst = " << std::endl;
    //epst.Print();
    //std::cout << "epsp = " << std::endl;
    //epsp.Print();
    //cout << "\n antes c  = " << fYC.GetCoes() << endl;
    if ( fhhatvel.size() != 0 ) {
        for ( Int ivar = 0; ivar < fhhatvel.size(); ivar++ ) {
            psist.Mult ( fhhatvel[ivar], hhat[ivar] );
        }
        fYC.updateatributes ( hhat );
    }
    //cout << "\n c  = " << fYC.GetCoes() <<endl;

    fYC.closestpointproj ( epst,epsp,projstress,projstrain,Dept,projgamma );

    //cout << "prj stress "<< endl;
    //projstress.Print();

    Dep.assign ( 3, 3, 0. );
    if ( Dept.nrows() >3 ) {
        Dep[0][0] = Dept[0][0];
        Dep[0][1] = Dept[0][1];
        Dep[0][2] = Dept[0][5];
        Dep[1][0] = Dept[1][0];
        Dep[1][1] = Dept[1][1];
        Dep[1][2] = Dept[1][5];
        Dep[2][0] = Dept[5][0];
        Dep[2][1] = Dept[5][1];
        Dep[2][2] = Dept[5][5];
    } else {
        Dep=Dept;
        //Dep.Print();
    }
    if ( fhhatvel.size() != 0 ) {
        fYC.restoreoriginalatributes();
    }
    //epspeint = epst - (epse);
    epspeint = epst;
    epspeint = epspeint - projstrain;

    fepspvec[fglobalcounter] = epspeint;
    //epspvecGLOBAL[globalcounter] = epspeint;

    fglobalcounter++;
    assembleBandN ( B, N, psis, GradPhi );
    N.Transpose ( NT );
    B.Transpose ( BT );

    //Dep.Print();

    //BT.Print();
    // BC=(BT*Dep)*B;
    BT.Mult ( Dep, BC );
    BC.Mult ( B, ek );




    stress[0][0] = projstress.XX();
    stress[1][0] = projstress.YY();
    stress[2][0] = projstress.XY();
    ek *= w*DetJ*fthickness;
    //std::cout << " stress "<< std::endl;
    //stress.Print();
    //ek.Print();
    //ef = (Transpose[BB].stress) weight DetJ;
    //ef2 = (Transpose[NShapes].{0, -bodyforce}) weight DetJ;
    BT.Mult ( stress, efint );
    NT.Mult ( fbodyforce, efbody );
    efint *= w*DetJ;
    efbody *= w*DetJ;
}
template <class YC>
void elastoplastic2D<YC>::CalcStiff ( NRmatrix<Doub>  &ek, NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, const NRmatrix<Doub>   &elcoords,NRmatrix<Doub>  eldisplace )
{
    MatDoub ptsweigths, ekt, eftint,eftbody;
    Doub xi, eta, w;
    Int nnodes = elcoords.nrows();
    ek.assign ( nnodes * 2, nnodes * 2, 0. );
    efint.assign ( nnodes * 2, 1, 0. );
    efbody.assign ( nnodes * 2, 1, 0. );
   // shapequad shape = shapequad ( fOrder, 1 );

    fshape->pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();

    NRvector<Doub> ptsw ( 3 );
    for ( Int ipt = 0; ipt < npts; ipt++ ) {

        ptsw[0]=ptsweigths[ipt][0];
        ptsw[1]= ptsweigths[ipt][1];
        ptsw[2]= ptsweigths[ipt][2];
        //std::cout << "integration point :" << " xi = " << xi << "  eta = "<<  eta << std::endl;
        Contribute ( ekt, eftint, eftbody, ptsw, elcoords, eldisplace );
        ek += ekt;
        efint += eftint;
        efbody += eftbody;
    }
}

template <class YC>
void elastoplastic2D<YC>::Contribute (  NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody,NRvector<Doub> ptsw, NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace )
{
    Doub xi,eta,w;
    xi=ptsw[0];
    eta=ptsw[1];
    w=ptsw[2];
    int type = 1;
    //shapequad objshapes ( fOrder, type );
	shape *objshapes = fshape;

    NRmatrix<Doub>  psis, GradPsi, elcoordst, xycoords, Jac, InvJac ( 2, 2 ), GradPhi, B, BT, N, NT, psist, C, BC, BCS, stress ( 3, 1, 0. ), temp, CS, KSt;
    objshapes->shapes ( psis, GradPsi, xi, eta );
    psis.Transpose ( psist );
    psist.Mult ( elcoords, xycoords );
    NRvector<NRmatrix<Doub> > hhat ( fhhatvel.size() );

    if ( fhhatvel.size() != 0 ) {
        for ( Int ivar = 0; ivar < fhhatvel.size(); ivar++ ) {
            psist.Mult ( fhhatvel[ivar], hhat[ivar] );
        }

    }
    GradPsi.Mult ( elcoords, Jac );
    Int nnodes = psis.nrows();
    Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
    if ( DetJ <= 0 ) {

        std::cout << "\n DetJ < 0 " << endl;
        std::cout << "\n xi " << xi << endl;
        std::cout << "\n eta " << eta << endl;

        std::cout << "\n w " << w << endl;
        Jac.Print();
        GradPsi.Print();
        elcoords.Print();
        GradPsi.Mult ( elcoords, Jac );
        xycoords.Print();

        psis.Print();

        objshapes->shapes ( psis, GradPsi, xi, eta );
        std::cout << "\n Det < 0 " <<std::endl;
        DebugStop();
    }
    InvJac[0][0] = Jac[1][1] / DetJ;
    InvJac[0][1] = -Jac[0][1] / DetJ;
    InvJac[1][0] = -Jac[1][0] / DetJ;
    InvJac[1][1] = Jac[0][0] / DetJ;
    InvJac.Mult ( GradPsi, GradPhi );
    NRmatrix<Doub>  gradu;
    GradPhi.Mult ( eldisplace, gradu );

    //{ {dudx, dudy}, { dvdx, dvdy }} = gradprevsol;
    Doub ex = gradu[0][0];// dudx;
    Doub ey = gradu[1][1];
    Doub exy = ( gradu[0][1] + gradu[1][0] );


    NRtensor<Doub>   epst ( 0. ), epsp ( 0. ), projstress ( 0. ), projstrain ( 0. ), epspeint ( 0. );

    NRmatrix<Doub>  Dep,Dept;
    Doub  projgamma=0.;
    epst.XX() = ex;
    epst.YY() = ey;
    epst.XY() = exy;
    //epsp = epspsoliternGLOBAL[globalcounter];
    epsp = fepspsolitern[fglobalcounter];

    //std::cout << "epst = " << std::endl;
    //epst.Print();
    //std::cout << "epsp = " << std::endl;
    //epsp.Print();
    //cout << "\n antes c  = " << fYC.GetCoes() << endl;
    if ( fhhatvel.size() != 0 ) {
        for ( Int ivar = 0; ivar < fhhatvel.size(); ivar++ ) {
            psist.Mult ( fhhatvel[ivar], hhat[ivar] );
        }
        fYC.updateatributes ( hhat );
    }
    //cout << "\n c  = " << fYC.GetCoes() <<endl;

    fYC.closestpointproj ( epst,epsp,projstress,projstrain,Dept,projgamma );

    //cout << "prj stress "<< endl;
    //projstress.Print();

    Dep.assign ( 3, 3, 0. );
    if ( Dept.nrows() >3 ) {
        Dep[0][0] = Dept[0][0];
        Dep[0][1] = Dept[0][1];
        Dep[0][2] = Dept[0][5];
        Dep[1][0] = Dept[1][0];
        Dep[1][1] = Dept[1][1];
        Dep[1][2] = Dept[1][5];
        Dep[2][0] = Dept[5][0];
        Dep[2][1] = Dept[5][1];
        Dep[2][2] = Dept[5][5];
    } else {
        Dep=Dept;
        //Dep.Print();
    }
    if ( fhhatvel.size() != 0 ) {
        fYC.restoreoriginalatributes();
    }
    //epspeint = epst - (epse);
    epspeint = epst;
    epspeint = epspeint - projstrain;

    fepspvec[fglobalcounter] = epspeint;
    //epspvecGLOBAL[globalcounter] = epspeint;

    fglobalcounter++;
    assembleBandN ( B, N, psis, GradPhi );
    N.Transpose ( NT );
    B.Transpose ( BT );


    stress[0][0] = projstress.XX();
    stress[1][0] = projstress.YY();
    stress[2][0] = projstress.XY();

	
    BT.Mult ( stress, efint );
    NT.Mult ( fbodyforce, efbody );
	
    efint *= w*DetJ;
    efbody *= w*DetJ;
}

template <class YC>
void elastoplastic2D<YC>::CalcStiff ( NRmatrix<Doub>  &efint, NRmatrix<Doub>  &efbody, const NRmatrix<Doub>   &elcoords,NRmatrix<Doub>  eldisplace )
{
    MatDoub ptsweigths,  eftint,eftbody;
    Doub xi, eta, w;
    Int nnodes = elcoords.nrows();
    efint.assign ( nnodes * 2, 1, 0. );
    efbody.assign ( nnodes * 2, 1, 0. );

    fshape->pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();

    NRvector<Doub> ptsw ( 3 );
    for ( Int ipt = 0; ipt < npts; ipt++ ) {
        ptsw[0]=ptsweigths[ipt][0];
        ptsw[1]= ptsweigths[ipt][1];
        ptsw[2]= ptsweigths[ipt][2];
        //std::cout << "integration point :" << " xi = " << xi << "  eta = "<<  eta << std::endl;
        Contribute (eftint, eftbody, ptsw, elcoords, eldisplace );
        efint += eftint;
        efbody += eftbody;
    }
}

template <class YC>
void elastoplastic2D<YC>::assembleBandN ( MatrixXd &B, MatrixXd  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi )
{
    B.resize ( 3,psis.nrows() * 2 );
    N.resize ( 2,psis.nrows() * 2 );


    Int j = 0, k = 0;
    for ( Int i = 0; i < psis.nrows() * 2; i++ ) {
        if ( i % 2 == 0 ) {
            B ( 0,i ) = GradPhi[0][j];
            B ( 1,i ) = 0;
            B ( 2,i ) = GradPhi[1][j];

            N ( 0,i ) = psis[j][0];
            N ( 1,i ) = 0;
            j++;
        } else {
            B ( 0,i ) = 0;
            B ( 1,i ) = GradPhi[1][k];
            B ( 2,i ) = GradPhi[0][k];

            N ( 0,i ) = 0;
            N ( 1,i ) = psis[k][0];

            k++;
        }

    }

}


template <class YC>
void elastoplastic2D<YC>::assembleBandN ( NRmatrix<Doub>  &B, NRmatrix<Doub>  &N, const NRmatrix<Doub>  &psis, const NRmatrix<Doub>  &GradPhi )
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


template <class YC>
void elastoplastic2D<YC>::GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords )
{
    elcoords.assign ( allcoords[el].size(), 2, 0. );
    for ( Int j = 0; j < allcoords[el].size(); j++ ) {
        Doub x = allcoords[el][j][0];
        Doub y = allcoords[el][j][1];
        elcoords[j][0] = x;
        elcoords[j][1] = y;
    }
}
template <class YC>
void elastoplastic2D<YC>::DirichletBC ( NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val )
{
    Int nodes = ids.size();
    Int sz = KG.nrows();

    for ( Int i = 0; i < nodes; i++ ) {
        Int pso = ids[i];
        if ( dir == 0 ) {
            for ( Int j = 0; j < sz; j++ ) {
                KG[2 * pso][j] = 0.;
                KG[j][2 * pso] = 0.;
            }
            KG[2 * pso][2 * pso] = 1.;
            FG[2 * pso][0] = val;
        } else {
            for ( Int j = 0; j < sz; j++ ) {
                KG[2 * pso + 1][j] = 0.;
                KG[j][2 * pso + 1] = 0.;
            }
            KG[2 * pso + 1][2 * pso + 1] = 1.;
            FG[2 * pso + 1][0] = val;
        }


    }
}

template <class YC>
void elastoplastic2D<YC>::DirichletBC ( SparseMatrix<double> & KG, VectorXd& FG, std::vector<int> ids, Int  dir, Int val )
{
    Int nodes = ids.size();
    Int sz = KG.rows();
    for ( Int i = 0; i < nodes; i++ ) {
        Int pso = ids[i];
        if ( dir == 0 ) {
            for ( Int j = 0; j < sz; j++ ) {
                KG.coeffRef ( 2 * pso,j ) = 0;
                KG.coeffRef ( j,2 * pso ) = 0;
            }
            KG.coeffRef ( 2 * pso,2 * pso ) = 1;
            FG ( 2 * pso ) = val;
        } else {
            for ( Int j = 0; j < sz; j++ ) {
                KG.coeffRef ( 2 * pso + 1,j ) = 0;
                KG.coeffRef ( j,2 * pso + 1 ) = 0;
            }
            KG.coeffRef ( 2 * pso + 1,2 * pso + 1 ) = 1;
            FG ( 2 * pso + 1 ) = val;
        }


    }
}



template <class YC>
void elastoplastic2D<YC>::ContributeLineNewan ( NRmatrix<Doub>  &KG, NRmatrix<Doub>  & FG, std::vector<int> ids, Int  dir, Int val )
{
    NRmatrix<Doub>  psis, gradpsis, ptsws;

    int type = 1;
    //shapequad objshapes ( fOrder, type );

	shape * objshapes = fshape;
    objshapes->pointsandweigths1D ( ptsws );
    Int npts = ptsws.nrows();

    Doub xi, w;
    NRmatrix<Doub>  DetJ;

    for ( Int ipt = 0; ipt < npts; ipt++ ) {
        xi = ptsws[ipt][0];
        w = ptsws[ipt][2];
        objshapes->shapes1D ( psis, gradpsis, xi );
        psis.Mult ( gradpsis, DetJ );

    }

}

template <class YC>
void elastoplastic2D<YC>::ContributeCurvedLine ( NRmatrix<Doub>  &KG, NRmatrix<Doub>  &FG, NRmatrix<Doub>  meshnodes, MatInt linetopology,Doub force )
{
    int type = 1;
    //shapequad objshapes ( fOrder, type );
	shape * objshapes = fshape;
    Int sz = 2*meshnodes.nrows();
    FG.assign ( sz, 1, 0. );
    //std::cout << "sz = " << sz << std::endl;
    MatDoub psis, gradpsis, ptsws, DetJ, psist;

    objshapes->pointsandweigths1D ( ptsws );

    Int npts = 1000,els=linetopology.nrows(),nodes= fOrder+1;
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
        Doub delta = 2. / npts;
        xi = -1.;
        for ( Int ipt = 0; ipt <= npts+1; ipt++ ) {
            objshapes->shapes1D ( psis, gradpsis, xi );
            psis.Transpose ( psist );
            temp = xycoords;
            psist.Mult ( elcoords, xycoords );
            if ( ipt > 0 ) {
                temp2 = xycoords;
                temp2 -= temp;
                diff += temp2;
            }
            xi += delta;
        }
        Doub jac=sqrt ( diff[0][0]* diff[0][0] + diff[0][1] * diff[0][1] ) /2.;
        Int npts = ptsws.nrows();
        MatDoub integral ( fOrder+1, 1, 0. ) ;
        for ( Int inode = 0; inode < fOrder + 1; inode++ ) {
            for ( Int ipt = 0; ipt < npts; ipt++ ) {
                xi = ptsws[ipt][0];
                w = ptsws[ipt][1];
                objshapes->shapes1D ( psis, gradpsis, xi );
                integral[inode][0] += psis[inode][0] * force*jac*w;

            }

            VecDoub normal ( 2,0. );
            Doub norm = sqrt ( elcoords[inode][0] * elcoords[inode][0] + elcoords[inode][1] * elcoords[inode][1] );
            normal[0]= elcoords[inode][0] / norm;
            normal[1] = elcoords[inode][1] / norm;
            //normal.Print();
            //std::cout << "iel " << iel << std::endl;
            //std::cout << "inode " << inode << std::endl;
            //std::cout << " linetopology[iel][inode] * 2 " <<linetopology[iel][inode] * 2 << std::endl;
            FG[ linetopology[iel][inode] * 2 ][0] += integral[inode][0]*normal[0];
            FG[ linetopology[iel][inode] * 2 + 1 ][0] += integral[inode][0]*normal[1];
        }

    }

    //FG.Print();
}

template <class YC>
void elastoplastic2D<YC>::SolPt ( mesh * inmesh, const Int &el, const  NRmatrix<Doub>  &nodalsol, const Doub &xi, const Doub &eta, NRmatrix<Doub>  &xycoords, NRmatrix<Doub>  &sol )
{
    const std::vector<std::vector< std::vector<Doub > > > allcoords = inmesh->GetAllCoords();
    const MatInt meshtopology  = inmesh->GetMeshTopology();

    MatDoub psis, GradPsi, elcoords, psist, solel;
    GetElCoords ( allcoords, el, elcoords );

    int type = 1;
    //shapequad objshapes ( fOrder, type );

    fshape->shapes ( psis, GradPsi, xi, eta );
    Int nodes = psis.nrows();
    Int nstatevars = 1;
    solel.assign ( nodes, 1, 0 );
    psis.Transpose ( psist );
    psist.Mult ( elcoords, xycoords );

    for ( Int inode = 0; inode < nodes; inode++ ) {
        solel[inode][0] = nodalsol[meshtopology[el][inode]][0];
    }
    psist.Mult ( solel, sol );
}

template <class YC>
//void elastoplastic2D<YC>::PostProcess(mesh & inmesh, const MatDoub & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly)
void elastoplastic2D<YC>::PostProcess ( std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnode, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly )
{
    int type = 1;
    //shapequad objshapes ( fOrder, type );
    //std::vector<std::vector< std::vector<Doub > > > allcoords = inmesh.GetAllCoords();
    //MatDoub meshnodes = inmesh.GetMeshNodes();
    //MatInt meshtopology = inmesh.GetMeshTopology();
    MatDoub elcoords, eltopology, psis, gradpsis, xycoords, psist;
    GetElCoords ( allcoords, 0, elcoords );
    Int rows = elcoords.nrows();
    Int cols = rows;
    Int nels = allcoords.size();
    Doub refine = 0.1;

    for ( Int iel = 0; iel < nels; iel++ ) {
        GetElCoords ( allcoords, iel, elcoords );
        for ( Doub xi = -1.; xi < 1 - refine; xi += refine ) {
            std::vector<double> sol ( 3 );
            for ( Doub eta = -1.; eta < 1 - refine; eta += refine ) {
                Doub approx = 0., approy = 0.;
               fshape->shapes ( psis, gradpsis, xi, eta );
                psis.Transpose ( psist );
                psist.Mult ( elcoords, xycoords );
                sol[0] = xycoords[0][0];
                sol[1] = xycoords[0][1];
                for ( Int inode = 0; inode < elcoords.nrows(); inode++ ) {
                    approx += psis[inode][0] * nodalsol[meshtopology[iel][inode] * 2][0];
                    approy += psis[inode][0] * nodalsol[meshtopology[iel][inode] * 2 + 1][0];
                }
                sol[2] = approx;
                solx.push_back ( sol );
                sol[2] = approy;
                soly.push_back ( sol );

            }
        }
    }
}
template <class YC>
//void elastoplastic2D<YC>::PostProcess(mesh & inmesh,Int var,const MatDoub & nodalsol, std::vector<std::vector<double>> &sol)
void elastoplastic2D<YC>::PostProcess ( std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnode, MatInt meshtopology, Int var, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol )
{
    int type = 1;
   // shapequad objshapes ( fOrder, type );
    //std::vector<std::vector< std::vector<Doub > > > allcoords = inmesh.GetAllCoords();
    //MatDoub meshnodes = inmesh.GetMeshNodes();
    //MatInt meshtopology = inmesh.GetMeshTopology();
    NRmatrix<Doub>  elcoords, eltopology, psis, gradpsis, xycoords, psist;
    GetElCoords ( allcoords, 0, elcoords );
    Int rows = elcoords.nrows();
    Int cols = rows;
    Int nels = allcoords.size();
    Doub refine = 0.05;

    for ( Int iel = 0; iel < nels; iel++ ) {
        GetElCoords ( allcoords, iel, elcoords );
        for ( Doub xi = -1.; xi < 1 - refine; xi += refine ) {
            std::vector<double> soli ( 3 );
            for ( Doub eta = -1.; eta < 1 - refine; eta += refine ) {
                Doub approx = 0., approy = 0.;
                fshape->shapes ( psis, gradpsis, xi, eta );
                psis.Transpose ( psist );
                psist.Mult ( elcoords, xycoords );
                soli[0] = xycoords[0][0];
                soli[1] = xycoords[0][1];
                for ( Int inode = 0; inode < elcoords.nrows(); inode++ ) {
                    approx += psis[inode][0] * nodalsol[meshtopology[iel][inode]][var];
                }
                soli[2] = approx;
                sol.push_back ( soli );

            }
        }
    }
}


template <class YC>
//void elastoplastic2D<YC>::PostProcessIntegrationPointVar(mesh & inmesh, const MatDoub & nodalsol, std::vector<std::vector<double>> &sol)
void elastoplastic2D<YC>::PostProcessIntegrationPointVar ( std::vector<std::vector< std::vector<Doub > > > allcoords, NRmatrix<Doub>  meshnode, MatInt meshtopology, const NRmatrix<Doub>  & nodalsol, std::vector<std::vector<double>> &sol )
{
//	std::vector<std::vector< std::vector<Doub > > > allcoords = inmesh.GetAllCoords();
//	MatDoub meshnodes = inmesh.GetMeshNodes();
//	MatInt meshtopology = inmesh.GetMeshTopology();
    int type = 1;
    //shapequad objshapes ( fOrder, type );

    NRmatrix<Doub>  ptsweigths,psis,psist,GradPsi,xycoords;
    Doub xi, eta;
   // shapequad shape = shapequad ( fOrder, 1 );
    fshape->pointsandweigths ( ptsweigths );
    Int npts = ptsweigths.nrows();
    Int nels = meshtopology.nrows();
    NRmatrix<Doub>  elcoords;
    Int counter=0;
    for ( Int iel = 0; iel < nels; iel++ ) {
        GetElCoords ( allcoords, iel, elcoords );
        for ( Int ipt = 0; ipt < npts; ipt++ ) {
            std::vector<double> soli ( 3 );
            xi = ptsweigths[ipt][0];
            eta = ptsweigths[ipt][1];

            fshape->shapes ( psis, GradPsi, xi, eta );
            psis.Transpose ( psist );
            psist.Mult ( elcoords, xycoords );
            Doub plasticj2 = fepspvec[counter].J2();
            soli[0] = xycoords[0][0];
            soli[1] = xycoords[0][1];
            soli[2] =sqrt ( plasticj2 );
            sol.push_back ( soli );
            counter++;

        }
    }
}

template <class YC>
void elastoplastic2D<YC>::ComputeSolAndDSol2 ( mesh * inmesh,NRmatrix<Doub>&sol,NRmatrix<Doub>&dsol )
{
    DebugStop();
}


template <class YC>
void elastoplastic2D<YC>::ComputeSolAndDSol ( mesh * inmesh,NRvector<NRmatrix<Doub>>&sol,NRvector<NRmatrix<Doub>>&dsol )
{
    DebugStop();
}


template <class YC>
void elastoplastic2D<YC>::ComputeSolAndDSol ( mesh * inmesh,NRmatrix<Doub>&sol,NRmatrix<Doub>&dsol )
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



template <class YC>
void elastoplastic2D<YC>::ComputeSolution ( mesh * inmesh,NRmatrix<Doub>  elcoords,NRmatrix<Doub>  eldisplace,NRvector<Doub> ptsw,NRmatrix<Doub> &sol,NRmatrix<Doub> &dsol )
{
    sol.resize ( 1,2 );
    dsol.resize ( 1,4 );
    Doub xi,eta,w;
    xi=ptsw[0];
    eta=ptsw[1];
    w=ptsw[2];
    int type = 1;
    //shapequad objshapes ( fOrder, type );

    NRmatrix<Doub>  psis,psist, GradPsi,xycoords,Jac,InvJac ( 2,2 ),GradPhi;
    fshape->shapes ( psis, GradPsi, xi, eta );
    psis.Transpose ( psist );
    psist.Mult ( elcoords, xycoords );

    GradPsi.Mult ( elcoords, Jac );
    Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
    if ( DetJ<=0 ) {
        std::cout<< "DetJ<0 "<<std::endl;
        DebugStop();
    }

    InvJac[0][0] = Jac[1][1] / DetJ;
    InvJac[0][1] = -Jac[0][1] / DetJ;
    InvJac[1][0] = -Jac[1][0] / DetJ;
    InvJac[1][1] = Jac[0][0] / DetJ;
    InvJac.Mult ( GradPsi, GradPhi );
    NRmatrix<Doub>  gradu;
    GradPhi.Mult ( eldisplace, gradu );

    Doub ex = gradu[0][0];// dudx;
    Doub ey = gradu[1][1];
    Doub exy = ( gradu[0][1] + gradu[1][0] );
}

template <class YC>
NRmatrix<Doub>  elastoplastic2D<YC>::GetSolution()
{
    return fdisplace;
}


template class elastoplastic2D<vonmises>;
template class elastoplastic2D<druckerprager>;
template class elastoplastic2D<mohrcoulomb>;
