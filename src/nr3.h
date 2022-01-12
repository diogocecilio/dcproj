#ifndef _NR3_H_
#define _NR3_H_
#include <math.h>
#include <cmath>
//#define _CHECKBOUNDS_ 1
//#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
//#define _TURNONFPES_ 1

// all the system #include's we'll ever need
//#include "boost/numeric/ublas/matrix.hpp"
//#include "boost/numeric/ublas/io.hpp"
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include "error.h"

#define _XX_ 0
#define _YY_ 1
#define _ZZ_ 2
#define _XZ_ 3
#define _YZ_ 4
#define _XY_ 5

using namespace std;


#define _CHECKBOUNDS_
#include <Eigen/Core>
//#include<Eigen/SparseCholesky>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>

//#include <mkl.h>


//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE
//#define EIGEN_NO_DEBUG
//#define EIGEN_DONT_PARALLELIZE
//#include <cblas.h>


using namespace Eigen;


#define _USE_MATH_DEFINES

#include <cmath>






// macro-like inline functions

template<class T>
inline T SQR ( const T a )
{
    return a*a;
}

template<class T>
inline const T &MAX ( const T &a, const T &b )
{
    return b > a ? ( b ) : ( a );
}

inline float MAX ( const double &a, const float &b )
{
    return b > a ? ( b ) : float ( a );
}

inline float MAX ( const float &a, const double &b )
{
    return b > a ? float ( b ) : ( a );
}

template<class T>
inline const T &MIN ( const T &a, const T &b )
{
    return b < a ? ( b ) : ( a );
}

inline float MIN ( const double &a, const float &b )
{
    return b < a ? ( b ) : float ( a );
}

inline float MIN ( const float &a, const double &b )
{
    return b < a ? float ( b ) : ( a );
}

template<class T>
inline T SIGN ( const T &a, const T &b )
{
    return b >= 0 ? ( a >= 0 ? a : -a ) : ( a >= 0 ? -a : a );
}

inline float SIGN ( const float &a, const double &b )
{
    return b >= 0 ? ( a >= 0 ? a : -a ) : ( a >= 0 ? -a : a );
}

inline float SIGN ( const double &a, const float &b )
{
    return ( float ) ( b >= 0 ? ( a >= 0 ? a : -a ) : ( a >= 0 ? -a : a ) );
}

template<class T>
inline void SWAP ( T &a, T &b )
{
    T dum = a;
    a = b;
    b = dum;
}

// exception handling

#ifndef _USENRERRORCLASS_
//#define throw(message) \
//{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#else
struct NRerror {
    char *message;
    char *file;
    int line;
    NRerror ( char *m, char *f, int l ) : message ( m ), file ( f ), line ( l ) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch ( NRerror err )
{
    printf ( "ERROR: %s\n     in file %s at line %d\n",
             err.message, err.file, err.line );
    exit ( 1 );
}
#endif

#ifdef _USESTDVECTOR_
#define NRvector vector
#else

template <class T>
class NRvector
{
private:
    int nn;	// size of array. upper index is nn-1
    T *v;
public:
    NRvector();
    explicit NRvector ( int n );		// Zero-based array
    NRvector ( int n, const T &a );	//initialize to constant value
    NRvector ( int n, const T *a );	// Initialize to array
    NRvector ( const NRvector &rhs );	// Copy constructor
    NRvector & operator= ( const NRvector &rhs );	//assignment
    typedef T value_type; // make T available externally
    inline T & operator[] ( const int i );	//i'th element
    inline const T & operator[] ( const int i ) const;
    inline int size() const;
    void resize ( int newn ); // resize (contents not preserved)
    void assign ( int newn, const T &a ); // resize and assign a constant val
    void Print();
    //void  Toboost(boost::numeric::ublas::vector<double> &boostvec);
    //void CopyFromMatrix(NRmatrix<T> & source);

    ~NRvector();
};




template <class T>
NRvector<T>::NRvector() : nn ( 0 ), v ( NULL ) {}

template <class T>
NRvector<T>::NRvector ( int n ) : nn ( n ), v ( n>0 ? new T[n] : NULL ) {}

template <class T>
NRvector<T>::NRvector ( int n, const T& a ) : nn ( n ), v ( n>0 ? new T[n] : NULL )
{
    for ( int i = 0; i<n; i++ ) v[i] = a;
}

template <class T>
NRvector<T>::NRvector ( int n, const T *a ) : nn ( n ), v ( n>0 ? new T[n] : NULL )
{
    for ( int i = 0; i<n; i++ ) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector ( const NRvector<T> &rhs ) : nn ( rhs.nn ), v ( nn>0 ? new T[nn] : NULL )
{
    for ( int i = 0; i<nn; i++ ) v[i] = rhs[i];
}

template <class T>
NRvector<T> & NRvector<T>::operator= ( const NRvector<T> &rhs )
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
    if ( this != &rhs ) {
        if ( nn != rhs.nn ) {
            if ( v != NULL ) delete[] ( v );
            nn = rhs.nn;
            v = nn>0 ? new T[nn] : NULL;
        }
        for ( int i = 0; i<nn; i++ )
            v[i] = rhs[i];
    }
    return *this;
}

template <class T>
inline T & NRvector<T>::operator[] ( const int i )	//subscripting
{
#ifdef _CHECKBOUNDS_
    if ( i<0 || i >= nn ) {
        throw ( "NRvector subscript out of bounds" );
    }
#endif
    return v[i];
}

template <class T>
inline const T & NRvector<T>::operator[] ( const int i ) const	//subscripting
{
#ifdef _CHECKBOUNDS_
    if ( i<0 || i >= nn ) {
        throw ( "NRvector subscript out of bounds" );
    }
#endif
    return v[i];
}

template <class T>
inline int NRvector<T>::size() const
{
    return nn;
}

template <class T>
void NRvector<T>::resize ( int newn )
{
    if ( newn != nn ) {
        if ( v != NULL ) delete[] ( v );
        nn = newn;
        v = nn > 0 ? new T[nn] : NULL;
    }
}

template <class T>
void NRvector<T>::assign ( int newn, const T& a )
{
    if ( newn != nn ) {
        if ( v != NULL ) delete[] ( v );
        nn = newn;
        v = nn > 0 ? new T[nn] : NULL;
    }
    for ( int i = 0; i<nn; i++ ) v[i] = a;
}

template <class T>
NRvector<T>::~NRvector()
{
    if ( v != NULL ) delete[] ( v );
}


template <class T>
void NRvector<T>::Print()
{
    for ( int i = 0; i<this->size(); i++ ) {
        std::cout << std::fixed << this->v[i] << " \n";
    }
}

// end of NRvector definitions

#endif //ifdef _USESTDVECTOR_

//using std::vector;
//using boost::numeric::ublas::matrix;
//using namespace boost::numeric::ublas;

template <class T>
class NRmatrix
{
private:
    int nn;
    int mm;
    T **v;
public:
    NRmatrix();
    NRmatrix ( int n, int m );			// Zero-based array
    NRmatrix ( int n, int m, const T &a );	//Initialize to constant
    NRmatrix ( int n, int m, const T *a );	// Initialize to array
    NRmatrix ( const NRmatrix &rhs );		// Copy constructor
    //NRmatrix(const std::vector<std::vector<T>>&input);		// Converto to Matdouble ou Matint
    NRmatrix & operator= ( const NRmatrix &rhs );	//assignment
    typedef T value_type; // make T available externally
    inline T* operator[] ( const int i );	//subscripting: pointer to row i
    inline const T* operator[] ( const int i ) const;
    NRmatrix & operator*= ( const T &multipl );
    // NRmatrix<T>  operator*();
    NRmatrix & operator+= ( NRmatrix &source );
    NRmatrix & operator-= ( NRmatrix &source );
    inline int nrows() const;
    inline int ncols() const;
    void resize ( int newn, int newm ); // resize (contents not preserved)
    void assign ( int newn, int newm, const T &a ); // resize and assign a constant value
    void Print();
    void Print2();
    void Transpose ( NRmatrix &transpose );
    void Mult ( const NRmatrix &A, NRmatrix &B );
    void Dot ( const NRmatrix &A, NRmatrix &B );
    void  CopyFromVector ( const std::vector<std::vector<T>>&input );
    void  FromFullToVoigt ( NRmatrix<T> & full );
    NRmatrix & FromEigen ( VectorXd data );
	
	void IdentityMatrix(int n);
	
    inline void FromEigen ( MatrixXd  data )
    {
        this->assign ( data.rows(), data.cols(), 0. );
        for ( int i = 0; i < this->nrows(); i++ ) {
            for ( int j = 0; j < this->ncols(); j++ ) {
                v[i][j]=data ( i,j );
            }
        }
    }
    inline void ToEigen ( MatrixXd & data )
    {
        data.resize ( this->nrows(), this->nrows() );
        for ( int i = 0; i < this->nrows(); i++ ) {
            for ( int j = 0; j < this->ncols(); j++ ) {
                data ( i, j ) = v[i][j];
            }
        }
    }

    T DetEig()
    {
        MatrixXd A ( this->nrows(), this->nrows() );
        for ( int i = 0; i < this->nrows(); i++ ) {
            for ( int j = 0; j < this->ncols(); j++ ) {
                A ( i, j ) = v[i][j];
            }
        }
        T det = A.determinant();
        return det;
    }

    void ComputeInverse ( NRmatrix<T> &InvOut )
    {
        MatrixXd A ( this->nrows(), this->nrows() ), Inv;
        InvOut.assign ( this->nrows(), this->nrows(),0. );
        for ( int i = 0; i < this->nrows(); i++ ) {
            for ( int j = 0; j < this->ncols(); j++ ) {
                A ( i, j ) = v[i][j];
            }
        }
        Inv = A.inverse();
        for ( int i = 0; i < this->nrows(); i++ ) {
            for ( int j = 0; j < this->ncols(); j++ ) {
                InvOut[i][j] = Inv ( i, j );
            }
        }
    }

    T NRmatrixNorm();
    T Det();

    void  CopyFromVecdouble ( const  NRvector<T>&input );
    //void Toboost( matrix<double> &boostmat);
    //void Fromboost(matrix<double> &boostmat);


    //inline void CopyFromNRtensor(NRtensor<T> source);

    ~NRmatrix();
};

template <class T>
void  NRmatrix<T>::IdentityMatrix(int n)
{

	this->assign ( n, n, 0. );
    for ( int i = 0; i < n; i++ ) {
        v[i][i] = 1.;
    }
	
}

template <class T>
void  NRmatrix<T>::CopyFromVecdouble ( const  NRvector<T>&input )
{
    this->assign ( input.nrows(), 1, 0. );
    for ( int i = 0; i < this->nrows(); i++ ) {
        v[i][0] = input[i];
    }
}

template <class T>
void  NRmatrix<T>::FromFullToVoigt ( NRmatrix<T> & voigt )
{
    voigt.assign ( 6, 1, 0. );
    voigt[0][0] = v[0][0];
    voigt[1][0] = v[1][1];
    voigt[2][0] = v[2][2];
    voigt[3][0] = v[0][2];
    voigt[4][0] = v[1][2];
    voigt[5][0] = v[0][1];
}

template <class T>
T NRmatrix<T>::NRmatrixNorm()
{
    double val = 0.;
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < this->ncols(); j++ ) {
            val += v[i][j] * v[i][j];
        }
    }
    return sqrt ( val );
}

template <class T>
T NRmatrix<T>::Det()
{
    double SMALL = 1.0E-30;
    T det = 1;
    NRmatrix<T> A = *this;
    int n = A.nrows();
    // Row operations for i = 0, ,,,, n - 2 (n-1 not needed)
    for ( int i = 0; i < n - 1; i++ ) {
        // Partial pivot: find row r below with largest element in column i
        int r = i;
        double maxA = abs ( A[i][i] );
        for ( int k = i + 1; k < n; k++ ) {
            double val = abs ( A[k][i] );
            if ( val > maxA ) {
                r = k;
                maxA = val;
            }
        }
        if ( r != i ) {
            for ( int j = i; j < n; j++ ) swap ( A[i][j], A[r][j] );
            det = -det;
        }

        // Row operations to make upper-triangular
        double pivot = A[i][i];
        if ( abs ( pivot ) < SMALL ) return 0.0;         // Singular matrix

        for ( int r = i + 1; r < n; r++ ) {                // On lower rows
            double multiple = A[r][i] / pivot;                // Multiple of row i to clear element in ith column
            for ( int j = i; j < n; j++ ) A[r][j] -= multiple * A[i][j];
        }
        det *= pivot;                                        // Determinant is product of diagonal
    }

    det *= A[n - 1][n - 1];

    return det;
}


template < class T >
NRmatrix<T> & NRmatrix<T>::operator-= ( NRmatrix &source )
{
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < this->ncols(); j++ ) {
            v[i][j] -= source[i][j];
        }
    }
    return *this;
}

template < class T >
NRmatrix<T> & NRmatrix<T>::operator+= ( NRmatrix &source )
{
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < this->ncols(); j++ ) {
            v[i][j] += source[i][j];
        }
    }
    return *this;
}




template <class T>
void NRmatrix<T>::Print()
{
    std::cout << " \n";
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < this->ncols(); j++ ) {
            std::cout << std::setprecision ( 8 ) << std::fixed << this->v[i][j] << " ";
        }
        std::cout << " \n";
    }
    std::cout << " \n";
}

template <class T>
void NRmatrix<T>::Print2()
{
    std::cout << " \n";
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < this->ncols(); j++ ) {
            std::cout << " { " <<this->v[i][j][0] << " "<< this->v[i][j][0] << " } " <<" ";
        }
        std::cout << " \n";
    }
    std::cout << " \n";
}




template <class T>
void NRmatrix<T>::Transpose ( NRmatrix &transpose )
{
    transpose.resize ( this->ncols(), this->nrows() );
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < this->ncols(); j++ ) {
            transpose[j][i] = v[i][j];
        }
    }
}

template <class T>
void NRmatrix<T>::Mult ( const NRmatrix &A, NRmatrix &B )
{
    if ( this->ncols() != A.nrows() ) {
        std::cout<< " Incompatible size multiplication: " << std::endl;
        std::cout<< " this->ncols() =  "<< this->ncols()  << "|  A.nrows() = "<< A.nrows() << std::endl;
        DebugStop();
    }
    B.assign ( this->nrows(), A.ncols(), 0. );
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < A.ncols(); j++ ) {
            for ( int k = 0; k < A.nrows(); k++ ) {
                B[i][j] += v[i][k] * A[k][j];
                //cout << "\n v[i][k]  = " << v[i][k] <<  " A[j][k] = " << A[k][j];
            }
        }
    }
}


template <class T>
void NRmatrix<T>::Dot ( const NRmatrix &A, NRmatrix &B )
{
    if ( this->ncols() != A.nrows() ) {
        cout << "\n incompatible size multiplication " << endl;
#ifdef _USENRERRORCLASS_
        throw ( message ) \
        {printf ( "ERROR: %s\n     in file %s at line %d\n", message, __FILE__, __LINE__ ); throw ( 1 );}
#endif
    }
    B.assign ( this->nrows(), A.ncols(), 0. );
    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < A.ncols(); j++ ) {
            for ( int k = 0; k < A.nrows(); k++ ) {
                B[i][j] += v[i][k] * A[k][j];
                //cout << "\n v[i][k]  = " << v[i][k] <<  " A[j][k] = " << A[k][j];
            }
        }
    }
}

template <class T>
void NRmatrix<T>::CopyFromVector ( const std::vector<std::vector<T>>&input )
{
    int newn = input.size();
    int newm = input[0].size();
    //v.resize(rows, cols);


    int i, nel;
    if ( newn != nn || newm != mm ) {
        if ( v != NULL ) {
            delete[] ( v[0] );
            delete[] ( v );
        }
        nn = newn;
        mm = newm;
        v = nn>0 ? new T*[nn] : NULL;
        nel = mm*nn;
        if ( v ) v[0] = nel>0 ? new T[nel] : NULL;
        for ( i = 1; i< nn; i++ ) v[i] = v[i - 1] + mm;
    }


    for ( int i = 0; i < newn; i++ ) {
        for ( int j = 0; j < newm; j++ ) {
            v[i][j] = input[i][j];
        }
    }
}


template <class T>
NRmatrix<T>::NRmatrix() : nn ( 0 ), mm ( 0 ), v ( NULL ) {}

template <class T>
NRmatrix<T>::NRmatrix ( int n, int m ) : nn ( n ), mm ( m ), v ( n>0 ? new T*[n] : NULL )
{
    int i, nel = m*n;
    if ( v ) v[0] = nel>0 ? new T[nel] : NULL;
    for ( i = 1; i<n; i++ ) v[i] = v[i - 1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix ( int n, int m, const T &a ) : nn ( n ), mm ( m ), v ( n>0 ? new T*[n] : NULL )
{
    int i, j, nel = m*n;
    if ( v ) v[0] = nel>0 ? new T[nel] : NULL;
    for ( i = 1; i< n; i++ ) v[i] = v[i - 1] + m;
    for ( i = 0; i< n; i++ ) for ( j = 0; j<m; j++ ) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix ( int n, int m, const T *a ) : nn ( n ), mm ( m ), v ( n>0 ? new T*[n] : NULL )
{
    int i, j, nel = m*n;
    if ( v ) v[0] = nel>0 ? new T[nel] : NULL;
    for ( i = 1; i< n; i++ ) v[i] = v[i - 1] + m;
    for ( i = 0; i< n; i++ ) for ( j = 0; j<m; j++ ) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix ( const NRmatrix &rhs ) : nn ( rhs.nn ), mm ( rhs.mm ), v ( nn>0 ? new T*[nn] : NULL )
{
    int i, j, nel = mm*nn;
    if ( v ) v[0] = nel>0 ? new T[nel] : NULL;
    for ( i = 1; i< nn; i++ ) v[i] = v[i - 1] + mm;
    for ( i = 0; i< nn; i++ ) for ( j = 0; j<mm; j++ ) v[i][j] = rhs[i][j];
}



template <class T>
NRmatrix<T> & NRmatrix<T>::operator= ( const NRmatrix<T> &rhs )
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
    if ( this != &rhs ) {
        int i, j, nel;
        if ( nn != rhs.nn || mm != rhs.mm ) {
            if ( v != NULL ) {
                delete[] ( v[0] );
                delete[] ( v );
            }
            nn = rhs.nn;
            mm = rhs.mm;
            v = nn>0 ? new T*[nn] : NULL;
            nel = mm*nn;
            if ( v ) {
                if ( nel>0 ) {
                    //if(v[0] != NULL)  delete v[0];

                    v[0] = new T[nel];
                } else {
                    v[0] = NULL;
                }
                //v[0] = nel>0 ? new T[nel] : NULL;
            }
            for ( i = 1; i< nn; i++ ) v[i] = v[i - 1] + mm;
        }
        for ( i = 0; i< nn; i++ ) for ( j = 0; j<mm; j++ ) v[i][j] = rhs[i][j];
    }
    return *this;
}



template < class T >
NRmatrix<T> & NRmatrix<T>::operator*= ( const T &multipl )
{

    for ( int i = 0; i < this->nrows(); i++ ) {
        for ( int j = 0; j < this->ncols(); j++ ) {
            v[i][j] *= multipl;
        }
    }
    return *this;
}

template <class T>
inline T* NRmatrix<T>::operator[] ( const int i )	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
    if ( i<0 || i >= nn ) {
        std::cout<< "NRmatrix subscript out of bounds" << std::endl;
        DebugStop();
    }
#endif
    return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[] ( const int i ) const
{
#ifdef _CHECKBOUNDS_
    if ( i<0 || i >= nn ) {
        std::cout<< "NRmatrix subscript out of bounds" << std::endl;
        DebugStop();
    }
#endif
    return v[i];
}

template <class T>
inline int NRmatrix<T>::nrows() const
{
    return nn;
}

template <class T>
inline int NRmatrix<T>::ncols() const
{
    return mm;
}

template <class T>
void NRmatrix<T>::resize ( int newn, int newm )
{
    int i, nel;
    if ( newn != nn || newm != mm ) {
        if ( v != NULL ) {
            delete[] ( v[0] );
            delete[] ( v );
        }
        nn = newn;
        mm = newm;
        v = nn>0 ? new T*[nn] : NULL;
        nel = mm*nn;
        if ( v ) v[0] = nel>0 ? new T[nel] : NULL;
        for ( i = 1; i< nn; i++ ) v[i] = v[i - 1] + mm;
    }
}

template <class T>
void NRmatrix<T>::assign ( int newn, int newm, const T& a )
{
    int i, j, nel;
    if ( newn != nn || newm != mm ) {
        if ( v != NULL ) {
            delete[] ( v[0] );
            delete[] ( v );
        }
        nn = newn;
        mm = newm;
        v = nn>0 ? new T*[nn] : NULL;
        nel = mm*nn;
        if ( v ) v[0] = nel>0 ? new T[nel] : NULL;
        for ( i = 1; i< nn; i++ ) v[i] = v[i - 1] + mm;
    }
    for ( i = 0; i< nn; i++ ) for ( j = 0; j<mm; j++ ) v[i][j] = a;
}

//template <class T>
//NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
//{
//	int i, nel = m*n;
//	if (v) v[0] = nel>0 ? new T[nel] : NULL;
//	for (i = 1;i<n;i++) v[i] = v[i - 1] + m;
//}


template <class T>
NRmatrix<T>::~NRmatrix()
{
    if ( v != NULL ) {
        delete[] ( v[0] );
        delete[] ( v );
    }
}

template <class T>
class NRMat3d
{
private:
    int nn;
    int mm;
    int kk;
    T ***v;
public:
    NRMat3d();
    NRMat3d ( int n, int m, int k );
    inline T** operator[] ( const int i );	//subscripting: pointer to row i
    inline const T* const * operator[] ( const int i ) const;
    inline int dim1() const;
    inline int dim2() const;
    inline int dim3() const;
    void assign ( int n, int m, int k, const T& a );
    void Print();
    ~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d() : nn ( 0 ), mm ( 0 ), kk ( 0 ), v ( NULL ) {}

template <class T>
NRMat3d<T>::NRMat3d ( int n, int m, int k ) : nn ( n ), mm ( m ), kk ( k ), v ( new T**[n] )
{
    int i, j;
    v[0] = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for ( j = 1; j<m; j++ ) v[0][j] = v[0][j - 1] + k;
    for ( i = 1; i<n; i++ ) {
        v[i] = v[i - 1] + m;
        v[i][0] = v[i - 1][0] + m*k;
        for ( j = 1; j<m; j++ ) v[i][j] = v[i][j - 1] + k;
    }
}

template <class T>
inline T** NRMat3d<T>::operator[] ( const int i ) //subscripting: pointer to row i
{
    return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[] ( const int i ) const
{
    return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
    return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
    return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
    return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
    if ( v != NULL ) {
        delete[] ( v[0][0] );
        delete[] ( v[0] );
        delete[] ( v );
    }
}


template <class T>
void NRMat3d<T>::assign ( int n, int m, int k, const T& a )
{
    int i, j,w;
    v[0] = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for ( j = 1; j<m; j++ ) v[0][j] = v[0][j - 1] + k;
    for ( i = 1; i<n; i++ ) {
        v[i] = v[i - 1] + m;
        v[i][0] = v[i - 1][0] + m*k;
        for ( j = 1; j<m; j++ ) v[i][j] = v[i][j - 1] + k;
    }
    for ( i = 0; i< nn; i++ ) for ( j = 0; j<mm; j++ ) for ( w = 0; w<kk; w++ ) v[i][j][w] = a;
}

template <class T>
void NRMat3d<T>::Print()
{
    std::cout << " \n";
    for ( int i = 0; i < nn; i++ ) {
        for ( int j = 0; j < mm; j++ ) {
            for ( int k = 0; k < kk; k++ ) {
                std::cout << std::setprecision ( 8 ) << std::fixed << this->v[i][j][k] << " ";
            }
            std::cout << " \n";
        }
        std::cout << " \n";
    }
    std::cout << " \n";
}

template <class T>
class NRtensor
{
    bool fail = false;
private:

    int nn = 6;	// size of array. upper index is nn-1
    T *v = new T[6];
public:
    NRtensor();
    NRtensor ( const T &a );	//initialize to constant value
    NRtensor ( const T *a );	// Initialize to array
    NRtensor ( const NRtensor &rhs );	// Copy constructor
//	NRtensor(const NRmatrix &input);	// Copy constructor
    NRtensor & operator= ( const NRtensor &rhs );	//assignment
    NRtensor & operator+ ( const NRtensor &rhs );	//assignment
    NRtensor & operator- ( const NRtensor &rhs );	//assignment
    NRtensor & operator*= ( const T &multipl );
    NRtensor& operator+= ( const NRtensor &rhs );
    NRtensor& operator-= ( const NRtensor &rhs );
    typedef T value_type; // make T available externally
    inline T & operator[] ( const int i );	//i'th element
    inline const T & operator[] ( const int i ) const;
    void assign ( int newn, const T &a ); // resize and assign a constant val
    void Print();
    ~NRtensor();
    bool isfailed()
    {
        return fail;
    }

    inline T & XX() const
    {
        return v[_XX_];
    }

    inline T & XY() const
    {
        return v[_XY_];
    }

    inline T & XZ() const
    {
        return v[_XZ_];
    }

    inline T & YY() const
    {
        return v[_YY_];
    }

    inline T & YZ() const
    {
        return v[_YZ_];
    }

    inline T & ZZ() const
    {
        return v[_ZZ_];
    }

    void Multiply ( const NRtensor<T> tensor, NRtensor<T> & resp ) const;

    void  DeviatoricDiagonal ( NRvector<T> & vec ) const;

    T J2() const;
	T J3() const;

    T I1() const;
    T I2() const;
    T I3() const;
	
	NRtensor<T> dI1() const;
	NRtensor<T> dI2() const;
    NRtensor<T> dI3() const;
    NRtensor<T> dJ2() const;
	NRtensor<T> dJ3() const;
	

    void EigenValue ( NRvector<T> &eigenval ) const;


    void  FromTensorToNRmatrix ( NRmatrix<T> & resp );

    void  FromTensorVoigtToFullTensor ( NRmatrix<T> & full );

    void  FromNRmatrixToTensor ( NRmatrix<T> & input, NRtensor<T> & resp );

    //void  FromNRmatrixToNRtensor(NRtensor<T> & resp);

    void  EigenSystem ( NRmatrix<T> & eigenvalues, NRmatrix<T> & eigenvectors );

    void ComputeS ( NRtensor<T> & out );

    inline T Norm()
    {

        T norm = 0.;
        for ( int i = 0; i < 6; i++ ) {
            norm += v[i] * v[i];
        }
        norm += v[_XY_] * v[_XY_];
        norm += v[_XZ_] * v[_XZ_];
        norm += v[_YZ_] * v[_YZ_];
        return sqrt ( norm );
    }

    void CopyFromNRmatrix ( NRmatrix<T> source );
};




template <class T>
void NRtensor<T>::CopyFromNRmatrix ( NRmatrix<T> source )
{

    XX() = source[0][0];
    YY() = source[1][0];
    ZZ() = source[2][0];
    XZ() = source[3][0];
    YZ() = source[4][0];
    XY() = source[5][0];
}

// NRtensor definitions

template <class T>
NRtensor<T>::NRtensor() //: nn(6), v(nn>0 ? new T[nn] : NULL)
{
    for ( int i = 0; i < 6; i++ ) v[i] = 0.;
}

//template <class T>
//NRtensor<T>::NRtensor(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
NRtensor<T>::NRtensor ( const T& a ) : nn ( 6 ), v ( nn>0 ? new T[nn] : NULL )
{
    for ( int i = 0; i<nn; i++ ) v[i] = a;
}

template <class T>
NRtensor<T>::NRtensor ( const T *a ) : nn ( 6 ), v ( nn>0 ? new T[nn] : NULL )
{
    for ( int i = 0; i<nn; i++ ) v[i] = *a++;
}

template <class T>
NRtensor<T>::NRtensor ( const NRtensor<T> &rhs ) : nn ( rhs.nn ), v ( nn>0 ? new T[nn] : NULL )
{
    for ( int i = 0; i<nn; i++ ) v[i] = rhs[i];
}


template <class T>
NRtensor<T> & NRtensor<T>::operator= ( const NRtensor<T> &rhs )
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
    if ( this != &rhs ) {
        for ( int i = 0; i<nn; i++ )
            v[i] = rhs[i];
    }
    return *this;
}

template <class T>
NRtensor<T> & NRtensor<T>::operator+ ( const NRtensor<T> &rhs )
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
    if ( this != &rhs ) {
        if ( nn != rhs.nn ) {
            if ( v != NULL ) delete[] ( v );
            nn = rhs.nn;
            v = nn>0 ? new T[nn] : NULL;
        }
        for ( int i = 0; i<nn; i++ )
            v[i] += rhs[i];
    }
    return *this;
}

template <class T>
NRtensor<T> & NRtensor<T>::operator- ( const NRtensor<T> &rhs )
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
    if ( this != &rhs ) {
        if ( nn != rhs.nn ) {
            if ( v != NULL ) delete[] ( v );
            nn = rhs.nn;
            v = nn>0 ? new T[nn] : NULL;
        }
        for ( int i = 0; i<nn; i++ )
            v[i] -= rhs[i];
    }
    return *this;
}


template <class T>
NRtensor<T> & NRtensor<T>::operator*= ( const T &multipl )
{
    int i;
    for ( i = 0; i < 6; i++ ) v[i] *= multipl;
    return *this;
}


template <class T>
NRtensor<T> & NRtensor<T>::operator+= ( const NRtensor<T> &rhs )
{
    int i;
    for ( i = 0; i < 6; i++ ) v[i] += rhs[i];
    return *this;
}

template <class T>
NRtensor<T> & NRtensor<T>::operator-= ( const NRtensor<T> &rhs )
{
    int i;
    for ( i = 0; i < 6; i++ ) v[i] -= rhs[i];
    return *this;
}

template <class T>
inline T & NRtensor<T>::operator[] ( const int i )	//subscripting
{
#ifdef _CHECKBOUNDS_
    if ( i<0 || i >= nn ) {
        std::cout<< "NRvector subscript out of bounds" << std::endl;
        DebugStop();
    }
#endif
    return v[i];
}

template <class T>
inline const T & NRtensor<T>::operator[] ( const int i ) const	//subscripting
{
#ifdef _CHECKBOUNDS_
    if ( i<0 || i >= nn ) {
        std::cout<< "NRvector subscript out of bounds" << std::endl;
        DebugStop();
    }
#endif
    return v[i];
}


template <class T>
void NRtensor<T>::assign ( int newn, const T& a )
{
    if ( newn != nn ) {
        if ( v != NULL ) delete[] ( v );
        nn = newn;
        v = nn > 0 ? new T[nn] : NULL;
    }
    for ( int i = 0; i<nn; i++ ) v[i] = a;
}

template <class T>
NRtensor<T>::~NRtensor()
{
    if ( v != NULL ) delete[] ( v );
}


template <class T>
void NRtensor<T>::Print()
{

    std::cout << " \n";
    std::cout << v[_XX_] << std::endl;
    std::cout << v[_YY_] << std::endl;
    std::cout << v[_ZZ_] << std::endl;
    std::cout << v[_XZ_] << std::endl;
    std::cout << v[_YZ_] << std::endl;
    std::cout << v[_XY_] << std::endl;
    std::cout << " \n";
}

template < class T >
void NRtensor<T>::DeviatoricDiagonal ( NRvector<T> & vec ) const
{
    T p = I1() / T ( 3. );
    vec[0] = v[_XX_] - p;
    vec[1] = v[_YY_] - p;
    vec[2] = v[_ZZ_] - p;
}

template < class T >
T NRtensor<T>::J2() const
{
    //NRvector<T> s(3);
  //  DeviatoricDiagonal(s);
  //  T value = -s[0] * s[1] - s[0] * s[2] - s[1] * s[2] + v[_XY_] * v[_XY_] + v[_XZ_] * v[_XZ_] + v[_YZ_] * v[_YZ_];

    T value = ( pow ( v[_XX_], 2 ) + 3 * pow ( v[_XY_], 2 ) + 3 * pow ( v[_XZ_], 2 ) + pow ( v[_YY_], 2 ) + 3 * pow ( v[_YZ_], 2 ) - v[_YY_] * v[_ZZ_] +
                pow ( v[_ZZ_], 2 ) - v[_XX_] * ( v[_YY_] + v[_ZZ_] ) ) / 3.;
	
  
// 	T sxx = v[_XX_],syy=v[_YY_],szz=v[_ZZ_],sxy=v[_XY_],sxz=v[_XZ_],syz=v[_YZ_];
// 	T value =(pow(sxx,2.) + 3.*pow(sxy,2.) + 3.*pow(sxz,2.) + pow(syy,2.) + 3.*pow(syz,2.) - 
//      syy*szz + pow(szz,2.) - sxx*(syy + szz))/3.;

    if ( value < 0 ) {
		cout<< "J2<0" <<endl;
        DebugStop();
    }
    return fabs ( value );
}

template < class T >
T NRtensor<T>::J3() const
{
   T sxx = v[_XX_],syy=v[_YY_],szz=v[_ZZ_],sxy=v[_XY_],sxz=v[_XZ_],syz=v[_YZ_];
    return (2.*pow(sxx,3.) - 18.*pow(sxz,2.)*syy + 2.*pow(syy,3.) + 54.*sxy*sxz*syz + 
     9.*syy*pow(syz,2.) + 9.*pow(sxy,2)*(syy - 2.*szz) + 9.*pow(sxz,2.)*szz - 
     3.*pow(syy,2.)*szz + 9.*pow(syz,2.)*szz - 3.*syy*pow(szz,2.) + 
     2.*pow(szz,3.) - 3.*pow(sxx,2.)*(syy + szz) + 
     3.*sxx*(3.*pow(sxy,2.) + 3.*pow(sxz,2.) - pow(syy,2.) - 6.*pow(syz,2.) + 
        4.*syy*szz - pow(szz,2.)))/27.;
}

template < class T >
T NRtensor<T>::I1() const
{
    return v[_XX_] + v[_YY_] + v[_ZZ_];
}


template < class T >
T NRtensor<T>::I2() const
{
    return - ( v[_XY_] * v[_XY_] +v[_XZ_] * v[_XZ_] +v[_YZ_] * v[_YZ_] )+ ( v[_XX_] * v[_YY_] +v[_YY_] * v[_ZZ_] +v[_XX_] * v[_ZZ_] );
}

template < class T >
T NRtensor<T>::I3() const
{
    return v[_XX_] * v[_YY_] * v[_ZZ_]
           + ( v[_XY_] * v[_XZ_] * v[_YZ_] ) * 2.
           - ( v[_XX_] * v[_YZ_] * v[_YZ_] +
               v[_YY_] * v[_XZ_] * v[_XZ_] +
               v[_ZZ_] * v[_XY_] * v[_XY_] );
}


template < class T >
NRtensor<T> NRtensor<T>::dI1() const{
	NRtensor<T> dI1;
	dI1.XX()=1.;dI1.YY()=1.;dI1.ZZ()=1.;
	return dI1;
}
template < class T >
NRtensor<T> NRtensor<T>::dI2() const{
	DebugStop();
}
template < class T >
NRtensor<T> NRtensor<T>::dI3() const{
	DebugStop();
}
template < class T >
NRtensor<T> NRtensor<T>::dJ2() const{
	NRtensor<T> s;
	T sxx = v[_XX_],syy=v[_YY_],szz=v[_ZZ_],sxy=v[_XY_],sxz=v[_XZ_],syz=v[_YZ_];
    s.XX() = ( 2*sxx - syy - szz ) /3.;
    s.YY() =  ( -sxx + 2*syy - szz ) /3.;
    s.ZZ() = ( -sxx - syy + 2*szz ) /3.;
    s.XZ() = 2*sxz;
    s.YZ() = 2*syz;
    s.XY() = 2*sxy;
	
		return s;
}
template < class T >
NRtensor<T> NRtensor<T>::dJ3() const
{
    NRtensor<T> dj3,s;
    T J2 = this->J2();
    T sxx = v[_XX_],syy=v[_YY_],szz=v[_ZZ_],sxy=v[_XY_],sxz=v[_XZ_],syz=v[_YZ_];
    s.XX() = ( 2*sxx - syy - szz ) /3.;
    s.YY() = ( -sxx + 2*syy - szz ) /3.;
    s.ZZ() = ( -sxx - syy + 2*szz ) /3.;
    s.XZ() = sxz;
    s.YZ() = syz;
    s.XY() = sxy;
    dj3.XX() = s[_YY_]*s[_ZZ_]-s[_YZ_]*s[_YZ_]+J2/3.;
    dj3.YY() = s[_XX_]*s[_ZZ_]-s[_XZ_]*s[_XZ_]+J2/3.;
    dj3.ZZ() = s[_XX_]*s[_YY_]-s[_XY_]*s[_XY_]+J2/3.;

    dj3.XY() = 2* ( s[_YZ_]*s[_XZ_]-s[_ZZ_]*s[_XY_] );
    dj3.YZ() = 2* ( s[_XZ_]*s[_XY_]-s[_XX_]*s[_YZ_] );
    dj3.XZ() = 2* ( s[_XY_]*s[_YZ_]-s[_YY_]*s[_XZ_] );

    return dj3;

}
	


template <class T>
void NRtensor<T>::EigenValue ( NRvector<T> &eigenval ) const
{
    const T I1 ( this->I1() ),
          I2 ( this->I2() ),
          I3 ( this->I3() );


    eigenval.resize ( 3 );

    const T R ( ( T ( -2. ) * ( I1 * I1 * I1 ) + T ( 9. ) * I1 * I2 - T ( 27. ) * I3 ) / T ( 54. ) );
    const T Q ( ( ( I1 * I1 ) - ( T ( 3. ) * I2 ) ) / T ( 9. ) );

    T theta ( 0. );
    if ( fabs ( Q ) < 1e-6 ) {
        if ( fabs ( R ) < 1e-6 ) {
            theta = M_PI_2;
        } else {
            DebugStop();
        }
    } else {
        if ( Q < 0. ) DebugStop();
        T val ( R * pow ( Q, T ( -1.5 ) ) ); // i.e. R/sqrt(Q*Q*Q);
        if ( val > +1. ) val = +1.;
        if ( val < -1. ) val = -1.;
        theta = acos ( val );
    }
    T sqrtQ ( 0. );
    if ( Q < 0. ) {
        if ( fabs ( Q ) <1.e-6 ) {
            sqrtQ = T ( 0. );
        } else {
            DebugStop();
        }
    } else {
        sqrtQ = sqrt ( Q );
    }

    eigenval[0] = T ( -2. ) * sqrtQ * cos ( theta / T ( 3. ) ) + I1 / T ( 3. );
    eigenval[1] = T ( -2. ) * sqrtQ * cos ( ( theta + T ( 2. * M_PI ) ) / T ( 3. ) ) + I1 / T ( 3. );
    eigenval[2] = T ( -2. ) * sqrtQ * cos ( ( theta - T ( 2. * M_PI ) ) / T ( 3. ) ) + I1 / T ( 3. );
}



template < class T >
void NRtensor<T>::Multiply ( const NRtensor<T> tensor, NRtensor<T> & resp ) const
{
    const T XX = this->v[_XX_];
    const T YY = this->v[_YY_];
    const T ZZ = this->v[_ZZ_];
    const T XY = this->v[_XY_];
    const T XZ = this->v[_XZ_];
    const T YZ = this->v[_YZ_];
    resp.v[_XX_] = XX * tensor.XX() + XY * tensor.XY() + XZ * tensor.XZ();
    resp.v[_YY_] = XY * tensor.XY() + YY * tensor.YY() + YZ * tensor.YZ();
    resp.v[_ZZ_] = XZ * tensor.XZ() + YZ * tensor.YZ() + ZZ * tensor.ZZ();
    resp.v[_XY_] = XX * tensor.XY() + XY * tensor.YY() + XZ * tensor.YZ();
    resp.v[_XZ_] = XX * tensor.XZ() + XY * tensor.YZ() + XZ * tensor.ZZ();
    resp.v[_YZ_] = XY * tensor.XZ() + YY * tensor.YZ() + YZ * tensor.ZZ();
}

template < class T >
void NRtensor<T>::FromTensorToNRmatrix ( NRmatrix<T> & resp )
{
    resp.assign ( 6, 1, 0. );
    resp[0][0] = this->XX();
    resp[1][0] = this->YY();
    resp[2][0] = this->ZZ();
    resp[3][0] = this->XZ();
    resp[4][0] = this->YZ();
    resp[5][0] = this->XY();

}


template < class T >
void NRtensor<T>::FromNRmatrixToTensor ( NRmatrix<T> & input, NRtensor<T> & resp )
{
    resp.XX() = input[0][0];
    resp.YY() = input[1][0];
    resp.ZZ() = input[2][0];
    resp.XZ() = input[3][0];
    resp.YZ() = input[4][0];
    resp.XY() = input[5][0];
}

template < class T >
void  NRtensor<T>::FromTensorVoigtToFullTensor ( NRmatrix<T> & full )
{
    full.assign ( 3, 3, 0. );
    full[0][0] = this->XX();
    full[0][1] = this->XY();
    full[0][2] = this->XZ();
    full[1][0] = this->XY();
    full[1][1] = this->YY();
    full[1][2] = this->YZ();
    full[2][0] = this->XZ();
    full[2][1] = this->YZ();
    full[2][2] = this->ZZ();

    //full[0][0] = v[_XX_];full[0][1] = v[_XY_];full[0][1] = v[_XZ_];
    //full[1][0] = v[_XY_];full[1][1] = v[_YY_];full[1][1] = v[_YZ_];
    //full[2][0] = v[_XZ_];full[2][1] = v[_YZ_];full[2][1] = v[_ZZ_];

}





template < class T >
void NRtensor<T>::ComputeS ( NRtensor<T> & out )
{
    NRvector<T> s ( 3 );
    DeviatoricDiagonal ( s );
    out.XX() = s[0];
    out.YY() = s[1];
    out.ZZ() = s[2];
    out.XZ() = 2*v[_XZ_];//PERIGO VERIFICAR ESSE 2*
    out.YZ() = 2*v[_YZ_];//PERIGO VERIFICAR ESSE 2*
    out.XY() = 2*v[_XY_];//PERIGO VERIFICAR ESSE 2*
}

template <class T>
class NRmatrix3D
{
private:
    NRmatrix<NRvector<T>>  vv;
    T ***v;
public:
    //NRtensor();
    //NRtensor(const T &a);	//initialize to constant value
    //NRtensor(const T *a);	// Initialize to array
    //NRtensor(const NRtensor &rhs);	// Copy constructor
    NRmatrix3D();
    NRmatrix3D ( int n, int m, int k );			// Zero-based array
    typedef NRmatrix<NRvector<T>> value_type; // make T available externally
    ~NRmatrix3D();


    void Print();
};

template <class T>
NRmatrix3D<T>::NRmatrix3D() {}

template <class T>
NRmatrix3D<T>::NRmatrix3D ( int n, int m, int k )
{
    value_type teste ( n, m );
    for ( int i = 0; i < teste.nrows(); i++ ) {
        for ( int j = 0; j < teste.ncols(); j++ ) {
            teste[i][j].assign ( k, 0. );
        }
    }
    vv = teste;
}
//
//template <class T>
//void NRmatrix3D<T>::Print()
//{
//	std::cout << " \n";
//	for (int i = 0;i < nn;i++)
//	{
//		for (int j = 0;j < mm;j++)
//		{
//			for (int k = 0;k < kk;k++)
//			{
//				std::cout << std::setprecision(8) << std::fixed << this->v[i][j][k] << " ";
//			}
//			std::cout << " \n";
//		}
//		std::cout << " \n";
//	}
//	std::cout << " \n";
//}








// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef double Doub; // default floating type
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;


typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;


// Tensors
//typedef const NRtensor<Doub> TensorDoub_I;
//typedef NRtensor<Doub> TensorDoub, TensorDoub_O, TensorDoub_IO;

typedef const NRmatrix3D<Doub> Matrix3DDoub_I;
typedef NRmatrix3D<Doub> Matrix3DDoub, Matrix3DDoub_O, Matrix3DDoub_IO;



//#include "eigen_unsym.h"
//#include "eigen_sym.h"




struct Jacobi {
    bool fail = false;
    const Int n;
    MatDoub a, v;
    VecDoub d;
    Int nrot;
    const Doub EPS;

    Jacobi ( MatDoub_I &aa ) : n ( aa.nrows() ), a ( aa ), v ( n, n ), d ( n ), nrot ( 0 ),
        EPS ( numeric_limits<Doub>::epsilon() )
    {
        Int i, j, ip, iq;
        Doub tresh, theta, tau, t, sm, s, h, g, c;
        VecDoub b ( n ), z ( n );
        for ( ip = 0; ip < n; ip++ ) {
            for ( iq = 0; iq < n; iq++ ) v[ip][iq] = 0.0;
            v[ip][ip] = 1.0;
        }
        for ( ip = 0; ip < n; ip++ ) {
            b[ip] = d[ip] = a[ip][ip];
            z[ip] = 0.0;
        }
        for ( i = 1; i <= 50; i++ ) {
            sm = 0.0;
            for ( ip = 0; ip < n - 1; ip++ ) {
                for ( iq = ip + 1; iq < n; iq++ )
                    sm += abs ( a[ip][iq] );
            }
            if ( sm == 0.0 ) {
                eigsrt ( d, &v );
                return;
            }
            if ( i < 4 )
                tresh = 0.2*sm / ( n*n );
            else
                tresh = 0.0;
            for ( ip = 0; ip<n - 1; ip++ ) {
                for ( iq = ip + 1; iq < n; iq++ ) {
                    g = 100.0*abs ( a[ip][iq] );
                    if ( i > 4 && g <= EPS*abs ( d[ip] ) && g <= EPS*abs ( d[iq] ) )
                        a[ip][iq] = 0.0;
                    else if ( abs ( a[ip][iq] ) > tresh ) {
                        h = d[iq] - d[ip];
                        if ( g <= EPS*abs ( h ) )
                            t = ( a[ip][iq] ) / h;
                        else {
                            theta = 0.5*h / ( a[ip][iq] );
                            t = 1.0 / ( abs ( theta ) + sqrt ( 1.0 + theta*theta ) );
                            if ( theta < 0.0 ) t = -t;
                        }
                        c = 1.0 / sqrt ( 1 + t*t );
                        s = t*c;
                        tau = s / ( 1.0 + c );
                        h = t*a[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a[ip][iq] = 0.0;
                        for ( j = 0; j < ip; j++ )
                            rot ( a, s, tau, j, ip, j, iq );
                        for ( j = ip + 1; j < iq; j++ )
                            rot ( a, s, tau, ip, j, j, iq );
                        for ( j = iq + 1; j < n; j++ )
                            rot ( a, s, tau, ip, j, iq, j );
                        for ( j = 0; j < n; j++ )
                            rot ( v, s, tau, j, ip, j, iq );
                        ++nrot;
                    }
                }
            }
            for ( ip = 0; ip < n; ip++ ) {
                b[ip] += z[ip];
                d[ip] = b[ip];
                z[ip] = 0.0;
            }
        }
        //throw("Too many iterations in routine jacobi");
        //std::cout << "Too many iterations in routine jacobi! Critical error!" << std::endl;
        fail = true;
    }
    inline void rot ( MatDoub_IO &a, const Doub s, const Doub tau, const Int i,
                      const Int j, const Int k, const Int l )
    {
        Doub g = a[i][j];
        Doub h = a[k][l];
        a[i][j] = g - s* ( h + g*tau );
        a[k][l] = h + s* ( g - h*tau );
    }


    inline void eigsrt ( VecDoub_IO &d, MatDoub_IO *v = NULL )
    {
        Int k;
        Int n = d.size();
        for ( Int i = 0; i < n - 1; i++ ) {
            Doub p = d[k = i];
            for ( Int j = i; j < n; j++ )
                if ( d[j] >= p ) p = d[k = j];
            if ( k != i ) {
                d[k] = d[i];
                d[i] = p;
                if ( v != NULL )
                    for ( Int j = 0; j < n; j++ ) {
                        p = ( *v ) [j][i];
                        ( *v ) [j][i] = ( *v ) [j][k];
                        ( *v ) [j][k] = p;
                    }
            }
        }
    }

};

//autovetor em colunas
template < class T >
void  NRtensor<T>::EigenSystem ( NRmatrix<T> & eigenvalues, NRmatrix<T> & eigenvectors )
{
    NRmatrix<T> fulltensor,internaleigenvectors;
    this->FromTensorVoigtToFullTensor ( fulltensor );
    Jacobi* system = new Jacobi ( fulltensor );
    if ( system->fail ) {
        fail = true;
    }
    eigenvectors = system->v;
    VecDoub valtemp = system->d;
    eigenvalues.resize ( valtemp.size(), 1 );
    for ( int i = 0; i < valtemp.size(); i++ ) eigenvalues[i][0] = valtemp[i];
    delete system;
}


// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
    turn_on_floating_exceptions()
    {
        int cw = _controlfp ( 0, 0 );
        cw &= ~ ( EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
        _controlfp ( cw, MCW_EM );
    }
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */









#endif /* _NR3_H_ */

