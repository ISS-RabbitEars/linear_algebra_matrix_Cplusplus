#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
using namespace std;

class matrix
{

public:
    
    double **e;
    matrix();
    ~matrix();
    void init(int); //initialize square matrix of dim m
    void init(int,int);	//initialize matrix of mxn dim
    void init(string); //initialize matrix and read from file
    void ident(); //identity matrix
    friend istream& operator>>(istream&,matrix&); //read in matrix from input stream (matrix form)
    friend ostream& operator<<(ostream&,matrix&); //output matrix to output stream (matrix form)
    matrix& operator=(const double&); //assign scalar value to all matrix elements
    matrix& operator=(const matrix&); //matrix assignment (initializes target memory)
    matrix operator+(matrix&); //matrix addition
    matrix operator*(double); //scalar product
    matrix operator*(matrix&); //matrix product
    friend matrix operator*(double,matrix&); //scalar product (commutative)
    matrix operator-(matrix&); //matrix subtraction
    matrix operator/(double); //scalar division
    matrix& operator+=(matrix&); //matrix addition self-reference operator
    matrix& operator-=(matrix&); //matrix subtraction self-reference operator
    matrix& operator*=(double); //scalar product self-reference operator
    matrix& operator*=(matrix&); //matrix product self-reference operator
    matrix& operator/=(double); //scalar division self-reference operator
    bool operator==(matrix&); //matrix equality comparison
    bool operator!=(matrix&); //matrix inequality comparison
    bool operator==(double); //all matrix elements are equal to scalar value
    bool operator!=(double); //at least one matrix element is not equal to scalar value
    double Norm(char); //matrix norm { '1' : column norm , 'I' : row norm , '2' : Spectral norm , 'F' : Frobenius norm (Euclidean) }
    matrix T(); //transpose
    double Tr(); //trace
    matrix reduce(int,int); //reduce matrix by 1 in each dim (m,n) by eliminating a row and a column (indexing here begins at 0)
    double det(); //determinant
    void Pivot(int,int); //interchange rows i and j
    double fround(double,int); //round double to int sig. figs.
    void FactorLU(matrix&,matrix&); //find Upper and Lower matrix if square; pass according to order previously stated.
    void FactorLU(matrix&,matrix&,int); //overload with fround int sig. figs.
    void FactorLU(matrix&,matrix&,matrix&); //overload with P (permutation matrix) 
    void FactorLU(matrix&,matrix&,matrix&,int); //overload with fround int sig. figs. w/pivoting
    matrix GaussElimLU(matrix&,matrix&,matrix&);  //find solution to Ax=b; pass U,L,b and returns x 
    matrix GaussElimLU(matrix&,matrix&,matrix&,int); //overload with fround int sig. figs.
    matrix GaussElimLU(matrix&,matrix&,matrix&,matrix&); //overload with P (permutation matrix)
    matrix GaussElimLU(matrix&,matrix&,matrix&,matrix&,int); //overload with fround int sig. figs. w/pivoting
    bool Symmetric(); //returns true if symmetric (A==A.T()) , false otherwise
    bool PosDef(); //returns true if positive definite , false otherwise. 
    matrix CholeskyL(); //returns the Cholesky Lower Matrix
    bool Tridiagonal(); //returns true if Tridiagonal , false otherwise
    void FactorTLU(matrix&,matrix&); //find Upper and Lower Decomposition for Tridiagonal matrix; pass according to order previously stated.
    matrix TridiagonalSolve(matrix&,matrix&,matrix&); //find solution to Ax=b when A is Tridiagonal; pass U,L,b and return x

private:

    int m,n;
    
};

#include "matrix.cpp"