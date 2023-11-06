#ifndef MATRICE_HPP
#define MATRICE_HPP

#include "storage.hpp"
#include <tuple>

class matrice  : public denseStorage
{
public:
    ///Constructeurs-destructeurs
    matrice(int n=0, int m=0, double x=0):denseStorage(n,m,x) {} //construire une matrice a n lignes et m colonnes
    matrice(const matrice& M):denseStorage(M.nrow_, M.ncol_){mat_ = M.mat_;} // constructeur par copie
    ~matrice(){mat_.clear(); ncol_=0; nrow_=0;}

    ///Assignation
    matrice& operator=(const matrice& M){nrow_=M.nrow_; ncol_=M.ncol_; mat_=M.mat_; return(*this);}
    matrice& operator=(double x){matrice W(nrow_, ncol_, x); (*this)=W; return(*this);}

    ///Oprerations d'acces
    //Line access
    vecteur& operator[](int i) {return(mat_.at(i));}
    vecteur operator[](int i) const {return(mat_.at(i));}
    //Element acces
    double& operator()(int i, int j) {return (*this)[i][j];}
    double operator()(int i, int j)const {return (*this)(i,j);}
    //Access to attributes
    vector<vecteur> get_matrix() const {return(mat_);}

    ///Operateurs algebriques unaires
    matrice& operator+=(const matrice& V); // U+=V;
    matrice& operator+=(double x); // U+=x;
    matrice& operator-=(const matrice& V); // U-=V
    matrice& operator-=(double x); // U-=x
    matrice& operator*=(double x); // U*=x
    matrice& operator/=(double x); // U/=x

    ///Affichage
    void print(ostream& out=cout) const; // ecriture dans le flux out
    void print(ofstream& out) const;
    void save(string file_name) const;

};
inline void dim_test_Matrice(const matrice M, const matrice N){if (M.shape()!=N.shape()) stop("Dim incompatibles");} // Constrole de dimension
matrice eye(int n);
vector<matrice> PLU_(const matrice& A) ;
vecteur solve_lower(const matrice& T, const vecteur& b);
vecteur solve_lower(const matrice& T, const vecteur& b);
vecteur solve(const matrice& A, const vecteur& b);

/// Opreations binaires
matrice operator+(const matrice& M); // +M
matrice operator-(const matrice& M); // -M
matrice operator+(const matrice& N, const matrice& M); // N+M
matrice operator-(const matrice& M, const matrice& N); // M-N
matrice operator*(double x, const matrice& M); // x*M
matrice operator/(const matrice& M, double x); // M/x
matrice operator*(const matrice& M, double x); // M*x
matrice operator+(double x, const matrice& M); // x+M
matrice operator+(const matrice& M, double x); // M+x
matrice operator-(double x, const matrice& M); // x-M
matrice operator-(const matrice& M, double x); // M-x
vecteur operator*(const matrice& M, const vecteur& V); // M*V
matrice operator*(const matrice& M, const matrice& N); // M*N

ostream& operator<<(ostream& out,const matrice M);


class matrice_creuse : public sparseStorage
{
public:
    ///Constructeurs-destructeurs
    matrice_creuse(int n=0, int m=0):sparseStorage(n,m) {} //construire une matrice_creuse a n lignes et m colonnes
    matrice_creuse(const matrice_creuse& M):sparseStorage(M.nrow_, M.ncol_){data_=M.data_; indices_=M.indices_; indptr_=M.indptr_;} // constructeur par copie
    matrice_creuse(const vector<int >& P);// construire une matrice de permutation a partir d'un vecteur de permutation
    ~matrice_creuse(){data_.clear(); indices_.clear(); indptr_.clear(); ncol_=0; nrow_=0;}

    ///Assignation
    matrice_creuse& operator=(const matrice_creuse& M){nrow_=M.nrow_; ncol_=M.ncol_; data_=M.data_; indices_=M.indices_; indptr_=M.indptr_; return(*this);}

    ///Operateurs algebriques unaires
    matrice_creuse& operator+=(const matrice_creuse& V); // U+=V;
    matrice_creuse& operator-=(const matrice_creuse& V); // U-=V
    matrice_creuse& operator*=(double x); // U*=x
    matrice_creuse& operator/=(double x); // U/=x

    ///Affichage
    void print(ostream& out=cout) const; // ecriture dans le flux out
    void print(ofstream& out) const;
    void save(string file_name) const;

    ///Tools
    void swap_rows(int row1, int row2);
};
ostream& operator<<(ostream& out, const matrice_creuse& M);
bool operator==(const matrice_creuse& M, const matrice_creuse& N);

inline void dim_test_Matrice(const matrice_creuse M, const matrice_creuse N){if (M.shape()!=N.shape()) stop("Dim incompatibles");} // Constrole de dimension
matrice_creuse sparse_eye(int n);
vector<matrice_creuse> PLU_(const matrice_creuse& M);
vecteur solve_lower(const matrice_creuse& T, const vecteur& b);
vecteur solve_lower(const matrice_creuse& T, const vecteur& b);
vecteur solve_sparse(const matrice_creuse& A, const vecteur& b);

///Operations binaires
matrice_creuse operator+(const matrice_creuse& M); // +M
matrice_creuse operator-(const matrice_creuse& M); // -M
matrice_creuse operator+(const matrice_creuse& N, const matrice_creuse& M); // N+M
matrice_creuse operator-(const matrice_creuse& M, const matrice_creuse& N); // M-N
matrice_creuse operator*(double x, const matrice_creuse& M); // x*M
matrice_creuse operator/(const matrice_creuse& M, double x); // M/x
matrice_creuse operator*(const matrice_creuse& M, double x); // M*x
vecteur operator*(const matrice_creuse& M, const vecteur& V);
matrice_creuse operator*(const matrice_creuse& M, const matrice_creuse& N);
matrice_creuse mult(const matrice_creuse& M, const matrice_creuse& N);


#endif // MATRICE_HPP
