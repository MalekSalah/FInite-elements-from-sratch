#include <iostream>
#include <fstream>
#include "matrice.hpp"
#include <cmath>

using namespace std;
///==========================================================================================================
///                          MATRICE DENSE
///==========================================================================================================
//==============================================================================================
//                            Opreations algebriques internes
//==============================================================================================
matrice& matrice::operator+=(const matrice& M) // U+=V
{
    dim_test_Matrice(*this, M);
    vector<vecteur>::const_iterator iv=(M.mat_).begin();
    for (vector<vecteur>::iterator ii=(*this).mat_.begin(); ii!=(*this).mat_.end(); ii++,  iv++){(*ii)+=(*iv);}
    return *this;
}

matrice& matrice::operator+=(double x) // U+=x;
{
    for (vector<vecteur>::iterator ii=(*this).mat_.begin(); ii!=(*this).mat_.end(); ii++) {*ii+=x;}
    return *this;
}

matrice& matrice::operator-=(const matrice& M) // U-=V
{
    dim_test_Matrice(*this, M);
    vector<vecteur>::const_iterator iv=(M.mat_).begin();
    for (vector<vecteur>::iterator ii=(*this).mat_.begin(); ii!=(*this).mat_.end(); ii++,  iv++){(*ii)-=(*iv);}
    return *this;
}

matrice& matrice::operator-=(double x) // U-=x
{
    for (vector<vecteur>::iterator ii=(*this).mat_.begin(); ii!=(*this).mat_.end(); ii++) {*ii-=x;}
    return *this;
}

matrice& matrice::operator*=(double x) // U*=x
{
    for (vector<vecteur>::iterator ii=(*this).mat_.begin(); ii!=(*this).mat_.end(); ii++) {(*ii)*=x;}
    return *this;
}

matrice& matrice::operator/=(double x) // U/=x
{
    for (vector<vecteur>::iterator ii=(*this).mat_.begin(); ii!=(*this).mat_.end(); ii++) {(*ii)/=x;}
    return *this;
}

//=================================================================================================
//                                   Affichage de matrice
//=================================================================================================
void matrice::print(ostream& out) const
{
    if ((ncol_==0) & (nrow_==0)) { out<< "[]"; return;}
    vector<vecteur>::const_iterator ii=(*this).mat_.begin();
    out << "[" << *ii << endl; ii++;
    while (ii != (((*this).mat_).end()))
    {
        out <<" " <<*ii << endl;
        ii++;
    }
}

ostream& operator<<(ostream& out,const matrice M)
{
    M.print(out);
    return out;
}

void matrice::print(ofstream& out) const
{
    for (int i=0; i<nrow_; i++){
        for (int j=0; j<ncol_; j++) out << (*this)[i][j] << " ";
        out << endl;
    }
}

void matrice::save(string file_name) const { ofstream out(file_name); (*this).print(out); out.close();}

matrice eye(int n)
{
    matrice M(n, n);
    for (int i=0; i<n; i++) M[i][i] = 1;
    return M;
}

//=================================================================================================
//                                   Linear System Tools
//=================================================================================================
vector<matrice> PLU_(const matrice& A)
{
    //check if matrix is square
    int n= A.shape().first, m=A.shape().second;
    if (n!=m) {cerr << "Matrix is not square in PLU"; exit(-1);}
    matrice  P= eye(n), L(n,n,0), U(A);

    //Parcours des colonnes de la matrice
    for (int j=0; j<n; j++)
    {
        //On cherche le pivot de la colonne j
        int pivot_index = j;
        double pivot = abs(U[j][j]);
        for (int i = j + 1; i < n; i++) {
            if (fabs(U[i][j]) > fabs(pivot)) {
                pivot = U[i][j];
                pivot_index = i;
            }
        }

        //If the pivot is not on the diagonal place it
        if (pivot_index!=j){
            std::swap(U[j], U[pivot_index]);
            std::swap(P[j], P[pivot_index]);
            std::swap(L[j], L[pivot_index]);

        }
        for (int i=j+1; i<n; i++){
                double mult = U[i][j]/pivot;
                L[i][j] = mult;
                U[i]-= U[j]*mult;
        }
    }
    for (int j=0; j<n; j++){L[j][j]=1;}
    return vector<matrice> {P, L, U};
}


vecteur solve_upper(const matrice& T, const vecteur& b)
{
    int n = b.size();
    vecteur sol(n,0);
    sol[n-1] = b[n-1]/T[n-1][n-1];

    for (int i = n-2; i >= 0; i--) {
        double s=0;
        for (int j=i+1; j<n; j++) s+= T[i][j]*sol[j];
        sol[i] = (b[i]-s)/ T[i][i];
    }
    return sol;
}

vecteur solve_lower(const matrice& T, const vecteur& b)
{
    int n = b.size();
    vecteur sol(n,0);
    sol[0] = b[0] / T[0][0];

    for (int i=1; i<n; i++) {
        double s= 0;
        for (int j=0; j<i; j++) s+= T[i][j]*sol[j];
        sol[i] = (b[i]-s)/ T[i][i];
    }
    return sol;
}


vecteur solve(const matrice& A, const vecteur& b)
{
    //check for dimensions
    if ((A.shape().first!=(int)b.size())){cerr << "Dimension error in defining the system"; exit(-1);}

    vector<matrice> PLU = PLU_(A); matrice P=PLU[0], L=PLU[1], U=PLU[2];
    vecteur sol(b.size(),0), sol_tmp(sol);

    sol_tmp = solve_lower(L, P*b);
    sol=solve_upper(U, sol_tmp);
    return sol;
}

//=================================================================================================
//                                   Operations binaires
//=================================================================================================
matrice operator+(const matrice& U) {return (U);} // +U
matrice operator-(const matrice& U) {matrice W(U); return(W*=(-1));} // -U
matrice operator+(const matrice& V, const matrice& U){matrice W(U); return(W+=V);} // V+U
matrice operator+(double x, const matrice& V) {matrice W(V); return(W+=x);}// x+U
matrice operator+(const matrice& V, double x) {return(x+V);} // U+x
matrice operator-(const matrice& V, const matrice& U) {return(V+(-U));}// U-V
matrice operator-(double x, const matrice& V){return(x+(-V));} // x-U
matrice operator-(const matrice& V, double x){return(-x+V);}; // U-x
matrice operator*(double x, const matrice& U){matrice W(U); return(W*=x);} // x*U
matrice operator*(const matrice& U, double x){return(x*U);} // U*x
matrice operator/(const matrice& U, double x){return(U*(1/x));} // U/x

vecteur operator*(const matrice& M, const vecteur& V)
{
    pair<int, int> shape=M.shape();
    if (shape.second!=(int)V.size()) stop("operartion matrice*vecteur, dimensions incompatibles");

    vecteur S(shape.first, 0);
    vector<vecteur> mat=M.get_matrix();
    for (int i=0; i<shape.first; i++) S[i]+=mat[i]|V ;
    return (S);
}

matrice operator*(const matrice& M, const matrice& N)
{
    pair<int, int> shapeM=M.shape(), shapeN=N.shape();
    if ( shapeM.second!=shapeN.first) stop("Operation matrice*matrice, dimensions incompatibles");

    matrice P(shapeM.first, shapeN.second, 0);
    for (int i=0; i<shapeM.first; i++){
        for(int j=0; j<shapeN.second; j++){
            for(int k=0; k<shapeM.second; k++){
                    P[i][j]+= M[i][k]*N[k][j];
            }
        }
    }
    return(P);
}


///==========================================================================================================
///                          MATRICE CREUSE
///==========================================================================================================
//==============================================================================================
//                            Constructeur
//==============================================================================================

matrice_creuse::matrice_creuse(const vector<int>& P):sparseStorage(P.size(),P.size())
{
    for (int i=0; i<nrow_+1; i++) indptr_[i] = i;
    data_= vector<double>(nrow_,1);
    indices_= P;
}


//==============================================================================================
//                            Opreations algebriques internes
//==============================================================================================
matrice_creuse& matrice_creuse::operator+=(const matrice_creuse& M)
{
    dim_test_Matrice((*this), M);
    for (int i=0; i<M.nrow_; i++){
        for (int k = M.indptr_[i]; k<M.indptr_[i+1];  k++) (*this)(i, M.indices_[k]) += M.data_[k];
    }
    return (*this);
}

matrice_creuse& matrice_creuse::operator-=(const matrice_creuse& M)
{
    dim_test_Matrice((*this), M);
    for (int i=0; i<nrow_; i++){
        for (int k = indptr_[i]; k<indptr_[i+1];  k++) (*this)(i, indices_[k]) -= data_[k];
    }
    return (*this);
}

matrice_creuse& matrice_creuse::operator*=(double x)
{
    for (auto& dp:data_) dp*=x;
    return (*this);
}

matrice_creuse& matrice_creuse::operator/=(double x)
{
    if (x==0) {cerr<< "ERREUR division par zero dans matrice_creuse/=x"; exit(-1);}
    for (auto& dp:data_) dp/=x;
    return (*this);
}

//=================================================================================================
//                                 Affichage des matrices creuses
//=================================================================================================
void matrice_creuse::print(ostream& out) const
{
    for (int i=0; i<nrow_;i++){
        out<< "[";
        for (int j=0; j<ncol_; j++){
            out << (*this)(i,j) << " ";
        }
        out << "}"<<endl;
    }
}

ostream& operator<<(ostream& out, const matrice_creuse& M)
{
    M.print(out);
    return (out);
}

//=================================================================================================
//                                 Operations
//=================================================================================================
matrice_creuse operator*(double x, const matrice_creuse& U){matrice_creuse W(U); return(W*=x);} // x*U
matrice_creuse operator*(const matrice_creuse& U, double x){return(x*U);} // U*x
matrice_creuse operator/(const matrice_creuse& U, double x){return(U*(1/x));} // U/x
matrice_creuse operator-(const matrice_creuse& V, const matrice_creuse& U) {return(V+(-U));}// U-V
matrice_creuse operator+(const matrice_creuse& U) {return (U);} // +U
matrice_creuse operator-(const matrice_creuse& U) {matrice_creuse W(U); return(W*=(-1));} // -U
matrice_creuse operator+(const matrice_creuse& V, const matrice_creuse& U){matrice_creuse W(U); return(W+=V);} // V+U

vecteur operator*(const matrice_creuse& M, const vecteur& V)
{
    pair<int, int> shape= M.shape();
    if (shape.second!=(int)V.size()) stop("operartion matrice*vecteur, dimensions incompatibles");

    auto [data, indptr, indices] = M.get_data();
    vecteur S(shape.first,0);
    for (int i=0; i<shape.first; i++){

        for (int k = indptr[i]; k<indptr[i+1];  k++){
                S[i]+= V[indices[k]]*data[k];
        }
    }
    return S;
}

matrice_creuse operator*(const matrice_creuse& M, const matrice_creuse& N)
{
    // Check dimensions of input matrices
    if (M.shape().second != N.shape().first) { cerr << "Dimensions incompatibles: matrice_creuse*matrice_creuse"; exit(-1);}

    // Initialize output matrix
    int n = M.shape().first, m = N.shape().second;
    matrice_creuse S(n, m);

    // Get data for input matrices
    auto [dataM, indptrM, indicesM] = M.get_data();
    auto [dataN, indptrN, indicesN] = N.get_data();

    // Perform matrix multiplication
    for (int i = 0; i < n; i++) {
        for (int j = indptrM[i]; j < indptrM[i+1]; j++) {
            int colM = indicesM[j];
            double valM_icolM = dataM[j];

            for (int k = indptrN[colM]; k < indptrN[colM+1]; k++) {
                int colN = indicesN[k];
                double valN_colMcolN = dataN[k];
                S(i, colN) += (valM_icolM * valN_colMcolN);
            }
        }
    }
    return S;
}


//=================================================================================================
//                                   Linear System Tools
//=================================================================================================
template<typename T>
void vect_swap(vector<T>& V, int i, int seq_i, int j, int seq_j)
{
    vector<T> tmp1(make_move_iterator(V.begin() + i), make_move_iterator(V.begin() + i + seq_i));
    vector<T> tmp2(make_move_iterator(V.begin() + j), make_move_iterator(V.begin() + j + seq_j));
    V.erase(V.begin() + j, V.begin() + j + seq_j);
    V.insert(V.begin() + j, make_move_iterator(tmp1.begin()), make_move_iterator(tmp1.end()));
    V.erase(V.begin() + i, V.begin() + i + seq_i);
    V.insert(V.begin() + i, make_move_iterator(tmp2.begin()), make_move_iterator(tmp2.end()));

}

void  matrice_creuse::swap_rows(int i, int j)
{
    int seq_i = indptr_[i+1]-indptr_[i];
    int seq_j = indptr_[j+1]-indptr_[j];
    if (i>j) swap(i,j);
    vect_swap(indices_, indptr_[i], seq_i, indptr_[j], seq_j);
    vect_swap(data_, indptr_[i], seq_i, indptr_[j], seq_j);
    for (int k=i+1; k<j+1; k++) indptr_[k]+=(seq_j-seq_i);
}



vector<matrice_creuse> PLU_(const matrice_creuse& A)
{
    double epsilon=10e-15;

    //check if matrix is square
    int n= A.shape().first, m=A.shape().second;
    if (n!=m) {cerr << "Matrix is not square in LUP"; exit(-1);}
    matrice_creuse  L(n,n), U(A);
    auto [data, indptr, indices] = U.get_data_ref();

    // initialize permutation matrix P and L;
    vector<int> p(n);
    for (int i = 0; i < n; i++) {
        L(i,i)=0;
        p[i] = i;
    }

    //Parcours des colonnes de la matrice
    for (int k=0; k<n-1; k++)
    {
        int pivot_index = k;
        double pivot_value = fabs(U(k,k));

        // Swap rows if necessary;
        for (int i=k+1; i<n; i++) {
            int t = indptr[i]; //  le k-eme element non nul de la ligne i
            if ((fabs(data[t]) > fabs(pivot_value)) && ((indices[t])==k) && (t< indptr[i+1])) {
                pivot_index = i;
                pivot_value = data[t];
            }
        }

        //Si le pivot est nul la matrice n'est pas inversible
        if (pivot_value==0){
            ofstream error_file("error_file.txt");
            cout<< "\nERROR: PIVOT IS NULL -> MATRIX IS NOT INVERTIBLE\n";
            error_file << "colonne=" << k << " ligne pivot="<< pivot_index << " pivot="<< pivot_value<<"\n" << U;
            exit(-1);
        }

        if (pivot_index != k) {
            U.swap_rows(k, pivot_index);
            L.swap_rows(k, pivot_index);
            swap(p[k], p[pivot_index]);
        }

        // k is now the position of the pivot !!!!!!!
        pivot_value = U(k,k);

        // Update L and U
        for ( int i=k+1; i<n; i++){//parcours des lignes de k+1 jusqu'a la derniere
            if (indices[indptr[i]] == k){ // verify if current_column(k) == column of the first nnz element in line i
                double mult = data[indptr[i]]/pivot_value;
                L(i,k) = mult;
                for (int c = indptr[k]; c<indptr[k+1]; c++){ //Parcours des elements non nuls de la ligne i
                        if (fabs((U(i,indices[c]) - (mult* U(k,indices[c])))) < epsilon)
                            {U.remove(i,indices[c]); }
                        else
                            { U(i,indices[c]) -= (mult* U(k,indices[c]));}
                }
            }
        }
    }

    //construction des matrices L et U
    matrice_creuse P(p);
    for (int i=0; i<n; i++) L(i,i)= 1;
    return  vector<matrice_creuse> {P,L,U};
}

vecteur solve_upper(const matrice_creuse& T, const vecteur& b)
{
    auto [data, indptr, indices] = T.get_data();
    int n = b.size();

    vecteur sol(n,0);
    sol[n-1] = b[n-1]/T(n-1,n-1);

    for (int i = n-2; i >= 0; i--) {
        double s=0;
        for (int j=indptr[i]+1; j<indptr[i+1]; j++) s+= data[j]*sol[indices[j]];
        sol[i] = (b[i]-s)/ T(i,i);
    }
    return sol;
}

vecteur solve_lower(const matrice_creuse& T, const vecteur& b)
{
    auto [data, indptr, indices] = T.get_data();
    int n = b.size();
    vecteur sol(n,0);
    sol[0] = b[0] / T(0,0);

    for (int i=1; i<n; i++) {
        double s= 0;
        for (int j=indptr[i]; j<indptr[i+1]; j++){
                s+= data[j]*sol[indices[j]];
        }
        sol[i] = (b[i]-s)/ T(i,i);
    }
    return sol;
}

vecteur solve_sparse(const matrice_creuse& A, const vecteur& b)
{
    //check for dimensions
    if ((A.shape().first!=(int)b.size())){cerr << "Dimension error in defining the system"; exit(-1);}
    vector<matrice_creuse> PLU = PLU_(A); matrice_creuse P=PLU[0], L=PLU[1], U=PLU[2];
    vecteur sol(b.size(),0), sol_tmp(sol);
    sol_tmp = solve_lower(L, P*b);
    sol=solve_upper(U, sol_tmp);
    return sol;
}
