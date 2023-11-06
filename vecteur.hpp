#ifndef VECTEUR
#define VECTEUR

#include <iostream>
#include <vector>
using namespace std;



class vecteur : public vector<double>
{
public:
    ///Constructeurs-destructeurs
    vecteur(int n=0, double x=0):vector<double>(n, x){} //constructeur dim valeur initiale
    vecteur(std::initializer_list<double> lst) : vector<double>(lst) {}
    vecteur(const vecteur& v):vector<double>(v){} //Constructeur par copie
    ~vecteur(){(*this).clear();} // Destructeur

    ///Assignation
    vecteur& operator=(std::initializer_list<double> lst) { assign(lst); return *this; }
    vecteur& operator=(const vecteur& v){(*this).assign(v.begin(), v.end()); return(*this);} // U = V
    vecteur& operator=(double x) {vecteur W((*this).size(), x); (*this)=W; return(*this);} // U=x(double)

    ///Opearteurs d'acces
    double& operator[](int i){return((*this).at(i));} // Acces a l'element i=0...n
    double operator[](int i) const {return((*this).at(i));}// Acces a l'element i=0...n

    ///Operateurs algebriques unaires
    vecteur& operator+=(const vecteur& V); // U+=V;
    vecteur& operator+=(double x); // U+=x;
    vecteur& operator-=(const vecteur& V); // U-=V
    vecteur& operator-=(double x); // U-=x
    vecteur& operator*=(double x); // U*=x
    vecteur& operator/=(double x); // U/=x

    ///Affichage
    void print(ostream& out=cout) const; // ecriture dans le flux out
};

///Tools
inline void stop(const char* arg){ cout << "ERREUR : " << arg<< endl; exit(-1);} // Arret
inline void dim_test(const vecteur v, const vecteur u){if (v.size()!=u.size()) stop("Dim incompatibles");} // Constrole de dimension


/// Opreations binaires
vecteur operator+(const vecteur& U); // +U
vecteur operator-(const vecteur& U); // -U
vecteur operator+(const vecteur& V, const vecteur& U); // V+U
vecteur operator+(double x, const vecteur& V); // x+U
vecteur operator+(const vecteur& V, double x); // U+x
vecteur operator-(const vecteur& V, const vecteur& U); // U-V
vecteur operator-(double x, const vecteur& V); // x-U
vecteur operator-(const vecteur& V, double x); // U-x
vecteur operator*(double x, const vecteur& U); // x*U
vecteur operator*(const vecteur& U, double x); // U*x
vecteur operator/(const vecteur& U, double x); // U/x
double operator|(const vecteur& U, const vecteur& V); // U|V

/// Surcharge de l'operateur flux <<
ostream& operator<<(ostream& out,const vecteur& v);

#endif // VECTEUR
