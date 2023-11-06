#include <iostream>
#include "vecteur.hpp"
using namespace std;


//==============================================================================================
//                            Opreations algebriques internes
//==============================================================================================
vecteur& vecteur::operator+=(const vecteur& V) // U+=V
{
    dim_test(*this, V);
    vecteur::const_iterator iv=V.begin();
    for (vecteur::iterator ii=(*this).begin(); ii!=(*this).end(); ii++,  iv++){(*ii)+=(*iv);}
    return *this;
}

vecteur& vecteur::operator+=(double x) // U+=x;
{
    for (vecteur::iterator ii=(*this).begin(); ii!=(*this).end(); ii++) {*ii+=x;}
    return *this;
}

vecteur& vecteur::operator-=(const vecteur& V) // U-=V
{
    dim_test(*this, V);
    vecteur::const_iterator iv=V.begin();
    for (vecteur::iterator ii=(*this).begin(); ii!=(*this).end(); ii++,  iv++) {*ii-=*iv;}
    return *this;
}

vecteur& vecteur::operator-=(double x) // U-=x
{
    for (vecteur::iterator ii=(*this).begin(); ii!=(*this).end(); ii++) {*ii-=x;}
    return *this;
}

vecteur& vecteur::operator*=(double x) // U*=x
{
    for (vecteur::iterator ii=(*this).begin(); ii!=(*this).end(); ii++) {*ii*=x;}
    return *this;
}

vecteur& vecteur::operator/=(double x) // U/=x
{
    for (vecteur::iterator ii=(*this).begin(); ii!=(*this).end(); ii++) {*ii/=x;}
    return *this;
}


//=================================================================================================
//                                   Affichage de vecteur
//=================================================================================================
void vecteur::print(ostream& out) const
{
    out << "[";
    for (vecteur::const_iterator ii=(*this).begin(); ii!=(*this).end(); ii++){out << " " << *ii;}
    out << " ]";
}

ostream& operator<<(ostream& out,const vecteur& v)
{
    v.print(out);
    return out;
}


//=================================================================================================
//                                   Opreations binaires externes
//=================================================================================================
vecteur operator+(const vecteur& U) {return (U);} // +U
vecteur operator-(const vecteur& U) {vecteur W(U); return(W*=(-1));} // -U
vecteur operator+(const vecteur& V, const vecteur& U){vecteur W(U); return(W+=V);} // V+U
vecteur operator+(double x, const vecteur& V) {vecteur W(V); return(W+=x);}// x+U
vecteur operator+(const vecteur& V, double x) {return(x+V);} // U+x
vecteur operator-(const vecteur& V, const vecteur& U) {return(V+(-U));}// U-V
vecteur operator-(double x, const vecteur& V){return(x+(-V));} // x-U
vecteur operator-(const vecteur& V, double x){return(-x+V);}; // U-x
vecteur operator*(double x, const vecteur& U){vecteur W(U); return(W*=x);} // x*U
vecteur operator*(const vecteur& U, double x){return(x*U);} // U*x
vecteur operator/(const vecteur& U, double x){return(U*(1/x));} // U/x

double operator|(const vecteur& U, const vecteur& V) // U|V
{
    dim_test(U, V);
    double S=0;
    vecteur::const_iterator iv=V.begin();
    for (vecteur::const_iterator iu=U.begin(); iu!=U.end(); iu++, iv++){S+=((*iv)*(*iu));}
    return(S);
}
