#ifndef maillage
#define maillage
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <cstdlib>
#include <ctime>

using namespace std;

///Classe POINT//
class Point
{
public :
    double x ;
    double y ;
    Point (double x0=0,double y0=0): x(x0),y(y0){} //constructeur commun
    Point& operator +=(const Point & P) {x+=P.x ; y+= P.y ; return *this;}
    Point& operator -=(const Point & P) {x-=P.x ; y-= P.y ; return *this;}
    Point& operator *=(double a)        {x*=a   ; y*=a    ; return *this;}
    Point& operator /=(double a)        {x/=a   ; y/=a    ; return *this;}
    Point& tf_affine(const vector<double>& , const vector <double> &);

    void print(ofstream& out) const;
};

Point operator + (const Point& P, const Point& Q);
Point operator - (const Point& P, const Point& Q);
Point operator * (const Point& P, double a);
Point operator * (double a, const Point& P);
Point operator / (const Point& P, double a);
bool operator == (const Point& P, const Point& Q);
bool operator != (const Point& P, const Point& Q);
bool operator < (const Point& P, const Point& Q);
ostream & operator <<(ostream &, const Point &) ;




///Classe Numeros
class Numeros : public vector<int>
{
public :
    Numeros(int i1 , int i2 , int i3)
    {
        resize(3) ;
        vector <int>::iterator itn=begin() ;
        *itn++=i1 ; *itn++ = i2; *itn=i3 ;
    }

    int operator()(int i) const {return (*this)[i-1];  }
    int& operator()(int i)      {return (*this)[i-1] ; }

    void print(ofstream& out)const ;

};
ostream &operator <<(ostream &,const Numeros &) ;


///Classe maillage
class Maillage
{
public :
    vector<Point> sommets ;
    vector<Numeros> numelts ;
    Maillage(int m, int n) {maille_carre_unite(m,n) ;}
    Maillage(int m, int n,double a ) {maille_carre(m,n,a) ;}
    void maille_carre_unite(int m, int n );
    void maille_carre(int m, int n ,double a);
    void affiche() const ;
    Maillage& tf_affine(const vector<double> &, const vector<double> &);
    Maillage& operator += (const Maillage &) ;
    void Save(const string) const ;
};
Maillage operator + (const Maillage &, const Maillage &);

#endif // maillage
