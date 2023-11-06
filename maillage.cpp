#include "maillage.hpp"


///Classe Point
Point& Point::tf_affine(const vector<double>& A , const vector<double> &t)
{
    double xx=x , yy=y ;
    x = A[0]*xx+A[1]*yy+t[0];
    y = A[2]*xx+A[3]*yy+t[1];
    return *this ;
}

void Point::print(ofstream& out) const{ out << x << " " << y << endl;}

Point operator+ (const Point & P, const Point & Q)
{Point R(P) ; return R+=Q ; }

Point operator- (const Point & P, const Point & Q)
{Point R(P) ; return R-=Q ; }

Point operator* (const Point & P, double a)
{Point R(P) ; return R*=a ; }

Point operator* (double a, const Point & P)
{Point R(P) ; return R*=a ; }


Point operator/ (const Point & P, double a)
{Point R(P) ; return R/=a ; }

bool operator == (const Point & P, const Point & Q)
{return (P.x==Q.x) && (P.y==Q.y) ;}

bool operator != (const Point & P, const Point & Q)
{return !(P==Q) ;}

bool operator < (const Point & P, const Point & Q)
{if (P.x<Q.x) return true    ;
if (P.x > Q.x) return false  ;
if (P.y <Q.y) return true     ;
return false                 ; }

ostream & operator << (ostream & os, const Point &P)
{os << "(" << P.x << "," << P.y << ")" ; return os;}




///Classe Numeros
ostream & operator << (ostream & os, const Numeros &ns)
{
    int i=1;
    for (; i<int(ns.size()) ; i++) {os<<ns(i)<<" " ; }
    os << ns(i) ;
    return os ;
}

void Numeros::print(ofstream& out) const { out << (*this)[0]+1 << " " << (*this)[1]+1 << " "<< (*this)[2]+1 << endl;}

///Classe maillage
void Maillage::maille_carre_unite(int m , int n)
{
    double dx = 1./m, dy = 1./n ;
    sommets.resize((m+1)*(n+1));
    vector<Point>::iterator its=sommets.begin() ;
    for (int j=0 ; j<n+1 ; j++){
        double y=j*dy ;
        for (int i =0 ; i<m+1 ; i++, its ++) *its = Point(i*dx,y) ;
    }

    for (int j=0; j<n ; j++)
    for (int i=0 ; i<m ; i++) {
        int q = j*(m+1)+i ;
        numelts.push_back(Numeros(q,q+1,q+m+2));
        numelts.push_back(Numeros(q,q+m+2,q+m+1)) ;
    }

}


void Maillage::maille_carre(int m , int n,double a)
{
    double dx = 1./m, dy = 1./n ;
    sommets.resize((m+1)*(n+1));
    vector<Point>::iterator its=sommets.begin() ;
    for (int j=0 ; j<n+1 ; j++){
        double y=a*j*dy;
        for (int i =0 ; i<m+1 ; i++, its ++) *its = Point(a*i*dx,y) ;
    }

    for (int j=0; j<n ; j++)
    for (int i=0 ; i<m ; i++) {
        int q = j*(m+1)+i ;
        numelts.push_back(Numeros(q,q+1,q+m+2));
        numelts.push_back(Numeros(q,q+m+2,q+m+1)) ;
    }

}




void Maillage::affiche() const
{

  cout << "Liste des sommets ("<<sommets.size()<<" points)\n";
  vector<Point>::const_iterator its=sommets.begin();
  int i=0;
  while(its!= sommets.end()){
    cout << "sommet "<<i<<" : "<<*its<< "\n";
    its++; i++;
  }

  cout << "Liste des Numeros ("<<numelts.size()<<" Numeros)\n";
  vector<Numeros>::const_iterator itn=numelts.begin() ;
  i=1 ;
  while (itn != numelts.end())
  {
      cout<<"triangle "<<i<<" : " <<*itn<<"\n";
      itn++; i++;
  }

}

Maillage& Maillage::tf_affine(const vector<double> &A , const vector<double> &t)
{
    vector<Point>::iterator its=sommets.begin() ;
    for (;its!=sommets.end();its++) its -> tf_affine(A,t);
    return *this ;
}

Maillage & Maillage::operator += (const Maillage &M)
{
    vector <Point>::const_iterator itp ;
    map<Point, int> ptrang ;
    map<Point , int>::iterator itm ;
    vector <int> num2 ;
    int k =0 ;
    for (itp=sommets.begin(); itp != sommets.end() ; ++itp,k++) ptrang.insert(pair<Point, int>(*itp,k)) ;
    int l = 0 ;
    num2.resize(M.sommets.size()) ;
    for (itp=M.sommets.begin() ; itp != M.sommets.end() ; ++itp , l++  ){
        map<Point, int> :: iterator itm=ptrang.find(*itp) ;
        if (itm != ptrang.end()) {num2[l] = itm->second ;   }
        else {ptrang.insert(pair<Point, int>(*itp, k)) ;
        num2[l]=k++;
        }
    }
    sommets.resize(ptrang.size()) ;
    for (itm=ptrang.begin(); itm != ptrang.end(); ++itm  ) sommets[itm->second] = itm->first ;

    //fusion des Numeros
    vector<Numeros>::const_iterator itn=M.numelts.begin() ;
    for (; itn != M.numelts.end(); itn++  ){
        Numeros nums=*itn ;
        for (int i=0; i < int (nums.size()) ; i++) nums[i] = num2[nums[i]] ;
        numelts.push_back(nums);
    }
    return *this ;
}

//concatenation de deux maillages
Maillage operator + (const Maillage & M1, const Maillage& M2) {
    Maillage M(M1) ; return M+= M2 ;}


//Export dans un fichier
void Maillage::Save (const string fn) const
{
    ofstream os(fn) ;
    os << "$Nodes" << endl;
    os << sommets.size() << endl;
    for (auto& p: sommets) p.print(os);
    os << "$Elements" << endl;
    os << numelts.size() << endl;
    for (auto& num: numelts) num.print(os);
    os.close() ;
}



