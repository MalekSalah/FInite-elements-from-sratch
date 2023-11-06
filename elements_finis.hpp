#ifndef ELEMENTS_FINIS_HPP
#define ELEMENTS_FINIS_HPP
#include "maillage.hpp"
#include "vecteur.hpp"
#include "matrice.hpp"

constexpr int m=20, n=20;
constexpr double a=2.0;
constexpr double r=0.05;

vector<double> Tspace(double t0, double tf, double dt);//Generate the time vector for t0 to tf with step dt
void output(string file_name,const matrice& M, const vector<double> T); // save the matrix M and the vector T to the file "file_name"

matrice mat_cov(); // Calculer la matrice de covariance
matrice A(const Point& p);// Calculer la matrice A en P
vecteur V(const Point& P);// Calculer le vecteur V en P

matrice W_ref(const vector<Point>& P); // Calculer la valeur des wi en un point P;
matrice grad_W(double& x1, double& x2, double& x3, double& y1, double& y2, double& y3);//Calculer les gradients Grad(wi) dans un triangle

pair<vector<Point>, vecteur> quad_2(); //Quadrature d'ordre 2
pair<vector<Point>, vecteur> quad_3(); //Quadrature d'ordre 3
pair<matrice, vecteur> trans_elt(double& x1, double& x2, double& x3, double& y1, double& y2, double& y3);//Calculer la matrice Bl et le vecteur Sl
vector<Point> ref_transform(const vector<Point>& P,pair<matrice, vecteur> p);// faire la transformation Bl*Point+Sl

matrice matM_el(const Point& S1, const Point& S2, const Point& S3);
matrice matK_el(const Point& S1, const Point& S2, const Point& S3);
matrice matB_el(const Point& S1, const Point& S2, const Point& S3);
#endif // ELEMENTS_FINIS_HPP

