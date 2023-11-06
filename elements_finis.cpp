#include "elements_finis.hpp"
#include "maillage.hpp"
#include <cmath>


///==============================================================================
///                      CALCUL DES MATRICES ELEMENTAIRES
///==============================================================================

//-----------------------------------------
//         Fonctions preliminaires
//-----------------------------------------
matrice mat_cov()
{
    matrice C(2,2);
    C[0]= {0.04, -0.024}; C[1]={-0.024, 0.04};
    return C;
}

matrice A(const Point& P)
{
    matrice mat_A(2,2), C=mat_cov();
    double x= P.x, y=P.y;
    mat_A[0] = { C[0][0]*x*x, C[0][1]*x*y};
    mat_A[1] = { C[1][0]*y*x, C[1][1]*y*y};
    return (mat_A*1/2);
}

vecteur V(const Point& P)
{
    vecteur vect_V(2); matrice C=mat_cov();
    double x=P.x, y=P.y;
    vect_V[0]= (C[0][0]+0.5*C[1][0]-r)*x;
    vect_V[1]= (C[1][1]+0.5*C[0][1]-r)*y;
    return vect_V;
}

matrice W_ref(const vector<Point>& P)
{
    matrice mat_W(4,3);
    int i=0;
    for (vector<Point>::const_iterator pp=P.begin(); pp!=P.end(); pp++, i++){
            double x=(*pp).x, y=(*pp).y;
            mat_W[i] = {1-x-y, x, y};
    }
    return mat_W;
}

matrice grad_W(double& x1, double& x2, double& x3, double& y1, double& y2, double& y3)
{
    matrice gard_W(3, 2);
    gard_W[0]= {y2-y3, -x2+x3};
    gard_W[1]= {y3-y1, -x3+x1};
    gard_W[2]= {y1-y2, -x1+x2};
    return gard_W;
}

pair<vector<Point>, vecteur> quad_2()
{
    vector<Point> S(3, Point());
    S[0] = Point(1./6, 1./6);
    S[1] = Point(2./3, 1./6);
    S[2] = Point(1./6, 2./3);
    vecteur C = {1./6, 1./6, 1./6};
    return std::make_pair( S, C);
}

pair<vector<Point>, vecteur> quad_3()
{
    vector<Point> S(4, Point());
    S[0] = Point(1.0/3, 1.0/3);
    S[1] = Point(0.6, 0.2);
    S[2] = Point(0.2, 0.6);
    S[3] = Point(0.2, 0.2);
    vecteur C;
    C= {-27.0/48, 25.0/48, 25.0/48, 25.0/48};

    return std::make_pair(S, C);
}

vector<Point> ref_transform(const vector<Point>& P,pair<matrice, vecteur> p)
{
    vector<Point> res;
    for (vector<Point>::const_iterator pp= P.begin(); pp!=P.end(); pp++){
            vecteur pt, new_pt;
            pt={(*pp).x, (*pp).y};
            new_pt = (p.first)*pt + p.second;
            res.push_back(Point(new_pt[0], new_pt[1]));
    }
    return (res);
}

pair<matrice, vecteur> trans_elt(double& x1, double& x2, double& x3, double& y1, double& y2, double& y3)
{
    matrice Bl(2,2,0); vecteur Sl;
    Bl[0] = {x2-x1, x3-x1};
    Bl[1] = {y2-y1, y3-y1};
    Sl = {x1, y1};
    return std::make_pair(Bl, Sl);
}

//-----------------------------------------
//         Calcul des matrices el
//-----------------------------------------

matrice matM_el(const Point& S1, const Point& S2, const Point& S3)
{
    //lecture des coordonnees de S1, S2, S3
    double x1= S1.x, y1= S1.y;
    double x2= S2.x, y2= S2.y;
    double x3= S3.x, y3= S3.y;

    double D = abs((x2-x3)*(y3-y1) - (x3-x1)*(y2-y3));

    //remplissage de la matrice M_el
    matrice M(3,3,D/12);
    for (int i=0; i<3; i++) M[i][i]= D/24;

    return(M);
}

matrice matK_el(const Point& S1, const Point& S2, const Point& S3)
{
    matrice K(3,3);

    //lecture des coordonnees de S1, S2, S3
    double x1= S1.x, y1= S1.y;
    double x2= S2.x, y2= S2.y;
    double x3= S3.x, y3= S3.y;

    //Calculs preliminaires
    double D = (x2-x3)*(y3-y1) - (x3-x1)*(y2-y3);
    double abs_det_Bl = abs((x2-x1)*(y3-y1) -(y2-y1)*(x3-x1));

    //Matrice Grad_W: Grad_W[i] = gard(wi)
    matrice Grad_W = grad_W(x1, x2, x3, y1, y2, y3);

    //Quadrature de gauss a 3 points d'ordre 2;
    pair<vector<Point>, vecteur> Q= quad_2();
    vector<Point> S_hat=Q.first; vecteur c_hat= Q.second;
    vector<Point> S_hat_trans= ref_transform(S_hat, trans_elt(x1, x2, x3, y1, y2, y3));

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            //Calcul de la fomule de quadrature
            for(int k=0; k<3; k++){
                K[i][j]+= (c_hat[k]* ((A(S_hat_trans[k])*Grad_W[i])|Grad_W[j])) * (abs_det_Bl/(D*D));
            }
        }
    }
    return K;
}

matrice matB_el(const Point& S1, const Point& S2, const Point& S3)
{
    matrice B(3,3);

    //lecture des coordonnees de S1, S2, S3
    double x1= S1.x, y1= S1.y;
    double x2= S2.x, y2= S2.y;
    double x3= S3.x, y3= S3.y;

    //Quadrature de gauss a 4 points d'ordre 3;
    pair<vector<Point>, vecteur> Q= quad_3();
    vector<Point> S_hat=Q.first; vecteur c_hat= Q.second;
    vector<Point> S_hat_trans= ref_transform(S_hat, trans_elt(x1, x2, x3, y1, y2, y3));

    //Calculs preliminaires
    double D = (x2-x3)*(y3-y1) - (x3-x1)*(y2-y3);
    double abs_det_Bl = abs((x2-x1)*(y3-y1) -(y2-y1)*(x3-x1));

    //Matrice Grad_W: Grad_W[i] = gard(wi)
    matrice Grad_W = grad_W(x1, x2, x3, y1, y2, y3);


    //Matrie W_ref: W_ref[i][j]= wj[Pi];
    matrice W = W_ref(S_hat);

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            //Calcul de la fomule de quadrature a 4 points
            for(int k=0; k<4; k++){
                B[i][j]+= c_hat[k]*(V(S_hat_trans[k])|Grad_W[j])* W[k][i]*(abs_det_Bl/D) ;
            }
        }
    }
    return B;
}


//==============================================================================
///                          Miscellanious tools
//==============================================================================
void output(string file_name,const matrice& M, const vector<double> T)
{
    ofstream out(file_name);
    out << "$Time" << endl;
    for (auto& t:T) out << t << " ";
    out << endl << "$Solution" << endl;
    M.print(out);
}

vector<double> Tspace(double t0, double tf, double dt)
{
    vector<double> Time;
    int steps = (int) (tf-t0)/dt;
    for (int i=0; i<=steps; i++) Time.push_back(t0+i*dt);
    if ((*Time.end())!= tf) Time.push_back(tf);
    return Time;
}
