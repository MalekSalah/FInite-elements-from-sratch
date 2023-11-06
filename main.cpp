//=======================================================================================
//                       Programme de resolution par elements P1
//=======================================================================================

//==========================================//
///            Preprocessing               ///
//==========================================//
#include <windows.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "elements_finis.hpp"
#include "maillage.hpp"
#include "matrice.hpp"

using namespace std;

///Fichiers de sortie
string fMaillage("maillage.txt");
string fSolution("Solution.txt");

///Definition des constentes du Probleme
//Parametrres du maillage

constexpr double K_ = 1; //Prix
constexpr double T = 365*2; //Periode du temps

//Discretisation temporelle
constexpr double dt = 3; //Pas du temps
constexpr double nbr_pas = T/dt; // Nobmre des pas de temps


int main()
{
    HANDLE console_handle = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_CURSOR_INFO cursor_info;
    GetConsoleCursorInfo(console_handle, &cursor_info);
    cursor_info.bVisible = false; // Hide the cursor
    SetConsoleCursorInfo(console_handle, &cursor_info);

    cout << "Creation du maillage:";
    Maillage mc(m, n, a);
    mc.Save(fMaillage);

    long int Nbpts=mc.sommets.size(),  Nbtri=mc.numelts.size();
    vector<Point> Coorneu = mc.sommets;
    vector<Numeros> Numtri = mc.numelts;
    cout << " 100%" << "\n";

    //Initialisation des matrices M, K, B
    matrice_creuse M(Nbpts, Nbpts);//
    matrice_creuse K(Nbpts, Nbpts);//
    matrice_creuse B(Nbpts, Nbpts);//

    ///Calcul des matrices elements finis
    for (int l=0; l<Nbtri; l++)
    {
        double r=(double(l)/Nbtri)*100;
        cout << "Calcul des matrices elements fini:" << round(r) << "%\r";

        Point S1 = Coorneu[Numtri[l][0]];
        Point S2 = Coorneu[Numtri[l][1]];
        Point S3 = Coorneu[Numtri[l][2]];

        //Calcul des matrices elementaires M_el, K_el, B_el
        matrice M_el, K_el, B_el;
        M_el = matM_el(S1, S2, S3);
        K_el = matK_el(S1, S2, S3);
        B_el = matB_el(S1, S2, S3);


        //Remplissage des matrices M, K, B
        for (int i=0; i<M_el.shape().first; i++){
            int I = Numtri[l][i];
            for (int j=0; j<M_el.shape().second; j++){
                int J = Numtri[l][j];
                M(I,J) += M_el[i][j];
                K(I,J) += K_el[i][j];
                B(I,J) += B_el[i][j];
            }
        }
    }
    cout << "\n";

    ///Shema numerique
    //Calcul du vecteur Q0
    vecteur Q_0(Nbpts, 0);
    for (int l=0; l<Nbpts; l++)
    {
        Point P = Coorneu[l];
        double x1=P.x, x2=P.y;
        Q_0[l] = (x1+x2-K_>0) ? x1+x2-K_ : 0;
    }

    //Calcul de la matrice D
    matrice_creuse D = M*r+ K + B;

    vector<double> Time = Tspace(0, T, dt);
    int steps = Time.size();
    matrice P(steps,Nbpts);
    P[0] = Q_0;
    for (int t=1; t<steps; t++) {
        double r=(double(t+1)/steps)*100;
        cout << "Calcul du schema numerique:" << round(r) << "%\r";
        //P[t] = solve_sparse((M+(dt/2)*D), (M-(dt/2)*D)*P[t-1]);
        P[t] = solve_sparse((M+dt*D), M*P[t-1]);
    }
    cout << "\n";
    cout << "==> " << fMaillage << "\n";
    output(fSolution, P, Time);
    cout << "==> " << fSolution << "\n";

    return 0;
}
