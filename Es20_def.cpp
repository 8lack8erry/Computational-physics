#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <ctime>
#include <functional>
#include <complex>
#include "zeri_di_funzione.h"

using namespace std;

complex<double> f (complex<double> z, vector<double> parameter){
    complex<double> u (1., 0);
    return (pow(z, 3) - u);
}

complex<double> f_dot (complex<double> z, vector<double> parameter){ // derivata di f
    complex<double> u (1., 0);
    return 3. * pow(z, 2);
}

int main(int argc, char **argv) {

    vector<double> par = {};
    double step = 0.001;
    double min = -2.;
    double max = 2.;
    pair<complex<double>, int> cont;

    ofstream file("data_Es20.dat"); // creazione file di destinazione
    file << "Re(z)\tIm(z)\tn\ttheta\n"; // intestazione file di destinazione

    // definizione griglia
    for(double i = min; i < max + step; i += step){
        for(double j = min; j < max + step; j += step){
            complex<double> z (i, j);
            cont = NewtonRaphson(function(f), function(f_dot), par, z, 0.00001);
            /*
                riempimento file di destinazione

                dato che gli zeri di f giacciono sulla circonferenza unitaria, per distinguere
                i diversi zeri si stampa l'angolo del valore in coordinate polari
            */
            file << i << "\t" << j << "\t" << cont.second + 1 << "\t" << arg(cont.first) << "\n"; 
        }
    }

    return 0;
}