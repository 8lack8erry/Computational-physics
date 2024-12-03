#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <iomanip>
#include "Metodi_integrazione_2.0.h"

using namespace std;
/*
    questa funzione scompone un numero intero in un vettore in cui ad ogni entrata corrisponde 
    una cifra del numero di partenza
    es. 3649 ---> v = {9, 4, 6, 3}
*/
vector<int> Scomposition(int n) {
    vector<int> coo;

    if(n == 0){
        coo.push_back(0);
        return coo;
    }

    // per scomporre il numero lo si divide ripetutamente per 10
    while(n > 0) {
        int last = n % 10;
        coo.push_back(last);
        n = (n - last) / 10;
    }
    return coo;
}

/*
    questa funzione ha lo scopo di aggiungere degli zeri se il numero delle cifre
    del numero in ingresso è minore di K
    es. K = 4 & n = 56 ---> v = {6, 5, 0, 0}
*/
vector<int> Scomposition8(int n, int K){
    vector<int> v = Scomposition(n);
    while(v.size() < K)   v.push_back(0);
    return v;    
}

/*
    questa funzione definisce la norma quadra di un vettore
*/
template<typename T>
T square_norm(vector<T> data){
    T square_sum = 0;
    for(auto n : data)  square_sum += pow(n, 2);
    return square_sum;
}

double f(vector<double> x){
    double n = square_norm(x);   
    if(n <= 1)  return sqrt(1 - n);
    else    return 0;
}

/*
    questa funzione esegue l'integrazione di midpoint su una griglia di dimensione 10 x 10 x 10 x ...
    il grosso limite di questa funzione è appunto che non si può definire il numero di punti della griglia
*/
double Integrazione_midpoint_10(int M, double f(vector<double>)){
    int K = 10;
    vector<int> coord(M - 1, 0); // vettore delle coordinate date dalla scomposizione di un numero intero
    vector<double> coord_prime; // vettore delle coordinate dei punti medi
    double h = 1. / K; // fattore di scala
    double I = 0;

    for (int i = 0; i < pow(K, M - 1); i++) {
        coord = Scomposition8(i, M - 1); 
        coord_prime.clear();
        for (auto n : coord) coord_prime.push_back(h * (static_cast<double>(n) + 1. / 2));

        I += f(coord_prime);    
    }
    return I * pow(h, M - 1);
}

// funzione per trovare i valori veri
double f_true(int M){
    return pow(M_PI, static_cast<double> (M) / 2) / tgamma(1 + static_cast<double> (M) / 2); // tgamma è l'implementazione della gamma function
}

/*
    algoritmo di Monte Carlo multidimensionale
*/
double Integrazione_Montecarlo(    
    double f(vector<double>),
    vector<double> min, // estremi inferiori dell' integrale
    vector<double> max, // estremi superiori dell' integrale
    int N // numero di volte in cui viene iterato l'algoritmo
){
    vector<double> r_v; // inizializzazione vettore di numeri random
    double sum = 0; 
    for(int i = 0; i < N; i++){
        r_v.clear();
        if(min.size() == max.size()){ // check 
            for(int j = 0; j < min.size(); j++){
                r_v.push_back(rand_range(min.at(j), max.at(j)));
            }
        }
        sum += f(r_v);      
    }  
    sum /= N;

    double V = 1; // inizializzazione ipervolume
    for(int i = 0; i < min.size(); i++) V *= (max.at(i) - min.at(i));

    return V * sum;
}

int main(int argc, char** argv) {

    int K = 10;
    int N = 1e5;
    int M_max = 10;

    ofstream file ("data_Es8.dat"); // creazione file di destinazione
    file << "M\tI(midpoint)\tt(midpoint)\tI(MonteCarlo)\tt(MonteCarlo)\tI(true)\n"; // intestazione file di destinazione

    for(int M = 2; M <= M_max; M++){
        clock_t start1, start2, end1, end2;
        double tempo_det;
        double tempo_MC;

        vector<double> min (M - 1, 0);
        vector<double> max (M - 1, 1);
        //inizio monitoraggio tempi di calcolo
        start1 = clock();
        // calcolo dell' integrale
        double I_det = Integrazione_midpoint_10(M, f);
        // fine monitoraggio tempi di calcolo
        end1 = clock();
        tempo_det = ((double)(end1 - start1))/CLOCKS_PER_SEC; // tempo in secondi
        
        //inizio monitoraggio tempi di calcolo
        start2 = clock();
        // calcolo dell' integrale
        double I_MC = Integrazione_Montecarlo(f, min, max, N);
        end2 = clock();
        // fine monitoraggio tempi di calcolo
        tempo_MC = ((double)(end2 - start2))/CLOCKS_PER_SEC; // tempo in secondi

        // riempimento file di destinazione  
        file << M << "\t" << I_det * pow(2, M) << "\t" << tempo_det << "\t" << I_MC * pow(2, M) << "\t" << tempo_MC << "\t" << f_true(M) << "\n";
    }

    return 0;
}