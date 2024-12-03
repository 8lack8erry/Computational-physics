#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <ctime>
#include "Ode_Resolution.h"

using namespace std;

template <typename T, typename Q = T>
T f_x(T t, vector<T> x, vector<Q> parameter){
    return x.at(1);
}

// oscillatore armonico semplice in approssimazione di piccole oscillazioni
template <typename T, typename Q = T>
T f_semplice(T t, vector<T> x, vector<Q> parameter){
    return -x.at(0);
}

// oscillatore armonico semplice
template <typename T, typename Q = T>
T f_reale(T t, vector<T> x, vector<Q> parameter){
    return -sin(x.at(0));
}

// oscillatore armonico smorzato
template <typename T, typename Q = T>
T f_smorzato(T t, vector<T> x, vector<Q> parameter){
    return -sin(x.at(0)) - parameter.at(0) * x.at(1);
}

// oscillatore armonico smorzato soggetto ad una forzante
template <typename T, typename Q = T>
T f_smorzatoForzato(T t, vector<T> x, vector<Q> parameter){
    return -sin(x.at(0)) - parameter.at(0) * x.at(1) + parameter.at(1) * sin((2./3) * t);
}

int main(int argc, char ** argv){

    vector<double> x_0 = {0, 1};
    vector<double> par = {1, 1};
    vector<string> h = {"x_n", "y_n"};
        
    double dt = 0.1;
    
    //definizione vettore di funzioni, sistema di equazioni differenziali
    vector<function<double(double, vector<double>, vector<double>)>> f = {f_x<double>, f_semplice<double>};

    // per lo studio richiesto è utile avere vari set di dati con differenti passi
    for(int i = 50; i < 1000 ; i += 50){
        ode_Solution(f, "E", 0., x_0, par, "data_semplice_E_" + to_string(i) + ".dat", h, M_PI / 2, i, 15);
        ode_Solution(f, "RK2", 0., x_0, par, "data_semplice_RK2_" + to_string(i) + ".dat", h, M_PI / 2,i, 15);
        ode_Solution(f, "RK4", 0., x_0, par, "data_semplice_RK4_" + to_string(i) + ".dat", h, M_PI / 2,i, 15);
    }

    double t_max = 20.;

    f = {f_x<double>, f_reale<double>};
    ode_Solution(f, "E", 0., x_0, par, "data_menoSemplice_E_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK2", 0., x_0, par, "data_menoSemplice_RK2_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK4", 0., x_0, par, "data_menoSemplice_RK4_" + to_string(dt) + ".dat", h, dt, t_max, 15);

    f = {f_x<double>, f_smorzato<double>};
    ode_Solution(f, "E", 0., x_0, par, "data_smorzato_E_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK2", 0., x_0, par, "data_smorzato_RK2_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK4", 0., x_0, par, "data_smorzato_RK4_" + to_string(dt) + ".dat", h, dt, t_max, 15);

    f = {f_x<double>, f_smorzatoForzato<double>};
    ode_Solution(f, "E", 0., x_0, par, "data_forzatoSmorzato_E_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK2", 0., x_0, par, "data_forzatoSmorzato_RK2_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK4", 0., x_0, par, "data_forzatoSmorzato_RK4_" + to_string(dt) + ".dat", h, dt, t_max, 15);

    return 0;
}