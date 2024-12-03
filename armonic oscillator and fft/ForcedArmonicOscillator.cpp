#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <iomanip>
#include "ode_resolution.h"

using namespace std;

template <typename T, typename Q = T>
T f_x(T t, vector<T> x, vector<Q> parameter){
    return x.at(1);
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
    vector<function<double(double, vector<double>, vector<double>)>> f = {f_x<double>, f_smorzatoForzato<double>};
    
    ode_Solution(f, "E", 0., x_0, par, "data_forzatoSmorzato_E_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK2", 0., x_0, par, "data_forzatoSmorzato_RK2_" + to_string(dt) + ".dat", h, dt, t_max, 15);
    ode_Solution(f, "RK4", 0., x_0, par, "data_forzatoSmorzato_RK4_" + to_string(dt) + ".dat", h, dt, t_max, 15);

    return 0;
}
