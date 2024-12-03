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
T fx(
    T t,
    vector<T> v,
    vector<Q> parameter
){
    return parameter.at(0) * (v.at(1) - v.at(0));
}

template <typename T, typename Q = T>
T fy(
    T t,
    vector<T> v,
    vector<Q> parameter
){
    return parameter.at(1) * v.at(0) - v.at(0) * v.at(2) - v.at(1);
}

template <typename T, typename Q = T>
T fz(
    T t,
    vector<T> v,
    vector<Q> parameter
){
    return v.at(0) * v.at(1) - parameter.at(2) * v.at(2);
}

int main(int argc, char ** argv){

    vector<double> x_0 = {2, 1, 1};
    vector<double> par = {10, 28, 8. / 3.};
    vector<string> h = {"x_n", "y_n", "z_n"};
        
    double dt = 0.01;
    double T = 25.;

    //definizione vettore di funzioni, sistema di equazioni differenziali
    vector<function<double(double, vector<double>, vector<double>)>> f = {fx<double>, fy<double>, fz<double>};

    ode_Solution(f, "E", 0., x_0, par, "data_Attrattore_E_" + to_string(dt) + ".dat", h, dt, T, 10);
    ode_Solution(f, "RK2", 0., x_0, par, "data_Attrattore_RK2_" + to_string(dt) + ".dat", h, dt, T, 10);
    ode_Solution(f, "RK4", 0., x_0, par, "data_Attrattore_RK4_" + to_string(dt) + ".dat", h, dt, T, 10);

    // set di dati per analizzare un moto caotico
    double epsilon = 1e-3;
    vector<double> x_epsilon = {2. + epsilon, 1. + epsilon, 1. + epsilon};
    ode_Solution(f, "E", 0., x_epsilon, par, "data_Attrattore_epsilon_E_" + to_string(dt) + ".dat", h, dt, T, 10);
    ode_Solution(f, "RK2", 0., x_epsilon, par, "data_Attrattore_epsilon_RK2_" + to_string(dt) + ".dat", h, dt, T, 10);
    ode_Solution(f, "RK4", 0., x_epsilon, par, "data_Attrattore_epsilon_RK4_" + to_string(dt) + ".dat", h, dt, T, 10);   

    return 0;
}