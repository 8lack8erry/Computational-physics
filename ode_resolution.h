#ifndef Ode_Resolution_h
#define Ode_Resolution_h

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <functional>

using namespace std;

template<typename T, typename Q = T>
vector<T> Eulero(
    vector<function<T(T, vector<T>, vector<Q>)>> functions,
    T dt,
    T t,
    vector<T> x_0,
    vector<Q> parameter
)
{
    vector<T> x_n = x_0;

    for(int i = 0; i < x_0.size(); i++) x_n.at(i) += dt * functions.at(i) (t, x_0, parameter);
    return x_n;
};

template<typename T, typename Q = T>
vector<T> RK2(
    vector<function<T(T, vector<T>, vector<Q>)>> functions,
    T dt,
    T t,
    vector<T> x_0,
    vector<Q> parameter
)
{
    vector<T> x_n = x_0;
    vector<T> k1(x_0.size(), 0);

    for(int i = 0; i < x_0.size(); i++) k1.at(i) += dt * functions.at(i) (t, x_0, parameter);

    vector<T> x_prev = x_0;
    for(int i = 0; i < x_0.size(); i++) x_prev.at(i) += (k1.at(i) / 2.);

    vector<T> k2(x_0.size(), 0);
    for(int i = 0; i < x_0.size(); i++){
        k2.at(i) += dt * functions.at(i) (t + (dt / 2), x_prev, parameter);  
        x_n.at(i) += k2.at(i);
    } 

    return x_n;
};

template<typename T, typename Q = T>
vector<T> RK4(
    vector<function<T(T, vector<T>, vector<Q>)>> functions,
    T dt,
    T t,
    vector<T> x_0,
    vector<Q> parameter
)
{
    vector<T> x_n = x_0;

    vector<T> k1(x_0.size(), 0);
    for(int i = 0; i < x_0.size(); i++) k1.at(i) += dt * functions.at(i) (t, x_0, parameter);

    vector<T> x_prev = x_0;
    for(int i = 0; i < x_0.size(); i++) x_prev.at(i) += (k1.at(i) / 2.);

    vector<T> k2(x_0.size(), 0);
    for(int i = 0; i < x_0.size(); i++) k2.at(i) += dt * functions.at(i) (t + (dt / 2.), x_prev, parameter);  
    
    x_prev = x_0;
    for(int i = 0; i < x_0.size(); i++) x_prev.at(i) += (k2.at(i) / 2.);

    vector<T> k3(x_0.size(), 0);
    for(int i = 0; i < x_0.size(); i++) k3.at(i) += dt * functions.at(i) (t + (dt / 2.), x_prev, parameter);  

    x_prev = x_0;
    for(int i = 0; i < x_0.size(); i++) x_prev.at(i) += k3.at(i);

    vector<T> k4(x_0.size(), 0);
    for(int i = 0; i < x_0.size(); i++) k4.at(i) += dt * functions.at(i) (t + dt, x_prev, parameter);  

    for(int i = 0; i < x_0.size(); i++) x_n.at(i) += (1. / 6.) * (k1.at(i) + 2 * k2.at(i) + 2 * k3.at(i) + k4.at(i));

    return x_n;
};

template<typename T, typename Q = T>
void ode_Solution(
    vector<function<T(T, vector<T>, vector<Q>)>> functions,
    const string & method,
    T t_min,
    vector<T> x_0,
    vector<Q> parameter,
    const string & fileNameDat,
    vector<string> header, //without time (t)
    T dt,
    int N_t, 
    int precision
)
{   
    function<vector<T>(vector<function<T(T, vector<T>, vector<Q>)>>, T, T, vector<T>, vector<Q>)> f; 
    T t = t_min;
    
    vector<T> x_n = x_0;

    if (method == "Eulero" or method == "E")    f = Eulero<T, Q>;
    if (method == "RungeKutta2" or method == "RK2")    f = RK2<T, Q>;
    if (method == "RungeKutta4" or method == "RK4")    f = RK4<T, Q>;

    string h = "t";
    for(auto s : header)    h += ("\t" + s);
    h += "\n";

    ofstream file (fileNameDat.c_str());
    file << h;
    for(int i = 0; i < N_t; i++){
        file << setprecision(precision) << t;
        for(int j = 0; j < x_n.size(); j++)  file << "\t" << x_n.at(j);
        file << "\n";

        vector<T> x_prev = x_n;
        x_n = f(functions, dt, t, x_prev, parameter);
        t += dt;
    }    
};

template<typename T, typename Q = T>
void ode_Solution(
    vector<function<T(T, vector<T>, vector<Q>)>> functions,
    const string & method,
    T t_min,
    vector<T> x_0,
    vector<Q> parameter,
    const string & fileNameDat,
    vector<string> header, //without time (t)
    int N_t,
    T t_max, 
    int precision
)
{   
    T dt = (t_max - t_min) / N_t;
    ode_Solution(functions, method, t_min, x_0, parameter, fileNameDat, header, dt, N_t + 1, precision);
};

template<typename T, typename Q = T>
void ode_Solution(
    vector<function<T(T, vector<T>, vector<Q>)>> functions,
    const string & method,
    T t_min,
    vector<T> x_0,
    vector<Q> parameter,
    const string & fileNameDat,
    vector<string> header, //without time (t)
    T dt,
    T t_max, 
    int precision
)
{   
    int N_t = static_cast<int>((t_max - t_min) / dt);
    ode_Solution(functions, method, t_min, x_0, parameter, fileNameDat, header, dt, N_t + 1, precision);
};

#endif
