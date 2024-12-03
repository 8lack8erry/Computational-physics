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
#include "ode_resolution.h"

using namespace std;

// vector <T> x = {x1, y1, z1, x2, y2, z2, x3, y3, z3, x_dot1, y_dot1, z_dot1, x_dot2, y_dot2, z_dot2, x_dot3, y_dot3, z_dot3}
//                 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9     , 10    , 11    , 12    , 13    , 14    , 15    , 16    , 17  
// vector <Q> parameter = {m1, m2, m3}
//                         0 , 1 , 2  

template <typename T, typename Q = T>
T f(
    int i,  // pianeta, i = 0, 1, 2
    int ax, // assi, ax = 0, 1, 2
    T t,
    vector<T> v,
    vector<Q> parameter
){
    return v.at(ax + i * 3 + 9);
}

template <typename T, typename Q = T> 
T f_dot(
    int i,  // pianeta, i = 0, 1, 2
    int ax, // assi, ax = 0, 1, 2
    T t,
    vector<T> v,
    vector<Q> parameter
){ 
    T x = 0; 
    for (int j = 0; j < 3; j++){
        if(i == j)  continue; 

        T r = 0;
        for(int k = 0; k < 3; k++)  r += pow(v.at(k + i * 3) - v.at(k + j * 3), 2);

        x += (-1) * parameter.at(j) * (v.at(ax + i * 3) - v.at(ax + j * 3)) / pow(r, 3. / 2);
    }
    return x;
}

/*
    per ridurre i parametri delle funzioni, in modo da poter applicare le funzioni 
    presentate nella libreria Ode_Resolution.h, si usa una lambda function
*/
template <typename T, typename Q = T>
auto make_function_vector(
    int i, // pianeta
    int ax, // assi
    T (*fun)(int, int, T, vector<T>, vector<Q>)
){
    return [i, ax, fun](T t, vector<T> v, vector<Q> parameter){
        return fun(i, ax, t, v, parameter);
    };
}

int main(int argc, char **argv) {
    //definizione vettore di funzioni, sistema di equazioni differenziali
    vector<function<double(double, vector<double>, vector<double>)>> g;
    //riempimento del vettore di funzioni
    for (int planet = 0; planet < 3; planet++){
        for (int axis = 0; axis < 3; axis++) {
            g.push_back(make_function_vector(planet, axis, f<double>));
        }
    }
    for (int planet = 0; planet < 3; planet++){
        for (int axis = 0; axis < 3; axis++) {
            g.push_back(make_function_vector(planet, axis, f_dot<double>));
        }
    }

    vector<string> h ={"x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "vx1", "vy1", "vz1", "vx2", "vy2", "vz2", "vx3", "vy3", "vz3"};

    vector<double> x_0 = {1, 0, 0,  // r1
          -1, 0, 0,  // r2
           0, 0, 0,  // r3
           0, 0.4, 0,   //v1
           0, -0.8, 0.7,    //v2
           0, -0.8, -0.7};     //v3

    vector<double> par = {1.6, 0.4, 0.4};

    ode_Solution(g, "E", 0., x_0, par, "3Corpi_E_2.dat", h, 0.0001, 20., 10);
    ode_Solution(g, "RK2", 0., x_0, par, "3Corpi_RK2_2.dat", h, 0.0001, 20., 10);
    ode_Solution(g, "RK4", 0., x_0, par, "3Corpi_RK4_2.dat", h, 0.0001, 20., 10);

    return 0;
}
