#ifndef Metodi_integrazione_2_h
#define Metodi_integrazione_2_h

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <limits>

#include "file_manipulation.h"

#define inf(T) numeric_limits<T>::infinity()

using namespace std;

template<typename T>
T rand_range(
	T min,
	T max
)
{
	return min + (max - min) * rand() / static_cast<double> (RAND_MAX);
};

// function (T x, vector<Q> parameter)

template<typename T, typename Q = T>
T Integrazione_Trapezio(
    T f(T, vector<Q>), 
    vector<Q> parameter,
    T min, 
    T max, 
    int N
)
{
    T sum = 0;
    T x = min;
    for(int i = 0; i < N; i++){
        if(i == 0 or i == N - 1)  sum += (f(x, parameter) / 2);
        else    sum += f(x, parameter);
        x += ((max - min) / (N - 1));
    }
    return (max - min) * sum / (N - 1);
};

template<typename T, typename Q = T>
T Integrazione_Simpson(
    T f(T, vector<Q>), 
    vector<Q> parameter,
    T min, 
    T max, 
    int N
)
{
    T sum = 0;
    T x = min;
    if(N % 2 == 0)  N++;
    for(int i = 0; i < N; i++){
        if(i == 0 or i == N - 1)  sum += f(x, parameter);
        else if(i % 2 == 0)   sum += 2 * f(x, parameter);
        else    sum += 4 * f(x, parameter);
        x += ((max - min) / (N - 1));
    }
    return (max - min) * sum / (3 * (N - 1));
};

template<typename T, typename Q = T>
T Integrazione_Romberg(
    T f(T, vector<Q>), 
    vector<Q> parameter, 
    T min, 
    T max,
    int j,
    int k 
)
{

    int N = pow(2, j) + 1;
    if(k == 0)    return Integrazione_Trapezio(f, parameter, min, max, N);
    
    T num = pow(4, k) * Integrazione_Romberg(f, parameter, min, max, j, k - 1) - Integrazione_Romberg(f, parameter, min, max, j - 1, k - 1);
    T den = pow(4, k) - 1;  
    return num / den;
};

template<typename T, typename Q = T>
T Integrazione_Legendre(
    T f(T, vector<Q>), 
    vector<Q> parameter,
    int n,
    T min,
    T max
)
{
    vector<T> weight;
    vector<T> x_i;

    T z1=-1;
    T z2=1;

    T alpha=(min - max) / (z1 - z2);
    T beta=(min - alpha * z1);

    Lettura_file(x_i, weight, "/mnt/c/Users/User/Desktop/info/Pesi_gaussiani/Legendre_p" + to_string(n) + ".txt");
    
    T I = 0;
    for(int i = 0; i < weight.size(); i++)  I += weight.at(i) * f(alpha * x_i.at(i) + beta, parameter);

    return alpha * I;
};

template<typename T, typename Q = T>
T Integrazione_Laguerre(
    T f(T, vector<Q>), 
    vector<Q> parameter,
    int n, 
    T min,
    T max = inf(T)
)
{
    vector<T> weight; 
    vector<T> x_i; 

    Lettura_file(x_i, weight, "/mnt/c/Users/User/Desktop/info/Pesi_gaussiani/Legendre_p" + to_string(n) + ".txt");

    //f(x + a) e^(x + a) e^(-x - a) = f(x + a) e^(x) * w
    T I = 0;
    if(max == inf(T)){
        for(int i = 0; i < weight.size(); i++)  I += weight.at(i) * f(x_i.at(i) + min, parameter) * exp(x_i.at(i));
    }
    else {
        T I1 = 0;
        T I2 = 0;

        for(int i = 0; i < weight.size(); i++){
            I1 += weight.at(i) * f(x_i.at(i) + min, parameter) * exp(x_i.at(i));
            I2 += weight.at(i) * f(x_i.at(i) + max, parameter) * exp(x_i.at(i)); 
        }  
        I = I1 - I2;
    }

    return I;
};

template<typename T, typename Q = T>
T Integrazione_Hermite(
    T f(T, vector<Q>), 
    vector<Q> parameter,
    int n,
    T min = (-1) * inf(T)
)
{
    vector<T> weight; 
    vector<T> x_i;

    Lettura_file(x_i, weight, "/mnt/c/Users/User/Desktop/info/Pesi_gaussiani/Legendre_p" + to_string(n) + ".txt");
    
    T I = 0;
    if(min == (-1) * inf(T)){
        for(int i = 0; i < weight.size(); i++)  I += weight.at(i) * f(x_i.at(i), parameter) * exp(pow(x_i.at(i), 2));
    }
    else{
        for(int i = 0; i < weight.size(); i++)  I += weight.at(i) * f(abs(x_i.at(i)) + min, parameter) * exp(pow(x_i.at(i), 2));
        I *= (1. / 2);
    }

    return I;
};

template<typename T, typename Q = T>
T Integrazione_MonteCarlo(
    T f(T, vector<Q>), 
    vector<Q> parameter,
    T min, 
    T max, 
    int N
)
{
    T sum = 0; 
    for(int i = 0; i < N; i++){
        T r_v=rand_range(min, max);
        sum += f(r_v, parameter); 
    }  

    sum /= N;

    return (max - min) * sum;
};

template<typename T, typename Q = T>
T Integrazione_HitOrMiss(
    T f(T, vector<Q>), 
    vector<Q> parameter,
    T x_min,
    T x_max,
    T y_min,
    T y_max,
    int N
)
{
    T area = (y_max - y_min) * (x_max - x_min);
    int cont = 0;
    T x;
    T y;

    for(int i = 0; i < N; i++){
        x = rand_range(x_min, x_max);
        y = rand_range(y_min, y_max);
        if(y <= f(x, parameter)) cont++;     
    }

    return area * cont / N;
};

#endif
