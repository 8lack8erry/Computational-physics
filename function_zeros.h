#ifndef zeri_di_funzione_h
#define zeri_di_funzione_h

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <limits>
#include <functional>

#define nan(T) std::numeric_limits<T>::quiet_NaN()

using namespace std;

template<typename T, typename Q = T>
vector<pair<T, T>> intervalli_bisezione(
    function<T(T, vector<Q>)> f,
    vector<Q> parameter,
    T min, 
    T max,
    int N_range
){
    vector<pair<T, T>> interval;
    T step = (max - min) / N_range;

    for(T i = min; i < max + step; i += step){
        if(f(i, parameter) * f(i + step, parameter) > 0) continue;
        else if(f(i, parameter) * f(i + step, parameter) == 0){
            if(f(i, parameter) == 0) interval.push_back(make_pair(i, i + step));
            else continue;    
        }
        else    interval.push_back(make_pair(i, i + step));
    }
    return interval;
};

template<typename T, typename Q = T>
T bisezione(    // Su singolo intervallo
    function<T(T, vector<Q>)> f,
    vector<Q> parameter,
    T min, 
    T max,
    double precision
){
    if(f(min, parameter) == 0)  return min;
    if(f(max, parameter) == 0)  return max;
    T x_star = min;
	
	while((max - min) > precision){
		x_star = 0.5 * (max + min);
		
		if(f(x_star, parameter) * f(min, parameter)>0.)	min = x_star;
		else	max = x_star;
	}
	return x_star;
};

template<typename T, typename Q = T>
vector<T> bisezione(
    function<T(T, vector<Q>)> f,
    vector<Q> parameter,
    T min, 
    T max,
    int N_range,
    double precision
){
    vector<pair<T, T>> interval = intervalli_bisezione(f, parameter, min, max, N_range);
    vector<T> zeros;
    
    for(int i = 0; i < interval.size(); i++)    zeros.push_back(bisezione(f, parameter, interval.at(i).first, interval.at(i).second, precision));

    return zeros;
};

template<typename T, typename Q = T>
pair<T, int> NewtonRaphson( //Su singolo punto
    function<T(T, vector<Q>)> f,
    function<T(T, vector<Q>)> f_dot,
    vector<Q> parameter,
    T x0,
    double precision, 
    int Max_iteration = 10000
){
    T x = x0;
    T d = 0;
    int n = 0;
    while(n < Max_iteration){
        if(abs(f(x, parameter)) < precision) return {x, n};
        d = - f(x, parameter) / f_dot(x, parameter);
        x += d;
        n++;
    }
    return {nan(T), Max_iteration};
}

template<typename T, typename Q = T>
vector<pair<T, int>> NewtonRaphson(
    function<T(T, vector<Q>)> f,
    function<T(T, vector<Q>)> f_dot,
    vector<Q> parameter,
    T min,
    T max, 
    int N_points,
    double precision, 
    int Max_iteration = 10000
){
    vector<pair<T, int>> zeros_tm;
    T step = (max - min) / N_points;
    pair<T, int> nr;
    for(double x = min; x < max + step; x += step){
        nr = NewtonRaphson(f, f_dot, parameter, x, precision);
        if(!isnan(nr.first))  zeros_tm.push_back(nr);
    } 

    vector<pair<T, int>> zeros;
    bool check;
    for (auto tm : zeros_tm){
        check = true;
        for(auto z : zeros){
            if(fabs(tm.first - z.first) < precision) check = false;
        }

        if(check)   zeros.push_back(tm);
    }
    return zeros;
}

#endif

