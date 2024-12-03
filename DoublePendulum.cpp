// g++ -o DoublePendulum DoublePendulum.cpp -lsfml-graphics -lsfml-window -lsfml-system -lm -std=c++17

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include "SFML/Graphics.hpp"
#include "ode_resolution.h"

using namespace std;
using namespace sf;

double f1 (
    double t,
    vector<double> x, 
    vector<double> parameter
){
    return x.at(2);
};

double f2(
    double t,
    vector<double> x, 
    vector<double> parameter
){
    return x.at(3);
}

double f1_dot (
    double t,
    vector<double> x, 
    vector<double> parameter
){
    return (- parameter.at(4) * (2 * parameter.at(2) + parameter.at(3)) * sin(x.at(0)) - parameter.at(3) * parameter.at(4) * sin(x.at(0) - 2 * x.at(1)) - 2 * sin(x.at(0) - x.at(1)) * parameter.at(3) * (pow(x.at(3), 2) * parameter.at(1) + pow(x.at(3), 2) * parameter.at(0) * cos(x.at(0) - x.at(1)))) / (parameter.at(0) * (2 * parameter.at(2) + parameter.at(3) - parameter.at(3) * cos(2 * (x.at(0) - x.at(1)))));
};

double f2_dot (
    double t,
    vector<double> x, 
    vector<double> parameter
){
    return (2 * sin(x.at(0) - x.at(1)) * (pow(x.at(2), 2) * parameter.at(0) * (parameter.at(2) + parameter.at(3)) + parameter.at(4) * (parameter.at(2) + parameter.at(3)) * cos(x.at(0)) + pow(x.at(3), 2) * parameter.at(1) * parameter.at(3) * cos(x.at(0) - x.at(1)))) / (parameter.at(1) * (2 * parameter.at(2) + parameter.at(3) - parameter.at(3) * cos(2 * (x.at(0) - x.at(1)))));
};

//  vector<double> theta = {theta1, theta2, theta1_dot, theta2_dot}
//  vector<double> x = {x1, y1, x2, y2}
//  vector<double> parameter = {l1, l2, m1, m2, g}

int main(int argc, char** argv)
{    
//Dati iniziali
    vector<double> par = {1., 1., 1., 1., 9.81};
    double theta_1_0;
    double theta_2_0;
    double omega_1_0;
    double omega_2_0;

    cout << "\ntheta 1 a t = 0 (deg):\t";
    cin >> theta_1_0;
    cout << "\ntheta 2 a t = 0 (deg):\t";
    cin >> theta_2_0;
    cout << "\nomega 1 a t = 0:\t";
    cin >> omega_1_0;
    cout << "\nomega 2 a t = 0:\t";
    cin >> omega_2_0;

    vector<double> theta_0 = {theta_1_0 * M_PI / 180, theta_2_0 * M_PI / 180, omega_1_0, omega_2_0};
    vector<double> x_0 = {par.at(0) * sin(theta_0.at(0)), par.at(0) * cos(theta_0.at(0)), par.at(0) * sin(theta_0.at(0)) + par.at(1) * sin(theta_0.at(0)), par.at(0) * cos(theta_0.at(0)) + par.at(1) * cos(theta_0.at(0))};
        
    double dt = 0.001;
    vector<function<double(double, vector<double>, vector<double>)>> f = {f1, f2, f1_dot, f2_dot};

    // apertura finestra grafica
    RenderWindow window(VideoMode(800, 800), "");
    window.setVerticalSyncEnabled(true);

    //dati visualizzazione
    double dt_display = 500; //tempo con cui si aggiorna la simulazione
    int scale = 130; //fattore di scala con cui viene visualizzata la simulazione
    int traccia_len = 3000; //numero di iterazioni massime visualizzate

    double t=0;
    vector<Vertex> L1;
    vector<Vertex> L2;
    vector<Vertex> Tr;
    vector<double> x;
    vector<double> theta;
    vector<pair<double, double>> traccia;

    // definizione assi
    Vertex X_top[2]; // definizione di un array di due vertici
    X_top[0].position = Vector2f(100, 100); // imposta la posizione del primo vertice
    X_top[0].color = Color::White; // imposta il colore del primo vertice
    X_top[1].position = Vector2f(700, 100); // imposta la posizione del secondo vertice
    X_top[1].color = Color::White; // imposta il colore del secondo vertice

    
    Vertex X_bottom[2]; // definizione di un array di due vertici
    X_bottom[0].position = Vector2f(100, 700); // imposta la posizione del primo vertice
    X_bottom[0].color  = Color::White; // imposta il colore del primo vertice
    X_bottom[1].position = Vector2f(700, 700); // imposta la posizione del secondo vertice
    X_bottom[1].color = Color::White; // imposta il colore del secondo vertice

    Vertex Y_left[2]; // definizione di un array di due vertici 
    Y_left[0].position = Vector2f(100, 100); // imposta la posizione del primo vertice
    Y_left[0].color  = Color::White; // imposta il colore del primo vertice
    Y_left[1].position = Vector2f(100, 700); // imposta la posizione del secondo vertice
    Y_left[1].color = Color::White; // imposta il colore del secondo vertice

    Vertex Y_right[2]; // definizione di un array di due vertici
    Y_right[0].position = Vector2f(700, 100); // imposta la posizione del primo vertice
    Y_right[0].color  = Color::White; // imposta il colore del primo vertice
    Y_right[1].position = Vector2f(700, 700); // imposta la posizione del secondo vertice
    Y_right[1].color = Color::White; // imposta il colore del secondo vertice
    
    // definizione di un vettore di vertici
    vector<Vertex> Axis;
    Axis.push_back(X_top[0]); // aggiunge il primo vertice superiore
    Axis.push_back(X_top[1]); // aggiunge il secondo vertice superiore
    Axis.push_back(X_bottom[0]); // aggiunge il primo vertice inferiore
    Axis.push_back(X_bottom[1]); // aggiunge il secondo vertice inferiore
    Axis.push_back(Y_left[0]); // aggiunge il primo vertice sinistro
    Axis.push_back(Y_left[1]); // aggiunge il secondo vertice sinistro
    Axis.push_back(Y_right[0]); // aggiunge il primo vertice destro
    Axis.push_back(Y_right[1]); // aggiunge il secondo vertice destro

    // ciclo principale del programma che continua fino a quando la finestra è aperta
    while (window.isOpen())
    {
        // pulisce il contenuto della finestra
        window.clear(); 
        
        // gestisce gli eventi all'interno del loop degli eventi
        Event event;
        while (window.pollEvent(event))
        {
            // se l'evento è di chiusura della finestra, chiude la finestra
            if (event.type == Event::Closed)
                window.close();
        }  
        L1.clear();
        L2.clear();
        Tr.clear();

        // aggiornamento della simulazione
        theta = RK4(f, dt, t, theta_0, par);
        x = {par.at(0) * sin(theta.at(0)), par.at(0) * cos(theta.at(0)), par.at(0) * sin(theta.at(0)) + par.at(1) * sin(theta.at(1)), par.at(0) * cos(theta.at(0)) + par.at(1) * cos(theta.at(1))};
        
        // definizione traccia
        if(traccia.size() < traccia_len)    traccia.push_back(make_pair(x.at(2), x.at(3)));      
        else if(traccia.size() == traccia_len){
            traccia.erase(traccia.begin());
            traccia.push_back(make_pair(x.at(2), x.at(3)));
        }

        for(int i = 0; i < traccia.size() - 1; i++){
            Vertex tr[2];
            tr[0].position = Vector2f(400 + scale * traccia.at(i).first, 400 + scale * traccia.at(i).second);
            tr[0].color = Color::Red;
            tr[1].position = Vector2f(400 + scale * traccia.at(i + 1).first, 400 + scale * traccia.at(i + 1).second);
            tr[1].color = Color::Red;

            Tr.push_back(tr[0]);
            Tr.push_back(tr[1]);
        }

        // definizione di L1 e L2
        Vertex l1[2];
        l1[0].position = Vector2f(400, 400);
        l1[0].color = Color::White;
        l1[1].position = Vector2f(400 + scale * x.at(0), 400 + scale * x.at(1));
        l1[1].color = Color::White;

        L1.push_back(l1[0]);
        L1.push_back(l1[1]);

        Vertex l2[2];
        l2[0].position = Vector2f(400 + scale * x.at(0), 400 + scale * x.at(1));
        l2[0].color = Color::White;
        l2[1].position = Vector2f(400 + scale * x.at(2), 400 + scale * x.at(3));
        l2[1].color = Color::White;

        L2.push_back(l2[0]);
        L2.push_back(l2[1]);

        // visualizzazione dei dati a schermo
        window.draw(&L1.at(0), L1.size(), Lines);
        window.draw(&L2.at(0), L2.size(), Lines);
        if(Tr.size() != 0)  window.draw(&Tr.at(0), Tr.size(), Lines);
        window.draw(&Axis.at(0), Axis.size(), Lines);

        //aggiornamento dati iniziali
        t += dt_display;
        theta_0 = theta;

        window.display();
    }

    return 0;
}
