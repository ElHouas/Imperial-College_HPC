#ifndef CLASS_BURGERS
#define CLASS_BURGERS
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <math.h>
#include "Model.h"

//Inheritance in order to allow Burger class access the values and functions of Model

class Burgers : public Model{
private:
        
        Model m;
        double **u;
        double **v;
        double **u_temp;
        double **u_new;//Array to store the computed the velocity field in 
        double **u_buf;//
        double **u_masterbuf;//Master array that receives all the data to 
        double **u_finalbuf;//Final array that gathers all the data
        double **v_temp;
        double **v_new;//Array to store the computed the velocity field in 
        double **v_buf;//
        double **v_masterbuf;//Master array that receives all the data to 
        double **v_finalbuf;//Final array that gathers all the data
        double r;
        double E;
    
    public:
         Burgers(Model m);
        ~Burgers();// Model class destructor#
        void iVelocity();
        void vIntegral(int rank, int size);
        void vEnergy();
        void vWrite(int rank, int size);
    };

#endif