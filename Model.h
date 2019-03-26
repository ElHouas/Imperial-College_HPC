#ifndef CLASS_MODEL
#define CLASS_MODEL
#include <stdexcept>
#include <iostream>

using namespace std;

//#include "mpi.h"

class Model {
	public:
		Model(int argc, char* argv[]);  // Model class constructor
        Model();  // Default Model class constructor
        ~Model();// Model class destructor
        
        void PrintParameters();

        bool IsValid();

        // Getters
		bool   IsVerbose() const { return verbose; }
        bool   IsHelp()    const { return help; }
        double GetX0()     const { return x0; }
        double GetY0()     const { return y0; }
        double GetLx()     const { return Lx; }
        double GetLy()     const { return Ly; }
        double GetT()      const { return T; }
        int    GetNx()     const { return Nx; }
        int    GetNy()     const { return Ny; }
        int    GetNt()     const { return Nt; }
        double GetDx()     const { return dx; }
        double GetDy()     const { return dy; }
        double GetDt()     const { return dt; }
        double GetAx()     const { return ax; }
        double GetAy()     const { return ay; }
        double GetB()      const { return b; }
        double GetC()      const { return c; }
        
        // Add any other getters here...
        
        
        void SetX0(double x0){ 
                    if (x0 < 0){ 
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Initial X is negative");
                    }else{
                        this->x0 = x0; 
                    }
        }
        
        void SetY0(double y0){
                    if (y0 < 0){ 
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Initial Y is negative");
                    }else{
                        this->y0 = y0; 
                    }
        }
            

        void SetLx(double Lx){ 
                    if (Lx < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Longitude X is negative");
                    }else{
                        this->Lx = Lx; 
                    }
        }    
        
        void SetLy(double Ly){ 
                    if (Ly < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Longitude Y is negative");
                    }else{
                        this->Ly = Ly; 
                    }
        }
        
        void SetT(double T){
                   if (T < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Time is negative");
                    }else{
                        this->T = T; 
                    } 
        }
          
        void SetNx(int Nx){
                    if (Nx < 0){
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Number of X points is negative");
                    }else{
                        this->Nx = Nx; 
                    }
        }
 
        void SetNy(int Ny){ 
                    if (Ny < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Number of Y points is negative");
                    }else{
                        this->Ny = Ny; 
                    }
        }
        
        void SetNt(int Nt){ 
                    if (Nt < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Number of time steps is negative");
                    }else{
                        this->Nt = Nt; 
                    }
        }
        
        void SetDx(double dx){ 
                    if (dx < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Delta X is negative");
                    }else{
                        this->dx = dx; 
                    }
        }
        
        void SetDy(double dy){ 
                    if (dy < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Delta Y is negative");
                    }else{
                        this->dy = dy; 
                    }
        }
        
        void SetDt(double dt){ 
                    if (dt < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Delta T is negative");
                    }else{
                        this->dt = dt; 
                    }
        }
        
        void SetAx(double ax){ 
                    if (ax < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Ax is negative");
                    }else{
                        this->ax = ax; 
                    } 
        }
        
        void SetAy(double ay){
                    if (ay < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("Ay is negative");
                    }else{
                        this->ay = ay; 
                    } 
        }
        
        void SetB(double b){
                   if (b < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("B is negative");
                    }else{
                        this->b = b; 
                    }
        }
        
        void SetC(double c){
                   if (c < 0){  
                    // Triggers an exception. We can catch this in the calling program
                    // using a 'try' and 'catch' block.
                        throw std::invalid_argument("C is negative");
                    }else{
                        this->c = c; 
                    }
        }
        
        
	private:
        void ParseParameters(int argc, char* argv[]);
        void ValidateParameters();

	    bool verbose;
        bool help;

	    // Numerics
	    double x0;
	    double y0;
	    double Lx;
	    double Ly;
	    double T;
	    int    Nx;
	    int    Ny;
	    int    Nt;
	    double dx;
	    double dy;
	    double dt;

	    // Physics
	    double ax;
	    double ay;
	    double b;
	    double c;

        // Add any additional parameters here...
};

#endif