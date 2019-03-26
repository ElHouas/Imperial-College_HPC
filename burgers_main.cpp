#include <chrono>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>

#include "Model.h"
#include "Burgers.h"
#include "MPI_functions.h"

using namespace std;

//Linux commnad to compile and run the code
//g++ burgers.cpp -o a
//

int main(int argc, char* argv[]) {

    Model m(argc, argv);
    
    // Call code to initialise the problem here
    
    Burgers B(m);   
    
    //typedef std::chrono::high_resolution_clock hrc;
    //typedef std::chrono::milliseconds ms;
    //hrc::time_point start = hrc::now();

    // Call code to perform time integration here
    
    B.iVelocity();
    
    //B.vEnergy();
    
        
    /* MPI variables */
        // Initialize the MPI environment
    MPI_Init(NULL, NULL);
                    
        // Find out rank, size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //this is the id of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size); //this is the num of the process 
    
    B.vIntegral(rank, size);
    B.vEnergy();
    B.vWrite(rank, size);
    
    MPI_Finalize();
    //hrc::time_point end = hrc::now();

    // Calculate final energy and write output
    
    
    
    //

    
    
    return 0;
            }