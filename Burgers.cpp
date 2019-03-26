#include <stdexcept>
#include <iostream>
#include "Model.h"
#include "Burgers.h"
#include "MPI_functions.h"
#include <cstdlib>
#include <iomanip>      
#include <mpi.h>
#include <math.h>


Burgers::Burgers(Model m) {
    
    this-> m = m;
}

Burgers::~Burgers(){ // Model class destructor#
    //cout << "Burguers class destroyed" << endl;
}


void Burgers::iVelocity(){

            int Nx = m.GetNx();
            int Ny = m.GetNy();
            double dx = m.GetDx();
            double dy = m.GetDy();
            
            u = alloc_2d_int(Nx,Ny);
            v = alloc_2d_int(Nx,Ny);
            
             for(unsigned i = 0; i < Nx; ++i){
                    
                    for(unsigned j = 0; j < Ny; ++j){
                        
                        r = sqrt( (pow(i*dx- m.GetX0(),2) )  + (pow(j*dy- m.GetY0(),2) ) );
                        
                        if(r<=1){                         
                         
                            u[i][j] = 2.0*pow((1-r),4)*((4*r)+1.0);
                           
                            v[i][j] = 2.0*pow((1-r),4)*((4*r)+1.0);
                        } else {
                            
                            u[i][j] = 0;
                         
                            v[i][j] = 0;
                        }
                                                
                        cout.precision(4);
                        cout << u[i][j] << " ";
                    }
                    cout << endl;
                    
                }
        }
        

void Burgers::vIntegral(int rank, int size){
    
            //Declare derivivatives
            
            double Udx;
            double Udy;
            double Udx_2;
            double Udy_2;
            
            double Vdx;
            double Vdy;
            double Vdx_2;
            double Vdy_2;
            
            //Variables from Model to be used
            
            int Nt = m.GetNt();
            int Nx = m.GetNx();
            int Ny = m.GetNy();
            double dt = m.GetDt();  
            double dx = m.GetDx();
            double dy = m.GetDy();
            double ax = m.GetAx(); 
            double ay = m.GetAy();
            double b = m.GetB();
            double c = m.GetC(); //Kinematic viscosity 
            
            
        /* MPI variables for initializing grid and MPI message request*/
        MPI_Comm mygrid;
        MPI_Request request1, request2, request3, request4, request5, request6, request7, request8;
        MPI_Request request9, request10, request11, request12, request13, request14, request15, request16;
        
        /* =================================================================================================
           *	START: Define our MPI Virtual Topology
           * ================================================================================================= */

          /* Set sizei and sizej to factors of MPI size */

        int P = 2;
        int dims[P] = {0,0};
        MPI_Dims_create(size, P, dims);
        int sizes[P] = {dims[0], dims[1]}; 
        
        
        /* Create virtual grid topology in 2D, 
           *  set to communicator 'grid_comm' and set coords */
            int periods[P] = {0, 0};
            int reorder = 0; //Runtime system allowed to reorder processes
            int coords[P];
            int ndims = 2; // Number of processes in the dimension i

            MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &mygrid);
            MPI_Cart_coords(mygrid, rank, ndims, coords);

        /* Now we set ranki and rankj to their grid topology */
          int ranki = coords[0];
          
          //cout <<"ranki"<<ranki<< endl;
          
          int rankj = coords[1];
          
          //cout <<"rankj"<<rankj<< endl;


        /* We use MPI Cart shifts to define four nearest neighbouring
           * processes to each process in our virtual MPI Topology   */
        
        int disp = 1; //Define the displacement, the distance to the neighbor
        int directionx = 0;//X-direction associated to 0
        int directiony = 1;//Y-direction associated to 0
        int rankleft, rankright, rankdown, rankup;
        MPI_Cart_shift(mygrid, directionx, disp, &rankleft, &rankright);
        MPI_Cart_shift(mygrid, directiony, disp, &rankdown, &rankup);

      /* =================================================================================================
       *	END: Define our cartesian virtual Topology
       * ================================================================================================= */
            
      /* =================================================================================================
       *    START: Define Number of points on each process and where they start/end
       * ===========================================================================
       * In the event that we have number of points not divisible by sizei (sizej). We shall 
       * increase the number of points on the first sizei-1 (or sizej-1) processes such that
       * all these arrays are of equal size. The last process shall have the remainder.
       * ================================================================================================= */

      int Nx_extended;
      int Ny_extended;
      
      /* This function sets Nx_extended (or Ny_extended) to a number divisible by p */
      padding (&Nx_extended, &Ny_extended, Nx, Ny, sizes[0], sizes[1]);
      
      /* Nx_p, Ny_p shall be chunksize per process in the i and j direction respectively */
      
      int Nx_P = Nx_extended/sizes[0];
      
      int Ny_P = Ny_extended/sizes[1];
      
      /* We define the point at which each process starts/ends in direction i and j 
       * This will not work if Nx_extended - Nx >= Nx_p (Ny_extended - Ny >= Ny_p) which should
       * not be the case - only if a high and prime number of processes are chosen */
      
      int starti = ranki*Nx_P;
      
      int endi = (ranki+1)*Nx_P;
      
    
      if (endi > Nx)
      {
        Nx_P = Nx_P - (endi-Nx);
        endi = Nx;
      }
      if (endi < starti) 
      {
        cout << "Bad number of processors in i direction!"<< endl;
        exit(0);
      }
      
      int startj = (rankj*Ny_extended)/sizes[1];
      int endj = ((rankj+1)*Ny_extended)/sizes[1];
      if (endj > Ny)
      {
        Ny_P = Ny_P - (endj-Ny);
        endj = Ny;
      }
      if (endj < startj)
      {
        cout <<"Bad number of processors in j direction"<< endl; 
        exit(0);
      }
      
      if (rank==0){
          cout <<"rank" <<rank<< endl;
          cout <<"Nx_extended" <<Nx_extended<< endl;
          cout <<"Nx_P" <<Nx_P<< endl;
          cout <<"sizes[0]"<<sizes[0]<< endl;
          cout <<"sizes[1]"<<sizes[1]<< endl;
          cout <<"Ny_extended" <<Ny_extended<< endl;
          cout <<"Ny_P" <<Ny_P<< endl;
          cout <<"starti"<<starti<< endl;
          cout <<"endi"<<endi<< endl;
      }
      
      
      
      /* =================================================================================================
       *    END: Define Number of pixels on each process. 
       * ================================================================================================= */

        /* We must define a new datatype to send non-contiguous data*/
          MPI_Datatype noncontiguous;
          MPI_Type_vector (Nx_P, 1, (Ny_P+2), MPI_DOUBLE,&noncontiguous);
          MPI_Type_commit (&noncontiguous); 
        
        
         /*Initialize 2d array to create a extended domain with ghost shells*/
        u_temp = alloc_2d_int(Nx_P+2,Ny_P+2);
        v_temp = alloc_2d_int(Nx_P+2,Ny_P+2);

        /*Initialize 2d array to compute the new velocity field each time step*/
        u_new = alloc_2d_int(Nx_P+2,Ny_P+2);  
        v_new = alloc_2d_int(Nx_P+2,Ny_P+2);

        /*Initialize 2d array to store the data for each rank*/
        u_buf = alloc_2d_int(Nx_P,Ny_P);
        v_buf = alloc_2d_int(Nx_P,Ny_P);
         
        /*Initialize 2d array to store the data of each buf*/
        u_masterbuf = alloc_2d_int(Nx,Ny);    
        v_masterbuf = alloc_2d_int(Nx,Ny);
        
        /*Initialize 2d array to store the gather all final he data */
        u_finalbuf = alloc_2d_int(Nx,Ny);       
        v_finalbuf = alloc_2d_int(Nx,Ny);
          
          
          /* =================================================================================================
           *	START: Read in data to masterbuf, this is done on all processes so they all have a copy.
           * ================================================================================================= */
                if (rank == 0)
                {
                    cout <<"sizei:"<< sizes[0] << "sizej:"<< sizes[1] << endl;
                    cout <<"Processing " << Nx << " x " << Ny << " velocity field"<< endl;

               }
               
 
                for (unsigned i=0; i < Nx; i++)
                {
                    for (unsigned j=0; j < Ny; j++)
                    {
                      u_masterbuf[i][j]=u[i][j];
                      v_masterbuf[i][j]=v[i][j];
                    }
                }     
               
            /* =================================================================================================
             *	END: Read in data to masterbuf.
             * ================================================================================================= */
          
          /* Set buf to appropriate values of masterbuf for each process */
          pack_buf (Nx_P,Ny_P, Nx, Ny, u_buf, u_masterbuf, starti, endi, startj, endj);
          pack_buf (Nx_P,Ny_P, Nx, Ny, v_buf, v_masterbuf, starti, endi, startj, endj);
          
          
          if(rank==0){
          cout << "RANK " << rank << endl;
          for (unsigned i=0; i<Nx_P; i++)
          {
              for (unsigned j=0; j<Ny_P; j++)
            {
               cout << u_buf[i][j] << " ";
               cout << v_buf[i][j] << " ";
            }
            cout << endl;
          }
          cout << endl;
          }
          

            
            
          /*Assign the velocity field to the temporal array with the ghost shells*/
          for (unsigned i=1; i < Nx_P+1; i++)
        {
            for (unsigned j=1; j < Ny_P+1; j++)
            {
               u_temp[i][j] = u[i-1][j-1];
               v_temp[i][j] = v[i-1][j-1];
            }
        }
         
         
        /*Main iteration for each time step*/
          for(unsigned n = 0; n < Nt; ++n){
               
               /* =================================================================================================
               * =================================================================================================
               *	START: Main Loop. This shall be the computationally expensive part of the code. START TIME.
               * ================================================================================================= 
               * ================================================================================================= */
                
            /* =================================================================================================
             *	START: Halo Exchange.
             * ================================================================================================= */
            
             /*U velocity field Halo Exchange.*/
            
            
            /* Halo swapping horizontally (using derived datatype for noncontiguous data)   */
            MPI_Issend(&(u_temp[1][Ny_P]), 1, noncontiguous, rankright, 1, MPI_COMM_WORLD, &request1);
            MPI_Irecv(&(u_temp[1][0]), Nx_P, noncontiguous, rankleft, 1, MPI_COMM_WORLD, &request2);
            MPI_Issend(&(u_temp[1][1]), 1, noncontiguous, rankleft, 0, MPI_COMM_WORLD, &request3);
            MPI_Irecv(&(u_temp[1][Ny_P+1]), Nx_P, noncontiguous, rankright, 0, MPI_COMM_WORLD, &request4);
            
            /* Halo swapping vertically (contiguous data)*/
            MPI_Issend(&(u_temp[1][1]), Ny_P, MPI_DOUBLE, rankup, 2, MPI_COMM_WORLD, &request5);
            MPI_Irecv(&(u_temp[Nx_P+1][1]), Ny_P, MPI_DOUBLE, rankdown, 2, MPI_COMM_WORLD, &request6);
            MPI_Issend(&(u_temp[Nx_P][1]), Ny_P, MPI_DOUBLE, rankdown, 3, MPI_COMM_WORLD, &request7);
            MPI_Irecv(&(u_temp[0][1]), Ny_P, MPI_DOUBLE, rankup, 3, MPI_COMM_WORLD, &request8);
           
            
            /* All Halo communication must be completed
            before applying next iteration MPI_Wait(&request5, MPI_STATUS_IGNORE);*/
            MPI_Wait(&request5, MPI_STATUS_IGNORE);
            MPI_Wait(&request6, MPI_STATUS_IGNORE);
            MPI_Wait(&request7, MPI_STATUS_IGNORE);
            MPI_Wait(&request8, MPI_STATUS_IGNORE); 
            MPI_Wait(&request1, MPI_STATUS_IGNORE);
            MPI_Wait(&request2, MPI_STATUS_IGNORE);
            MPI_Wait(&request3, MPI_STATUS_IGNORE);
            MPI_Wait(&request4, MPI_STATUS_IGNORE);
            
             

            
            /*V velocity field Halo Exchange.*/
            
            /* Halo swapping horizontally(using derived datatype for noncontiguous data)  */
            MPI_Issend(&(v_temp[1][Ny_P]), 1, noncontiguous, rankright, 1, MPI_COMM_WORLD, &request9);
            MPI_Irecv(&(v_temp[1][0]), Nx_P, noncontiguous, rankleft, 1, MPI_COMM_WORLD, &request10);
            MPI_Issend(&(v_temp[1][1]), 1, noncontiguous, rankleft, 0, MPI_COMM_WORLD, &request11);
            MPI_Irecv(&(v_temp[1][Ny_P+1]), Nx_P, noncontiguous, rankright, 0, MPI_COMM_WORLD, &request12);
            
            /* Halo swapping vertically (contiguous data)*/
            MPI_Issend(&(v_temp[1][1]), Ny_P, MPI_DOUBLE, rankup, 6, MPI_COMM_WORLD, &request13);
            MPI_Irecv(&(v_temp[Nx_P+1][1]), Ny_P, MPI_DOUBLE, rankdown, 6, MPI_COMM_WORLD, &request14);
            MPI_Issend(&(v_temp[Nx_P][1]), Ny_P, MPI_DOUBLE, rankdown, 7, MPI_COMM_WORLD, &request15);
            MPI_Irecv(&(v_temp[0][1]), Ny_P, MPI_DOUBLE, rankup, 7, MPI_COMM_WORLD, &request16);
            
            /* All Halo communication must be completed
            before applying next iteration  */
            MPI_Wait(&request9, MPI_STATUS_IGNORE);
            MPI_Wait(&request10, MPI_STATUS_IGNORE);
            MPI_Wait(&request11, MPI_STATUS_IGNORE);
            MPI_Wait(&request12, MPI_STATUS_IGNORE);
            MPI_Wait(&request13, MPI_STATUS_IGNORE);
            MPI_Wait(&request14, MPI_STATUS_IGNORE);
            MPI_Wait(&request15, MPI_STATUS_IGNORE);
            MPI_Wait(&request16, MPI_STATUS_IGNORE);  
            
            
           if(rank==1){
              
                cout << "Temporal AFTER RECEIVE Rank: " << rank <<endl ;
              
                for (unsigned i=0; i < Nx_P+2; i++)
            {
                for (unsigned j=0; j < Ny_P+2; j++)
                {
                  cout<< u_temp[i][j] << " ";
                }
                cout << endl;
            }
            }
            cout << endl;
         
                
            /* =================================================================================================
            *	END: Halo Exchange.
            * ================================================================================================= */
                
                
                /* Perform main iteration */
                
                for(unsigned i = 1; i < Nx_P+1; ++i){
                    
                    
                    for(unsigned j = 1; j < Ny_P+1; ++j){
                       
                        /*Main iteration for u velocity field*/
                        
                        Udx = (u_temp[i][j]-u_temp[i-1][j])/dx;
                        Udy = (u_temp[i][j]-u_temp[i][j-1])/dy;
                      

                        Udx_2 = (u_temp[i+1][j]-2*u_temp[i][j]+u_temp[i-1][j])/(pow(dx,2));
                        Udy_2 = (u_temp[i][j+1]-2*u_temp[i][j]+u_temp[i][j-1])/(pow(dy,2));
                        
                         
                         u_new[i][j] = u_temp[i][j]+ dt*(c*(Udx_2 + Udy_2)- (ax+b*u_temp[i][j])*Udx - (ay+b*v_temp[i][j])*Udy);
                        
                         /*Main iteration for v velocity field*/
                         
                         Vdx = (v_temp[i][j]-v_temp[i-1][j])/dx;
                         
                         Vdy = (v_temp[i][j]-v_temp[i][j-1])/dy;
            
                         Vdx_2 = (v_temp[i+1][j]-2*v_temp[i][j]+v_temp[i-1][j])/(pow(dx,2));
                         Vdy_2 = (v_temp[i][j+1]-2*v_temp[i][j]+v_temp[i][j-1])/(pow(dy,2));
                         
                         v_new[i][j] =  v_temp[i][j] + dt*(c*(Vdx_2 + Vdy_2) - (ax+b*u_temp[i][j])*Vdx - (ay+b*v_temp[i][j])*Vdy);
                         
                         /*Assign values to temporal to loop another time*/
                         u_temp[i][j] = u_new[i][j];
                         v_temp[i][j] = v_new[i][j];
 

                    }
                   
                }
                
                         
               
            }         

            /* =================================================================================================
             * =================================================================================================
             *	END: Main Loop. END TIME
             * ================================================================================================= 
             * ================================================================================================= */
                
                for (unsigned i=1;i<Nx_P+1;i++)
                {
                    for (unsigned j=1;j<Ny_P+1;j++)
                    {
                      u_buf[(i-1)][(j-1)]=u_temp[i][j];
                      v_buf[(i-1)][(j-1)]=v_temp[i][j];
   
                    }
                
                }
                  
                          
                        
              /* Set appropriate value of masterbuf to values of buf for each process */
              unpack_buf (Nx_P, Ny_P, Nx, Ny, u_buf, u_masterbuf, starti, endi, startj, endj);
              unpack_buf (Nx_P, Ny_P, Nx, Ny, v_buf, v_masterbuf, starti, endi, startj, endj);
              
              /*cout << "HIIII Masterbuf" <<rank<<endl ;*/
               
                for (unsigned i=0; i < Nx; i++)
                {
                    for (unsigned j=0; j < Ny; j++)
                    {
                      cout.precision(5);
                    }
                    cout << endl;
                }
                cout << endl;
              
              /* rank 0 will now be given the entire velocity by this reductions */
              MPI_Reduce(&u_masterbuf[0][0], &u_finalbuf[0][0], (Nx*Ny), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              MPI_Reduce(&v_masterbuf[0][0], &v_finalbuf[0][0], (Nx*Ny), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

              
}
        

void Burgers::vEnergy(){
    
                double Area = m.GetLx()*m.GetLy();
                
                for(unsigned i = 0; i < m.GetNx(); ++i){
                    
                    for(unsigned j = 0; j < m.GetNy(); ++j){
                            
                         E += 0.5*(pow(u_finalbuf[i][j],2) + pow(v_finalbuf[i][j],2));
                    }
                }
                E = E*Area;
                
                cout << "Energy: " << E << endl;
        
    }




void Burgers::vWrite(int rank, int size){
    
                if (rank==0){
                  
                cout << "HIIII final" <<rank<<endl ;
               
                for (unsigned i=0; i < m.GetNy(); i++)
                {
                    for (unsigned j=0; j < m.GetNx(); j++)
                    {
                      cout.precision(5);
                      
                      cout << " " << u_finalbuf[i][j] << " ";
                      //cout << " " << v_finalbuf[i][j] << " ";
                    }
                    cout << endl;
                }
                cout << endl;
              
                        
            
                ofstream vOut("MPI_velocity_field.txt", ios::out | ios::trunc);
                vOut.precision(4);
                
                for(unsigned i = 0; i < m.GetNy(); ++i){
                    
                    for(unsigned j = 0; j < m.GetNx(); ++j){

                        vOut << setw(15) << i*m.GetDy(); 
                        vOut << setw(20) << j*m.GetDx() ;
                        vOut << setw(20) << u_finalbuf[i][j];
                        vOut << setw(20) << v_finalbuf[i][j] << endl; 
                    }
                }
           
               vOut.close();
                
              }

        }