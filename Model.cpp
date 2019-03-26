//Implentation of Model class

#include <stdexcept>
#include <iostream>
#include "Model.h"
#include <cstdlib>
//#include <boost/program_options.hpp>


//NOT GOOD PRACTICE TO CALL METHOD INSIDE CONSTRUCTOR
Model:: Model(int argc, char* argv[]) {
    ParseParameters(argc, argv);
}

Model:: Model(){ 
    //cout << "Class Model constructed" << endl;
}

Model:: ~Model(){ 
    //cout << "Class Model destroyed" << endl;
}


void Model::ParseParameters(int argc, char* argv[]){
       
    
    //cout << "You have entered " << argc 
    //     << " arguments:" << "\n"; 
  
    //for (int i = 0; i < argc; ++i){
    //    cout << argv[i] << "\n"; 
    //}
    // An alias to reduce typing
    /*namespace po = boost::program_options;
    //atoi stdlib function to convert string to int
    po::options_description opts("Compute Burgers equation.");
    opts.add_options()
        ("Test case", po::value<int>()->default_value(1),
				 "Test case .")
        ("Nx",  po::value<int>()->default_value(15),
				 "Minimum value of numbers to generate.")
        ("Ny",  po::value<int>()->default_value(15),
				 "Maximum value of numbers to generate.")
		("Nt",  po::value<int>()->default_value(100),
				 "Maximum value of numbers to generate.");
        //("help",       "Print help message.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);
    

    if (vm.count("help")) {
        cout << "Performs Burgers equation computation "
             << "Plot in Python" << endl;
        cout << opts << endl;
    }*/

    
    //int Test = vm["Test case"].as<int>();
    
    int Test = atoi(argv[1]);
    
    cout << "Selected Test case:" << Test << endl;
    switch (Test) {
        case 1:
            {
            cout << "Selected Test case 1 with:" << endl;
            cout << "ax = 0, ay = 0, b = 0, c = 1" << endl;
                this->SetAx(0);
                this->SetAy(0);
                this->SetB(0);
                this->SetC(1); 
            break;
            }
        case 2:
            {
            cout << "Selected menu item 2 with:" ;
            cout << "ax = 1, ay = 0, b = 0, c = 0" << endl;
                this->SetAx(1);
                this->SetAy(0);
                this->SetB(0);
                this->SetC(0); 
            break;
            }
        case 3:
            {
            cout << "Selected menu item 3 with:" ;
            cout << "ax = 0, ay = 1, b = 0, c = 0" << endl;
                this->SetAx(0);
                this->SetAy(1);
                this->SetB(0);
                this->SetC(0);
            break;
            }
        case 4:
            {
            cout << "Selected menu item 4 with:" ;
            cout << "ax = 1.0, ay = 0.5, b = 1.0 and c = 0.02" << endl;
                this->SetAx(1.0);
                this->SetAy(0.5);
                this->SetB(1.0);
                this->SetC(0.02);
            break;
            }  
        default:
            {
            cout << "Unknown menu item!" << endl;
            }
        
    }
    
    // The user enters initial conditions of the problem

    this->SetNx(atoi(argv[2]));
    //this->SetNx(vm["Nx"].as<int>());
    cout << "Entered positive number of points in x direction: "<< this->GetNx() <<endl;
    
    this->SetNy(atoi(argv[3]));
    //this->SetNy(vm["Ny"].as<int>());
    cout << "Entered positive number of points in y direction: "<< this->GetNy() <<endl;
    
    this->SetNt(atoi(argv[4]));
    //this->SetNt(vm["Nt"].as<int>());
    cout << "Entered positive number of time steps: "<< this->GetNt() << endl;    
    

    cout << "The problem has established parameters L=10 and T=1: "<< endl;
    this->SetLx(10);
    this->SetLy(10);
    this->SetT(1);
    
    this->SetX0(this->GetLx()/2);
    cout << "Entered positive initial X: "<< this->GetX0() <<endl;
    
    
    this->SetY0(this->GetLy()/2);
    cout << "Entered positive initial Y: "<< this->GetY0() <<endl;
    
    this->SetDt(this->GetT()/this->GetNt());
    this->SetDx(this->GetLx()/this->GetNx());
    this->SetDy(this->GetLy()/this->GetNy());
  
}