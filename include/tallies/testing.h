#ifndef TESTING_H
#define TESTING_H

#include <memory>
#include <array>
#include <vector>

#include <tallies/box_position_filter.h>
#include <tallies/energy_filter.h>
#include <vector>
#include <simulation/tracker.hpp>
#include <tallies/cylinder_position_filter.h>
#include <tallies/cartesian_filter.h>
#include <tallies/mesh_position_filter.h>
#include <tallies/ITally.h>
#include <tallies/general_tally.h>
#include <ndarray.hpp>


void test_function(){
    std::cout<<"*********---- Starting of Testing Mode -----*********\n";
    Position LOW1(0. ,0., 0.);
    Position HIGH1(1., 1.0, 1.0);
    std::shared_ptr<PositionFilter> pos_f = 
        std::make_shared<BoxPositionFilter>(LOW1, HIGH1);
    
    std::vector<double> ebound {1., 2., 3.};
    std::shared_ptr<EnergyFilter> ef = std::make_shared<EnergyFilter>(ebound);

    std::shared_ptr<ITally> ta = 
        std::make_shared<GeneralTally>(ITally::Quantity::Absorption, 
                                ITally::Estimator::Collision,
                                pos_f, ef);
    //ta->score_tally();
    std::cout<<"*********---- End of Testing Mode -----*********\n\n";
    
}



void test_function2(){
    std::cout<<"*********---- Starting of Testing Mode -----*********\n";
    Position LOW1(0. ,0., 0.);
    Position HIGH1(1., 1.0, 1.0);
    
    std::shared_ptr<CylinderPositionFilter> B1 = 
        std::make_shared<CylinderPositionFilter>(LOW1, 0.1, 0.5, 0.50, 0.5, 
                                            100, 10, 10, 
                                            CylinderPositionFilter::Orientation::X);
        //std::make_shared<MeshPositionFilter>(LOW1, HIGH1, 10, 10, 10);
    std::cout<<"-------\n";
    
    
    Direction D1( 1.0,1.0, 1.0);
    Position R(0.1, 0.1, 0.2);//<<<<<
    Tracker tktr(R, D1, true);

    std::array<int, 3> index;
    std::cout<<", "<<B1->Nx()<<"\n";

    if ( B1->get_indices(tktr, index)){
        std::cout<<"Found: Nx = "<<index[0]<<", Ny = "<<index[1]<<", Nz = "<<index[2]<<"\n";
    }else{
        std::cout<<"Not found\n";
    }
    std::cout<<"-------\n";

    
    std::cout<<"*********---- End of Testing Mode -----*********\n\n";
}


/*void test_function(){
    std::cout<<"*********---- Starting of Testing Mode -----*********\n";
    

    std::cout<<"*********---- End of Testing Mode -----*********\n\n";
    
}*/

#endif