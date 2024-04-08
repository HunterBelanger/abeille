#ifndef TESTING_H
#define TESTING_H

#include <memory>
#include <array>
#include <vector>

#include <tallies/box_position_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <vector>
#include <simulation/tracker.hpp>
#include <tallies/cylinder_position_filter.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/mesh_position_filter.hpp>
#include <tallies/itally.hpp>
#include <tallies/general_tally.hpp>
#include <ndarray.hpp>


void test_function(){
    std::cout<<"*********---- Starting of Testing Mode -----*********\n";
    Position LOW1(0. ,0., 0.);
    Position HIGH1(1., 1.0, 1.0);
    
    std::shared_ptr<BoxPositionFilter> B1 = std::make_shared<BoxPositionFilter>(LOW1, HIGH1);
        //std::make_shared<CylinderPositionFilter>(LOW1, 0.1, 0.5, 0.50, 0.5, 
          //                                  100, 10, 10, 
        //                                    CylinderPositionFilter::Orientation::X);
        //std::make_shared<MeshPositionFilter>(LOW1, HIGH1, 10, 10, 10);
        
    std::cout<<"-------\n";
    
    
    Direction D1( 1.0, 1.0, 0.0);
    Position R(-0.1, 0.9, 0.2);//<<<<<
    Tracker tktr(R, D1, true);

    auto temp__ = B1->get_indices_tracklength(tktr, 1.0);
    
    std::cout<<"*********---- End of Testing Mode -----*********\n\n";
}


#endif