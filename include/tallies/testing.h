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
#include <tallies/legendre_fet.hpp>
#include <ndarray.hpp>


void test_function(){
    std::cout<<"*********---- Starting of Testing Mode -----*********\n";
    Position LOW1(0. ,0., 0.);
    Position HIGH1(1., 1.0, 1.0);
    
    std::shared_ptr<CartesianFilter> B1 = std::make_shared<BoxPositionFilter>(LOW1, HIGH1);
        //std::make_shared<CylinderPositionFilter>(LOW1, 0.1, 0.5, 0.50, 0.5, 
          //                                  100, 10, 10, 
        //                                    CylinderPositionFilter::Orientation::X);
        //std::make_shared<MeshPositionFilter>(LOW1, HIGH1, 10, 10, 10);
        
    std::cout<<"-------\n";
    
    
    Direction D1( 1.0, 1.0, 0.0);
    Position R(-0.1, 0.9, 0.2);//<<<<<
    Tracker tktr(R, D1, true);

    std::vector<double> energy_f{0.0, 1.0, 2.0};
    std::shared_ptr<EnergyFilter> const_ef = std::make_shared<EnergyFilter>(energy_f);

    size_t FET_order = 6;
    std::vector<LegendreFET::Axis> vec_axes { LegendreFET::Axis::X };
    std::shared_ptr<ITally> LFET = 
                std::make_shared<LegendreFET>(FET_order, LegendreFET::Quantity::Flux, 
                                          LegendreFET::Estimator::Collision, "temp_FET", 
                                          vec_axes, B1, const_ef );

    //auto temp__ = B1->get_indices_tracklength(tktr, 1.0);
    
    std::cout<<"*********---- End of Testing Mode -----*********\n\n";
}


#endif