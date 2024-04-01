#ifndef TESTING_H
#define TESTING_H

#include <tallies/box_position_filter.h>
#include <memory>
#include <array>
#include <tallies/energy_filter.h>
#include <vector>
#include <simulation/tracker.hpp>
#include <tallies/cylinder_position_filter.h>
#include <tallies/cartesian_filter.h>
#include <tallies/mesh_position_filter.h>

#include <ndarray.hpp>

class parent{
    protected:
        NDArray<int> FF;
    public:
    parent() = default;
    parent( size_t i) : FF(){
        FF.reallocate({i});
        FF.fill(0);
    }

    ~parent()= default;

    virtual void info(){
        std::cout<<"+++\n";
        for (size_t i = 0; i< FF.size(); i++)
            std::cout<<FF[i]<<"\t";
                std::cout<<"+++\n";
    }
};

class child : public parent{
    public:
        child( size_t i) {
            FF.reallocate({i});
            FF.fill(2);
        }
        ~child() = default;

};

void test_function(){
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

    std::cout<<"Upper class nd array reallocation";
    std::shared_ptr<parent> ch = std::make_shared<child>(2);
    ch->info();
    std::cout<<"*********---- End of Testing Mode -----*********\n\n";
}

#endif