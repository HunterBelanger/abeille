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

#include <boost/container/static_vector.hpp>

#include <ndarray.hpp>

namespace bc = boost::container;

using StaticVector = bc::static_vector<size_t, 9>;

void test_function(){
  std::cout<<"*********---- Starting of Testing Mode -----*********\n";
  boost::container::static_vector<size_t, 3> parth{1, 2, 3};
  
  std::cout<<"Made the static vector.\n";
  StaticVector part{11, 2, 3};
  std::cout<<"-->"<<part[0]<<"\n";
  for ( size_t x : parth){
    std::cout<<x<<"\t";
  }
  std::cout<<std::endl;
  //parth.push_back(4);
  NDArray<double> temp_nd_parth;
  std::vector<size_t> id{1, 2,3};
  temp_nd_parth.reallocate(id);
  std::cout<<"We got it here.\n";
  NDArray<double> temp_nd;
  temp_nd.reallocate(parth);
  temp_nd.fill(1.0);
  temp_nd({0,0,0}) = 10;

  boost::container::static_vector<size_t, 3> num{0, 0, 0};
  std::cout<<"Value = "<<temp_nd(num)<<"\n";

  std::cout<<"*********---- End of Testing Mode -----*********\n";

}



#endif