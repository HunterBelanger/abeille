#ifndef ENREGY_FILTER_H
#define ENREGY_FILTER_H

#include <iostream>
#include <vector>
#include <string>

#include "filter.h"


class EnergyFilter : public Filter{
    public:
        EnergyFilter(const std::vector<double>& energy_bounds_) : energy_bounds(energy_bounds_) {
            // Check the size of energy_bounds should be multiple of 2
            // since any bounds has 2 end-point-boundary-elements to define the range. 
            if ( energy_bounds.size() < 2){
                std::cout<<"Fatal error, Size of energy bounds is less than 2.\n";
            }
        }
        
        ~EnergyFilter() = default;

        bool get_index(const double& E, std::size_t& index_E)const {
            //index_E = -1;
            for ( std::size_t i = 0; i<energy_bounds.size() - 1; i++){
                if (E >= energy_bounds[i] && E<= energy_bounds[i+1]){
                    index_E = i;
                    return true;
                    break;
                }
            }
            return false;
        }

        unsigned int get_size() { return static_cast<unsigned int> (energy_bounds.size()-1); }

        FilterType type()const override { return FilterType::Energy_Filter; }

        std::string type_str() const override { return "energyfilter"; }

        std::vector<double> get_energy_bounds() { return energy_bounds; }

    private:
        std::vector<double>energy_bounds;
};

#endif