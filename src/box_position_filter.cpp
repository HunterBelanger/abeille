#include <tallies/box_position_filter.hpp>

std::vector<TracklengthDistance> BoxPositionFilter::get_indices_tracklength(const Tracker& tktr, double d_flight){

        std::vector<TracklengthDistance> indexes_tracklength;

        Position r = tktr.r();
        const Direction u_ = tktr.u();
        const Position final_loc = r + d_flight * u_;
        bool inside_bin = false;

        // Check if particle inside the bin
        if ( (r_low.x() <= r.x() && r_high.x() >= r.x())
        && (r_low.y() <= r.y() && r_high.y() >= r.y())
        && (r_low.z() <= r.z() && r_high.z() >= r.z())){

            inside_bin = true;
        }
            
        // it is possible that particle is not inside the box, but can intersect
        if ( inside_bin == false){
            bool cross_in_x = ( r.x() < r_low.x() && final_loc.x() > r_low.x() ) || ( r.x() > r_high.x() && final_loc.x() < r_high.x() );
            bool cross_in_y = ( r.y() < r_low.y() && final_loc.y() > r_low.y() ) || ( r.y() > r_high.y() && final_loc.y() < r_high.y() );
            bool cross_in_z = ( r.z() < r_low.z() && final_loc.z() > r_low.z() ) || ( r.z() > r_high.z() && final_loc.z() < r_high.z() );

            inside_bin = cross_in_x || cross_in_y || cross_in_z ;

        }

        if (inside_bin){

            

        }

        return indexes_tracklength;
    }