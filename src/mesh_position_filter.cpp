#include <tallies/position_filter.hpp>
#include <tallies/mesh_position_filter.hpp>
#include <utils/output.hpp>

std::vector<size_t> MeshPositionFilter::get_indices(const Tracker& tktr){
    std::vector<size_t> indexes;
    const Position r = tktr.r();
    int index_x = static_cast<int> (std::floor( (r.x() - r_low.x()) * dx_inv) );
    int index_y = static_cast<int> (std::floor( (r.y() - r_low.y()) * dy_inv) );
    int index_z = static_cast<int> (std::floor( (r.z() - r_low.z()) * dz_inv) );

    if ( (index_x >= 0 && index_x < static_cast<int>(Nx_) ) 
    && ( index_y >= 0 && index_y < static_cast<int>(Ny_) ) 
    && (index_z >= 0 && index_z < static_cast<int>(Nz_) )){

        indexes.push_back(index_x);
        indexes.push_back(index_y);
        indexes.push_back(index_z);
        xmin = index_x * dx;
        ymin = index_y * dy;
        zmin = index_z * dz;
        reduce_dimension(indexes);
    }

    return indexes;
}        


//std::vector<TracklengthDistance> MeshPositionFilter::get_indices_tracklength(Position r, const Direction& u_, double d_flight){
std::vector<TracklengthDistance> MeshPositionFilter::get_indices_tracklength(const Tracker& trkr, double d_flight){

    std::vector<TracklengthDistance> indexes_tracklength;
    TracklengthDistance trlen_d;

    Position r = trkr.r();
    const Direction u_ = trkr.u();
    bool inside_bin = false;

    int i =0, j=0, k=0;
    std::array<int, 3> on;
    on.fill(0.0);
    initialize_indices(r, u_, i, j, k, on);

    // Check if particle inside the bin
    if ( i >= 0 && i < static_cast<int>(Nx_)
    && j >= 0 && j < static_cast<int>(Ny_)
    && k >= 0 && k < static_cast<int>(Nz_)){

        inside_bin = true;
    }
        
    // it is possible that particle is not inside the box, but can intersect
    if ( inside_bin == false){
        if ( find_entry_point(r, u_, d_flight) == false ) {
            return indexes_tracklength;
        }
        initialize_indices(r, u_, i, j, k, on);
        if ( i >= 0 && i < static_cast<int>(Nx_)
        && j >= 0 && j < static_cast<int>(Ny_)
        && k >= 0 && k < static_cast<int>(Nz_) ){

            inside_bin = true;
        } else {
        // This is a problem, in theory, we should now be inside the tally
        // region. We will therefore spew a warning here.
        warning("Could not locate tile after fast forward to mesh entry.\n");
        }
    } 
 
    // Distance remaining to tally
    double distance_remaining = d_flight;

    while (distance_remaining > 0.) {
        // Distance we will travel in this cell
        //indexes_tracklength.clear();
        auto next_tile = distance_to_next_index(r, u_, on, i, j, k );

        if (next_tile.first == INF) {
        // Something went wrong.... Don't score.
        Output::instance().save_warning("Problem encountered with mesh tally in box_tally.\n");
        return indexes_tracklength;
       //break;
        } else if (next_tile.first < 0.) {
        // Something went wrong.... Don't score.
        warning("Negative distance encountered with mesh tally.");
        }

        double d_tile = std::min(next_tile.first, distance_remaining);

        // Make the score if we are in a valid cell
        if (i >= 0 && i < static_cast<int>(Nx_)
        && j >= 0 && j < static_cast<int>(Ny_)
        && k >= 0 && k < static_cast<int>(Nz_)) {

        size_t ui = static_cast<size_t>(i);
        size_t uj = static_cast<size_t>(j);
        size_t uk = static_cast<size_t>(k);
        
        std::vector<size_t> u_indexes{ ui, uj, uk};
        
        reduce_dimension(u_indexes);
        trlen_d.indexes_ = u_indexes;
        trlen_d.distance_in_bin = d_tile;
        indexes_tracklength.push_back(trlen_d);


        } else {
        // If we arrive here, it means that we have left the tally region
        // when were we initially inside it. We can return here, as it's
        // impossible to go back in.
            return indexes_tracklength;
        }

        // Remove the traveled distance
        distance_remaining -= d_tile;

       if (distance_remaining <= 0.) break;

        // Update the position and cell indices
        r = r + d_tile * u_;
        update_indices(next_tile.second, i, j, k, on);
        
    
    }  // While we still have to travel

    return indexes_tracklength;
}
