#include <tallies/mesh_position_filter.hpp>

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
    }

    return indexes;
}        
 