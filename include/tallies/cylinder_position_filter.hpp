#ifndef CYLINDER_POSITION_FILTER_H
#define CYLINDER_POSITION_FILTER_H

#include <array>

#include <tallies/position_filter.hpp>
#include <utils/position.hpp>
#include <utils/error.hpp>


class CylinderPositionFilter : public PositionFilter{
    public:
        enum class Orientation {X, Y, Z}; // update this

        CylinderPositionFilter(Position origin_, double radius_,  double dx, double dy, double dz_, 
            size_t nx_ = 1, size_t ny_ = 1, size_t nz_ = 1, Orientation z_ = Orientation::Z)
            : origin(origin_),  new_origin_(),
            Nx_(nx_), Ny_(ny_), Nz_(nz_), axial_dir(z_),
            radius(radius_), pitch_x(dx), pitch_y(dy), dz(dz_),
            inv_pitch_x(), inv_pitch_y(), inv_dz()
            {   
                // Mapping the orientation to make the length direction to
                origin = map_corrdinate(origin_);
                if (axial_dir == Orientation::X){
                    pitch_x = dz;
                    dz = dx;
                    Nz_ = nx_;
                    Nx_ = nz_;
                }

                if (axial_dir == Orientation::Y){
                    pitch_y = dz;
                    dz = dy;
                    Nz_ = ny_;
                    Ny_ = nz_;
                }
                
                // Check for the valid inputs
                if (radius <= 0.0)      
                    fatal_error("Incorrect value of radius is provided.");

                if ( (pitch_x == 0.0) &&  (Nx_ != 1) ){
                    warning("x-direction pitch is zero, Only one will be considered in x-direction.");
                    Nx_ = static_cast<size_t> (1);
                    inv_pitch_x = 0.0;
                }
                else{
                    if ( pitch_x < (2*radius) )
                        fatal_error("Pitch is less than the radius");
                    inv_pitch_x = 1.0 / pitch_x;}
                
                if ((pitch_y == 0.0) && (Ny_ != 1)){
                    warning("y-direction pitch is zero, Only one will be considered in y-direction.");
                    Ny_ = static_cast<size_t> (1);
                    inv_pitch_y = 0.0;
                }
                else{
                    if ( pitch_y < (2*radius) )
                        fatal_error("Pitch is less than the radius ");
                    inv_pitch_y = 1.0 / pitch_y;}
                
                if (dz == 0.0)  
                    fatal_error("Length of the cylinder shouldn't zero.");
                else
                    inv_dz = 1.0 / dz;
        }

        ~CylinderPositionFilter() = default;

        std::size_t Nx()const override{
            if ( axial_dir == Orientation::X)
                return Nz_;   

            return Nx_;
        }

        std::size_t Ny()const override{
            if ( axial_dir == Orientation::Y)
                return Nz_;  

            return Ny_; 
        }

        std::size_t Nz()const override{
            if ( axial_dir == Orientation::X)
                return Nx_;

            if ( axial_dir == Orientation::Y)
                return Ny_;   

            return Nz_;  
        }

        bool get_indices(const Tracker& tktr, std::array<int, 3>& indices);

        double get_radius()const { return radius; }
        
        double get_dx()const { return pitch_x; }
        double get_dy()const { return pitch_y; }
        double get_dz()const { return dz;}

        
        Position new_origin()const { return new_origin_; } // to get the new_origin for zernike; z() will zmin
        double axial_min(){ return new_origin_.z(); }
        double axial_max(){ return new_origin_.z() + dz; }

        //perhaps not needed
        FilterType type()const override { return FilterType::Cylinder_Position_Filter; }
        std::string type_str()const override { return "cylinderpostionfilter"; }


    private:
        Position origin, new_origin_;
        size_t Nx_, Ny_, Nz_;
        Orientation axial_dir;
        double radius, pitch_x, pitch_y, dz, inv_pitch_x, inv_pitch_y, inv_dz; // dz is not required to store, neither length_z, only inv_dz is sufficient
        

        // To map the indexes to either converting into class co-ordinate or into original
        void map_indexes(std::array<int, 3>& indexes)const {
            if (axial_dir == Orientation::Z)
                return;
            
            if (axial_dir == Orientation::Y){
                int nz = indexes[2];
                indexes[2] = indexes[1];
                indexes[1] = nz;
            }

            if (axial_dir == Orientation::X){
                int nz = indexes[2];
                indexes[2] = indexes[0];
                indexes[0] = nz;
            }
        }
        // To map the positions to either converting into class co-ordinate or into original
        Position map_corrdinate(const Position& point) const{
            if ( axial_dir == Orientation::Z ) { return point; }
            if ( axial_dir == Orientation::Y ) { return Position( point.x(), point.z(), point.y() ); }
            if ( axial_dir == Orientation::X ){ return Position( point.z(), point.y(), point.x() ); }
            return point;
        }

};

bool CylinderPositionFilter::get_indices(const Tracker& tktr, std::array<int, 3>& indices){

    const Position r = map_corrdinate(tktr.r());
    int nx = 0, ny = 0;

    if (pitch_x != 0)
        nx = static_cast<int> (std::floor( (r.x() - origin.x()) * inv_pitch_x ));
        
    if (pitch_y != 0)  
        ny = static_cast<int> (std::floor( (r.y() - origin.y()) * inv_pitch_y ));

    int nz = static_cast<int> (std::floor( (r.z() - origin.z()) * inv_dz ));

    double new_origin_x = origin.x() + pitch_x * static_cast<double> (nx);
    double new_origin_y = origin.y() + pitch_y * static_cast<double> (ny);

    //std::cout<<">>>>> "<<nx<<", "<<ny<<", "<<nz<<"\n";
    if (  nx >= 0 && nx < static_cast<int>(Nx_)
    && ny >= 0 && ny < static_cast<int>(Ny_)
    && nz >= 0 && nz < static_cast<int>(Nz_) )
    {
        if ( sqrt( (new_origin_x - r.x()) * (new_origin_x - r.x()) 
        + (new_origin_y - r.y()) * (new_origin_y - r.y()) ) <= ( radius + 1E-15) ){
            
            indices[0] = static_cast<int> (nx);
            indices[1] = static_cast<int> (ny);
            indices[2] = static_cast<int> (nz);
            map_indexes(indices);

            new_origin_ = Position(new_origin_x, new_origin_y, (origin.z() + nz * dz) );

            return true;
        }
        else{

        indices.fill(static_cast<int> (-1));
        return false;
        }
    }
    else{
        indices.fill(static_cast<int> (-1));
        return false;
    }

    return false;
}



#endif