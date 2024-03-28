#ifndef CYLINDER_POSITION_FILTER_H
#define CYLINDER_POSITION_FILTER_H

#include "position_filter.h"
#include <utils/position.hpp>

#include <array>


class CylinderPositionFilter : public PositionFilter{
    public:
        enum class Orientation {X, Y, Z};

        CylinderPositionFilter(Position origin_, double radius_,  double dx, double dy, double dz_, 
            int Nx_ = 1, int Ny_ = 1, int Nz_ = 1, Orientation z_ = Orientation::Z)
            : origin(), radius(radius_), pitch_x(dx), pitch_y(dy), dz(dz_)
            inv_pitch_x(), inv_pitch_y(), inv_dz(),
            Nx(Nx_), Ny(Ny_), Nz(Nz_), z(z_),
            new_origin()
            {   
                // Mapping the orientation to make the length direction to
                origin = map_corrdinate(origin_);
                if (z == Orientation::X){
                    pitch_x = dz;
                    dz = dx;
                    Nz = Nx_;
                    Nx = Nz_;
                }

                if (z == Orientation::Y){
                    pitch_y = dz;
                    dz = dy;
                    Nz = Ny_;
                    Ny = Nz_;
                    std::cout<<"\n";
                }
                
                // Check for the valid inputs
                if (radius <= 0.0)      std::cout<<"Fatal Error, Incorrect value of radius is provided.\n";
                if (pitch_x == 0.0){
                    std::cout<<"Warning, x-direction pitch is zero, Only one will be considered in x-direction.\n";
                    Nx = static_cast<int> (1);
                    inv_pitch_x = 0.0;
                }
                else
                    inv_pitch_x = 1.0 / pitch_x;
                
                if (pitch_y == 0.0){
                    std::cout<<"Warning, y-direction pitch is zero, Only one will be considered in y-direction.\n";
                    Ny = static_cast<int> (1);
                    inv_pitch_y = 0.0;
                }
                else    
                    inv_pitch_y = 1.0 / pitch_y;
                
                if (dz == 0.0)  std::cout<<"Fatal Error, length of the cylinder is zero.\n";

                //const double dz = length/ static_cast<double> (Nz);
                // inv_dz = 1.0 / dz;
                inv_dz = static_cast<double> (Nz) / dz;
                std::cout<<"N "<<Nx<<", "<<Ny<<", "<<Nz<<", "<<inv_dz<<"\n";
        }

        ~CylinderPositionFilter() = default;

        FilterType type()const override { return FilterType::Cylinder_Position_Filter; }
        std::string type_str()const override { return "cylinderpostionfilter"; }


        bool get_index(onst Tracker& tktr, std::array<int, 3>& indices)override final{

            const Position r = map_corrdinate(tktr.r());
            std::cout<<"Position --> "<<r.x()<<", "<<r.y()<<", "<<r.z()<<"\n";
            int nx = 0, ny = 0;

            if (pitch_x != 0)
                nx = static_cast<int> (std::floor( (r.x() - origin.x()) * inv_pitch_x ));
                
            if (pitch_y != 0)  
                ny = static_cast<int> (std::floor( (r.y() - origin.y()) * inv_pitch_y ));

            int nz = static_cast<int> (std::floor( (r.z() - origin.z()) * inv_dz ));

            std::cout<<"Nx = "<<nx<<", "<<ny<<", "<<nz<<".."<<std::floor( (r.z() - origin.z()) * inv_dz )<<" .\n";

            double new_origin_x = origin.x() + pitch_x * static_cast<double> (nx);
            double new_origin_y = origin.y() + pitch_y * static_cast<double> (ny);

            if (  nx >= 0 && nx < Nx
            && ny >= 0 && ny < Ny
            && nz >= 0 && nz < Nz )
            {
                std::cout<<"<<< "<<sqrt( (new_origin_x - r.x()) * (new_origin_x - r.x()) 
                + (new_origin_y - r.y()) * (new_origin_y - r.y()) )<<"\n";

                if ( sqrt( (new_origin_x - r.x()) * (new_origin_x - r.x()) 
                + (new_origin_y - r.y()) * (new_origin_y - r.y()) ) <= ( radius + 1E-15) ){
                    
                    indices[0] = static_cast<int> (nx);
                    indices[1] = static_cast<int> (ny);
                    indices[2] = static_cast<int> (nz);
                    map_indexes(indices);

                    new_origin = Position(new_orign_x, new_origin_y, (origin.z() + nz*dz) );

                    return true;
                }
                else{
                std::cout<<"<<< "<<indices[0]<<", "<<indices[1]<<", "<<indices[2]<<"\n";
                indices.fill(static_cast<int> (-1));
                std::cout<<"<<< "<<indices[0]<<", "<<indices[1]<<", "<<indices[2]<<"\n";
                return false;
                }
            }
            else{
                indices.fill(static_cast<int> (-1));
                return false;
            }

            return false;
        }

        int get_Nx()const override{ return Nx; }
        int get_Ny()const override { return Ny; }
        int get_Nz()const override { return Nz; }

        double get_radius()const { return radius; }
        
        double get_dx()const { return pitch_x; }
        double get_dy()const { return pitch_y; }
        double get_dz()const { return 1.0 / inv_dz;}

        double inv_dz()const { return inv_dz;}

        Position new_origin()const { return new_origin; } // to get the new_rogin for zernike; z() will zmin
        double z_min()const override { return new_origin.z(); }



    private:
        Position origin, new_origin;

        double pitch_x, pitch_y, inv_pitch_x, inv_pitch_y, inv_dz, dz; // dz is not required to store, neither length_z, only inv_dz is sufficient
        double  radius; //,length_z;
        
        int Nx, Ny, Nz;
        Orientation z;

        // To map the indexes to either converting into class co-ordinate or into original
        void map_indexes(std::array<int, 3>& indexes)const {
            if (z == Orientation::Z)
                return;
            
            if (z == Orientation::Y){
                int nz = indexes[2];
                indexes[2] = indexes[1];
                indexes[2] = nz;
            }

            if (z == Orientation::X){
                int nz = indexes[2];
                indexes[2] = indexes[0];
                indexes[0] = nz;
            }
        }
        // To map the positions to either converting into class co-ordinate or into original
        Position map_corrdinate(const Position& point) const{
            if ( z == Orientation::Z ) { return point; }
            if ( z == Orientation::Y ) { return Position( point.x(), point.z(), point.y() ); }
            if ( z== Orientation::X ){ return Position( point.z(), point.y(), point.x() ); }
            return point;
        }

};


#endif