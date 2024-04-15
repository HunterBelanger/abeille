#include <tallies/box_position_filter.hpp>
#include <utils/constants.hpp>
#include <utils/output.hpp>

//std::vector<TracklengthDistance> BoxPositionFilter::get_indices_tracklength(Position r, const Direction& u_, double d_flight)
std::vector<TracklengthDistance> BoxPositionFilter::get_indices_tracklength(const Tracker& trkr, double d_flight){

    std::vector<TracklengthDistance> indexes_tracklength;
    TracklengthDistance trlen_d;

    //const Position final_loc = r + d_flight * u_;

    Position r = trkr.r();
    const Direction u_ = trkr.u();
    bool inside_bin = false;

    int i =0, j=0, k=0;
    std::array<int, 3> on;
    on.fill(0.0);

    // Check if particle inside the bin
    if ( (r_low.x() <= r.x() && r_high.x() >= r.x())
    && (r_low.y() <= r.y() && r_high.y() >= r.y())
    && (r_low.z() <= r.z() && r_high.z() >= r.z())){

        inside_bin = true;
    }
        
    // it is possible that particle is not inside the box, but can intersect
    if ( inside_bin == false){
        if ( find_entry_point(r, u_, d_flight) == false ) {
            return indexes_tracklength;

        //check_on_boundary(tktr, on);
        }
        initialize_indices(r, u_, i, j, k, on);
        if ( i == 0 && j == 0 && k ==0 ){
            inside_bin = true;
        } else {
        // This is a problem, in theory, we should now be inside the tally
        // region. We will therefore spew a warning here.
        warning("Could not locate tile after fast forward to mesh entry.\n");
        }
    } 
 
    // Distance remaining to tally
    double distance_remaining = d_flight;

    auto next_tile = distance_to_next_index(r, u_, on, i, j, k );

    if (next_tile.first == INF) {
    // Something went wrong.... Don't score.
        Output::instance().save_warning("Problem encountered with mesh tally in box_tally.\n");
        return indexes_tracklength;
        
    } else if (next_tile.first < 0.) {
    // Something went wrong.... Don't score.
    warning("Negative distance encountered with mesh tally.");
    }

    double d_tile = std::min(next_tile.first, distance_remaining);

    // Make the score if we are in a valid cell
    if ( i == 0 && j == 0 && k ==0  ) {
        trlen_d.indexes_ = std::vector<size_t>{0};
        trlen_d.distance_in_bin = d_tile;
        indexes_tracklength.push_back(trlen_d);
    }
    return indexes_tracklength;
}

void BoxPositionFilter::check_on_boundary(const Tracker tktr, std::array<int, 3> on){
        const Position r = tktr.r();
        const Direction u = tktr.u();
        if (std::abs(r_low.x() - r.x()) < SURFACE_COINCIDENT) {
            if (u.x() < 0.) {
                on[0] = 1;
            } else {
                on[0] = -1;
            }
        } else if (std::abs(r_high.x()- r.x()) < SURFACE_COINCIDENT) {
            if (u.x() < 0.) {
                on[0] = 1;
            } else {
                on[0] = -1;
            }
        }

        if (std::abs(r_low.y() - r.y()) < SURFACE_COINCIDENT) {
            if (u.y() < 0.) {
                on[1] = 1;
            } else {
                on[1] = -1;
            }
        } else if (std::abs(r_high.y() - r.y()) < SURFACE_COINCIDENT) {
            if (u.y() < 0.) {
                on[1] = 1;
            } else {
                on[1] = -1;
            }
        }

        if (std::abs(r_low.z() - r.z()) < SURFACE_COINCIDENT) {
            if (u.z() < 0.) {
                on[2] = 1;
            } else {
               on[2] = -1;
            }
        } else if (std::abs(r_high.z() - r.z()) < SURFACE_COINCIDENT) {
            if (u.z() < 0.) {
                on[2] = 1;
            } else {
                on[2] = -1;
            }
        }
}


//template<typename T>
std::shared_ptr<PositionFilter> make_box_position_filter(const YAML::Node &node){
    if (!node["low"])
        fatal_error("For box position-filter \"low\" co-ordinates is not provided.");
    if (!node["high"])
        fatal_error("For box position-filter \"high\" co-ordinates is not provided.");

    std::vector<double> low_point = node["low"].as<std::vector<double>>();
    std::vector<double> high_point = node["high"].as<std::vector<double>>();

    Position r_low_(low_point[0], low_point[1], low_point[2]);
    Position r_high_(high_point[0], high_point[1] , high_point[2]);

    std::shared_ptr<PositionFilter> box_type_filter 
                    = std::make_shared<BoxPositionFilter>(r_low_, r_high_);

    return box_type_filter;
}