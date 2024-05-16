/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Abeille Monte Carlo code (Abeille).
 *
 * Abeille is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Abeille is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Abeille. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#include <geometry/surfaces/cross.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>

Cross::Cross(Position origin, const std::vector<double>& dists,
             BoundaryType bound, uint32_t id, const std::string& name)
    : Surface(bound, id, name), origin_(), dists_(dists) {
  // Get the origin. We set z=0, as we are parallel to the z-axis
  origin_ = Position(origin.x(), origin.y(), 0.);

  if (dists_.size() < 2) {
    fatal_error("Must have at least 2 distances in cross surface.");
  }

  // Make sure distances are sorted
  if (std::is_sorted(dists_.begin(), dists_.end()) == false) {
    fatal_error("Cross surface distances are not sorted.");
  }

  if (dists_.front() <= 0.) {
    fatal_error("Cross surface distances must be > 0.");
  }

  if (boundary_ == BoundaryType::Reflective) {
    fatal_error("Cannot have reflective boundary condition on cross surface.");
  }
}

int Cross::sign(const Position& r, const Direction& u) const {
  const double dxo = r.x() - origin_.x();
  const double dyo = r.y() - origin_.y();
  const Position r_loc(std::abs(dxo), std::abs(dyo), r.z());
  const Direction u_loc(dxo < 0. ? -u.x() : u.x(), dyo < 0. ? -u.y() : u.y(),
                        u.z());

  if (r_loc.x() - dists_.back() > SURFACE_COINCIDENT ||
      r_loc.y() - dists_.back() > SURFACE_COINCIDENT) {
    // Completely outside of the cross
    return 1;
  } else if (std::abs(r_loc.x() - dists_.back()) < SURFACE_COINCIDENT ||
             std::abs(r_loc.y() - dists_.back()) < SURFACE_COINCIDENT) {
    if (std::abs(r_loc.x() - dists_.back()) < SURFACE_COINCIDENT &&
        std::abs(r_loc.y() - dists_.back()) < SURFACE_COINCIDENT) {
      return 1;
    } else if (std::abs(r_loc.x() - dists_.back()) < SURFACE_COINCIDENT) {
      if (u.x() > 0.)
        return 1;
      else
        return -1;
    } else {
      if (u.y() > 0.)
        return 1;
      else
        return -1;
    }
  }

  // Iterate over all possible distance combinations
  int i = 0;
  int j = static_cast<int>(dists_.size()) - 1;
  while (i < static_cast<int>(dists_.size()) && j >= 0) {
    const std::size_t i_indx = static_cast<std::size_t>(i);
    const std::size_t j_indx = static_cast<std::size_t>(j);

    if (r_loc.x() - dists_[i_indx] > SURFACE_COINCIDENT &&
        r_loc.y() - dists_[j_indx] > SURFACE_COINCIDENT) {
      return 1;
    } else if (std::abs(r_loc.x() - dists_[i_indx]) < SURFACE_COINCIDENT ||
               std::abs(r_loc.x() - dists_[i_indx]) < SURFACE_COINCIDENT) {
      if (std::abs(r_loc.x() - dists_[i_indx]) < SURFACE_COINCIDENT &&
          std::abs(r_loc.x() - dists_[i_indx]) < SURFACE_COINCIDENT) {
        if (u.x() > 0. && u.y() > 0.) {
          return 1;
        } else if (u.x() < 0. && u.y() < 0.) {
          return -1;
        } else if (u.x() < 0. && u.y() > 0.) {
          // Random choice
          return 1;
        } else {
          // Complementary random choice
          return -1;
        }
      } else if (std::abs(r_loc.x() - dists_[i_indx]) < SURFACE_COINCIDENT) {
        if (u.x() > 0.)
          return 1;
        else
          return -1;
      } else {
        if (u.y() > 0.)
          return 1;
        else
          return -1;
      }
    }

    i++;
    j--;
  }

  return -1;
}

double Cross::distance(const Position& r, const Direction& u,
                       bool on_surf) const {
  const Position r_loc(r.x() - origin_.x(), r.y() - origin_.y(), r.z());

  double d_min = INF;

  for (const auto p : dists_) {
    // Check positive x distance
    {
      const double diff_xp = p - r_loc.x();
      const double dxp = diff_xp / u.x();
      if ((on_surf && std::abs(diff_xp) < SURFACE_COINCIDENT) || u.x() == 0. ||
          dxp < 0.) {
        // Skip this surface, we cannot intersect it
      } else if (dxp < d_min) {
        d_min = dxp;
      }
    }

    // Check negative x distance
    {
      const double diff_xn = -p - r_loc.x();
      const double dxn = diff_xn / u.x();
      if ((on_surf && std::abs(diff_xn) < SURFACE_COINCIDENT) || u.x() == 0. ||
          dxn < 0.) {
        // Skip this surface, we cannot intersect it
      } else if (dxn < d_min) {
        d_min = dxn;
      }
    }

    // Check positive y distance
    {
      const double diff_yp = p - r_loc.y();
      const double dyp = diff_yp / u.y();
      if ((on_surf && std::abs(diff_yp) < SURFACE_COINCIDENT) || u.y() == 0. ||
          dyp < 0.) {
        // Skip this surface, we cannot intersect it
      } else if (dyp < d_min) {
        d_min = dyp;
      }
    }

    // Check negative y distance
    {
      const double diff_yn = -p - r_loc.y();
      const double dyn = diff_yn / u.y();
      if ((on_surf && std::abs(diff_yn) < SURFACE_COINCIDENT) || u.y() == 0. ||
          dyn < 0.) {
        // Skip this surface, we cannot intersect it
      } else if (dyn < d_min) {
        d_min = dyn;
      }
    }
  }

  return d_min;
}

Direction Cross::norm(const Position& /*r*/) const {
  // Cannot get normal for this surface
  fatal_error("Cannot compute normal vector for cross surface.");

  // NEVER GETS HERE
  return Direction(1., 0., 0.);
}

std::shared_ptr<Cross> make_cross(const YAML::Node& surface_node) {
  // Get id
  uint32_t id = 1;
  if (surface_node["id"])
    id = surface_node["id"].as<uint32_t>();
  else {
    fatal_error(
        "Surface must have an id attribute with a unique positive integer.");
  }

  // Get name
  std::string name;
  if (surface_node["name"])
    name = surface_node["name"].as<std::string>();
  else
    name = "";

  // Get list of distances
  if (!surface_node["distances"] ||
      surface_node["distances"].IsSequence() == false) {
    std::stringstream mssg;
    mssg << "Cross surface with id " << id
         << " does not have a valid distances entry.";
    fatal_error(mssg.str());
  } else if (surface_node["distances"].size() < 2) {
    std::stringstream mssg;
    mssg << "Cross surface with id " << id
         << " must have at least 2 distances.";
    fatal_error(mssg.str());
  }
  std::vector<double> dists =
      surface_node["distances"].as<std::vector<double>>();

  if (std::is_sorted(dists.begin(), dists.end()) == false) {
    std::stringstream mssg;
    mssg << "Cross surface with id " << id << " has unsorted distances.";
    fatal_error(mssg.str());
  }

  if (dists.front() <= 0.) {
    std::stringstream mssg;
    mssg << "Distances on cross surface with id " << id << " must all be > 0.";
    fatal_error(mssg.str());
  }

  // Get the origin
  Position origin(0., 0., 0.);
  if (surface_node["origin"] && surface_node["origin"].IsSequence() == false) {
  } else if (surface_node["origin"] && surface_node["origin"].size() != 2) {
  } else if (surface_node["origin"]) {
    std::vector<double> origin_vec =
        surface_node["origin"].as<std::vector<double>>();
    origin = Position(origin_vec[0], origin_vec[1], 0.);
  }

  // Get boundary type
  BoundaryType boundary = BoundaryType::Normal;
  if (surface_node["boundary"]) {
    std::string boundary_string = surface_node["boundary"].as<std::string>();
    if (boundary_string == "vacuum") {
      boundary = BoundaryType::Vacuum;
    } else if (boundary_string == "reflective") {
      std::stringstream mssg;
      mssg << "Cross surface with id " << id
           << " has a reflective boundary condition. Reflective boundaries are "
              "not allowed on cross surfaces.";
      fatal_error(mssg.str());
    } else if (boundary_string == "normal") {
      boundary = BoundaryType::Normal;
    } else {
      fatal_error("Unknown boundary type \"" + boundary_string + "\".");
    }
  } else {
    boundary = BoundaryType::Normal;
  }

  return std::make_shared<Cross>(origin, dists, boundary, id, name);
}
