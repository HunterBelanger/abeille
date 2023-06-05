#include <simulation/entropy.hpp>
#include <source_location>
#include <sstream>
#include <utils/error.hpp>
#include <utils/mpi.hpp>

void Entropy::add_point(const Position& r, const double& w) {
  // Get bin idecies
  int32_t nx =
      static_cast<std::int32_t>(std::floor((r.x() - lower_corner.x()) / dx));
  int32_t ny =
      static_cast<std::int32_t>(std::floor((r.y() - lower_corner.y()) / dy));
  int32_t nz =
      static_cast<std::int32_t>(std::floor((r.z() - lower_corner.z()) / dz));

  // Make sure position is in a valid bin
  if (nx >= 0 && nx < static_cast<std::int32_t>(shape[0]) && ny >= 0 &&
      ny < static_cast<std::int32_t>(shape[1]) && nz >= 0 &&
      nz < static_cast<std::int32_t>(shape[2])) {
    // Add to total weight, no matter sign being tallied
    this->total_weight += w;

    // Get sign of particle
    Sign p_sign;
    if (sign == Sign::Total)
      p_sign = Sign::Total;
    else if (w > 0.)
      p_sign = Sign::Positive;
    else
      p_sign = Sign::Negative;

    // Add weight to bin if sign agrees
    if (p_sign == sign) {
      this->bins[static_cast<std::size_t>(
          (shape[1] * shape[2]) * static_cast<std::uint32_t>(nx) +
          (shape[2]) * static_cast<std::uint32_t>(ny) +
          static_cast<std::uint32_t>(nz))] += w;
    }
  }
}

void Entropy::synchronize_entropy_across_nodes() {
  // All worker threads must send their generation score to the master.
  // Master must recieve all generations scores from workers and add
  // them to it's own generation score. We do this with MPI_Reduce.
  mpi::Reduce_sum(bins, 0);
  mpi::Reduce_sum(total_weight, 0);
}

double Entropy::calculate_entropy() const {
  // Set sum to zero
  double sum = 0.;

  for (const auto& b : bins) {
    double p = std::fabs(b) / this->total_weight;

    if (p > 1.0) {
      std::stringstream mssg;
      mssg << " Negative entropy: p = " << p << ", bin = " << b
           << ", total_weight = " << this->total_weight << "\n";
      warning(mssg.str());
    } else if (p != 0.) {
      sum -= p * std::log2(p);
    }
  }

  return sum;
}

double Entropy::calculate_empty_fraction() const {
  double num_empty_bins = 0.;

  for (const auto& b : bins) {
    if (b == 0.) {
      num_empty_bins += 1.;
    }
  }

  return num_empty_bins / static_cast<double>(bins.size());
}

void Entropy::zero() {
  // Zero all bins
  for (auto& b : this->bins) b = 0.;

  this->total_weight = 0.;
}
