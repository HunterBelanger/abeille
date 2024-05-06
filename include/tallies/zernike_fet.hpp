#ifndef ZERNIKE_FET_H
#define ZERNIKE_FET_H

#include <tallies/cylinder_position_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <tallies/itally.hpp>
#include <tallies/zernike_polynomial.hpp>

#include <yaml-cpp/yaml.h>

class ZernikeFET : public ITally {
 public:
  ZernikeFET(std::shared_ptr<CylinderPositionFilter> cylinder_filter,
             std::shared_ptr<EnergyFilter> energy_filter,
             std::size_t zernike_order, std::size_t lengendre_order,
             ZernikeFET::Quantity quantity, ZernikeFET::Estimator estimator,
             std::string name);

  void score_collision(const Particle& p, const Tracker& tktr,
                       MaterialHelper& mat) override final;

  void score_flight(const Particle& /*p*/, const Tracker& /*trkr*/,
                    double /*d_flight*/,
                    MaterialHelper& /*mat*/) override final {}

  void write_tally() override final;

  std::size_t get_zernike_order() { return zr_order_; }
  std::size_t get_legendre_order() { return legen_order_; }

 private:
  std::shared_ptr<CylinderPositionFilter> cylinder_filter_;
  std::shared_ptr<EnergyFilter> energy_filter_;

  // Zernike Polynomials can hold the polynomials upto that order
  // can return the std::vector<double> calculated for each order
  ZernikePolynomials zr_polynomial_;
  std::size_t zr_order_, legen_order_;

  CylinderPositionFilter::Orientation axial_direction_;
};

std::shared_ptr<ZernikeFET> make_zernike_fet(const YAML::Node& node);

#endif