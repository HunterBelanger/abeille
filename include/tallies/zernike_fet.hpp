#ifndef ZERNIKE_FET_H
#define ZERNIKE_FET_H

#include <tallies/cylinder_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <tallies/itally.hpp>
#include <tallies/legendre_polynomial.hpp>
#include <tallies/zernike_polynomial.hpp>
#include <utils/error.hpp>

#include <yaml-cpp/yaml.h>

#include <span>

class ZernikeFET : public ITally {
 public:
  // following constructor will be called when zernike and legendre both needs
  // to evaluated
  ZernikeFET(std::shared_ptr<CylinderFilter> cylinder_filter,
             std::shared_ptr<EnergyFilter> energy_filter,
             std::size_t zernike_order, std::size_t legendre_order,
             Quantity quantity, Estimator estimator, std::string name);

  // following constructor will be called when only zernike fet needs to be
  // evaluated
  ZernikeFET(std::shared_ptr<CylinderFilter> cylinder_filter,
             std::shared_ptr<EnergyFilter> energy_filter,
             std::size_t zernike_order, Quantity quantity, Estimator estimator,
             std::string name);

  void score_collision(const Particle& p, const Tracker& trkr,
                       MaterialHelper& mat) override final;

  void score_flight(const Particle& /*p*/, const Tracker& /*trkr*/,
                    double /*d_flight*/,
                    MaterialHelper& /*mat*/) override final {
    fatal_error("the track-length for the zernike-fet is not supoorted yet.");
  }

  void score_source(const BankedParticle& p) override final;

  double evaluate(const Position& r, const double& E) const override final;
  std::vector<double> evaluate(const std::vector<std::pair<Position, double>> r_E) const override final;

  void write_tally() override final;

  std::size_t get_zernike_order() { return zr_order_; }
  std::size_t get_legendre_order() { return legen_order_; }

 private:
  std::shared_ptr<CylinderFilter> cylinder_filter_;
  std::shared_ptr<EnergyFilter> energy_filter_;

  // Zernike and Legendre Polynomials can hold the polynomials upto that order
  // can return the std::vector<double> calculated for each order
  ZernikePolynomials zr_polynomial_;
  LegendrePolynomials legendre_polynomial_;
  std::size_t zr_order_, legen_order_;

  CylinderFilter::Orientation axial_direction_;

  bool check_for_legendre = true;
};

std::shared_ptr<ZernikeFET> make_zernike_fet(const YAML::Node& node);

#endif