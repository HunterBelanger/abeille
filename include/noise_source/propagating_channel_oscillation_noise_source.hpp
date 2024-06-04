#ifndef PROPAGATING_CHANNEL_OSCILLATION_NOISE_SOURCE_H
#define PROPAGATING_CHANNEL_OSCILLATION_NOISE_SOURCE_H

#include <noise_source/oscillation_noise_source.hpp>

#include <yaml-cpp/yaml.h>

class PropagatingChannelOscillationNoiseSource : public OscillationNoiseSource {
 public:
  PropagatingChannelOscillationNoiseSource(Position low, Position hi,
                                           double radius, char axis,
                                           double eps_tot, double eps_fis,
                                           double eps_sct,
                                           double angular_frequency,
                                           double velocity, double phase);

  bool is_inside(const Position& r) const override final;
  std::complex<double> dEt(const Position& r, double E,
                           double w) const override final;

  std::complex<double> dEt_Et(const Position& r, double E,
                              double w) const override final;
  std::complex<double> dEf_Ef(const Position& r, double E,
                              double w) const override final;
  std::complex<double> dEelastic_Eelastic(const Position& r, double E,
                                          double w) const override final;
  std::complex<double> dEmt_Emt(uint32_t mt, const Position& r, double E,
                                double w) const override final;

 private:
  Position low_, hi_;
  Position r_origin_;
  double radius_;
  char axis_;
  double w0_, phase_, velocity_;
  double eps_t_;
  double eps_f_;
  double eps_s_;

  double calc_center_line_radius(const Position& r) const;
  double calc_total_phase(const Position& r) const;
};

std::shared_ptr<OscillationNoiseSource>
make_propagating_channel_oscillation_noise_source(const YAML::Node& snode);

#endif