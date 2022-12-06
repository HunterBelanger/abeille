#ifndef CE_NUCLIDE_H
#define CE_NUCLIDE_H

#include <PapillonNDL/st_neutron.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>
#include <materials/nuclide.hpp>
#include <memory>

class CENuclide : public Nuclide {
 public:
  CENuclide(const std::shared_ptr<pndl::STNeutron>& ce,
            const std::shared_ptr<pndl::STThermalScatteringLaw>& tsl);

  bool fissile() const override final;
  double total_xs(double E_in, std::size_t i) const override final;
  double disappearance_xs(double E_in, std::size_t i) const override final;
  double fission_xs(double E_in, std::size_t i) const override final;
  double nu_total(double E_in, std::size_t i) const override final;
  double nu_prompt(double E_in, std::size_t i) const override final;
  double nu_delayed(double E_in, std::size_t i) const override final;
  double reaction_xs(uint32_t mt, double E_in, size_t i) const override final;
  double elastic_xs(double E_in, std::size_t i) const override final;
  std::size_t energy_grid_index(double E) const override final;
  std::size_t num_delayed_groups() const override final;
  double delayed_group_constant(std::size_t g) const override final;
  double delayed_group_probability(std::size_t g,
                                   double E) const override final;
  ScatterInfo sample_scatter(double Ein, const Direction& u, std::size_t i,
                             pcg32& rng) const override final;
  ScatterInfo sample_scatter_mt(uint32_t mt, double Ein, const Direction& u,
                                std::size_t i, pcg32& rng) const override final;
  FissionInfo sample_fission(double Ein, const Direction& u, std::size_t i,
                             double Pdelayed, pcg32& rng) const override final;
  FissionInfo sample_prompt_fission(double Ein, const Direction& u,
                                    std::size_t i,
                                    pcg32& rng) const override final;
  FissionInfo sample_delayed_fission(double Ein, const Direction& u,
                                     std::size_t g,
                                     pcg32& rng) const override final;

  double max_energy() const override final;
  double min_energy() const override final;
  double speed(double E, std::size_t i) const override final;

  // CENuclide specific options
  const std::shared_ptr<pndl::STNeutron>& cedata() const { return cedata_; }
  const std::shared_ptr<pndl::STThermalScatteringLaw>& tsl() const {
    return tsl_;
  }

 private:
  std::shared_ptr<pndl::STNeutron> cedata_;
  std::shared_ptr<pndl::STThermalScatteringLaw> tsl_;

  void elastic_scatter(double Ein, const Direction& uin, double& Eout,
                       Direction& uout, pcg32& rng) const;

  void thermal_scatter(double Ein, const Direction& uin, double& Eout,
                       Direction& uout, pcg32& rng) const;
};

#endif
