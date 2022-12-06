#include <array>
#include <iomanip>
#include <ios>
#include <materials/ce_nuclide.hpp>
#include <sstream>
#include <stdexcept>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

CENuclide::CENuclide(const std::shared_ptr<pndl::STNeutron>& ce,
                     const std::shared_ptr<pndl::STThermalScatteringLaw>& tsl)
    : cedata_(ce), tsl_(tsl) {
  if (!cedata_) {
    std::string mssg =
        "CENuclide instance must have a valid pndl::STNeutron instance.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

bool CENuclide::fissile() const { return cedata_->fissile(); }

double CENuclide::total_xs(double E_in, std::size_t i) const {
  if (tsl_ && E_in < tsl_->max_energy()) {
    const double thermal = tsl_->xs(E_in);
    const double total = cedata_->total_xs()(E_in, i);
    const double elastic = cedata_->elastic_xs()(E_in, i);
    return total - elastic + thermal;
  }

  return cedata_->total_xs()(E_in, i);
}

double CENuclide::disappearance_xs(double E_in, std::size_t i) const {
  return cedata_->disappearance_xs()(E_in, i);
}

double CENuclide::fission_xs(double E_in, std::size_t i) const {
  return cedata_->fission_xs()(E_in, i);
}

double CENuclide::nu_total(double E_in, std::size_t /*i*/) const {
  return cedata_->fission().nu_total()(E_in);
}

double CENuclide::nu_prompt(double E_in, std::size_t /*i*/) const {
  return cedata_->fission().nu_prompt()(E_in);
}

double CENuclide::nu_delayed(double E_in, std::size_t /*i*/) const {
  return cedata_->fission().nu_delayed()(E_in);
}

double CENuclide::reaction_xs(uint32_t mt, double E_in, std::size_t i) const {
  if (!cedata_->has_reaction(mt)) return 0.;

  return cedata_->reaction(mt).xs()(E_in, i);
}

double CENuclide::elastic_xs(double E_in, std::size_t i) const {
  if (tsl_ && E_in < tsl_->max_energy()) return tsl_->xs(E_in);

  return cedata_->elastic_xs()(E_in, i);
}

std::size_t CENuclide::energy_grid_index(double E) const {
  return cedata_->energy_grid().get_lower_index(E);
}

std::size_t CENuclide::num_delayed_groups() const {
  return cedata_->fission().n_delayed_families();
}

double CENuclide::delayed_group_constant(std::size_t g) const {
  return cedata_->fission().delayed_family(g).decay_constant();
}

double CENuclide::delayed_group_probability(std::size_t g, double E) const {
  return cedata_->fission().delayed_family(g).probability()(E);
}

double CENuclide::max_energy() const {
  return cedata_->energy_grid().max_energy();
}

double CENuclide::min_energy() const {
  return cedata_->energy_grid().min_energy();
}

double CENuclide::speed(double E, std::size_t /*i*/) const {
  constexpr double inv_mass_t2 =
      2. / (N_MASS_EV / (C_CM_S * C_CM_S));  // 2/Mass in [cm ^ 2 / (eV * s^2)]
  return std::sqrt(E * inv_mass_t2);         // Speed in [cm / s]
}

// These are all possible MTs which result in an exit neutron
// (other than fission) from the ENDF manual.
static std::array<uint32_t, 106> MT_LIST{
    2,   5,   11,  16,  17,  22,  23,  24,  25,  28,  29,  30,  32,  33,
    34,  35,  36,  37,  41,  42,  44,  45,  51,  52,  53,  54,  55,  56,
    57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
    71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,
    85,  86,  87,  88,  89,  90,  91,  152, 153, 154, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173,
    174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188,
    189, 190, 194, 195, 196, 198, 199, 200};

#include <utils/output.hpp>
ScatterInfo CENuclide::sample_scatter(double Ein, const Direction& u,
                                      std::size_t i, pcg32& rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));
  std::array<double, 106> XS;
  XS.fill(0.);

  // Get XS value for all of the possible scattering reactions
  const double xi = RNG::rand(rng);
  const double scatter_xs = cedata_->total_xs()(Ein, i) -
                            cedata_->disappearance_xs()(Ein, i) -
                            cedata_->fission_xs()(Ein, i);
  double xs_sum = 0;
  uint32_t MT = 0;
  for (size_t j = 0; j < MT_LIST.size(); j++) {
    if (MT_LIST[j] != 2 && cedata_->has_reaction(MT_LIST[j]) &&
        Ein > cedata_->reaction(MT_LIST[j]).threshold()) {
      XS[j] = cedata_->reaction(MT_LIST[j]).xs()(Ein, i);
    } else if (MT_LIST[j] == 2) {
      XS[j] = cedata_->elastic_xs()(Ein, i);
    }

    xs_sum += XS[j];

    if (xi < xs_sum / scatter_xs) {
      MT = MT_LIST[j];
      break;
    }
  }

  // If for some reason we didn't sample the mt before,
  // we do it now with a discrete distribution
  if (MT == 0) {
    MT = MT_LIST[static_cast<std::size_t>(
        RNG::discrete(rng, &XS[0], &XS[XS.size() - 1]))];
  }

  double yield = 1.;
  double Eout = 0.;
  Direction uout = Direction(1., 0., 0.);

  // Check for elastic scattering
  if (MT == 2) {
    elastic_scatter(Ein, u, Eout, uout, rng);
  } else {
    // Get yield for the reaction
    yield = cedata_->reaction(MT).yield()(Ein);

    // Sample outgoing info
    pndl::AngleEnergyPacket ae_out =
        cedata_->reaction(MT).sample_neutron_angle_energy(Ein, rngfunc);

    Eout = ae_out.energy;
    uout = rotate_direction(u, ae_out.cosine_angle, 2. * PI * rngfunc());
  }

  ScatterInfo info;
  info.mt = MT;
  info.yield = yield;
  info.energy = Eout;
  info.direction = uout;
  return info;
}

ScatterInfo CENuclide::sample_scatter_mt(uint32_t mt, double Ein,
                                         const Direction& u, std::size_t /*i*/,
                                         pcg32& rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));

  if (cedata_->has_reaction(mt) == false) {
    std::stringstream mssg;
    mssg << "Nuclide " << this->id() << ", ZAID " << cedata_->zaid();
    mssg << ", has no reaction data for MT " << mt << ".";
    fatal_error(mssg.str(), __FILE__, __LINE__);
  }

  double yield = 1.;
  double E_out;
  Direction u_out = Direction(1., 0., 0.);

  // Check for elastic scattering
  if (mt == 2) {
    elastic_scatter(Ein, u, E_out, u_out, rng);
  } else {
    // Get yield for the reaction
    yield = cedata_->reaction(mt).yield()(Ein);

    // Sample outgoing info
    pndl::AngleEnergyPacket ae_out =
        cedata_->reaction(mt).sample_neutron_angle_energy(Ein, rngfunc);

    E_out = ae_out.energy;
    u_out = rotate_direction(u, ae_out.cosine_angle, 2. * PI * rngfunc());
  }

  ScatterInfo info;
  info.mt = mt;
  info.yield = yield;
  info.energy = E_out;
  info.direction = u_out;
  return info;
}

FissionInfo CENuclide::sample_fission(double Ein, const Direction& u,
                                      std::size_t i, double Pdelayed,
                                      pcg32& rng) const {
  if (RNG::rand(rng) < Pdelayed) {
    // Make delayed neutron
    // Must first sample the delayed family
    const std::size_t ngroups = cedata_->fission().n_delayed_families();
    std::vector<double> group_probs(ngroups, 0.);
    for (size_t g = 0; g < ngroups; g++) {
      group_probs[g] = cedata_->fission().delayed_family(g).probability()(Ein);
    }
    std::size_t group =
        static_cast<std::size_t>(RNG::discrete(rng, group_probs));

    return this->sample_delayed_fission(Ein, u, group, rng);
  } else {
    // Sample prompt neutron
    return this->sample_prompt_fission(Ein, u, i, rng);
  }
}

FissionInfo CENuclide::sample_prompt_fission(double Ein, const Direction& u,
                                             std::size_t /*i*/,
                                             pcg32& rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));
  pndl::AngleEnergyPacket ae_out;

  bool sampled = false;
  while (!sampled) {
    ae_out =
        cedata_->fission().prompt_spectrum().sample_angle_energy(Ein, rngfunc);

    // Check outgoing energy
    if (ae_out.energy < settings::max_energy) sampled = true;
  }

  FissionInfo info;
  info.energy = ae_out.energy;
  double phi = 2. * PI * rngfunc();
  info.direction = rotate_direction(u, ae_out.cosine_angle, phi);
  info.delayed = false;
  info.delayed_family = 0;
  info.precursor_decay_constant = 0.;
  return info;
}

FissionInfo CENuclide::sample_delayed_fission(double Ein, const Direction& u,
                                              std::size_t g, pcg32& rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));
  double Eout = 0.;
  bool sampled = false;

  while (!sampled) {
    Eout = cedata_->fission().delayed_family(g).sample_energy(Ein, rngfunc);

    if (Eout < settings::max_energy) sampled = true;
  }

  double mu = 2 * rngfunc() - 1.;
  double phi = 2. * PI * rngfunc();

  FissionInfo info;
  info.energy = Eout;
  info.direction = rotate_direction(u, mu, phi);
  info.delayed = true;
  info.delayed_family = static_cast<uint32_t>(g);
  info.precursor_decay_constant =
      cedata_->fission().delayed_family(g).decay_constant();
  return info;
}

void CENuclide::elastic_scatter(double Ein, const Direction& uin, double& Eout,
                                Direction& uout, pcg32& rng) const {
  // Check if we need to use the thermal scattering law
  if (tsl_ && Ein < tsl_->max_energy()) {
    thermal_scatter(Ein, uin, Eout, uout, rng);
    return;
  }

  // Not using thermal scattering law, use whatever the STNeutron instance has
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));
  auto ae = cedata_->elastic().sample_angle_energy(Ein, rngfunc);
  Eout = ae.energy;
  double mu = ae.cosine_angle;
  uout = rotate_direction(uin, mu, 2. * PI * rngfunc());
}

void CENuclide::thermal_scatter(double Ein, const Direction& uin, double& Eout,
                                Direction& uout, pcg32& rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));

  // Must first grab all xs values to compute probabilities.
  // This algo should be mixed-elastic ready !
  const double xi = RNG::rand(rng);
  const double ii_xs = tsl_->incoherent_inelastic().xs(Ein);
  const double ie_xs = tsl_->incoherent_elastic().xs(Ein);
  const double ce_xs = tsl_->coherent_elastic().xs(Ein);
  const double inv_tot_xs = 1. / (ii_xs + ie_xs + ce_xs);

  pndl::AngleEnergyPacket ae_out;

  if (xi < (ii_xs * inv_tot_xs) && ii_xs > 0.) {
    // Incoherent Inelastic
    ae_out = tsl_->incoherent_inelastic().sample_angle_energy(Ein, rngfunc);
  } else if (xi < ((ii_xs + ie_xs) * inv_tot_xs) && ie_xs > 0.) {
    // Incoherent Elastic
    ae_out = tsl_->incoherent_elastic().sample_angle_energy(Ein, rngfunc);
  } else if (xi < ((ii_xs + ie_xs + ce_xs) * inv_tot_xs) && ie_xs > 0.) {
    // Coherent Elastic
    ae_out = tsl_->coherent_elastic().sample_angle_energy(Ein, rngfunc);
  } else {
    // This is bad. For some reason, we couldn't sample something....
    // Let's just write a warning, and then sample II ?
    std::string mssg =
        " Could not sample a TSL reaction. Using Incoherent Inelastic.\n";
    warning(mssg, __FILE__, __LINE__);
    ae_out = tsl_->incoherent_inelastic().sample_angle_energy(Ein, rngfunc);
  }

  uout = rotate_direction(uin, ae_out.cosine_angle, 2. * PI * rngfunc());
  Eout = ae_out.energy;
}
