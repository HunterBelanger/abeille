/*=============================================================================*
 * Copyright (C) 2021, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est un programme informatique servant à faire des comparaisons
 * entre les méthodes de transport qui sont capable de traiter les milieux
 * continus avec la méthode Monte Carlo. Il résoud l'équation de Boltzmann
 * pour les particules neutres, à une vitesse et dans une dimension.
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *============================================================================*/
#include <materials/material_helper.hpp>
#include <simulation/oscillation_noise_source.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

OscillationNoiseSource::OscillationNoiseSource(Position low, Position hi,
                                               double eps_tot, double eps_fis,
                                               double eps_sct,
                                               double angular_frequency)
    : low_(low),
      hi_(hi),
      w0_(angular_frequency),
      eps_t_(eps_tot),
      eps_f_(eps_fis),
      eps_s_(eps_sct) {
  // Check low and high
  if (low_.x() >= hi_.x() || low_.y() >= hi_.y() || low_.z() >= hi_.z()) {
    std::string mssg =
        "Low is greater than or equal to hi in OscillationNoiseSource.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make sure frequency is positive
  if (w0_ <= 0.) {
    std::string mssg =
        "Negative or zero frequency provided to OscillationNoiseSource.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make sure epsilons are positive
  if (eps_t_ <= 0.) {
    std::string mssg =
        "Negative or zero epsilon total provided to OscillationNoiseSource.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (eps_f_ <= 0.) {
    std::string mssg =
        "Negative or zero epsilon fission provided to OscillationNoiseSource.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (eps_s_ <= 0.) {
    std::string mssg =
        "Negative or zero epsilon scatter provided to OscillationNoiseSource.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

bool OscillationNoiseSource::is_inside(const Position& r,
                                       const Direction& /*u*/) const {
  if (r.x() > low_.x() && r.y() > low_.y() && r.z() > low_.z() &&
      r.x() < hi_.x() && r.y() < hi_.y() && r.z() < hi_.z())
    return true;

  return false;
}

std::complex<double> OscillationNoiseSource::dEt(const Position& r,
                                                 const Direction& u, double E,
                                                 double w) const {
  // Get the material
  Tracker trkr(r, u);
  auto material = trkr.material();

  // Make sure material is valid
  if (!material) {
    std::string mssg = "";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make material helper and evaluate XS
  MaterialHelper mat(material, E);
  double xs = mat.Et(E);

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_t_ * xs * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::complex<double> OscillationNoiseSource::dEt_Et(const Position& /*r*/,
                                                    const Direction& /*u*/,
                                                    double /*E*/,
                                                    double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_t_ * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::complex<double> OscillationNoiseSource::dEf(const Position& r,
                                                 const Direction& u, double E,
                                                 double w) const {
  // Get the material
  Tracker trkr(r, u);
  auto material = trkr.material();

  // Make sure material is valid
  if (!material) {
    std::string mssg = "";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make material helper and evaluate XS
  MaterialHelper mat(material, E);
  double xs = mat.Ef(E);

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_f_ * xs * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::complex<double> OscillationNoiseSource::dEf_Ef(const Position& /*r*/,
                                                    const Direction& /*u*/,
                                                    double /*E*/,
                                                    double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_f_ * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::complex<double> OscillationNoiseSource::dEelastic(const Position& r,
                                                       const Direction& u,
                                                       double E,
                                                       double w) const {
  // Get the material
  Tracker trkr(r, u);
  auto material = trkr.material();

  // Make sure material is valid
  if (!material) {
    std::string mssg = "";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make material helper and evaluate XS
  MaterialHelper mat(material, E);
  double xs = mat.Eelastic(E);

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_s_ * xs * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::complex<double> OscillationNoiseSource::dEelastic_Eelastic(
    const Position& /*r*/, const Direction& /*u*/, double /*E*/,
    double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_s_ * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::complex<double> OscillationNoiseSource::dEmt(uint32_t mt,
                                                  const Position& r,
                                                  const Direction& u, double E,
                                                  double w) const {
  // Get the material
  Tracker trkr(r, u);
  auto material = trkr.material();

  // Make sure material is valid
  if (!material) {
    std::string mssg = "";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make material helper and evaluate XS
  MaterialHelper mat(material, E);
  double xs = mat.Emt(mt, E);

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_s_ * xs * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::complex<double> OscillationNoiseSource::dEmt_Emt(uint32_t /*mt*/,
                                                      const Position& /*r*/,
                                                      const Direction& /*u*/,
                                                      double /*E*/,
                                                      double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if ((n == 1 || n == -1) && std::abs(err) < 0.01) {
    return {eps_s_ * PI, 0.};
  } else {
    return {0., 0.};
  }
}

std::shared_ptr<OscillationNoiseSource> make_oscillation_noise_source(
    const YAML::Node& snode) {
  // Get low
  if (!snode["low"] || !snode["low"].IsSequence() ||
      !(snode["low"].size() == 3)) {
    std::string mssg = "No valid low entry for oscillation noise source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double xl = snode["low"][0].as<double>();
  double yl = snode["low"][1].as<double>();
  double zl = snode["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  // Get hi
  if (!snode["hi"] || !snode["hi"].IsSequence() || !(snode["hi"].size() == 3)) {
    std::string mssg = "No valid hi entry for oscillation noise source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double xh = snode["hi"][0].as<double>();
  double yh = snode["hi"][1].as<double>();
  double zh = snode["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  // Get frequency
  if (!snode["angular-frequency"] || !snode["angular-frequency"].IsScalar()) {
    std::string mssg =
        "No valid angular-frequency entry for oscillation noise source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double w0 = snode["angular-frequency"].as<double>();

  if (w0 <= 0.) {
    std::string mssg =
        "Angular frequency can not be negative or zero for oscillation noise "
        "source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get epsilons
  if (!snode["epsilon-total"] || !snode["epsilon-total"].IsScalar()) {
    std::string mssg =
        "No valid epsilon-total entry for oscillation noise source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double eps_t = snode["epsilon-total"].as<double>();

  if (eps_t <= 0.) {
    std::string mssg =
        "Epsilon total can not be negative or zero for oscillation noise "
        "source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (!snode["epsilon-fission"] || !snode["epsilon-fission"].IsScalar()) {
    std::string mssg =
        "No valid epsilon-fission entry for oscillation noise source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double eps_f = snode["epsilon-fission"].as<double>();

  if (eps_f <= 0.) {
    std::string mssg =
        "Epsilon fission can not be negative or zero for oscillation noise "
        "source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (!snode["epsilon-scatter"] || !snode["epsilon-scatter"].IsScalar()) {
    std::string mssg =
        "No valid epsilon-scatter entry for oscillation noise source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double eps_s = snode["epsilon-scatter"].as<double>();

  if (eps_s <= 0.) {
    std::string mssg =
        "Epsilon scatter can not be negative or zero for oscillation noise "
        "source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  return std::make_shared<OscillationNoiseSource>(r_low, r_hi, eps_t, eps_f,
                                                  eps_s, w0);
}
