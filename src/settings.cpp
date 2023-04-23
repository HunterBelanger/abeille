/*=============================================================================*
 * Copyright (C) 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
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
#include <memory>
#include <utils/constants.hpp>
#include <utils/settings.hpp>

namespace settings {
int nparticles = 100000;  // Number of particles per batch / generation
int ngenerations = 120;
int nignored = 20;
int nskip =
    10;  // Number of power iteration generations between each noise batch

uint32_t ngroups = 0;
int n_cancel_noise_gens = INF_INT;

Timer alpha_omega_timer;
double max_time = INF;

double min_energy = 0.;       // [MeV]
double max_energy = 100000.;  // [MeV]

double target_at_rest_threshold = 400.;

SimulationMode mode = SimulationMode::K_EIGENVALUE;
TrackingMode tracking = TrackingMode::SURFACE_TRACKING;
EnergyMode energy_mode = EnergyMode::CE;

uint64_t rng_seed = 19073486328125;
uint64_t rng_stride = 152917;
pcg32 rng;

// double wgt_cutoff = 0.8;  // Needs to be high for noise !
double wgt_cutoff = 0.25;
double wgt_survival = 1.0;
double wgt_split = 2.0;

double w_noise = -1.;
double eta = 1.;   // Used in noise transport
double keff = 1.;  // Used in noise transport

bool use_urr_ptables = true;

bool converged = false;

bool pair_distance_sqrd = false;
bool families = false;

bool regional_cancellation = false;
bool regional_cancellation_noise = false;

bool inner_generations = true;
bool normalize_noise_source = true;
bool rng_stride_warnings = false;

bool save_source = false;
bool load_source_file = false;

bool branchless_splitting = false;
bool branchless_combing = true;
bool branchless_material = true;

// Energy bounds for multi-group mode
std::vector<double> energy_bounds{};
bool chi_matrix = false;
bool use_virtual_collisions = true;
std::size_t group(double E) {
  if (energy_bounds.size() <= 1) return 0;

  for (std::size_t g = 0; g < energy_bounds.size() - 1; g++) {
    if (energy_bounds[g] <= E && E <= energy_bounds[g + 1]) return g;
  }

  // Should never get here
  return 0;
}
std::vector<double> sample_xs_ratio{};

std::string output_file_name = "output.txt";
std::string source_file_name = "source.txt";
std::string in_source_file_name = "";

std::unique_ptr<NDDirectory> nd_directory = nullptr;
std::string nd_directory_fname;
bool use_dbrc = true;
TempInterpolation temp_interpolation = TempInterpolation::Linear;

void initialize_global_rng() {
  rng.seed(rng_seed);
  rng.set_stream(2);
}
}  // namespace settings
