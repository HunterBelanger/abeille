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
#ifndef PARTICLE_H
#define PARTICLE_H

#include <pcg_random.hpp>
#include <utils/direction.hpp>
#include <utils/position.hpp>
#include <vector>

struct BankedParticle {
  // Particle info
  Position r;
  Direction u;
  double E;
  double wgt;
  double wgt2;

  // Info about how it was made
  uint64_t parent_history_id;
  uint64_t parent_daughter_id;

  // For use in performing cancellation
  bool parents_previous_was_virtual = false;
  Position parents_previous_position = Position();     // R1
  Direction parents_previous_direction = Direction();  // U1
  double parents_previous_previous_energy = 0;         // E1
  double parents_previous_energy = 0;                  // E3
  double Esmp_parent = 0.;

  bool operator<(const BankedParticle& rhs) const {
    if (parent_history_id < rhs.parent_history_id) return true;
    if (parent_history_id == rhs.parent_history_id &&
        parent_daughter_id < rhs.parent_daughter_id)
      return true;
    return false;
  }
};

class Particle {
 private:
  struct ParticleState {
    Position position;
    Direction direction;
    double energy;
    double weight;
    double weight2;
  };

 public:
  Particle(Position r, Direction u, double engy, double wgt, uint64_t id = 0);
  Particle(Position r, Direction u, double engy, double wgt, double wgt2,
           uint64_t id = 0);

  Position& r() { return state.position; }
  const Position& r() const { return state.position; }
  Direction& u() { return state.direction; }
  const Direction& u() const { return state.direction; }
  double wgt() const { return state.weight; }
  double wgt2() const { return state.weight2; }
  bool is_alive() const { return alive; }
  bool is_reflected() const { return reflected; }
  double E() const { return state.energy; }
  Position previous_r() const { return previous_position; }
  Direction previous_u() const { return previous_direction; }
  double previous_E() const { return previous_energy; }
  uint64_t history_id() const { return history_id_; }
  uint64_t secondary_id() const { return secondary_id_; }
  uint64_t daughter_counter() { return daughter_counter_++; }
  double Esmp() const { return Esmp_; }
  const Position& r_birth() const { return r_birth_; }

  std::size_t num_secondaries() const { return secondaries.size(); }

  void set_position(Position r) {
    previous_position = state.position;
    state.position = r;
  }
  void set_direction(Direction u) {
    previous_direction = state.direction;
    state.direction = u;
  }
  void set_weight(double w) { state.weight = w; }
  void set_weight2(double w) { state.weight2 = w; }
  void set_energy(double E) {
    previous_energy = state.energy;
    state.energy = E;
  }
  void set_previous_r(Position r) { previous_position = r; }
  void set_reflected(bool ref) { reflected = ref; }
  void set_history_id(uint64_t i) { history_id_ = i; }
  void set_secondary_id(uint64_t i) { secondary_id_ = i; }
  void set_Esmp(double new_Esmp) { Esmp_ = new_Esmp; }

  void move(double dist) {
    if (reflected == false) {
      previous_position = state.position;
      state.position = state.position + dist * state.direction;
    } else {
      state.position = state.position + dist * state.direction;
      reflected = false;
    }
  }

  void kill() { alive = false; }

  void add_fission_particle(const BankedParticle& fiss_particle) {
    history_fission_bank.push_back(fiss_particle);
  }

  void add_noise_particle(const BankedParticle& noise_particle) {
    history_noise_bank.push_back(noise_particle);
  }

  void empty_fission_bank(std::vector<BankedParticle>& bank) {
    bank.insert(std::end(bank), std::begin(history_fission_bank),
                std::end(history_fission_bank));

    history_fission_bank.clear();
    history_fission_bank.shrink_to_fit();
  }

  void empty_noise_bank(std::vector<BankedParticle>& bank) {
    bank.insert(std::end(bank), std::begin(history_noise_bank),
                std::end(history_noise_bank));

    history_noise_bank.clear();
    history_noise_bank.shrink_to_fit();
  }

  void make_secondary(Direction u, double E, double wgt, double wgt2 = 0.) {
    secondaries.push_back({state.position, u, E, wgt, wgt2});
  }

  void split(int n_new) {
    if (n_new > 1) {
      this->set_weight(this->wgt() / static_cast<double>(n_new));
      this->set_weight2(this->wgt2() / static_cast<double>(n_new));
      for (int np = 0; np < n_new - 1; np++) {
        this->make_secondary(this->u(), this->E(), this->wgt(), this->wgt2());
      }
    }
  }

  void resurect() {
    if (!alive && !secondaries.empty()) {
      alive = true;
      state = secondaries.back();
      secondaries.pop_back();
      // r_birth_ = state.position;

      // The particle is now the next particle in the history
      // so we advance the secondary id.
      secondary_id_++;
    }
  }

  void initialize_rng(uint64_t seed, uint64_t stride) {
    rng.seed(seed);
    uint64_t n_advance = stride * history_id_;
    rng.advance(n_advance);
    histories_initial_rng = rng;
  }

  void set_initial_rng(pcg32 initial_rng) {
    histories_initial_rng = initial_rng;
  }

  bool previous_collision_virtual() const {
    return this->previous_collision_virtual_;
  }

  void set_previous_collision_virtual() {
    this->previous_collision_virtual_ = true;
  }

  void set_previous_collision_real() {
    this->previous_collision_virtual_ = false;
  }

  uint64_t number_of_rng_calls() const { return rng - histories_initial_rng; }

  pcg32 rng;

 private:
  ParticleState state;

  uint64_t history_id_;

  std::vector<ParticleState> secondaries;
  std::vector<BankedParticle> history_fission_bank;
  std::vector<BankedParticle> history_noise_bank;

  uint64_t secondary_id_ = 0;
  uint64_t daughter_counter_ = 0;

  Position previous_position = Position();
  Direction previous_direction = Direction();
  double previous_energy = 0;

  double Esmp_ = 0.;

  bool alive = true;
  bool reflected = false;
  bool previous_collision_virtual_ = false;

  Position r_birth_ = Position();

  pcg32 histories_initial_rng;
};  // Particle

#endif  // MG_PARTICLE_H
