/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
 * Copyright 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
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
#ifndef PARTICLE_H
#define PARTICLE_H

#include <utils/direction.hpp>
#include <utils/position.hpp>

#include <pcg_random.hpp>

#include <vector>

// !!! IMPORTANT !!!
// If you change BankedParticle you must also change
// register_banked_particle_type in mpi.cpp !!!!
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
  uint64_t family_id;

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
  uint64_t family_id() const { return family_id_; }
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
  void set_family_id(uint64_t i) { family_id_ = i; }
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

  uint64_t family_id_ = 0;

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
