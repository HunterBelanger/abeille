/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
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
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/nd_directory.hpp>
#include <utils/output.hpp>

#include <PapillonNDL/elastic_dbrc.hpp>

#include <algorithm>
#include <memory>
#include <sstream>

const std::shared_ptr<pndl::STNeutron>& NDDirectory::NeutronACEList::get_temp(
    const std::string& key, double T) {
  std::size_t closest_tmp_indx = neutron_ace_files.size();
  double closest_diff = 1000000.;

  for (std::size_t i = 0; i < neutron_ace_files.size(); i++) {
    if (std::abs(T - neutron_ace_files[i].temperature()) < closest_diff) {
      closest_tmp_indx = i;
      closest_diff = std::abs(T - neutron_ace_files[i].temperature());
    }
  }

  if (closest_diff > 0.1) {
    // We couldn't find the desired temeprature... This is bad. Shouldn't
    // happen.
    std::stringstream mssg;
    mssg << "Could not find required temperature " << T << " for " << key
         << '.';
    fatal_error(mssg.str());
  }

  if (!neutron_ace_files[closest_tmp_indx].loaded()) {
    std::stringstream mssg;
    mssg << " Reading Free-Gas Neutron data for " << key << " at "
         << neutron_ace_files[closest_tmp_indx].temperature() << " K.\n";
    Output::instance().write(mssg.str());

    std::unique_ptr<pndl::ACE> ace{nullptr};

    // Only the master reads the file from disk.
    if (mpi::rank == 0) {
      const std::string ace_fname =
          neutron_ace_files[closest_tmp_indx].ace_entry.fname.string();
      const pndl::ACE::Type ace_type =
          neutron_ace_files[closest_tmp_indx].ace_entry.type;
      ace = std::make_unique<pndl::ACE>(ace_fname, ace_type);
    }

    if (mpi::size > 1) {
      std::vector<std::stringstream::char_type> ace_data;
      // Master loads a vector with the binary ACE data
      if (mpi::rank == 0) {
        std::stringstream ace_strm;
        ace->save_binary(ace_strm);
        auto ace_strm_view = ace_strm.view();
        ace_data = std::vector<std::stringstream::char_type>(
            ace_strm_view.begin(), ace_strm_view.end());
      }

      // Binary ACE data sent to all other ranks
      mpi::Bcast(ace_data, 0);

      // If we are a rank just recieving the ACE data, we can now construct our
      // ACE pointer instance.
      if (mpi::rank != 0) {
        std::stringstream ace_strm;
        ace_strm.write(ace_data.data(), ace_data.size());
        ace_strm.seekg(0);
        ace = std::make_unique<pndl::ACE>(ace_strm, pndl::ACE::Type::BINARY);
      }
    }

    if (first_loaded) {
      neutron_ace_files[closest_tmp_indx].neutron_data =
          std::make_unique<pndl::STNeutron>(*ace, *first_loaded);
    } else {
      neutron_ace_files[closest_tmp_indx].neutron_data =
          std::make_unique<pndl::STNeutron>(*ace);
      first_loaded = neutron_ace_files[closest_tmp_indx].neutron_data;
    }

    // Turn off Target-At-Rest approximation for H1
    if (neutron_ace_files[closest_tmp_indx].neutron_data->awr() < 1.) {
      neutron_ace_files[closest_tmp_indx].neutron_data->elastic().set_use_tar(
          false);
    }
  }

  return neutron_ace_files[closest_tmp_indx].neutron_data;
}

const std::shared_ptr<pndl::STThermalScatteringLaw>&
NDDirectory::TSLACEList::get_temp(const std::string& key, double T) {
  std::size_t closest_tmp_indx = tsl_ace_files.size();
  double closest_diff = 1000000.;

  for (std::size_t i = 0; i < tsl_ace_files.size(); i++) {
    if (std::abs(T - tsl_ace_files[i].temperature()) < closest_diff) {
      closest_tmp_indx = i;
      closest_diff = std::abs(T - tsl_ace_files[i].temperature());
    }
  }

  if (closest_diff > 0.1) {
    // We couldn't find the desired temeprature... This is bad. Shouldn't
    // happen.
    std::stringstream mssg;
    mssg << "Could not find required temperature " << T << " for " << key
         << '.';
    fatal_error(mssg.str());
  }

  if (!tsl_ace_files[closest_tmp_indx].loaded()) {
    std::stringstream mssg;
    mssg << " Reading Thermal Scattering Law data for " << key << " at "
         << tsl_ace_files[closest_tmp_indx].temperature() << " K.\n";
    Output::instance().write(mssg.str());

    std::unique_ptr<pndl::ACE> ace{nullptr};

    // Only the master reads the file from disk.
    if (mpi::rank == 0) {
      const std::string ace_fname =
          tsl_ace_files[closest_tmp_indx].ace_entry.fname.string();
      const pndl::ACE::Type ace_type =
          tsl_ace_files[closest_tmp_indx].ace_entry.type;
      ace = std::make_unique<pndl::ACE>(ace_fname, ace_type);
    }

    if (mpi::size > 1) {
      std::vector<std::stringstream::char_type> ace_data;
      // Master loads a vector with the binary ACE data
      if (mpi::rank == 0) {
        std::stringstream ace_strm;
        ace->save_binary(ace_strm);
        auto ace_strm_view = ace_strm.view();
        ace_data = std::vector<std::stringstream::char_type>(
            ace_strm_view.begin(), ace_strm_view.end());
      }

      // Binary ACE data sent to all other ranks
      mpi::Bcast(ace_data, 0);

      // If we are a rank just recieving the ACE data, we can now construct our
      // ACE pointer instance.
      if (mpi::rank != 0) {
        std::stringstream ace_strm;
        ace_strm.write(ace_data.data(), ace_data.size());
        ace_strm.seekg(0);
        ace = std::make_unique<pndl::ACE>(ace_strm, pndl::ACE::Type::BINARY);
      }
    }

    tsl_ace_files[closest_tmp_indx].tsl_data =
        std::make_unique<pndl::STThermalScatteringLaw>(*ace);
  }

  return tsl_ace_files[closest_tmp_indx].tsl_data;
}

bool NDDirectory::has_nuclide_entry(const std::string& key) const {
  if (nuclides.find(key) == nuclides.end()) return false;
  return true;
}

bool NDDirectory::has_neutron_list(const std::string& key) const {
  if (neutron_dir.find(key) == neutron_dir.end()) return false;
  return true;
}

bool NDDirectory::has_tsl_list(const std::string& key) const {
  if (tsl_dir.find(key) == tsl_dir.end()) return false;
  return true;
}

NDDirectory::NuclideEntry& NDDirectory::get_nuclide_entry(
    const std::string& key) {
  if (!has_nuclide_entry(key)) {
    std::stringstream mssg;
    mssg << "No Nuclide entry for " << key << '\n';
    fatal_error(mssg.str());
  }

  return nuclides.at(key);
}

NDDirectory::NeutronACEList& NDDirectory::get_neutron_list(
    const std::string& key) {
  if (!has_neutron_list(key)) {
    std::stringstream mssg;
    mssg << "No NeutronACEList for " << key << '\n';
    fatal_error(mssg.str());
  }

  return neutron_dir.at(key);
}

NDDirectory::TSLACEList& NDDirectory::get_tsl_list(const std::string& key) {
  if (!has_tsl_list(key)) {
    std::stringstream mssg;
    mssg << "No TSLACEList for " << key << '\n';
    fatal_error(mssg.str());
  }

  return tsl_dir.at(key);
}

NDDirectory::NuclideEntry::NuclideEntry(const YAML::Node& node)
    : neutron(), tsl(std::nullopt), temps(), loaded(), awr(0.), dbrc(false) {
  // Make sure node is a map
  if (node.IsMap() == false) {
    fatal_error("Nuclide entry of nuclear data directory must be a map.");
  }

  // Get the mandatory neutron entry
  if (!node["neutron"] || !node["neutron"].IsScalar()) {
    fatal_error("A Nuclide entry is missing a valid \"neutron\" entry.");
  }
  neutron = node["neutron"].as<std::string>();

  // Now we get the AWR
  if (!node["awr"] || !node["awr"].IsScalar()) {
    fatal_error("No awr entry present in Nuclide entry " + neutron + ".");
  }
  awr = node["awr"].as<double>();
  if (awr < 0.) {
    fatal_error("The awr in a Nuclide entry must be > 0.");
  }

  // Now we get the TSL if present
  if (node["tsl"] && node["tsl"].IsScalar()) {
    tsl = node["tsl"].as<std::string>();
  } else if (node["tsl"]) {
    fatal_error("The tsl entry of a Nuclide entry must be a string.");
  }

  // Now we get the temperature list
  if (!node["temperatures"] || !node["temperatures"].IsSequence()) {
    fatal_error(
        "The temperatures entry of a Nuclide entry must be a sequence of "
        "positive floats.");
  }
  temps = node["temperatures"].as<std::vector<double>>();

  // Sort the temps
  std::sort(temps.begin(), temps.end());

  // Make sure all temps are positive
  if (temps.front() < 0.) {
    fatal_error("All temperatures in a Nuclide entry must be >= 0.");
  }

  // Fill loaded vector will nullptr
  loaded.resize(temps.size(), nullptr);
}

int NDDirectory::NuclideEntry::closest_temp(double T) const {
  int closest_tmp_indx = -1;
  double closest_temp = -10.;
  double closest_diff = 1000000.;

  for (std::size_t i = 0; i < temps.size(); i++) {
    if (std::abs(T - temps[i]) < closest_diff) {
      closest_tmp_indx = static_cast<int>(i);
      closest_temp = temps[i];
      closest_diff = std::abs(T - closest_temp);
    }
  }

  return closest_tmp_indx;
}

std::pair<int, int> NDDirectory::NuclideEntry::bounding_temps(double T) const {
  if (T < temps.front()) return {-1, 0};

  if (T > temps.back()) return {static_cast<int>(temps.size()) - 1, -1};

  for (int i = 0; i < static_cast<int>(temps.size()) - 1; i++) {
    if (temps[static_cast<std::size_t>(i)] <= T &&
        T <= temps[static_cast<std::size_t>(i + 1)]) {
      return {i, i + 1};
    }
  }

  // SHOULD NEVER GET HERE
  return {-1, -1};
}

NDDirectory::ACEEntry::ACEEntry(const std::filesystem::path& basename,
                                const YAML::Node& node)
    : fname(), type(pndl::ACE::Type::ASCII), temperature() {
  // Make sure node is a map
  if (node.IsMap() == false) {
    fatal_error("ACE entry of nuclear data directory must be a map.");
  }

  // Get the file path
  if (!node["file"] || !node["file"].IsScalar()) {
    fatal_error(
        "ACE entry of nuclear data directory must have a file entry with a "
        "single file name.");
  }
  fname = basename / node["file"].as<std::string>();

  // Get the ACE type (binary or ascii)
  if (node["binary"] && node["binary"].IsScalar()) {
    if (node["binary"].as<bool>()) type = pndl::ACE::Type::BINARY;
  } else if (node["binary"]) {
    fatal_error(
        "ACE entry of nuclear data directory can only have a \"binary\" entry "
        "that is boolean.");
  }

  // Get temperature
  if (!node["temperature"] || !node["temperature"].IsScalar()) {
    fatal_error(
        "ACE entry of nuclear data directory must have a \"temperature\" entry "
        "that is positive float.");
  }
  temperature = node["temperature"].as<double>();

  if (temperature < 0.) {
    fatal_error(
        "ACE entry of nuclear data directory must have a \"temperature\" entry "
        "that is positive float.");
  }
}

NDDirectory::NeutronACE::NeutronACE(const std::filesystem::path& basename,
                                    const YAML::Node& node)
    : ace_entry(basename, node), neutron_data(nullptr) {}

NDDirectory::TSLACE::TSLACE(const std::filesystem::path& basename,
                            const YAML::Node& node)
    : ace_entry(basename, node), tsl_data(nullptr) {}

NDDirectory::NeutronACEList::NeutronACEList(
    const std::filesystem::path& basename, const YAML::Node& node)
    : neutron_ace_files(), first_loaded(nullptr) {
  // Make sure the node is a sequence
  if (!node.IsSequence()) {
    fatal_error("Neutron entry list must be a sequence of ACEEntry items.");
  }

  if (node.size() == 0) {
    fatal_error(
        "Neutron entry list must be a non-empty sequence of ACEEntry items.");
  }

  // Iterate through all elements in the node
  for (std::size_t i = 0; i < node.size(); i++) {
    neutron_ace_files.emplace_back(basename, node[i]);
  }
}

NDDirectory::TSLACEList::TSLACEList(const std::filesystem::path& basename,
                                    const YAML::Node& node)
    : tsl_ace_files() {
  // Make sure the node is a sequence
  if (!node.IsSequence()) {
    fatal_error("TSL entry list must be a sequence of ACEEntry items.");
  }

  if (node.size() == 0) {
    fatal_error(
        "TSL entry list must be a non-empty sequence of ACEEntry items.");
  }

  // Iterate through all elements in the node
  for (std::size_t i = 0; i < node.size(); i++) {
    tsl_ace_files.emplace_back(basename, node[i]);
  }
}

NDDirectory::NDDirectory(const std::filesystem::path& fname,
                         TemperatureInterpolation interp, bool dbrc)
    : basename_(),
      neutron_dir(),
      tsl_dir(),
      nuclides(),
      name_(),
      date_(),
      code_(),
      notes_(),
      interp_(interp),
      use_dbrc_(dbrc) {
  // First, open the YAML file
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file " << fname << " does not exist.";
    fatal_error(mssg.str());
  }
  YAML::Node xsdir = YAML::LoadFile(fname.string());

  // Read the basename
  if (xsdir["basename"] && xsdir["basename"].IsScalar()) {
    basename_ = xsdir["basename"].as<std::string>();
  } else if (xsdir["basename"]) {
    fatal_error(
        "The basename entry of the nuclear data directory must be a valid "
        "system path.");
  }

  // Make sure the basename exists (if given)
  if (basename_.empty() == false) {
    if (std::filesystem::exists(basename_) == false) {
      std::stringstream mssg;
      mssg << "The basename \"" << basename_
           << "\" for the nuclear data directory does not exist.";
      fatal_error(mssg.str());
    }
  }

  // We now have the basename. We now iterate through the neutron-dir node
  if (!xsdir["neutron-dir"] || !xsdir["neutron-dir"].IsMap()) {
    fatal_error(
        "No valid \"neutron-dir\" entry in the nuclear data directory.");
  }
  for (YAML::const_iterator it = xsdir["neutron-dir"].begin();
       it != xsdir["neutron-dir"].end(); it++) {
    std::string key = it->first.as<std::string>();

    if (has_neutron_list(key)) {
      std::stringstream mssg;
      mssg << "The key \"" << key
           << "\" appears multiple times in neutron-dir.";
      fatal_error(mssg.str());
    }

    if (it->second.IsSequence() == false) {
      std::stringstream mssg;
      mssg << "The key \"" << key << "\" of neutron-dir is not a sequence.";
      fatal_error(mssg.str());
    }

    neutron_dir[key] = NeutronACEList(basename_, it->second);
  }

  // We now iterate through the tsl-dir node, which is actually optional
  if (xsdir["tsl-dir"] && !xsdir["tsl-dir"].IsMap()) {
    fatal_error("Invalid \"tsl-dir\" entry in the nuclear data directory.");
  } else if (xsdir["tsl-dir"]) {
    for (YAML::const_iterator it = xsdir["tsl-dir"].begin();
         it != xsdir["tsl-dir"].end(); it++) {
      std::string key = it->first.as<std::string>();

      if (has_tsl_list(key)) {
        std::stringstream mssg;
        mssg << "The key \"" << key << "\" appears multiple times in tsl-dir.";
        fatal_error(mssg.str());
      }

      if (it->second.IsSequence() == false) {
        std::stringstream mssg;
        mssg << "The key \"" << key << "\" of tsl-dir is not a sequence.";
        fatal_error(mssg.str());
      }

      tsl_dir[key] = TSLACEList(basename_, it->second);
    }
  }

  // Finally, we get all the nuclide entries
  if (!xsdir["nuclides"] || !xsdir["nuclides"].IsMap()) {
    fatal_error("Invalid \"nuclides\" entry in the nuclear data directory.");
  }
  for (YAML::const_iterator it = xsdir["nuclides"].begin();
       it != xsdir["nuclides"].end(); it++) {
    std::string key = it->first.as<std::string>();

    if (has_nuclide_entry(key)) {
      std::stringstream mssg;
      mssg << "The key \"" << key << "\" appears multiple times in nuclides.";
      fatal_error(mssg.str());
    }

    if (it->second.IsMap() == false) {
      std::stringstream mssg;
      mssg << "The key \"" << key << "\" of nuclides is not a map.";
      fatal_error(mssg.str());
    }

    nuclides[key] = NuclideEntry(it->second);
  }

  // Get the library info
  if (xsdir["library-info"] && xsdir["library-info"].IsMap()) {
    std::string key = "name";
    if (xsdir["library-info"][key] && xsdir["library-info"][key].IsScalar()) {
      name_ = xsdir["library-info"][key].as<std::string>();
    } else if (xsdir["library-info"][key]) {
      std::stringstream mssg;
      mssg << "No valid \"" << key << "\" in library-info.";
      fatal_error(mssg.str());
    }

    key = "date";
    if (xsdir["library-info"][key] && xsdir["library-info"][key].IsScalar()) {
      date_ = xsdir["library-info"][key].as<std::string>();
    } else if (xsdir["library-info"][key]) {
      std::stringstream mssg;
      mssg << "No valid \"" << key << "\" in library-info.";
      fatal_error(mssg.str());
    }

    key = "code";
    if (xsdir["library-info"][key] && xsdir["library-info"][key].IsScalar()) {
      code_ = xsdir["library-info"][key].as<std::string>();
    } else if (xsdir["library-info"][key]) {
      std::stringstream mssg;
      mssg << "No valid \"" << key << "\" in library-info.";
      fatal_error(mssg.str());
    }

    key = "notes";
    if (xsdir["library-info"][key] && xsdir["library-info"][key].IsScalar()) {
      notes_ = xsdir["library-info"][key].as<std::string>();
    } else if (xsdir["library-info"][key]) {
      std::stringstream mssg;
      mssg << "No valid \"" << key << "\" in library-info.";
      fatal_error(mssg.str());
    }
  } else if (xsdir["library-info"]) {
    fatal_error("No valid \"library-info\" entry in nuclear data directory.");
  }
}

std::shared_ptr<CENuclide> NDDirectory::get_cenuclide(const std::string& key,
                                                      NuclideEntry& nuclide,
                                                      double T) {
  std::size_t closest_T_indx =
      static_cast<std::size_t>(nuclide.closest_temp(T));
  double closest_T = nuclide.temps[closest_T_indx];

  if (std::abs(T - closest_T) > 0.1) {
    std::stringstream mssg;
    mssg << "The nuclide \"" << key << "\" has no temperature " << T << ".";
    fatal_error(mssg.str());
  }

  if (!nuclide.loaded[closest_T_indx]) {
    std::shared_ptr<pndl::STNeutron> cedata = nullptr;
    std::shared_ptr<pndl::STThermalScatteringLaw> tsl = nullptr;

    auto& neutron_list = get_neutron_list(nuclide.neutron);
    cedata = neutron_list.get_temp(nuclide.neutron, closest_T);

    if (nuclide.dbrc && this->use_dbrc_) {
      auto& cedata_0k = neutron_list.get_temp(nuclide.neutron, 0.);
      cedata->elastic().set_elastic_doppler_broadener(
          std::make_shared<pndl::ElasticDBRC>(cedata_0k->elastic_xs()));
    }

    if (nuclide.tsl) {
      auto& tsl_list = get_tsl_list(nuclide.tsl.value());
      tsl = tsl_list.get_temp(nuclide.tsl.value(), closest_T);
    }

    nuclide.loaded[closest_T_indx] = std::make_shared<CENuclide>(cedata, tsl);
  }

  return nuclide.loaded[closest_T_indx];
}

NDDirectory::CENuclidePacket NDDirectory::load_nuclide(const std::string& key,
                                                       double T) {
  // First check if we have the nuclide
  if (!has_nuclide_entry(key)) {
    std::stringstream mssg;
    mssg << "No entry in nuclear data directory for " << key << ".";
    fatal_error(mssg.str());
  }

  NuclideEntry& nuclide = get_nuclide_entry(key);

  CENuclidePacket out;
  out.nuclide_1 = std::nullopt;
  out.nuclide_2 = std::nullopt;

  double fraction1 = 0.;
  double fraction2 = 0.;

  std::shared_ptr<CENuclide> nuclide1 = nullptr;
  std::shared_ptr<CENuclide> nuclide2 = nullptr;

  // Solution strategy depends on the interpolation method
  if (interp_ == TemperatureInterpolation::Exact ||
      interp_ == TemperatureInterpolation::Nearest) {
    std::size_t closest_T_indx =
        static_cast<std::size_t>(nuclide.closest_temp(T));
    double closest_T = nuclide.temps[closest_T_indx];

    if (interp_ == TemperatureInterpolation::Exact &&
        std::abs(T - closest_T) > 0.1) {
      std::stringstream mssg;
      mssg << "The nuclide \"" << key << "\" has no temperature " << T << ".";
      fatal_error(mssg.str());
    }

    fraction1 = 1.;
    nuclide1 = get_cenuclide(key, nuclide, closest_T);

    out.nuclide_1 = CENuclideFraction(fraction1, nuclide1);
    out.nuclide_2 = std::nullopt;
  } else {
    int Tli, Thi;
    std::tie(Tli, Thi) = nuclide.bounding_temps(T);

    if (Tli < 0 && Thi < 0) {
      // This should never happen
      std::stringstream mssg;
      mssg << "Could not find any temperatures for nuclide " << key << ", at "
           << T << " Kelvin.";
      fatal_error(mssg.str());
    } else if (Tli < 0) {
      // Check if the upper temp is within 0.1 K of desired temp
      if (std::abs(nuclide.temps[static_cast<std::size_t>(Thi)] - T) < 0.1) {
        fraction1 = 1.;
        nuclide1 = get_cenuclide(key, nuclide,
                                 nuclide.temps[static_cast<std::size_t>(Thi)]);
        out.nuclide_1 = CENuclideFraction(fraction1, nuclide1);
      } else {
        std::stringstream mssg;
        mssg << "Could not find bouding temperatures for " << key << ", at "
             << T << " Kelvin.";
        fatal_error(mssg.str());
      }
    } else if (Thi < 0) {
      // Check if the lower temp is within 0.1 K of desired temp
      if (std::abs(nuclide.temps[static_cast<std::size_t>(Tli)] - T) < 0.1) {
        fraction1 = 1.;
        nuclide1 = get_cenuclide(key, nuclide,
                                 nuclide.temps[static_cast<std::size_t>(Tli)]);
        out.nuclide_1 = CENuclideFraction(fraction1, nuclide1);
      } else {
        std::stringstream mssg;
        mssg << "Could not find bouding temperatures for " << key << ", at "
             << T << " Kelvin.";
        fatal_error(mssg.str());
      }
    } else {
      const double Tl = nuclide.temps[static_cast<std::size_t>(Tli)];
      const double Th = nuclide.temps[static_cast<std::size_t>(Thi)];

      fraction1 = (Th - T) / (Th - Tl);
      fraction2 = 1. - fraction1;

      if (fraction1 < 0.) {
        std::stringstream mssg;
        mssg << "fraction1 = " << fraction1 << ", for " << key << " at " << T
             << " Kelvin.";
        fatal_error(mssg.str());
      }

      if (fraction2 < 0.) {
        std::stringstream mssg;
        mssg << "fraction2 = " << fraction2 << ", for " << key << " at " << T
             << " Kelvin.";
        fatal_error(mssg.str());
      }

      if (std::abs(Tl - T) < 0.1 || std::abs(Th - T) < 0.1) {
        // Only return 1 of the temps
        if (std::abs(Tl - T) < std::abs(Th - T)) {
          // Return low temp
          nuclide1 = get_cenuclide(key, nuclide, Tl);
          out.nuclide_1 = CENuclideFraction(1., nuclide1);
        } else {
          // Return high temp
          nuclide1 = get_cenuclide(key, nuclide, Th);
          out.nuclide_1 = CENuclideFraction(1., nuclide1);
        }
        return out;
      }

      nuclide1 = get_cenuclide(key, nuclide, Tl);
      nuclide2 = get_cenuclide(key, nuclide, Th);

      out.nuclide_1 = CENuclideFraction(fraction1, nuclide1);
      out.nuclide_2 = CENuclideFraction(fraction2, nuclide2);
    }
  }

  return out;
}

void NDDirectory::set_dbrc_nuclides(const std::vector<std::string>& dbrc_nucs) {
  // First, go through all entries in the neutron_dir and turn off dbrc
  for (auto& n : neutron_dir) {
    n.second.dbrc = false;
  }

  // Now set all those in list to true
  for (const auto& k : dbrc_nucs) {
    auto it = neutron_dir.find(k);

    if (it == neutron_dir.end()) {
      // We couldn't find this nuclide. Warn users.
      warning("Could not find neutron entry for DBRC nuclide key " + k + ".");
    } else {
      it->second.dbrc = true;
    }
  }

  // Now that these are all set, we go through all nuclides to set dbrc false
  for (auto& n : nuclides) {
    n.second.dbrc = neutron_dir.at(n.second.neutron).dbrc;
  }
}
