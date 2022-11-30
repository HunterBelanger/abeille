#ifndef ABEILLE_ND_DIRECTORY_H
#define ABEILLE_ND_DIRECTORY_H

#include <yaml-cpp/node/node.h>
#include <yaml-cpp/yaml.h>

#include <PapillonNDL/st_neutron.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

// Temporary struct to hold the place of future CENuclide class.
struct CENuclide {
  std::shared_ptr<pndl::STNeutron> cedata;
  std::shared_ptr<pndl::STThermalScatteringLaw> tsl;
};

class NDDirectory {
 public:
  enum class TemperatureInterpolation { Exact, Nearest, Linear };

  struct CENuclideFraction {
    double fraction;
    std::shared_ptr<CENuclide> nuclide;

    CENuclideFraction(double f, std::shared_ptr<CENuclide> n)
        : fraction(f), nuclide(n) {}
  };

  struct CENuclidePacket {
    std::optional<CENuclideFraction> nuclide_1;
    std::optional<CENuclideFraction> nuclide_2;
  };

  NDDirectory(const std::filesystem::path& fname,
              TemperatureInterpolation interp);

  CENuclidePacket load_nuclide(const std::string& key, double T);

  TemperatureInterpolation interpolation() const { return interp_; }

  const std::string& name() const { return name_; }

  const std::string& date() const { return date_; }

  const std::string& code() const { return code_; }

  const std::string& notes() const { return notes_; }

 private:
  struct ACEEntry {
    std::filesystem::path fname;
    pndl::ACE::Type type;
    double temperature;

    ACEEntry(const std::filesystem::path& basename, const YAML::Node& node);
  };

  struct NeutronACE {
    ACEEntry ace_entry;
    std::shared_ptr<pndl::STNeutron> neutron_data;

    NeutronACE(const std::filesystem::path& basename, const YAML::Node& node);
    bool loaded() const { return neutron_data != nullptr; }
    double temperature() const { return ace_entry.temperature; }
  };

  struct TSLACE {
    ACEEntry ace_entry;
    std::shared_ptr<pndl::STThermalScatteringLaw> tsl_data;

    TSLACE(const std::filesystem::path& basename, const YAML::Node& node);
    bool loaded() const { return tsl_data != nullptr; }
    double temperature() const { return ace_entry.temperature; }
  };

  struct NeutronACEList {
    std::vector<NeutronACE> neutron_ace_files;
    std::shared_ptr<pndl::STNeutron> first_loaded;

    NeutronACEList() {}
    NeutronACEList(const std::filesystem::path& basename,
                   const YAML::Node& node);
    const std::shared_ptr<pndl::STNeutron>& get_temp(const std::string& key,
                                                     double T);
  };

  struct TSLACEList {
    std::vector<TSLACE> tsl_ace_files;

    TSLACEList() {}
    TSLACEList(const std::filesystem::path& basename, const YAML::Node& node);
    const std::shared_ptr<pndl::STThermalScatteringLaw>& get_temp(
        const std::string& key, double T);
  };

  struct NuclideEntry {
    std::string neutron;
    bool dbrc;
    std::optional<std::string> tsl;
    std::vector<double> temps;
    std::vector<std::shared_ptr<CENuclide>> loaded;

    NuclideEntry() {}
    NuclideEntry(const YAML::Node& node);
    int closest_temp(double T) const;
    std::pair<int, int> bounding_temps(double T) const;
  };

  // Build in decending order !
  std::filesystem::path basename_;
  std::map<std::string, NeutronACEList> neutron_dir;
  std::map<std::string, TSLACEList> tsl_dir;
  std::map<std::string, NuclideEntry> nuclides;
  std::string name_;
  std::string date_;
  std::string code_;
  std::string notes_;
  TemperatureInterpolation interp_;

  bool has_nuclide_entry(const std::string& key) const;
  NuclideEntry& get_nuclide_entry(const std::string& key);

  bool has_neutron_list(const std::string& key) const;
  NeutronACEList& get_neutron_list(const std::string& key);

  bool has_tsl_list(const std::string& key) const;
  TSLACEList& get_tsl_list(const std::string& key);

  std::shared_ptr<CENuclide> get_cenuclide(const std::string& key,
                                           NuclideEntry& nuclide, double T);
};

#endif
