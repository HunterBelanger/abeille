#include <materials/material_helper.hpp>
#include <materials/nuclide.hpp>

MaterialHelper::MaterialHelper(Material* material, double E)
    : mat(material), E_(E), xs_(), zaid_to_urr_rand_() {
  // Add all nuclides to the xs_ map
  for (const auto& nuc : nuclides) {
    xs_[nuc.second.get()] = std::nullopt;
  }

  // For all ZAIDs we found with a URR, add them to the
  // zaid_to_urr_rand_ map
  for (const auto& za : zaids_with_urr) {
    zaid_to_urr_rand_[za] = std::nullopt;
  }
}