#pragma once

#include <Eigen/Dense>
#include "Resizer.hh"
#include <vector>
#include <string>

namespace rsz {

class MLGateSizer : public Resizer{
 public:
  MLGateSizer();
  
  void loadWeights(const std::string& weight_file);
  void addToken(const std::vector<float>& pin_data, const std::string& gate_type);
  void classifyAndResize();

  // Function to retrieve endpoints and critical paths (also used for debugging)
  void getEndpointsAndCriticalPaths();

 private:
  void extractWorstPaths();
  void transformPinDataToTokens();
  void applyTransformerModel();
  void feedforwardNetwork(Eigen::MatrixXf& input);
  void resizeGates();

  std::vector<std::vector<float>> pin_tokens_;
  std::vector<std::string> gate_types_;
  Eigen::MatrixXf transformer_weights_;
};

} // namespace rsz


