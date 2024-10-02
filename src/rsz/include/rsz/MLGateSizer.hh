#ifndef RSZ_ML_GATE_SIZER_HH
#define RSZ_ML_GATE_SIZER_HH

#include <Eigen/Dense>
#include <vector>
#include <string>

namespace rsz {

class MLGateSizer {
 public:
  MLGateSizer();
  
  void loadWeights(const std::string& weight_file);
  void addToken(const std::vector<float>& pin_data, const std::string& gate_type);
  void classifyAndResize();

 private:
  void applyTransformerModel();
  void feedforwardNetwork(Eigen::MatrixXf& input);
  void resizeGates();

  std::vector<std::vector<float>> pin_tokens_;
  std::vector<std::string> gate_types_;
  Eigen::MatrixXf transformer_weights_;
};

} // namespace rsz

#endif

