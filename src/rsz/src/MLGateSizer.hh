#pragma once

// #include <Eigen/Dense>  
#include "db_sta/dbSta.hh" // Includes taken from RecoverPower.hh
#include "sta/FuncExpr.hh"
#include "sta/Graph.hh"
#include "sta/MinMax.hh"
#include "sta/StaState.hh"
#include "utl/Logger.h"

#include <vector>
#include <string>

namespace sta {
class PathExpanded;
}

namespace rsz {

// Initial Declaration taken from RecoverPower.hh
class Resizer;

using std::vector;

using utl::Logger;

using sta::Corner;
using sta::dbNetwork;
using sta::dbSta;
using sta::DcalcAnalysisPt;
using sta::LibertyCell;
using sta::LibertyPort;
using sta::MinMax;
using sta::Net;
using sta::PathExpanded;
using sta::PathRef;
using sta::Pin;
using sta::Slack;
using sta::StaState;
using sta::TimingArc;
using sta::Vertex;


class MLGateSizer : public sta::dbStaState
{
 public:
  MLGateSizer(Resizer* resizer);
  
  void loadWeights(const std::string& weight_file);
  void addToken(const std::vector<float>& pin_data, const std::string& gate_type);
  // void classifyAndResize();

  // Function to retrieve endpoints and critical paths (also used for debugging)
  void getEndpointAndCriticalPaths();
  // void resizewithML();

 private:
  void init();
  // void extractWorstPaths();
  // void transformPinDataToTokens();
  // void applyTransformerModel();
  // void feedforwardNetwork(Eigen::MatrixXf& input);

  Logger* logger_ = nullptr;
  dbNetwork* network_ = nullptr;
  Resizer* resizer_;


  std::vector<std::vector<float>> pin_tokens_;
  std::vector<std::string> gate_types_;
  // Eigen::MatrixXf transformer_weights_;
};

} // namespace rsz


