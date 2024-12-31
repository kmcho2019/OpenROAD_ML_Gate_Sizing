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
using sta::PinSeq;
using sta::PinSet;
using sta::Slack;
using sta::StaState;
using sta::TimingArc;
using sta::Vertex;
using sta::VertexSet;

// Forward declaration of PinMetrics to avoid circular dependencies if Resizer needs to know about it.
struct PinMetrics;

// Type alias for a sequence of PinMetrics (for better readability)
using PinDataSequence = std::vector<PinMetrics>;


class MLGateSizer : public sta::dbStaState
{
 public:
  MLGateSizer(Resizer* resizer);

  // Structure to represent the metrics associated with a single pin
  struct PinMetrics {
      std::string pinName;
      double pin_tran;
      double slack;
      double input_pin_cap;
      double output_pin_cap;
      double load_cap;

      // Constructor for easier initialization
      PinMetrics(const std::string& name, double tran, double slack, double in_cap, double out_cap, double load_cap)
          : pinName(name), pin_tran(tran), slack(slack), input_pin_cap(in_cap), output_pin_cap(out_cap), load_cap(load_cap) {};
  };  
  
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
  dbNetwork* db_network_ = nullptr;
  Resizer* resizer_;


  std::vector<std::vector<float>> pin_tokens_;
  std::vector<std::string> gate_types_;
  // Eigen::MatrixXf transformer_weights_;
};

} // namespace rsz


