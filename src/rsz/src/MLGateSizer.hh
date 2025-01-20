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
#include <unordered_map>
#include <tuple>
#include <iostream>
#include <chrono>

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

// Forward declaration of PinMetrics and PinSequenceCollector to avoid circular dependencies if Resizer needs to know about it.
struct PinMetrics;
class PinSequenceCollector;
class SequenceArrayBuilder;


class MLGateSizer : public sta::dbStaState
{
 public:
  MLGateSizer(Resizer* resizer);

  // Structure to represent the metrics associated with a single pin
  struct PinMetrics {
      std::string pin_name;
      std::string cell_name;
      std::string cell_type;
      float x_loc{0.0};
      float y_loc{0.0};
      float p2p_dist{0.0};
      float hpwl{0.0};
      float input_pin_cap{-1.0};
      float wire_cap{0.0};
      float pin_cap{0.0};
      float total_cap{0.0};
      float fanout{0.0};
      float arc_delay{0.0};
      float reachable_endpoints{0.0};
      bool is_in_clock_nets{false};
      bool is_in_clock{false};
      bool is_sequential{false};
      bool is_macro{false};
      bool is_port{false};
      float max_cap{0.0};
      float max_slew{0.0};
      float rise_slew{0.0};
      float fall_slew{0.0};
      float slack{0.0}; // use maximum slack (setup time) as it is the one relevant for gate sizing, min/hold time slack have to be fixed with buffer insertion
      float rise_arrival_time{0.0};
      float fall_arrival_time{0.0};

      // Default Constructor
      PinMetrics() {}


      // Constructor for easier initialization
    // Constructor for easier initialization
    PinMetrics(const std::string& pin_name, const std::string& cell_name, const std::string& cell_type,
               float x_loc, float y_loc, float p2p_dist, float hpwl, float input_pin_cap, float wire_cap,
               float pin_cap, float total_cap, float fanout, float arc_delay, float reachable_endpoints,
               bool is_in_clock_nets, bool is_in_clock, bool is_sequential, bool is_macro, bool is_port,
               float max_cap, float max_slew, float rise_slew, float fall_slew, float slack,
               float rise_arrival_time, float fall_arrival_time)
        : pin_name(pin_name), cell_name(cell_name), cell_type(cell_type), x_loc(x_loc), y_loc(y_loc),
          p2p_dist(p2p_dist), hpwl(hpwl), input_pin_cap(input_pin_cap), wire_cap(wire_cap), pin_cap(pin_cap),
          total_cap(total_cap), fanout(fanout), arc_delay(arc_delay), reachable_endpoints(reachable_endpoints),
          is_in_clock_nets(is_in_clock_nets), is_in_clock(is_in_clock), is_sequential(is_sequential),
          is_macro(is_macro), is_port(is_port), max_cap(max_cap), max_slew(max_slew), rise_slew(rise_slew),
          fall_slew(fall_slew), slack(slack), rise_arrival_time(rise_arrival_time),
          fall_arrival_time(fall_arrival_time) {}

      // Print function for debugging
      void print() const {
          std::cout << "Pin Name: " << pin_name << std::endl;
          std::cout << "Cell Name: " << cell_name << std::endl;
          std::cout << "Cell Type: " << cell_type << std::endl;
          std::cout << "X: " << x_loc << std::endl;
          std::cout << "Y: " << y_loc << std::endl;
          std::cout << "Pin-to-Pin Distance: " << p2p_dist << std::endl;
          std::cout << "HPWL: " << hpwl << std::endl;
          std::cout << "Input Pin Cap: " << input_pin_cap << std::endl;
          std::cout << "Wire Cap: " << wire_cap << std::endl;
          std::cout << "Pin Cap: " << pin_cap << std::endl;
          std::cout << "Total Connected Cap: " << total_cap << std::endl;
          std::cout << "Fanout: " << fanout << std::endl;
          std::cout << "Arc Delay: " << arc_delay << std::endl;
          std::cout << "Reachable Endpoints: " << reachable_endpoints << std::endl;
          std::cout << "Is In Clock Nets: " << is_in_clock_nets << std::endl;
          std::cout << "Is In Clock: " << is_in_clock << std::endl;
          std::cout << "Is Port: " << is_port << std::endl;
          std::cout << "Is Sequential: " << is_sequential << std::endl;
          std::cout << "Is Macro: " << is_macro << std::endl;
          std::cout << "Max Cap: " << max_cap << std::endl;
          std::cout << "Max Slew: " << max_slew << std::endl;
          std::cout << "Rise/Fall Slew: " << rise_slew << "/" << fall_slew << std::endl;
          std::cout << "Slack: " << slack << std::endl;
          std::cout << "Rise Arrival Time: " << rise_arrival_time << std::endl;
          std::cout << "Fall Arrival Time: " << fall_arrival_time << std::endl;
      }
  };  

  class PinSequenceCollector {
  public:
      using PinSequence = std::vector<PinMetrics>;
      using SequenceCollection = std::vector<PinSequence>;
      
      void processPin(const PinMetrics& metrics) {
          if (metrics.is_port || metrics.is_sequential || metrics.is_macro) {
              // End current sequence if it has pins
              if (!current_sequence.empty()) {
                  if (current_sequence.size() > 1) { // Store only if sequence has multiple pins
                      sequences.push_back(current_sequence);
                  }
                  current_sequence.clear();
              }
          } else {
              // Add pin to current sequence if it's combinatorial
              current_sequence.push_back(metrics);
          }
      }
      
      void finalize() {
          if (!current_sequence.empty() && current_sequence.size() > 1) {
              sequences.push_back(current_sequence);
              current_sequence.clear();
          }
      }
      
      const SequenceCollection& getSequences() const { return sequences; }
      
  private:
      PinSequence current_sequence;
      SequenceCollection sequences;
  };
  
  class SequenceArrayBuilder {
  private:
      const PinSequenceCollector::SequenceCollection& sequences_;
      const std::unordered_map<std::string, int>& pin_name_to_id_;
      const std::unordered_map<std::string, int>& cell_name_to_id_;
      const std::unordered_map<std::string, int>& cell_type_to_id_;
      const std::unordered_map<int, std::vector<float>>& cell_type_embeddings_;
      
      size_t max_seq_len_;
      size_t embedding_dim_;
      size_t num_numerical_features_;
      
  public:
      SequenceArrayBuilder(
          const PinSequenceCollector::SequenceCollection& sequences,
          const std::unordered_map<std::string, int>& pin_id_map,
          const std::unordered_map<std::string, int>& cell_id_map,
          const std::unordered_map<std::string, int>& cell_type_map,
          const std::unordered_map<int, std::vector<float>>& embeddings)
          : sequences_(sequences),
            pin_name_to_id_(pin_id_map),
            cell_name_to_id_(cell_id_map),
            cell_type_to_id_(cell_type_map),
            cell_type_embeddings_(embeddings) {
          
          max_seq_len_ = findMaxSeqLen();
          embedding_dim_ = embeddings.begin()->second.size();
          num_numerical_features_ = 18; // Number of numerical features in PinMetrics
      }
      
      std::tuple<std::vector<std::vector<std::vector<float>>>,
                std::vector<std::vector<int>>,
                std::vector<std::vector<int>>,
                std::vector<std::vector<int>>> 
      build() {
          size_t N = sequences_.size();
          size_t L = max_seq_len_;
          size_t D = num_numerical_features_ + embedding_dim_;
          
          // Initialize arrays with padding
          std::vector<std::vector<std::vector<float>>> data_array(
              N, std::vector<std::vector<float>>(L, std::vector<float>(D, 0.0)));
          std::vector<std::vector<int>> pin_ids(N, std::vector<int>(L, -1));
          std::vector<std::vector<int>> cell_ids(N, std::vector<int>(L, -1));
          std::vector<std::vector<int>> cell_type_ids(N, std::vector<int>(L, -1));
          
          // Fill arrays
          for (size_t i = 0; i < N; i++) {
              const auto& sequence = sequences_[i];
              for (size_t j = 0; j < sequence.size(); j++) {
                  const auto& metrics = sequence[j];
                  
                  // Fill numerical features
                  std::vector<float> features = getNumericalFeatures(metrics);
                  std::copy(features.begin(), features.end(), data_array[i][j].begin());
                  
                  // Add embedding
                  int type_id = cell_type_to_id_.at(metrics.cell_type);
                  const auto& embedding = cell_type_embeddings_.at(type_id);
                  std::copy(embedding.begin(), embedding.end(), 
                          data_array[i][j].begin() + num_numerical_features_);
                  
                  // Fill lookup arrays
                  pin_ids[i][j] = pin_name_to_id_.at(metrics.pin_name);
                  cell_ids[i][j] = cell_name_to_id_.at(metrics.cell_name);
                  cell_type_ids[i][j] = type_id;
              }
          }
          
          return {data_array, pin_ids, cell_ids, cell_type_ids};
      }
      
  private:
      size_t findMaxSeqLen() {
          size_t max_len = 0;
          for (const auto& seq : sequences_) {
              max_len = std::max(max_len, seq.size());
          }
          return max_len;
      }
      
      std::vector<float> getNumericalFeatures(const PinMetrics& metrics) {
          return {
              metrics.x_loc,
              metrics.y_loc,
              metrics.p2p_dist,
              metrics.hpwl,
              metrics.input_pin_cap,
              metrics.wire_cap,
              metrics.pin_cap,
              metrics.total_cap,
              metrics.arc_delay,
              metrics.max_cap,
              metrics.max_slew,
              metrics.rise_slew,
              metrics.fall_slew,
              metrics.slack,
              metrics.rise_arrival_time,
              metrics.fall_arrival_time,
              static_cast<float>(metrics.fanout),
              static_cast<float>(metrics.reachable_endpoints)
          };
      }
  };

  // Naive encoder-like transformer forward pass (initial implementation)
  std::vector<std::vector<std::vector<float>>> runTransformer(
      const std::vector<std::vector<std::vector<float>>>& data_array,
      int num_heads,
      size_t N,
      size_t L,
      size_t D_in,
      size_t D_model,
      size_t FF_hidden_dim,
      int num_encoder_layers);

  // Eigen-based version, more efficient
  std::vector<std::vector<std::vector<float>>> runTransformerEigen(
      const std::vector<std::vector<std::vector<float>>>& data_array,
      int num_heads,
      size_t N,
      size_t L,
      size_t D_in,
      size_t D_model,
      size_t FF_hidden_dim,
      int num_encoder_layers);

  // Compare two transformer outputs for correctness
  bool compareOutputs(
      const std::vector<std::vector<std::vector<float>>>& A,
      const std::vector<std::vector<std::vector<float>>>& B,
      float tol = 1e-4f) const;

  // Simple timing utility that returns elapsed microseconds
  template <class Func>
  long long benchmark(Func&& func);

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


