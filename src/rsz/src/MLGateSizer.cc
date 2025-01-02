#include "MLGateSizer.hh"

#include <iostream>

#include "db_sta/dbNetwork.hh"
#include "rsz/Resizer.hh"
#include "sta/Corner.hh"
#include "sta/DcalcAnalysisPt.hh"
#include "sta/Fuzzy.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/InputDrive.hh"
#include "sta/Liberty.hh"
#include "sta/Parasitics.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/PathRef.hh"
#include "sta/PathVertex.hh"
#include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/Search.hh"
#include "sta/TimingArc.hh"
#include "sta/Units.hh"
#include "utl/Logger.h"

namespace rsz {

MLGateSizer::MLGateSizer(Resizer* resizer) : resizer_(resizer)
{
  // Initialize model, load weights, etc.
}

void MLGateSizer::init()
{
  logger_ = resizer_->logger_;
  dbStaState::init(resizer_->sta_);
  db_network_ = resizer_->db_network_;
}

void MLGateSizer::loadWeights(const std::string& weight_file)
{
  // Load transformer model weights from PyTorch model and store in Eigen matrix
  // Example: transformer_weights_ = LoadFromPyTorch(weight_file);
}

void MLGateSizer::addToken(const std::vector<float>& pin_data,
                           const std::string& gate_type)
{
  pin_tokens_.push_back(pin_data);
  gate_types_.push_back(gate_type);
}

// Eigen-related and sta_-related errors seem to exist in the following
// functions Commented out for now to avoid compilation errors, try to fix later
/*
  void MLGateSizer::classifyAndResize()
  {
    // Extract worst paths from the design
    extractWorstPaths();
    // Transform pin data into the format expected by the transformer model
    applyTransformerModel();

    // Once classification is done, resize gates based on predictions
    //resizeGates();
  }

  void MLGateSizer::extractWorstPaths()
  {
    // Example logic to extract worst timing paths using STA
    VertexSeq worst_paths;
    sta_->worstPaths(3, true, worst_paths); // per-endpoint

    for (auto& path : worst_paths) {
      // For each path, process the pins
      for (const Pin* pin : path->pins()) {
        // Extract pin details like slew, cap, slack, etc.
        float slew = sta_->slew(pin);
        float cap = sta_->cap(pin);
        float slack = sta_->slack(pin);
        std::string gate_type = sta_->gateType(pin);

        // Collect data for ML model
        transformPinDataToTokens(pin, slew, cap, slack, gate_type);
      }
    }
  }

  void MLGateSizer::transformPinDataToTokens(const Pin* pin, float slew, float
  cap, float slack, const std::string& gate_type)
  {
    // Turn the extracted pin data into tokens that can be used by the
  transformer model std::vector<float> tokenized_data = {slew, cap, slack};

    // Add token to the gateSizer for later classification
    gateSizer.addToken(tokenized_data, gate_type);
  }

  void MLGateSizer::applyTransformerModel()
  {
    // Convert pin tokens into an Eigen matrix
    Eigen::MatrixXf input_tokens(pin_tokens_.size(), pin_tokens_[0].size());
    for (size_t i = 0; i < pin_tokens_.size(); ++i) {
      for (size_t j = 0; j < pin_tokens_[i].size(); ++j) {
        input_tokens(i, j) = pin_tokens_[i][j];
      }
    }

    // Perform transformer model calculations
    feedforwardNetwork(input_tokens);
  }

  void MLGateSizer::feedforwardNetwork(Eigen::MatrixXf& input)
  {
    // Implement feedforward network using Eigen operations and loaded
  transformer weights
    // This could include positional encoding, multihead attention, and the like
    Eigen::MatrixXf output = input * transformer_weights_; // Example of matrix
  multiplication
    // Apply activation functions, attention mechanisms, etc.
  }

  void MLGateSizer::resizewithML()
  {
    // Based on the output from the transformer model, resize gates
    for (size_t i = 0; i < gate_types_.size(); ++i) {
      std::string gate_type = gate_types_[i];
      // Output the classification result (e.g., gate size)
      std::cout << "Resizing gate of type: " << gate_type << " to predicted
  size" << std::endl;
      // Resize the gate accordingly in the design
    }
  }

*/

void MLGateSizer::getEndpointAndCriticalPaths()
{
  // Print out statement to indicate the function is running
  std::cout << "Retrieving endpoints and critical paths..." << std::endl;

  // Retrieve endpoints and critical paths for debugging or further analysis
  init();
  // sta::Network* network = sta_->network();
  // sta::Graph* graph = sta_->graph();

  // Retrieve endpoints
  sta::VertexSet* endpoints = sta_->endpoints();

  sta::PinSet* pinset = new sta::PinSet(network_);

  // Print out the number of endpoints found
  std::cout << "Found " << endpoints->size() << " endpoints." << std::endl;

  // Iterate through first 10 endpoints terminate early if there aren't enough
  int count = 0;
  for (const Vertex* endpoint : *endpoints) {
    if (count >= 10) {
      break;
    }

    // Retrieve the name of the endpoint
    std::string endpoint_name = network_->name(endpoint->pin());
    std::cout << "Endpoint (Vertex): " << endpoint->name(network_) << std::endl;
    std::cout << "Endpoint (Pin): " << endpoint_name << std::endl;

    // Add the endpoint to the pinset
    pinset->insert(endpoint->pin());

    // Increment the count
    count++;
  }

  std::cout << "Debug Point 1" << std::endl;
  sta::ExceptionTo* exception_to = new sta::ExceptionTo(
      pinset, nullptr, nullptr, nullptr, nullptr, true, network_);
  std::cout << "Debug Point 2" << std::endl;

  // Check if exception_to is properly initialized
  if (exception_to == nullptr) {
    std::cout << "ExceptionTo is null" << std::endl;
  } else {
    std::cout << "ExceptionTo is not null" << std::endl;
  }
  // Check if pinset is properly initialized
  if (pinset == nullptr) {
    std::cout << "Pinset is null" << std::endl;
  } else {
    std::cout << "Pinset is not null" << std::endl;
    std::cout << "Pinset size: " << pinset->size() << std::endl;
  }

  // Retrieve the critical path for the endpoint
  sta::PathEndSeq path_ends = sta_->search()->findPathEnds(
      nullptr,
      nullptr,
      nullptr,  // exception_to, // test if exception_to is causing the issue
      false,
      sta_->cmdCorner(),
      sta::MinMaxAll::max(),  // using min leads to no critical paths found (max
                              // seems to be setup time slack which is more
                              // relevant for gatesizing min/holdtime slack have
                              // to be fixed with buffer insertion)
      10,                     // group_count
      10,                     // endpoint_count
      true,                   // unique_pins
      -sta::INF,
      sta::INF,  // slack_min, slack_max
      true,      // sort_by_slack
      nullptr,   // group_names
      true,
      false,
      false,
      false,
      false,
      false);
  std::cout << "Debug Point 3" << std::endl;
  // If no critical path is found, print a message
  if (path_ends.empty()) {
    std::cout << "No critical paths found " << std::endl;
  } else {
    // Print out the critical path
    // for (const sta::PathEnd* path_end : path_ends) {
    //  path_end->reportPath(std::cout, network_, graph, 2);
    //}
    int path_count = 0;
    // Declare tempoary vector to store the slack of each path
    std::vector<float> path_slacks;
    // Get all clk_nets to check if the pins are within sequential cells
    std::set<dbNet*> clk_nets = sta_->findClkNets();

    for (auto& path_end : path_ends) {  // similar usage found in TritonPart.cpp
                                        // BuildTimingPaths()
      std::cout << "Critical Path " << path_count << std::endl;
      auto* path = path_end->path();
      float slack = path_end->slack(sta_);
      path_slacks.push_back(slack);
      std::cout << "Slack: " << slack << std::endl;
      sta::PathExpanded expand(path, sta_);
      expand.path(expand.size() - 1);
      float p2p_dist = 0.0;
      float prev_x = 0.0;
      float prev_y = 0.0;
      bool is_port = false;        // check if pin is a port, port shouldn't be
                                   // included in transsizer data
      bool is_sequential = false;  // check if the cell is a sequential cell
      sta::Pin* prev_pin = nullptr; // used to store previous pin for arc delay
      const sta::DcalcAnalysisPt* dcalc_ap = path->dcalcAnalysisPt(sta_); // used to get arc delay
      for (size_t i = 0; i < expand.size(); i++) {
        // PathRef is reference to a path vertex
        sta::PathRef* ref = expand.path(i);
        sta::Pin* pin = ref->vertex(sta_)->pin();
        // Pin name
        std::cout << "Pin(" << path_count << "-" << i
                  << "): " << network_->name(pin) << std::endl;
        // Data to extract from pin:
        //[x, y, p2p_dist, hpwl, wire_cap, arc_delay, fanout, reach_end,
        //gate_type_id, mdelay, num_refs]
        // X, Y coordinates of the pin
        Point pin_loc = db_network_->location(pin);
        std::cout << "X: " << pin_loc.x() << std::endl;
        std::cout << "Y: " << pin_loc.y() << std::endl;

        // Pin-to-pin distance (p2p_dist)
        if (i > 0) {
          p2p_dist = std::sqrt(std::pow(pin_loc.x() - prev_x, 2)
                               + std::pow(pin_loc.y() - prev_y, 2));
        } else {
          p2p_dist = 0.0;
        }
        std::cout << "Pin-to-Pin Distance: " << p2p_dist << std::endl;

        // HPWL, first retrieve the net connected to the pin, then retrieve the
        // HPWL from the net
        if (i > 0) {
          // get max, min x and y coordinates from net and calculate HPWL
          // 0.5 * multiplier seems incorrect, should be 1.0 (compare with
          // TransSizer python code and fix if needed to unify)
          std::cout << "HPWL: "
                    << 0.5
                           * (std::abs(pin_loc.x() - prev_x)
                              + std::abs(pin_loc.y() - prev_y))
                    << std::endl;
        } else {
          std::cout << "HPWL: 0" << std::endl;
        }
        // sta::Net* net = network_->net(pin);
        //  check if pin is actually a pin or if it is a port or net
        if (network_->isTopLevelPort(
                pin))  // if pin is a port no need to check for net, also
                       // shouldn't be included in transsizer data
        {
          std::cout << "Current Pin is a port" << std::endl;
          is_port = true;
        } else {
          std::cout << "Current Pin is a pin" << std::endl;
          std::cout << "Net-Name: " << network_->name(network_->net(pin))
                    << std::endl;
          is_port = false;
        }

        // Wire capacitance

        // Arc delay
        // Calculates the arc delay between the current pin and the previous pin of path.
        // Skips the first element of path as there is no previous pin (first element is usually a port)
        // Method taken from repairPath() at RepairSetup.cc. 
        if (prev_pin != nullptr) {
          sta::TimingArc* prev_arc = expand.prevArc(i);
          sta::Edge* prev_edge = ref->prevEdge(prev_arc, sta_);
          sta::Delay delay = graph_->arcDelay(prev_edge, prev_arc, dcalc_ap->index());
          std::cout << "Arc Delay: " << delay << std::endl;
        }
        // Fanout

        // Reach end (reachable endpoints from the pin)
        int reachable_endpoints = 0;
        sta::Net* pin_net = network_->net(pin);
        if (pin_net) {
          // Use connectedPinIterator to get all pins connected to this net
          sta::NetConnectedPinIterator* pin_iter
              = network_->connectedPinIterator(pin_net);
          while (pin_iter->hasNext()) {
            const Pin* connected_pin = pin_iter->next();
            Vertex* connected_vertex = graph_->pinLoadVertex(connected_pin);
            // Check if the connected pin is an endpoint
            if (search_->isEndpoint(connected_vertex)) {
              reachable_endpoints++;
            }
          }
          delete pin_iter;
        }
        std::cout << "Reachable Endpoints: " << reachable_endpoints
                  << std::endl;

        // Cell type (gate type), retrieve the cell from the pin, then retrieve
        // the cell type
        // std::cout << "Pin's Cell Name: " <<
        // network_->name(network_->instance(pin)) << std::endl; std::cout <<
        // "Cell Type: " <<
        // network_->libertyCell(network_->instance(pin))->name() << std::endl;

        if (is_port == false) {
          std::cout << "Pin's Cell Name: "
                    << network_->name(network_->instance(pin)) << std::endl;
          std::cout << "Cell Type: "
                    << network_->libertyCell(network_->instance(pin))->name()
                    << std::endl;
        } else  // if current pin is a port, then the cell name or cell type is
                // not applicable
        {
          std::cout << "Cell Type: Port" << std::endl;
        }

        // Check if the cell is a sequential cell (as only combinational cells
        // are considered for gate sizing in TransSizer) Retrieve the clk_nets
        // and check if the pin_net is in the clk_nets If it is, then the cell
        // is a sequential cell If it is not, then the cell is a combinational
        // cell
        if (is_port == false) {
          sta::Net* sta_net = network_->net(pin);
          if (sta_net) {
            dbNet* db_net = db_network_->staToDb(sta_net);
            // Check if current pin's net is in clock nets
            if (clk_nets.find(db_net) != clk_nets.end()) {
              std::cout << "Is In Clock Nets: True" << std::endl;
            } else {
              std::cout << "Is In Clock Nets: False" << std::endl;
            }
          }
        }

        // Get the instance and cell containing the pin
        Instance* inst = network_->instance(pin);
        if (inst) {
          LibertyCell* cell = network_->libertyCell(inst);
          if (cell) {
            // Method 1: Check if cell has any sequential elements
            if (!cell->sequentials().empty()) {
              is_sequential = true;
            }
          }
        }
        // Is the cell in clock?
        bool is_in_clock = false;
        dbInst* db_inst = db_network_->staToDb(network_->instance(pin));
        if (db_inst) {
          for (odb::dbITerm* iterm : db_inst->getITerms()) {
            dbNet* net = iterm->getNet();
            if (net != nullptr && net->getSigType() == odb::dbSigType::CLOCK) {
              is_in_clock = true;
              break;
            } 
          }
        }

        std::cout << "Is In Clock: " << is_in_clock << std::endl;

        std::cout << "Is Sequential: " << is_sequential << std::endl;

        // maxcap
        // maxtran
        // tran


        // Slack of the pin
        std::cout << "Slack (max):" << sta_->pinSlack(pin, sta::MinMax::max())
                  << std::endl;
        std::cout << "Slack (min):" << sta_->pinSlack(pin, sta::MinMax::min())
                  << std::endl;

        std::cout << "Rise Slack (max):"
                  << sta_->pinSlack(
                         pin, sta::RiseFall::rise(), sta::MinMax::max())
                  << std::endl;
        std::cout << "Fall Slack (max):"
                  << sta_->pinSlack(
                         pin, sta::RiseFall::fall(), sta::MinMax::max())
                  << std::endl;
        std::cout << "Rise Slack (min):"
                  << sta_->pinSlack(
                         pin, sta::RiseFall::rise(), sta::MinMax::min())
                  << std::endl;
        std::cout << "Fall Slack (min):"
                  << sta_->pinSlack(
                         pin, sta::RiseFall::fall(), sta::MinMax::min())
                  << std::endl;

        //rise arrival time
        float rise_arrival_time = 0.0;
        rise_arrival_time = sta_->pinArrival(pin, sta::RiseFall::rise(), sta::MinMax::max());//Timing::getPinArrival(pin, sta::RiseFall::rise(), sta::MinMax::max());
        std::cout << "Rise Arrival Time: "
                  << rise_arrival_time << std::endl;
                  

        //fall arrival time
        float fall_arrival_time = 0.0;
        fall_arrival_time = sta_->pinArrival(pin, sta::RiseFall::rise(), sta::MinMax::max());//Timing::getPinArrival(pin, sta::RiseFall::fall(), sta::MinMax::max());
        std::cout << "Fall Arrival Time: "
                  << fall_arrival_time << std::endl;

        // input_pin_cap (if not an input pin, then -1)
        // Taken from getPortCap() in Timing.cc
        float input_pin_cap = -1.0;
        // also check if the "pin" is a port or a real pin
        if (is_port == false) {
          if (network_->direction(pin) == sta::PortDirection::input()) {
            sta::Corner* corner = sta_->cmdCorner();
            sta::LibertyPort* lib_port = network_->libertyPort(pin);
            input_pin_cap = sta_->capacitance(lib_port, corner, sta::MinMax::max());
          }
        }
        std::cout << "Input Pin Cap: " << input_pin_cap << std::endl;

        std::cout << std::endl;
        prev_x = pin_loc.x();
        prev_y = pin_loc.y();
        prev_pin = pin;
      }
      path_count++;
    }
    // Print out the slack of each path (for debugging), remove later
    for (size_t i = 0; i < path_slacks.size(); i++) {
      std::cout << "Path " << i << " Slack: " << path_slacks[i] << std::endl;
    }
  }

  // In addition or alternatively,
  // consider using vertexWorstSlackPath to find the critical path for each
  // endpoint

  // Next step is extracting the data from each pin
  // Then, transform the data into the format expected by the transformer model
  // Finally, apply the transformer model to classify and resize the gates

  // Use PinMetrics and PinDataSequence to store the extracted data

  // Tempoarily commented out to avoid compilation errors, try to fix later
  /*
  if (endpoints.empty()) {
    std::cout << "No endpoints found." << std::endl;
    return;
  } else {
    std::cout << "Found " << endpoints.size() << " endpoints." << std::endl;

    // Iterate through each endpoint and retrieve the critical path for each
  endpoint for (const Pin* endpoint : endpoints) { std::string endpoint_name =
  network_->name(endpoint); std::cout << "Endpoint: " << endpoint_name <<
  std::endl;

      PathEndSeq path_ends = sta->findPathEnds(
        nullptr, nullptr, new ExceptionTo(endpoint), false, nullptr,
        MinMaxAll::max(), 1, 1, false, -1e30, 1e30, true, nullptr,
        true, false, false, false, false, false);

      if (path_ends.empty()) {
        std::cout << "No critical path found for endpoint " << endpoint_name <<
  std::endl; } else { for (const PathEnd* path_end : path_ends) {
          path_end->reportPath(std::cout, network_, graph, 2);
        }
      }
    }
  }
  */
}

}  // namespace rsz
