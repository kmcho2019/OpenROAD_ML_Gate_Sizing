#include "MLGateSizer.hh"



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

  // Print out the number of endpoints found
  std::cout << "Found " << endpoints->size() << " endpoints." << std::endl;


  std::cout << "Debug Point 1" << std::endl;

  std::cout << "Debug Point 2" << std::endl;

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
    // Get all clk_nets to check if the pins are connected to a clock net
    std::set<dbNet*> clk_nets = sta_->findClkNets();
    // Used during data retrieval but stays uniform for all paths
    sta::LibertyLibrary* lib = network_->defaultLibertyLibrary();
    sta::Corner* corner = sta_->cmdCorner();

    // Get Library Cell to index mapping
    std::unordered_map<std::string, int> libcell_to_id;
    std::unordered_map<std::string, int> libcell_to_type_id; // equivalent libcells share the same, Uses EquivCells
    std::unordered_map<int, std::vector<float>> libcell_to_type_embedding;
    std::unordered_map<std::string, int> pin_name_to_id; // Initialize them as empty and build them up during pin retrieval process
    std::unordered_map<std::string, int> cell_name_to_id; // Initialize them as empty and build them up during pin retrieval process

    int libcell_id = 0;
    int libcell_type_id = 0;
    int pin_id = 0;
    int cell_id = 0;
    odb::dbDatabase* db = sta_->db();
    for (odb::dbLib* lib : db->getLibs()) {
      for (odb::dbMaster* master : lib->getMasters()) {
        std::string cell_name = master->getName();
        if (libcell_to_id.find(cell_name) == libcell_to_id.end()) {
          libcell_to_id[cell_name] = libcell_id++;
        }
      }
    }
    
    // Initialize libcell_to_type_id 
    // EquivCell part taken from Timing::EquivCells()
    for (odb::dbLib* lib : db->getLibs()) {
      for (odb::dbMaster* master : lib->getMasters()) {
        sta::Cell* sta_cell = db_network_->dbToSta(master);
        std::string cell_name = master->getName();
        std::vector<odb::dbMaster*> masterSeq;
        if (sta_cell) {
          sta::LibertyCell* lib_cell = network_->libertyCell(sta_cell);
          sta::LibertyCellSeq* equiv_cells = sta_->equivCells(lib_cell);
          if (equiv_cells) {
            for (sta::LibertyCell* equiv_cell : *equiv_cells) {
              odb::dbMaster* equiv_master = db_network_->staToDb(equiv_cell);
              masterSeq.emplace_back(equiv_master);
            } 
          } else {
            masterSeq.emplace_back(master);
          }
        }

        // Check if any of the equivalent cells already have a type ID
        bool type_id_found = false;
        for (odb::dbMaster* equiv_master : masterSeq) {
          std::string equiv_cell_name = equiv_master->getName();
          if (libcell_to_type_id.find(equiv_cell_name) != libcell_to_type_id.end()) {
            // Found an equivalent cell with an existing type ID, use it
            libcell_to_type_id[cell_name] = libcell_to_type_id[equiv_cell_name];
            type_id_found = true;
            break; // No need to check further
          }
        }

        // If none of the equivalent cells have a type ID, assign a new one
        if (!type_id_found) {
          libcell_to_type_id[cell_name] = libcell_type_id++;
          // Optionally generate and assign embeddings here, e.g.,
          // cell_type_embeddings[libcell_to_type_id[cell_name]] = generateEmbedding(master);
        }

      }
    }

    // Initialize libcell_to_type_embedding with placeholder values
    // Proper embedding comes from passing through libcell name through a sentence transformer and averaging the values for same libcell type
    for (const auto& type_pair : libcell_to_type_id) {
      int type_id = type_pair.second;

      // Create a 16-float embedding filled with 0.0, 1.0, 2.0, ...
      std::vector<float> embedding(16);
      for (int i = 0; i < embedding.size(); ++i) {
        embedding[i] = static_cast<float>(type_id); // Assign the type_id as the value for all elements in the embedding
      }

      libcell_to_type_embedding[type_id] = embedding;
    }


    // Print Cell Names and ID from libcell_to_id and Type ID
    for (auto& cell : libcell_to_id) {
      std::cout << "Cell Name: " << cell.first << " ID: " << cell.second << " Type ID: " << libcell_to_type_id[cell.first] << std::endl;
    }

  


    

    
    PinSequenceCollector collector;

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
      sta::Pin* prev_pin = nullptr; // used to store previous pin for arc delay or p2p_dist
      const sta::DcalcAnalysisPt* dcalc_ap = path->dcalcAnalysisPt(sta_); // used to get arc delay



      // Data to extract from pin:
      // [x, y, p2p_dist, hpwl, wire_cap, arc_delay, fanout, reach_end, 
      // gate_type_id, mdelay, num_refs, maxcap, maxtran, tran, slack, 
      // risearr, fallarr, cap]
      for (size_t i = 0; i < expand.size(); i++) {
          // PinMetrics object initialization
          PinMetrics pin_metrics;
          // Core object retrieval
          sta::PathRef* ref = expand.path(i); // PathRef is reference to a path vertex
          sta::Pin* pin = ref->vertex(sta_)->pin();
          sta::Net* net = network_->net(pin);
          sta::Instance* inst = network_->instance(pin);
          sta::LibertyPort* lib_port = network_->libertyPort(pin);

          sta::Vertex* vertex = graph_->pinLoadVertex(pin);

          // DB object retrieval
          dbInst* db_inst = db_network_->staToDb(inst);
          dbNet* db_net = net ? db_network_->staToDb(net) : nullptr;
          odb::dbITerm* iterm;
          odb::dbBTerm* bterm;
          odb::dbModITerm* moditerm;
          odb::dbModBTerm* modbterm;
          db_network_->staToDb(pin, iterm, bterm, moditerm, modbterm);

          // Pin state flags
          is_port = network_->isTopLevelPort(pin);
          bool is_supply_pin = (iterm && iterm->getSigType().isSupply()) || 
                              (bterm && bterm->getSigType().isSupply());
          bool is_in_clock_nets = false; // check if pin is in clock nets indicates /CLK pin
          bool is_in_clock = false; // check if pin's instance is connected to a clock net, either a sequential cell or a buffer/inverter cell connected to a clock net
          bool is_sequential = false; // check if the cell is a sequential cell
          bool is_macro = false; // check if cell is a macro/block cell
          if (is_port) {
              is_in_clock = false;
          } else {
              if (net && db_net && clk_nets.find(db_net) != clk_nets.end()) {
                  is_in_clock_nets = true;
              }
          }
          if (db_inst) {
              for (odb::dbITerm* term : db_inst->getITerms()) {
                  if (term->getNet() && term->getNet()->getSigType() == odb::dbSigType::CLOCK) {
                      is_in_clock = true;
                      break;
                  }
              }
          }
          if (lib_port) { // check if the cell is a sequential cell, referenced from Resizer::isRegister()
            sta::LibertyCell* lib_cell = lib_port->libertyCell();
            is_sequential = lib_cell && lib_cell->hasSequentials();
            is_macro = lib_cell && lib_cell->isMacro();
          }

          // Location calculations
          // Skips the first element of path as there is no previous pin (first element is usually a port)
          Point pin_loc = db_network_->location(pin);
          // Pin-to-pin distance (p2p_dist)
          p2p_dist = (i > 0) ? 
              std::sqrt(std::pow(pin_loc.x() - prev_x, 2) + std::pow(pin_loc.y() - prev_y, 2)) : 0.0;
          // HPWL calculation (hpwl)
          // -0.5 * multiplier is used in original TransSizer code, but it seems incorrect
          // -Unify implementation with Python code if needed
          float hpwl = (i > 0) ? 
              0.5 * (std::abs(pin_loc.x() - prev_x) + std::abs(pin_loc.y() - prev_y)) : 0.0;

          // Capacitance calculations
          // Wire capacitance (net_cap)
          // -Retrieves the capacitance of the net connected to the pin
          // -Total connected capacitance (total_cap = pin_cap + wire_cap)
          // -Corresponds to total_cap and net_cap found in CircuitOps's net_properties.csv
          float pin_cap = 0.0;
          float wire_cap = 0.0;
          if (net) {
              sta_->connectedCap(net, corner, sta::MinMax::max(), pin_cap, wire_cap);
          }

          // Fanout calculation (fanout) (calculates number of pins connected to the net of current pin)
          // -Referenced from connectedPins() in rsz/src/SteinerTree.cc
          int fanout = 0;
          if (net) {
              sta::NetConnectedPinIterator* connected_pins = network_->connectedPinIterator(net);
              while (connected_pins->hasNext()) {
                  connected_pins->next();
                  fanout++;
              }
              delete connected_pins;
              fanout--; // Subtract current pin as it is also included in the fanout
          }

          // Arc delay calculation (arc_delay)
          // -Calculates the arc delay between the current pin and the previous pin of path.
          // -Skips the first element of path as there is no previous pin (first element is usually a port)
          // -Referenced from repairPath() at RepairSetup.cc.
          sta::Delay arc_delay = 0.0;
          if (prev_pin != nullptr) {
              sta::TimingArc* prev_arc = expand.prevArc(i);
              sta::Edge* prev_edge = ref->prevEdge(prev_arc, sta_);
              arc_delay = graph_->arcDelay(prev_edge, prev_arc, dcalc_ap->index());
          }

          // Reachable endpoints (reach_end)
          // Calculates the number of reachable endpoints from the current pin
          int reachable_endpoints = 0;
          if (net) {
              sta::NetConnectedPinIterator* pin_iter = network_->connectedPinIterator(net);
              while (pin_iter->hasNext()) {
                  const Pin* connected_pin = pin_iter->next();
                  if (search_->isEndpoint(graph_->pinLoadVertex(connected_pin))) {
                      reachable_endpoints++;
                  }
              }
              delete pin_iter;
          }

          // Timing constraints
          // Max cap and max slew calculations (maxcap, maxtran)
          // -Referenced from getMaxCapLimit() and getMaxSlewLimit() in Timing.cc
          // -Check if the pin is a ground or power pin, if it is, then max_cap is 0.0
          // -If it is not, then get max_cap/max_slew from the liberty library
          // -If max_cap/max_slew is not found, then get default max_cap/max_slew from the liberty library
          float max_cap = 0.0, max_slew = 0.0;
          bool max_cap_exists = false, max_slew_exists = false;
          if (!is_supply_pin && !is_port) {
              lib_port->capacitanceLimit(sta::MinMax::max(), max_cap, max_cap_exists);
              if (!max_cap_exists) {
                  lib->defaultMaxCapacitance(max_cap, max_cap_exists);
              }
              lib_port->slewLimit(sta::MinMax::max(), max_slew, max_slew_exists);
              if (!max_slew_exists) {
                  lib->defaultMaxSlew(max_slew, max_slew_exists);
              }
          }

          // Timing measurements
          // Rise and fall slew calculations of pin (tran)
          // -CircuitOps used in TransSizer seems to check rise transition time only, check if this is true
          // -CircuitOps seems to use getPinSlew() which is also based on slewAllCorners() in Timing.cc
          // -slewAllCorners() uses sta->vertexSlew(vertex, sta::RiseFall::rise(), corner, minmax))
          // -This seems to show that rise slew is used in CircuitOps's (tran) calculation
          // Rise and fall arrival time calculations of pin (risearr, fallarr)
          // Input pin capacitance from CircuitOp's pin_properties.csv (cap) 
          // -If pin is an input pin, then input pin capacitance is calculated
          // -If pin is not an input pin, then input pin capacitance is -1.0
          // Slack of the pin (slack)
          // -Use MinMax::max() to get max slack as used in CircuitOps
          // -Circuits Ops used the min of rise and fall slack, but sta_->pinSlack(pin, sta::MinMax::max()) does seem to do this automatically
          float rise_slew = vertex ? sta_->vertexSlew(vertex, sta::RiseFall::rise(), sta::MinMax::max()) : 0.0;
          float fall_slew = vertex ? sta_->vertexSlew(vertex, sta::RiseFall::fall(), sta::MinMax::max()) : 0.0;
          float rise_arrival_time = sta_->pinArrival(pin, sta::RiseFall::rise(), sta::MinMax::max());
          float fall_arrival_time = sta_->pinArrival(pin, sta::RiseFall::fall(), sta::MinMax::max());
          float input_pin_cap = (!is_port && network_->direction(pin) == sta::PortDirection::input()) ?
              sta_->capacitance(lib_port, corner, sta::MinMax::max()) : -1.0;
          float slack = sta_->pinSlack(pin, sta::MinMax::max());
          float min_slack = sta_->pinSlack(pin, sta::MinMax::min());

          // Cell type (gate type), retrieve the cell from the pin, then retrieve the cell type
          std::string cell_type = is_port ? "Port" : network_->libertyCell(inst)->name();
          std::string cell_name = is_port ? "Port" : network_->name(inst);

          // Network name of the pin
          std::string net_name = net ? network_->name(net) : "None";

          // Pin Name
          std::string pin_name = network_->name(pin);


          // Assign pin_id and cell_id
          if (pin_name_to_id.find(pin_name) == pin_name_to_id.end()) {
            // Pin name doesn't exist, add it to the map with the current pin_id
            pin_name_to_id[pin_name] = pin_id;
            pin_id++;
          }
          if (cell_name_to_id.find(cell_name) == cell_name_to_id.end()) {
            // Pin name doesn't exist, add it to the map with the current pin_id
            cell_name_to_id[cell_name] = cell_id;
            cell_id++;
          }

          // Output all collected data (for debugging)
          std::cout << "Pin(" << path_count << "-" << i << "): " << pin_name << "\n"
                    << "X: " << pin_loc.x() << "\n"
                    << "Y: " << pin_loc.y() << "\n"
                    << "Pin-to-Pin Distance: " << p2p_dist << "\n"
                    << "HPWL: " << hpwl << "\n"
                    << "Wire Cap: " << wire_cap << "\n"
                    << "Pin Cap: " << pin_cap << "\n"
                    << "Total Connected Cap: " << wire_cap + pin_cap << "\n"
                    << "Fanout: " << fanout << "\n"
                    << "Arc Delay: " << arc_delay << "\n"
                    << "Reachable Endpoints: " << reachable_endpoints << "\n"
                    << "Pin's Net Name: " << net_name << "\n"
                    << "Pin's Cell Name: " << cell_name << "\n"
                    << "Cell Type: " << cell_type << "\n"
                    << "Is In Clock Nets: " << is_in_clock_nets << "\n"
                    << "Is In Clock: " << is_in_clock << "\n"
                    << "Is Sequential: " << is_sequential << "\n"
                    << "Is Macro: " << is_macro << "\n"
                    << "Is Supply Pin: " << is_supply_pin << "\n"
                    << "Max Cap: " << max_cap << "\n"
                    << "Max Slew: " << max_slew << "\n"
                    << "Rise/Fall Slew: " << rise_slew << "/" << fall_slew << "\n"
                    << "Slack (max/min): " << slack << "/"
                    << min_slack << "\n"
                    << "Rise/Fall Arrival Time: " << rise_arrival_time << "/" << fall_arrival_time << "\n"
                    << "Input Pin Cap: " << input_pin_cap << "\n\n";
          // Fill in the PinMetrics object
          pin_metrics.pin_name = pin_name;
          pin_metrics.cell_name = cell_name;
          pin_metrics.cell_type = cell_type;
          pin_metrics.x_loc = pin_loc.x();
          pin_metrics.y_loc = pin_loc.y();
          pin_metrics.p2p_dist = p2p_dist;
          pin_metrics.hpwl = hpwl;
          pin_metrics.input_pin_cap = input_pin_cap;
          pin_metrics.wire_cap = wire_cap;
          pin_metrics.pin_cap = pin_cap;
          pin_metrics.total_cap = wire_cap + pin_cap;
          pin_metrics.fanout = fanout;
          pin_metrics.arc_delay = arc_delay;
          pin_metrics.reachable_endpoints = reachable_endpoints;
          pin_metrics.is_in_clock_nets = is_in_clock_nets;
          pin_metrics.is_in_clock = is_in_clock;
          pin_metrics.is_port = is_port;
          pin_metrics.is_sequential = is_sequential;
          pin_metrics.is_macro = is_macro;
          pin_metrics.max_cap = max_cap;
          pin_metrics.max_slew = max_slew;
          pin_metrics.rise_slew = rise_slew;
          pin_metrics.fall_slew = fall_slew;
          pin_metrics.slack = slack;
          pin_metrics.rise_arrival_time = rise_arrival_time;
          pin_metrics.fall_arrival_time = fall_arrival_time;

          // Process the pin
          collector.processPin(pin_metrics);


          


          // Update previous values for next iteration
          prev_x = pin_loc.x();
          prev_y = pin_loc.y();
          prev_pin = pin;
      }
      path_count++;

      std::cout << "Debug Point 4" << std::endl;

      // Finalize the collection
      collector.finalize();





    }

    // Start processing the collected data

    std::cout << "Debug Point 5" << std::endl;

    const auto& sequences = collector.getSequences();

    // Need to convert the sequences into a format that can be used by the transformer model

    std::cout << "Debug Point 6" << std::endl;

    // Placeholder for string to id maps and id to embedding maps
    //const std::unordered_map<std::string, int> pin_name_to_id;
    //const std::unordered_map<std::string, int> cell_name_to_id;
    //const std::unordered_map<std::string, int> cell_type_to_id;
    //const std::unordered_map<int, std::vector<float>> cell_type_embeddings;


    std::cout << "Debug Point 7" << std::endl;

    
    
    auto builder = SequenceArrayBuilder(collector.getSequences(),
                                      pin_name_to_id,
                                      cell_name_to_id,
                                      libcell_to_type_id,
                                      libcell_to_type_embedding);

    auto [data_array, pin_ids, cell_ids, cell_type_ids] = builder.build();


    // Debugging print statements to check data_array and cell_type_ids

    // Print shape of the data_array, pin_ids, cell_ids, and cell_type_ids
    // N = number of sequences, L = max sequence length, D = number of features(token dimensions = pin_data + embedding)
    // data_array shape: (N, L, D)
    // pin_ids shape: (N, L): Pin IDs are used to lookup the pin name
    // cell_ids shape: (N, L): Cell IDs are used to lookup the cell name
    // cell_type_ids shape: (N, L): Cell Type IDs are used to lookup the cell type name
    std::cout << "Data Array Shape: (" << data_array.size() << ", " << data_array[0].size() << ", " << data_array[0][0].size() << ")" << std::endl;
    std::cout << "Pin IDs Shape: (" << pin_ids.size() << ", " << pin_ids[0].size() << ")" << std::endl;
    std::cout << "Cell IDs Shape: (" << cell_ids.size() << ", " << cell_ids[0].size() << ")" << std::endl;
    std::cout << "Cell Type IDs Shape: (" << cell_type_ids.size() << ", " << cell_type_ids[0].size() << ")" << std::endl;

    // Print example of the data_array (N, L, D), first 5 tokens of the first 2 sequence
    for (size_t o = 0; o < ((2 > data_array.size()) ? data_array.size() : 2); o++) {
      std::cout << "Sequence " << o << ": " << std::endl;
      for (size_t i = 0; i < ((5 > data_array[o].size()) ? data_array[o].size() : 5); i++) {
        std::cout << "Token " << i << ": ";
        for (size_t j = 0; j < data_array[0][0].size(); j++) {
          std::cout << data_array[o][i][j] << " ";
        }
        std::cout << std::endl;
      }
    }

    // Print cell_type_ids for the first 5 tokens of the first 2 sequences
    for (size_t o = 0; o < ((2 > data_array.size()) ? data_array.size() : 2); o++) {
      std::cout << "Sequence " << o << ": " << std::endl;
      for (size_t i = 0; i < ((5 > data_array[o].size()) ? data_array[o].size() : 5); i++) {
        std::cout << "Cell Type ID " << i << ": " << cell_type_ids[o][i] << std::endl;
      }
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
