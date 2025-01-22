#include "MLGateSizer.hh"
#include <cmath>  // for sqrt, pow
#include <cassert>
#include <iostream>
#include <random>


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

#include <Eigen/Dense>
#include <algorithm>
#include <fstream>

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


void MLGateSizer::getEndpointAndCriticalPaths()
{
  // Print out statement to indicate the function is running
  std::cout << "Retrieving endpoints and critical paths..." << std::endl;

  // Retrieve endpoints and critical paths for debugging or further analysis
  init();


  // Retrieve endpoints
  sta::VertexSet* endpoints = sta_->endpoints();

  // Print out the number of endpoints found
  std::cout << "Found " << endpoints->size() << " endpoints." << std::endl;


  //std::cout << "Debug Point 1" << std::endl;

  //std::cout << "Debug Point 2" << std::endl;

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
      2000, //10 * endpoints->size(), // group_count
      endpoints->size(),      // endpoint_count
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
  // In addition or alternatively,
  // consider using vertexWorstSlackPath to find the critical path for each
  // endpoint
  //std::cout << "Debug Point 3" << std::endl;

  // If no critical path is found, print a message
  if (path_ends.empty()) {
    std::cout << "No critical paths found " << std::endl;
  } else {

    int path_count = 0;
    // Declare tempoary vector to store the slack of each path
    std::vector<float> path_slacks;
    // Get all clk_nets to check if the pins are connected to a clock net
    std::set<dbNet*> clk_nets = sta_->findClkNets();
    // Used during data retrieval but stays uniform for all paths
    sta::LibertyLibrary* lib = network_->defaultLibertyLibrary();
    sta::Corner* corner = sta_->cmdCorner();

    // Get Library Cell to index mapping
    std::vector<std::string> all_libcell_names_temp; // Initialize them as empty and build them up, used for generating ordered_libcells_
    std::unordered_map<std::string, int> libcell_to_id_temp;
    std::unordered_map<std::string, int> libcell_to_type_id_temp; // equivalent libcells share the same, Uses EquivCells

    //std::unordered_map<int, std::vector<int>> libcell_type_id_to_libcell_id; // Maps libcell type id to libcell id
    //std::unordered_map<int, int> libcell_id_to_libcell_type_id; // Maps libcell id to libcell type id

    //std::unordered_map<int, std::vector<float>> libcell_to_embedding;
    //std::unordered_map<int, std::vector<float>> libcell_to_type_embedding;

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
        if (libcell_to_id_temp.find(cell_name) == libcell_to_id_temp.end()) {
          libcell_to_id_temp[cell_name] = libcell_id++;
          all_libcell_names_temp.push_back(cell_name);
        }
      }
    }

    // Use generateLibcellOrdering to generate ordered_libcells_
    generateLibcellOrdering(all_libcell_names_temp);

    // Assign the id and type id based on the ordered_libcells_
    // Convetion being that the first cell is assigned id 0
    // For type id, id increments from 0 for each unique cell type
    // The ordering is based on the smallest id of the libcells within the type
    // For example, if NAND2 libcell type has libcell id 0, 1, 2 and AND2 libcell type has libcell id 3, 4, 5
    // Then the type id for NAND2 is 0 and for AND2 is 1




    
    // Initialize libcell_to_type_id_temp to be used for updateLibcellTypeMap which updates internal data structures regarding libcell types 
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
          if (libcell_to_type_id_temp.find(equiv_cell_name) != libcell_to_type_id_temp.end()) {
            // Found an equivalent cell with an existing type ID, use it
            libcell_to_type_id_temp[cell_name] = libcell_to_type_id_temp[equiv_cell_name];
            type_id_found = true;
            break; // No need to check further
          }
        }

        // If none of the equivalent cells have a type ID, assign a new one
        if (!type_id_found) {
          libcell_to_type_id_temp[cell_name] = libcell_type_id++;
          // Optionally generate and assign embeddings here, e.g.,
          // cell_type_embeddings[libcell_to_type_id[cell_name]] = generateEmbedding(master);
        }

      }
    }

    // Initialize libcell_to_type_id_, libcell_type_id_to_libcell_ids_ and libcell_id_to_libcell_type_id_ based on libcell_to_id_ and libcell_to_type_id_temp
    updateLibcellTypeMap(libcell_to_type_id_temp);

    // Print the current tech_name
    //std::string tech_name = db->getTech()->getName();
    //std::cout << "Tech Name: " << tech_name << std::endl;
    // Print number of libcells and libcell types
    std::cout << "Total Libcells: " << libcell_to_id_.size() << std::endl;
    std::cout << "Total Libcell Types: " << libcell_type_id_to_libcell_ids_.size() << std::endl;

    std::cout << "All Libcell Names: " << std::endl;
    for (const auto& cell : ordered_libcells_) {
      std::cout << cell << std::endl;
    }

    std::cout << std::endl;

    // Print Type ID from 0, 1, 2, ... then print the libcell names and IDs associated with the type ID
    for (auto& type_pair : libcell_type_id_to_libcell_ids_) {
      int type_id = type_pair.first;
      std::cout << "Type ID: " << type_id << std::endl;
      for (int cell_id : type_pair.second) {
        std::string cell_name = ordered_libcells_[cell_id];
        std::cout << "Cell Name: " << cell_name << " ID: " << libcell_to_id_[cell_name] << std::endl;
      }
    }

    std::cout << std::endl;

    // Print Cell Names and ID from libcell_to_id_ and Type ID
    for (auto& cell : libcell_to_id_) {
      std::cout << "Cell Name: " << cell.first << " ID: " << cell.second << " Type ID: " << libcell_to_type_id_[cell.first] << std::endl;
    }

    // Attempt to load libcell embeddings from a binary file
    size_t embedding_size = 768; // Embedding size of deberta-v3-base model

    // Load the embeddings based on the number of libcells
    // Currently supports ASAP7 and Nangate45 temporary fix needs to be more robust and general
    if (libcell_to_id_.size() == 216) { // ASAP7
      loadEmbeddingsBinary("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/embedding_generation/ASAP7_libcell_embeddings.bin", embedding_size);
      updateLibcellTypeEmbeddings();
  }
    else if (libcell_to_id_.size() == 135) {  // Nangate45
      loadEmbeddingsBinary("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/embedding_generation/nangate45_libcell_embeddings.bin", embedding_size);
      updateLibcellTypeEmbeddings();
    }
    else {
      std::cout << "Unsupported, unable to load embeddings" << std::endl;

      std::cout << "Generate placeholder embeddings" << std::endl;
      libcell_id_to_embedding_.clear();
      embedding_size_ = embedding_size;
      // Initialize libcell_to_embeddings with placeholder values
      // Proper embedding comes from passing through libcell name through a sentence transformer and averaging the values for same libcell type
      for (const auto& type_pair : libcell_to_id_) {
        int libcell_id = type_pair.second;

        // Create a 16-float embedding filled with 0.0, 1.0, 2.0, ...
        std::vector<float> embedding(embedding_size);
        for (int i = 0; i < embedding.size(); ++i) {
          embedding[i] = static_cast<float>(libcell_id); // Assign the libcell_id as the value for all elements in the embedding
        }

        libcell_id_to_embedding_[libcell_id] = embedding;
      }

      updateLibcellTypeEmbeddings();

    }

    // Print the embeddings for each libcell id first 5 elements
    for (int i = 0; i < 5; i++) {
      std::cout << "Libcell ID: " << i << " Libcell Name: " << ordered_libcells_[i]<< " Embedding: ";
      for (int j = 0; j < 5; j++) {
        std::cout << libcell_id_to_embedding_[i][j] << " ";
      }
      std::cout << std::endl;
    }


    // Normalization constants for PinMetrics
    float pin_loc_x_max = 0.0;
    float pin_loc_x_min = 0.0;
    float pin_loc_x_mean = 0.0;
    float pin_loc_x_std = 0.0;
    float pin_loc_y_max = 0.0;
    float pin_loc_y_min = 0.0;
    float pin_loc_y_mean = 0.0;
    float pin_loc_y_std = 0.0;
    float p2p_dist_max = 0.0;
    float p2p_dist_min = 0.0;
    float p2p_dist_mean = 0.0;
    float p2p_dist_std = 0.0;
    float hpwl_max = 0.0;
    float hpwl_min = 0.0;
    float hpwl_mean = 0.0;
    float hpwl_std = 0.0;
    float input_pin_cap_max = 0.0;
    float input_pin_cap_min = 0.0;
    float input_pin_cap_mean = 0.0;
    float input_pin_cap_std = 0.0;
    float wire_cap_max = 0.0;
    float wire_cap_min = 0.0;
    float wire_cap_mean = 0.0;
    float wire_cap_std = 0.0;
    float pin_cap_max = 0.0;
    float pin_cap_min = 0.0;
    float pin_cap_mean = 0.0;
    float pin_cap_std = 0.0;
    float total_cap_max = 0.0;
    float total_cap_min = 0.0;
    float total_cap_mean = 0.0;
    float total_cap_std = 0.0;
    float fanout_max = 0.0;
    float fanout_min = 0.0;
    float fanout_mean = 0.0;
    float fanout_std = 0.0;
    float arc_delay_max = 0.0;
    float arc_delay_min = 0.0;
    float arc_delay_mean = 0.0;
    float arc_delay_std = 0.0;
    float reachable_endpoints_max = 0.0;
    float reachable_endpoints_min = 0.0;
    float reachable_endpoints_mean = 0.0;
    float reachable_endpoints_std = 0.0;
    float rise_slew_max = 0.0;
    float rise_slew_min = 0.0;
    float rise_slew_mean = 0.0;
    float rise_slew_std = 0.0;
    float fall_slew_max = 0.0;
    float fall_slew_min = 0.0;
    float fall_slew_mean = 0.0;
    float fall_slew_std = 0.0;
    float slack_max = 0.0;
    float slack_min = 0.0;
    float slack_mean = 0.0;
    float slack_std = 0.0;
    float rise_arrival_time_max = 0.0;
    float rise_arrival_time_min = 0.0;
    float rise_arrival_time_mean = 0.0;
    float rise_arrival_time_std = 0.0;
    float fall_arrival_time_max = 0.0;
    float fall_arrival_time_min = 0.0;
    float fall_arrival_time_mean = 0.0;
    float fall_arrival_time_std = 0.0;




  


    

    
    PinSequenceCollector collector;

    for (auto& path_end : path_ends) {  // similar usage found in TritonPart.cpp
                                        // BuildTimingPaths()
      //std::cout << "Critical Path " << path_count << std::endl;
      auto* path = path_end->path();
      float slack = path_end->slack(sta_);
      path_slacks.push_back(slack);
      //std::cout << "Slack: " << slack << std::endl;
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


          // Cell type (gate type), retrieve the cell from the pin, then retrieve the cell type
          std::string cell_type = is_port ? "Port" : network_->libertyCell(inst)->name();
          std::string cell_name = is_port ? "Port" : network_->name(inst);



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

          /*
          // Output all collected data (for debugging)
          float min_slack = sta_->pinSlack(pin, sta::MinMax::min());
          // Network name of the pin
          std::string net_name = net ? network_->name(net) : "None";
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
          */
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

      //std::cout << "Debug Point 4" << std::endl;

      // Finalize the collection
      collector.finalize();





    }

    // Start processing the collected data

    //std::cout << "Debug Point 5" << std::endl;

    // const auto& sequences = collector.getSequences();

    // Need to convert the sequences into a format that can be used by the transformer model

    //std::cout << "Debug Point 6" << std::endl;


    //std::cout << "Debug Point 7" << std::endl;

    
    
    auto builder = SequenceArrayBuilder(collector.getSequences(),
                                      pin_name_to_id,
                                      cell_name_to_id,
                                      libcell_to_type_id_,
                                      libcell_type_id_to_embedding_);

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



    
    /*
    // Print out the slack of each path (for debugging), remove later
    for (size_t i = 0; i < path_slacks.size(); i++) {
      std::cout << "Path " << i << " Slack: " << path_slacks[i] << std::endl;
    }
    */

    int num_heads = 8;//4;
    size_t N = data_array.size(); // Number of sequences (100)
    size_t L = data_array[0].size();  // Max sequence length (78)
    size_t D_in = data_array[0][0].size();  // Input feature dimensions (17)
    size_t D_model = 128;//64; // Transformer model hidden dimensions (64)
    size_t FF_hidden_dim = 4 * D_model; // Feedforward hidden dimensions (256)
    int num_encoder_layers = 6;//2;

    // Calculate the number of parameters in the model
    // Initial Projection weights: (D_in x D_model)
    // Transformer Encoder layers: num_encoder_layers
    // - Q, K, V, O weights: (D_model x D_model) x 4
    // - FF weights: (D_model x FF_hidden_dim) x 2 + FF_hidden_dim + D_model
    // Final Output Projection weights: (D_model x D_in)
    size_t model_params = ((D_model * D_model) * 4 + (D_model * FF_hidden_dim) * 2 + FF_hidden_dim + D_model) * num_encoder_layers + (D_in * D_model) * 2;

    auto outputs = runTransformer(data_array, num_heads, N, L, D_in, D_model, FF_hidden_dim, num_encoder_layers);  // Run the transformer model
    auto eigen_outputs = runTransformerEigen(data_array, num_heads, N, L, D_in, D_model, FF_hidden_dim, num_encoder_layers);  // Run the transformer model using Eigen

    // Debugging print statements to check outputs
    /*
    // 1) Benchmark NAIVE
    long long naive_us = benchmark([&](){
      auto out_naive = runTransformer(data_array, num_heads, N, L, D_in, D_model, FF_hidden_dim, num_encoder_layers);
    });
    // 2) Benchmark EIGEN
    long long eigen_us = benchmark([&](){
      auto out_eigen = runTransformerEigen(data_array, num_heads, N, L, D_in, D_model, FF_hidden_dim, num_encoder_layers);
    });
    */
    // 3) Actually store the outputs to compare correctness
    auto start = std::chrono::steady_clock::now();
    auto out_naive  = runTransformer(data_array, num_heads, N, L, D_in, D_model, FF_hidden_dim, num_encoder_layers);
    auto end = std::chrono::steady_clock::now();
    auto naive_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    start = std::chrono::steady_clock::now();
    auto out_eigen  = runTransformerEigen(data_array, num_heads, N, L, D_in, D_model, FF_hidden_dim, num_encoder_layers);
    end = std::chrono::steady_clock::now();
    auto eigen_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    bool ok = compareOutputs(out_naive, out_eigen, 1e-3f);

    // 4) Print the speed & correctness info
    // total tokens processed: N*L
    size_t total_tokens = data_array.size() * data_array[0].size();

    std::cout << "Model parameters: " << model_params << std::endl;

    std::cout << "Total tokens: " << total_tokens << ": " << data_array.size() << "X" << data_array[0].size() << std::endl;

    std::cout << "NAIVE  time: " << naive_us << " us   => "
              << (1e6 * double(total_tokens) / double(naive_us))
              << " tokens/sec\n";
    std::cout << "EIGEN  time: " << eigen_us << " us   => "
              << (1e6 * double(total_tokens) / double(eigen_us))
              << " tokens/sec\n";
    std::cout << "Outputs match? " << (ok ? "YES" : "NO") << std::endl;

    // If the outputs don't match, print the output first 5 tokens of the first 2 sequences or before the data ends
    // Print both NAIVE and EIGEN outputs for easy comparison and debugging
    if (!ok) {
      std::cout << "NAIVE:" << std::endl;
      for (size_t o = 0; o < ((2 > out_naive.size()) ? out_naive.size() : 2); o++) {
        std::cout << "Sequence " << o << ": " << std::endl;
        for (size_t i = 0; i < ((5 > out_naive[o].size()) ? out_naive[o].size() : 5); i++) {
          std::cout << "Token " << i << ": ";
          for (size_t j = 0; j < out_naive[0][0].size(); j++) {
            std::cout << out_naive[o][i][j] << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "EIGEN:" << std::endl;
      for (size_t o = 0; o < ((2 > out_eigen.size()) ? out_eigen.size() : 2); o++) {
        std::cout << "Sequence " << o << ": " << std::endl;
        for (size_t i = 0; i < ((5 > out_eigen[o].size()) ? out_eigen[o].size() : 5); i++) {
          std::cout << "Token " << i << ": ";
          for (size_t j = 0; j < out_eigen[0][0].size(); j++) {
            std::cout << out_eigen[o][i][j] << " ";
          }
          std::cout << std::endl;
        }
      }
    }



  }



  // Next step is extracting the data from each pin
  // Then, transform the data into the format expected by the transformer model
  // Finally, apply the transformer model to classify and resize the gates

  // Use PinMetrics and PinDataSequence to store the extracted data

}


void MLGateSizer::writeBinaryFile(const std::string& filename, 
                                 const std::vector<std::vector<std::vector<float>>>& data)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        logger_->error(utl::RSZ, 1001, "Cannot open file {} for writing", filename);
        return;
    }

    // Write dimensions
    size_t N = data.size();
    size_t L = N > 0 ? data[0].size() : 0;
    size_t D = L > 0 ? data[0][0].size() : 0;

    out.write(reinterpret_cast<const char*>(&N), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&L), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&D), sizeof(size_t));

    // Write data in flattened format
    for (const auto& seq : data) {
        for (const auto& token : seq) {
            out.write(reinterpret_cast<const char*>(token.data()), D * sizeof(float));
        }
    }

    if (!out) {
        logger_->error(utl::RSZ, 1002, "Error writing to file {}", filename);
    }
}

std::vector<std::vector<std::vector<float>>> MLGateSizer::readBinaryFile(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        logger_->error(utl::RSZ, 1003, "Cannot open file {} for reading", filename);
        return {};
    }

    // Read dimensions
    size_t N, L, D;
    in.read(reinterpret_cast<char*>(&N), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&L), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&D), sizeof(size_t));

    if (!in) {
        logger_->error(utl::RSZ, 1004, "Error reading dimensions from file {}", filename);
        return {};
    }

    // Initialize 3D vector
    std::vector<std::vector<std::vector<float>>> data(N, 
        std::vector<std::vector<float>>(L, std::vector<float>(D)));

    // Read data
    for (auto& seq : data) {
        for (auto& token : seq) {
            in.read(reinterpret_cast<char*>(token.data()), D * sizeof(float));
            if (!in) {
                logger_->error(utl::RSZ, 1005, "Error reading data from file {}", filename);
                return {};
            }
        }
    }

    return data;
}

void MLGateSizer::generateLibcellOrdering(const std::vector<std::string>& libcells)
{
    ordered_libcells_ = libcells;
    std::sort(ordered_libcells_.begin(), ordered_libcells_.end(), LibcellComparator());
    
    // Rebuild libcell_to_id map
    libcell_to_id_.clear();
    for (size_t i = 0; i < ordered_libcells_.size(); ++i) {
        libcell_to_id_[ordered_libcells_[i]] = i;
    }
}

 // update libcell_to_type_id_, libcell_type_id_to_libcell_ids_, and libcell_id_to_libcell_type_id_
void MLGateSizer::updateLibcellTypeMap(const std::unordered_map<std::string, int>& libcell_to_type_id)
{
    libcell_to_type_id_ = libcell_to_type_id;
    libcell_type_id_to_libcell_ids_.clear();
    libcell_id_to_libcell_type_id_.clear();

	for (const auto& libcell : libcell_to_type_id_) {
		std::string libcell_name = libcell.first;
		int libcell_id = libcell.second;
		int libcell_type_id = libcell_to_type_id_[libcell_name];
		libcell_type_id_to_libcell_ids_[libcell_type_id].push_back(libcell_id);
		libcell_id_to_libcell_type_id_[libcell_id] = libcell_type_id;
	}

}


void MLGateSizer::saveEmbeddingsBinary(const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        logger_->error(utl::RSZ, 3001, "Cannot open embeddings file {}", filename);
        return;
    }

    // Write embedding dimensions first
    const size_t num_items = ordered_libcells_.size();
    out.write(reinterpret_cast<const char*>(&num_items), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&embedding_size_), sizeof(size_t));

    // Write embeddings in alphabetical order
    for (const auto& libcell : ordered_libcells_) {
        const auto& emb = libcell_id_to_embedding_[libcell_to_id_[libcell]];
        out.write(reinterpret_cast<const char*>(emb.data()), emb.size() * sizeof(float));
    }
}

void MLGateSizer::loadEmbeddingsBinary(const std::string& filename, size_t embedding_size)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        logger_->error(utl::RSZ, 3002, "Cannot open embeddings file {}", filename);
        return;
    }

    // Read dimensions
    size_t num_items, file_embedding_size;
    in.read(reinterpret_cast<char*>(&num_items), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&file_embedding_size), sizeof(size_t));
    
    if (num_items != ordered_libcells_.size() || file_embedding_size != embedding_size) {
        logger_->error(utl::RSZ, 3003, 
            "Embedding file dimensions mismatch. Expected {} items of size {}, got {} items of size {}",
            ordered_libcells_.size(), embedding_size, num_items, file_embedding_size);
        return;
    }

    // Read embeddings in alphabetical order
    libcell_id_to_embedding_.clear();
    embedding_size_ = embedding_size;
    std::vector<float> embedding(embedding_size);
    
    for (size_t i = 0; i < num_items; ++i) {
        in.read(reinterpret_cast<char*>(embedding.data()), embedding_size * sizeof(float));
        if (!in) {
            logger_->error(utl::RSZ, 3004, "Error reading embedding {} from file", i);
            return;
        }
        libcell_id_to_embedding_[i] = embedding;
    }
}

// Generate libcell type embeddings from libcell embeddings by averaging over all libcells of the same type
void MLGateSizer::updateLibcellTypeEmbeddings()
{
	libcell_type_id_to_embedding_.clear();
	for (const auto& libcell_type : libcell_type_id_to_libcell_ids_) {
		int libcell_type_id = libcell_type.first;
		const auto& libcell_ids = libcell_type.second;
		std::vector<float> embedding(embedding_size_, 0.0f);
		for (const auto& libcell_id : libcell_ids) {
			const auto& libcell_emb = libcell_id_to_embedding_[libcell_id];
			for (size_t i = 0; i < embedding_size_; ++i) {
				embedding[i] += libcell_emb[i];
			}
		}
		for (size_t i = 0; i < embedding_size_; ++i) {
			embedding[i] /= libcell_ids.size();
		}
		libcell_type_id_to_embedding_[libcell_type_id] = embedding;
	}
}

// ------------------------------------------------------------
// Helper function: naive matrix multiply
// A: [M x K], B: [K x N] => C: [M x N]
// ------------------------------------------------------------
static std::vector<std::vector<float>> matMul(
    const std::vector<std::vector<float>>& A,
    const std::vector<std::vector<float>>& B)
{
  const size_t M = A.size();
  const size_t K = A[0].size();
  //const size_t K2 = B.size(); // should match K
  const size_t N = B[0].size();

  assert(K == B.size() && "Inner dimensions must match for matMul");

  std::vector<std::vector<float>> C(M, std::vector<float>(N, 0.0f));
  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < N; j++) {
      float sum = 0.0f;
      for (size_t k = 0; k < K; k++) {
        sum += A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
  return C;
}

// ------------------------------------------------------------
// Helper function: naive LayerNorm across "D" dimension
// For each token T, we have a vector of dimension D.
// We'll do LN(token) = (token - mean)/sqrt(var + eps) * gamma + beta
// This example uses no gamma/beta for simplicity, or sets them to 1/0.
// ------------------------------------------------------------
static void layerNorm(std::vector<std::vector<float>>& seq, float eps = 1e-5f)
{
  // seq is shape [L x D]
  // For each row (token) in seq, compute mean & var, then normalize in place
  const size_t L = seq.size();
  if (L == 0) return;
  const size_t D = seq[0].size();

  for (size_t i = 0; i < L; i++) {
    // 1) mean
    float mean = 0.0f;
    for (size_t d = 0; d < D; d++) {
      mean += seq[i][d];
    }
    mean /= static_cast<float>(D);

    // 2) variance
    float var = 0.0f;
    for (size_t d = 0; d < D; d++) {
      float diff = seq[i][d] - mean;
      var += diff * diff;
    }
    var /= static_cast<float>(D);

    // 3) normalize in place
    for (size_t d = 0; d < D; d++) {
      seq[i][d] = (seq[i][d] - mean) / std::sqrt(var + eps);
      // If you want to add gamma/beta:
      // seq[i][d] = gamma[d] * seq[i][d] + beta[d];
    }
  }
}

// ------------------------------------------------------------
// Helper function: naive multi-head self-attention for one sequence
// We'll do a single pass of multi-head attention for a shape [L x D].
//
// Heads: H
// We split D into H heads each of size Dh = D/H.
// For simplicity, we do random Wq, Wk, Wv. Real code would load from TransSizer.
// Model loading has not been implemented yet.
// We do the standard formula: softmax(Q*K^T / sqrt(Dh)) * V
// Then we re-concatenate the heads.
//
// This is extremely naive, with no parallelization, etc.
// Try to use this as a reference for more optimized code.
// ------------------------------------------------------------
static std::vector<std::vector<float>> multiHeadSelfAttention(
    const std::vector<std::vector<float>>& seq,
    int num_heads,
    size_t D,
    size_t L,
    std::vector<std::vector<float>>& Wq,
    std::vector<std::vector<float>>& Wk,
    std::vector<std::vector<float>>& Wv,
    std::vector<std::vector<float>>& Wo)
{
  // seq shape: [L x D]
  // Return shape: [L x D]
  //const size_t L = seq.size();
  if (L == 0) {
    return seq; // nothing to do
  }
  //const size_t D = seq[0].size();
  assert(D % num_heads == 0 && "D must be divisible by num_heads");
  size_t Dh = D / num_heads;

  // 1) Create random (or dummy) Wq, Wk, Wv of shape [D x D].
  // For demonstration, let's fill them with small random values.
  static std::mt19937 rng(42);
  std::uniform_real_distribution<float> dist(-0.1f, 0.1f);

  auto randomMatrix = [&](size_t rows, size_t cols) {
    std::vector<std::vector<float>> mat(rows, std::vector<float>(cols, 0.0f));
    for (size_t r = 0; r < rows; r++) {
      for (size_t c = 0; c < cols; c++) {
        mat[r][c] = dist(rng);
      }
    }
    return mat;
  };

  //std::vector<std::vector<float>> Wq = randomMatrix(D, D);
  //std::vector<std::vector<float>> Wk = randomMatrix(D, D);
  //std::vector<std::vector<float>> Wv = randomMatrix(D, D);

  // 2) Flatten seq from [L x D] to a form we can do matMul with:
  // For matMul, we treat it as [L x D].
  // Q = seq * Wq => shape [L x D]
  // K = seq * Wk => shape [L x D]
  // V = seq * Wv => shape [L x D]

  // We'll define a small helper to do matMul (LxD) * (DxD) => (LxD)
  auto Q = matMul(seq, Wq);
  auto K = matMul(seq, Wk);
  auto V = matMul(seq, Wv);

  // 3) Reshape Q,K,V into [L x H x Dh], do attention per head
  // We'll store them in 3D structures: shape: Q_3d[H][L][Dh]
  auto to3D = [&](const std::vector<std::vector<float>>& X) {
    // X is [L x D]; we want [H x L x Dh]
    std::vector<std::vector<std::vector<float>>> X3(
        num_heads, std::vector<std::vector<float>>(L, std::vector<float>(Dh, 0.0f)));
    for (size_t l = 0; l < L; l++) {
      for (size_t h = 0; h < (size_t)num_heads; h++) {
        for (size_t d = 0; d < Dh; d++) {
          X3[h][l][d] = X[l][h * Dh + d];
        }
      }
    }
    return X3;
  };

  auto Q3 = to3D(Q);
  auto K3 = to3D(K);
  auto V3 = to3D(V);

  // 4) Self-attention per head
  // outHead[h] shape [L x Dh]
  auto softmax = [&](std::vector<float>& logits) {
    float max_val = logits[0];
    for (auto v : logits) max_val = std::max(max_val, v);
    float sum = 0.0f;
    for (auto& v : logits) {
      v = std::exp(v - max_val);
      sum += v;
    }
    for (auto& v : logits) {
      v /= sum;
    }
  };

  std::vector<std::vector<std::vector<float>>> outHeads(
      num_heads, std::vector<std::vector<float>>(L, std::vector<float>(Dh, 0.0f)));

  for (size_t h = 0; h < (size_t)num_heads; h++) {
    // For each head, we do attention
    // For each query position lq in [0..L-1]:
    for (size_t lq = 0; lq < L; lq++) {
      // 1) Q3[h][lq]: shape [Dh]
      // compute attention logits vs all K positions => shape [L]
      std::vector<float> attn_logits(L, 0.0f);
      for (size_t lk = 0; lk < L; lk++) {
        // dot(Q3[h][lq], K3[h][lk]) / sqrt(Dh)
        float dot_val = 0.0f;
        for (size_t d = 0; d < Dh; d++) {
          dot_val += Q3[h][lq][d] * K3[h][lk][d];
        }
        attn_logits[lk] = dot_val / std::sqrt((float)Dh);
      }
      // 2) softmax over attn_logits
      softmax(attn_logits);

      // 3) Weighted sum of V3[h][lk]
      // outHeads[h][lq] = sum over lk( attn_logits[lk]*V3[h][lk] )
      for (size_t lk = 0; lk < L; lk++) {
        for (size_t d = 0; d < Dh; d++) {
          outHeads[h][lq][d] += attn_logits[lk] * V3[h][lk][d];
        }
      }
    }
  }

  // 5) Concatenate heads back => shape [L x D]
  std::vector<std::vector<float>> outSeq(L, std::vector<float>(D, 0.0f));
  for (size_t l_ = 0; l_ < L; l_++) {
    for (size_t h_ = 0; h_ < (size_t)num_heads; h_++) {
      for (size_t d_ = 0; d_ < Dh; d_++) {
        outSeq[l_][h_ * Dh + d_] = outHeads[h_][l_][d_];
      }
    }
  }

  // 6) Output projection
  // We'll define Wo: [D x D] for simplicity
  //std::vector<std::vector<float>> Wo = randomMatrix(D, D);

  // outSeq * Wo: [L x D] * [D x D] => [L x D]
  auto finalOut = matMul(outSeq, Wo);

  return finalOut;
}

// ------------------------------------------------------------
// Helper function: naive feed-forward network
// We'll do something like: FF(x) = ReLU(x*W1 + b1)*W2 + b2
// For demonstration, we just define random W1, b1, W2, b2
// ------------------------------------------------------------
static std::vector<std::vector<float>> feedForward(
    const std::vector<std::vector<float>>& seq,
    size_t L,
    size_t D,
    size_t H,
    std::vector<std::vector<float>>& W1,
    std::vector<std::vector<float>>& W2,
    std::vector<float>& b1,
    std::vector<float>& b2)
{
  // seq shape: [L x D]
  // Let hidden dim = 2D for example
  if (seq.empty()) return seq;

  //const size_t L = seq.size();
  //const size_t D = seq[0].size();
  //const size_t H = 2 * D; // arbitrary

  static std::mt19937 rng(123);
  std::uniform_real_distribution<float> dist(-0.1f, 0.1f);

  // We define W1: [D x H], b1: [H], W2: [H x D], b2: [D]
  auto randomMatrix = [&](size_t rows, size_t cols) {
    std::vector<std::vector<float>> mat(rows, std::vector<float>(cols, 0.0f));
    for (size_t r = 0; r < rows; r++) {
      for (size_t c = 0; c < cols; c++) {
        mat[r][c] = dist(rng);
      }
    }
    return mat;
  };

  auto randomVector = [&](size_t size) {
    std::vector<float> vec(size, 0.0f);
    for (size_t i = 0; i < size; i++) {
      vec[i] = dist(rng);
    }
    return vec;
  };

  //std::vector<std::vector<float>> W1 = randomMatrix(D, H);
  //std::vector<std::vector<float>> W2 = randomMatrix(H, D);
  //std::vector<float> b1 = randomVector(H);
  //std::vector<float> b2 = randomVector(D);

  // 1) Compute seq * W1 => shape [L x H]
  auto tmp = matMul(seq, W1);
  // add b1 and ReLU
  for (size_t l = 0; l < L; l++) {
    for (size_t h = 0; h < H; h++) {
      tmp[l][h] += b1[h];
      if (tmp[l][h] < 0.0f) {
        tmp[l][h] = 0.0f; // ReLU
      }
    }
  }
  // 2) Multiply by W2 => [L x D]
  auto out = matMul(tmp, W2);
  // add b2
  for (size_t l = 0; l < L; l++) {
    for (size_t d = 0; d < D; d++) {
      out[l][d] += b2[d];
    }
  }
  return out;
}


// --------------------------------------------------------------------
// runTransformer: Projects data_array from D_in -> D_model (divisible by H),
// runs M encoder layers, then projects back to D_in if desired.
//
// data_array: shape [N x L x D_in]
// returns shape [N x L x D_in] if we do the final projection back
//
// Main "encoder" demonstration (one or more layers).
// We'll do 2 encoder layers as a initial test example:
//   For each layer:
//     Y = LN( X + MultiHeadAttn(X) )
//     Z = LN( Y + FF(Y) )
// Return Z
//
// FF is a simple 2-layer feed-forward network: ReLU(x*W1 + b1)*W2 + b2
// - x: [L x D_model]
// - W1: [D_model x FF_hidden_dim]
// - W2: [FF_hidden_dim x D_model]
// - b1: [FF_hidden_dim]
// - b2: [D_model]
//   Returns: [L x D_model]
//
// Currently uses output projection where Z : [L x D_model] is projected back to [L x D_in]
// Z : [L x D_model] * W_out [D_model x D_in]=> [L x D_in]
// For Gate sizing need to take encoder output and pass it through a classification head.
// Each class corresponding to the libcell.
// --------------------------------------------------------------------
std::vector<std::vector<std::vector<float>>> MLGateSizer::runTransformer(
    const std::vector<std::vector<std::vector<float>>>& data_array,
    int num_heads,
    size_t N,
    size_t L,
    size_t D_in,
    size_t D_model,
    size_t FF_hidden_dim,
    int num_encoder_layers)
{
  //size_t N = data_array.size();
  if (N == 0) return {};

  //size_t L = data_array[0].size();
  if (L == 0) return {};

  // Original input dimension
  //size_t D_in = data_array[0][0].size();

  // Decide how many heads we want
  //const int num_heads = 4;

  // Decide a model dimension that is divisible by num_heads
  //size_t D_model = 64;

  // Number of encoder layers
  //const int num_encoder_layers = 2;

  // We will random-initialize projection weights:
  static std::mt19937 rng(999);
  std::uniform_real_distribution<float> dist(-0.05f, 0.05f);

  // W_in: [D_in x D_model]
  auto randomMatrix = [&](size_t rows, size_t cols) {
    std::vector<std::vector<float>> mat(rows, std::vector<float>(cols, 0.f));
    for (size_t r = 0; r < rows; r++) {
      for (size_t c = 0; c < cols; c++) {
        mat[r][c] = dist(rng);
      }
    }
    return mat;
  };
  auto randomVector = [&](size_t size) {
    std::vector<float> vec(size, 0.0f);
    for (size_t i = 0; i < size; i++) {
      vec[i] = dist(rng);
    }
    return vec;
  };

  // If you want to map back to D_in:
  // W_out: [D_model x D_in]
  auto W_in  = randomMatrix(D_in, D_model);
  auto W_out = randomMatrix(D_model, D_in);

  // Allocate Q, K, V, O projection matrices
  auto Wq = randomMatrix(D_model, D_model);
  auto Wk = randomMatrix(D_model, D_model);
  auto Wv = randomMatrix(D_model, D_model);
  auto Wo = randomMatrix(D_model, D_model);

  // Allocate FF weights
  // FF_W1: [D_model x FF_hidden_dim]
  // FF_b1: [FF_hidden_dim]
  // FF_W2: [FF_hidden_dim x D_model]
  // FF_b2: [D_model]
  auto FF_W1 = randomMatrix(D_model, FF_hidden_dim);
  auto FF_W2 = randomMatrix(FF_hidden_dim, D_model);
  auto FF_b1 = randomVector(FF_hidden_dim);
  auto FF_b2 = randomVector(D_model);

  // Allocate output
  std::vector<std::vector<std::vector<float>>> output(N,
      std::vector<std::vector<float>>(L, std::vector<float>(D_in, 0.f)));

  // For each sequence in the batch:
  for (size_t n = 0; n < N; n++) {
    // current seq: [L x D_in]
    std::vector<std::vector<float>> seq = data_array[n];

    // 1) Project seq => [L x D_model]
    auto seq_proj = matMul(seq, W_in);

    // 2) Run through encoder blocks
    std::vector<std::vector<float>> x = seq_proj; // shape [L x D_model]
    for (int layer = 0; layer < num_encoder_layers; layer++) {
      // (a) Multi-head self-attention
      auto attn_out = multiHeadSelfAttention(x, num_heads, D_model, L, Wq, Wk, Wv, Wo);
      // (b) Residual
      for (size_t l_ = 0; l_ < L; l_++) {
        for (size_t d_ = 0; d_ < D_model; d_++) {
          attn_out[l_][d_] += x[l_][d_];
        }
      }
      // (c) LayerNorm
      layerNorm(attn_out);

      // (d) FeedForward
      auto ff_out = feedForward(attn_out, L, D_model, FF_hidden_dim, FF_W1, FF_W2, FF_b1, FF_b2);
      // (e) Residual
      for (size_t l_ = 0; l_ < L; l_++) {
        for (size_t d_ = 0; d_ < D_model; d_++) {
          ff_out[l_][d_] += attn_out[l_][d_];
        }
      }
      // (f) LayerNorm
      layerNorm(ff_out);

      // next layer input
      x = ff_out;
    }

    // x is now [L x D_model]. Project back to [L x D_in]
    auto final_out = matMul(x, W_out);

    // Store final_out in output[n]
    output[n] = final_out; // shape [L x D_in]
  }

  // Return shape: [N x L x D_in]
  return output;
}

static Eigen::MatrixXf eigenSelfAttention(
  const Eigen::MatrixXf& seq, 
  int num_heads,
  size_t D,
  size_t L,
  Eigen::MatrixXf& Wq,
  Eigen::MatrixXf& Wk,
  Eigen::MatrixXf& Wv,
  Eigen::MatrixXf& Wo)
{
  // seq shape: [L x D]
  // Return shape: [L x D]
  //const int L = seq.rows();
  //const int D = seq.cols();

  assert(D % num_heads == 0 && "D must be divisible by num_heads");
  size_t Dh = D / num_heads;


  // Random Wq, Wk, Wv of shape [D x D].
  static std::mt19937 rng(42);
  std::uniform_real_distribution<float> dist(-0.1f, 0.1f);

  auto randomMatrix = [&](int r, int c){
    Eigen::MatrixXf m(r,c);
    for(int i=0; i<r; i++){
      for(int j=0; j<c; j++){
        m(i,j) = dist(rng);
      }
    }
    return m;
  };

  // Wq,Wk,Wv => [D x D]
  //Eigen::MatrixXf Wq = randomMatrix(D,D);
  //Eigen::MatrixXf Wk = randomMatrix(D,D);
  //Eigen::MatrixXf Wv = randomMatrix(D,D);
  // Output projection
  // We'll define Wo: [D x D] for simplicity
  //Eigen::MatrixXf Wo = randomMatrix(D,D);

  // Q,K,V => [L x D] * [D x D] => [L x D]
  Eigen::MatrixXf Q = seq * Wq;
  Eigen::MatrixXf K = seq * Wk;
  Eigen::MatrixXf V = seq * Wv;

  // attention output => L x D
  Eigen::MatrixXf out(L, D);

  // Store Q, K, V for each head in separate 2D matrices.
  // Qsplit[h] shape [L x Dh], vice versa for Ksplit, Vsplit
  std::vector<Eigen::MatrixXf> Qsplit(num_heads), Ksplit(num_heads), Vsplit(num_heads);
  for (int h = 0; h < num_heads; h++)
  {
    Qsplit[h].resize(L, Dh);
    Ksplit[h].resize(L, Dh);
    Vsplit[h].resize(L, Dh);

    // Fill each sub-matrix by slicing out the corresponding [Dh] chunk
    for (int l_ = 0; l_ < L; l_++) {
      for (int d_ = 0; d_ < Dh; d_++) {
        Qsplit[h](l_,d_) = Q(l_,h*Dh+d_);
        Ksplit[h](l_,d_) = K(l_,h*Dh+d_);
        Vsplit[h](l_,d_) = V(l_,h*Dh+d_);
      }
    }
  }

  // Accumulate per-head outputs (each [L x Dh]) into outHeads
  std::vector<Eigen::MatrixXf> outHeads(num_heads, Eigen::MatrixXf::Zero(L, Dh));

  // Self-attention per head
  // outHead[h] shape [L x Dh]

  // softmax function
  auto softmax = [&](Eigen::VectorXf& logits){
    float max_val = logits.maxCoeff();
    logits = (logits.array() - max_val).exp();
    float sum_ = logits.sum();
    logits /= sum_;
    return logits;
  };


  // Compute self-attention for each head
  for (size_t h = 0; h < (size_t) num_heads; h++) {
    // Qsplit[h], Ksplit[h], Vsplit[h] => [L x Dh]
    // For each head, we do attention
    // For each query position lq in [0..L-1]:
    for (int lq = 0; lq < L; lq++) {
      // 1) Q3[h][lq]: shape [Dh]
      // compute attention logits vs all K positions => shape [L]
      Eigen::VectorXf attn_logits(L);
      for (int lk = 0; lk < L; lk++) {
        // dot(Q3[h][lq], K3[h][lk]) / sqrt(Dh)
        float dot_val = Qsplit[h].row(lq).dot(Ksplit[h].row(lk));
        attn_logits(lk) = dot_val / std::sqrt(float(Dh));
      }
      // 2) softmax over attn_logits
      softmax(attn_logits);

      // 3) Weighted sum of V3[h][lk]
      // outHeads[h][lq] = sum over lk( attn_logits[lk]*V3[h][lk] )
      for (int lk = 0; lk < L; lk++) {
        outHeads[h].row(lq).noalias() += attn_logits(lk) * Vsplit[h].row(lk);
      }
    }
  }

  // Concatenate heads back => shape [L x D]
  for(int l_ = 0; l_< L; l_++){
    for(int h = 0; h < num_heads; h++){
      for(int d = 0; d < Dh; d++){
        out(l_,h*Dh+d) = outHeads[h](l_,d);
      }
    }
  }



  Eigen::MatrixXf finalOut = out * Wo;

  return finalOut;
}

static Eigen::MatrixXf eigenFF(
  const Eigen::MatrixXf& seq,
  size_t L,
  size_t D,
  size_t H,
  Eigen::MatrixXf& W1,
  Eigen::MatrixXf& W2,
  Eigen::VectorXf& b1,
  Eigen::VectorXf& b2)
{
  // seq: L x D, do 2D hidden layer
  //int L = seq.rows();
  //int D = seq.cols();
  //int H = 2*D;

  static std::mt19937 rng(123);
  std::uniform_real_distribution<float> dist(-0.1f, 0.1f);

  auto rMat = [&](int r, int c){
    Eigen::MatrixXf m(r,c);
    for(int i=0; i<r; i++){
      for(int j=0; j<c; j++){
        m(i,j) = dist(rng);
      }
    }
    return m;
  };
  auto rVec = [&](int len){
    Eigen::VectorXf v(len);
    for(int i=0; i<len; i++){
      v(i) = dist(rng);
    }
    return v;
  };

  //Eigen::MatrixXf W1 = rMat(D,H);
  //Eigen::VectorXf b1 = rVec(H);
  //Eigen::MatrixXf W2 = rMat(H,D);
  //Eigen::VectorXf b2 = rVec(D);

  // tmp = seq*W1 + b1 => [L x H]
  Eigen::MatrixXf tmp = seq * W1;
  for(int i=0; i<L; i++){
    tmp.row(i) += b1.transpose(); // broadcast add
  }
  // ReLU
  tmp = tmp.cwiseMax(0.f);

  // out = tmp*W2 + b2 => [L x D]
  Eigen::MatrixXf out = tmp * W2;
  for(int i=0; i<L; i++){
    out.row(i) += b2.transpose();
  }
  return out;
}

static void eigenLayerNorm(Eigen::MatrixXf& seq, float eps=1e-5f)
{
  // LN across columns for each row
  int L = seq.rows();
  int D = seq.cols();
  for(int i=0; i<L; i++){
    float mean = seq.row(i).mean();
    // variance
    float var = 0.f;
    for(int d=0; d<D; d++){
      float diff = seq(i,d) - mean;
      var += diff*diff;
    }
    var /= float(D);
    float denom = 1.f / std::sqrt(var + eps);
    for(int d=0; d<D; d++){
      seq(i,d) = (seq(i,d) - mean)*denom; // ignoring gamma/beta
    }
  }
}

// --------------------------------------------------------------------
// runTransformerEigen: Projects data_array from D_in -> D_model (divisible by H),
// runs M encoder layers, then projects back to D_in if desired.
//
// data_array: shape [N x L x D_in]
// returns shape [N x L x D_in] if we do the final projection back
//
// Main "encoder" demonstration (one or more layers).
// We'll do 2 encoder layers as a initial test example:
//   For each layer:
//     Y = LN( X + MultiHeadAttn(X) )
//     Z = LN( Y + FF(Y) )
// Return Z
//
// FF is a simple 2-layer feed-forward network: ReLU(x*W1 + b1)*W2 + b2
// - x: [L x D_model]
// - W1: [D_model x FF_hidden_dim]
// - W2: [FF_hidden_dim x D_model]
// - b1: [FF_hidden_dim]
// - b2: [D_model]
//   Returns: [L x D_model]
//   
// Uses Eigen for matrix operations.
//
// Currently uses output projection where Z : [L x D_model] is projected back to [L x D_in]
// Z : [L x D_model] * W_out [D_model x D_in]=> [L x D_in]
// For Gate sizing need to take encoder output and pass it through a classification head.
// Each class corresponding to the libcell.
// --------------------------------------------------------------------
std::vector<std::vector<std::vector<float>>> MLGateSizer::runTransformerEigen(
    const std::vector<std::vector<std::vector<float>>>& data_array,
    int num_heads,
    size_t N,
    size_t L,
    size_t D_in,
    size_t D_model,
    size_t FF_hidden_dim,
    int num_encoder_layers)
{
  //size_t N = data_array.size();
  if (N == 0) return {};

  //size_t L = data_array[0].size();
  if (L == 0) return {};

  // Original input dimension
  //size_t D_in = data_array[0][0].size();

  // Decide how many heads we want
  //const int num_heads = 4;

  // Decide a model dimension that is divisible by num_heads
  //size_t D_model = 64;

  // Number of encoder layers
  //const int num_encoder_layers = 2;

  // We will random-initialize projection weights:
  static std::mt19937 rng(999);
  std::uniform_real_distribution<float> dist(-0.05f, 0.05f);

  // W_in: [D_in x D_model]
  auto randomMatrix = [&](size_t r, size_t c){
    Eigen::MatrixXf m(r,c);
    for(int i=0; i<r; i++){
      for(int j=0; j<c; j++){
        m(i,j) = dist(rng);
      }
    }
    return m;
  };
  auto randomVector = [&](size_t len){
    Eigen::VectorXf v(len);
    for(int i=0; i<len; i++){
      v(i) = dist(rng);
    }
    return v;
  };

  // If you want to map back to D_in:
  // W_out: [D_model x D_in]
  auto W_in  = randomMatrix(D_in, D_model);
  auto W_out = randomMatrix(D_model, D_in);

  // Allocate Q, K, V, O projection matrices
  // Wq/Wk/Wv/Wo: [D_model x D_model]
  auto Wq = randomMatrix(D_model, D_model);
  auto Wk = randomMatrix(D_model, D_model);
  auto Wv = randomMatrix(D_model, D_model);
  auto Wo = randomMatrix(D_model, D_model);

  // Allocate FF weights
  // FF_W1: [D_model x FF_hidden_dim]
  // FF_W2: [FF_hidden_dim x D_model]
  // FF_b1: [FF_hidden_dim]
  // FF_b2: [D_model]
  auto FF_W1 = randomMatrix(D_model, FF_hidden_dim);
  auto FF_W2 = randomMatrix(FF_hidden_dim, D_model);
  auto FF_b1 = randomVector(FF_hidden_dim);
  auto FF_b2 = randomVector(D_model);


  // Allocate output
  std::vector<std::vector<std::vector<float>>> output(N,
      std::vector<std::vector<float>>(L, std::vector<float>(D_in, 0.f)));
  
  
  // For each sequence in the batch:
  for (size_t n = 0; n < N; n++) {
    // current seq: [L x D_in]
    Eigen::MatrixXf seq(L, D_in);
    for (size_t l = 0; l < L; l++) {
      for (size_t d = 0; d < D_in; d++) {
        seq(l, d) = data_array[n][l][d];
      }
    }

    // 1) Project seq => [L x D_model]
    Eigen::MatrixXf seq_proj = seq * W_in;

    // 2) Run through encoder blocks
    Eigen::MatrixXf x = seq_proj; // shape [L x D_model]
    for (int layer = 0; layer < num_encoder_layers; layer++) {
      // (a) Multi-head self-attention
      auto attn_out = eigenSelfAttention(x, num_heads, D_model, L, Wq, Wk, Wv, Wo);
      // (b) Residual
      x += attn_out;
      // (c) LayerNorm
      eigenLayerNorm(x);

      // (d) FeedForward
      auto ff_out = eigenFF(x, L, D_model, FF_hidden_dim, FF_W1, FF_W2, FF_b1, FF_b2);
      // (e) Residual
      x += ff_out;
      // (f) LayerNorm
      eigenLayerNorm(x);
    }

    // x is now [L x D_model]. Project back to [L x D_in]
    Eigen::MatrixXf final_out = x * W_out;

    // Store final_out in output[n]
    for (size_t l = 0; l < L; l++) {
      for (size_t d = 0; d < D_in; d++) {
        output[n][l][d] = final_out(l, d);
      }
    }
  }


  // Return shape: [N x L x D_in]
  return output;
}

template <class Func>
long long MLGateSizer::benchmark(Func&& func)
{
  auto start = std::chrono::steady_clock::now();
  func();  // Execute the callable
  auto end = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  return elapsed;
}

bool MLGateSizer::compareOutputs(
    const std::vector<std::vector<std::vector<float>>>& A,
    const std::vector<std::vector<std::vector<float>>>& B,
    float tol) const
{
  if (A.size() != B.size()) return false;
  for (size_t n = 0; n < A.size(); n++) {
    if (A[n].size() != B[n].size()) return false;
    for (size_t l = 0; l < A[n].size(); l++) {
      if (A[n][l].size() != B[n][l].size()) return false;
      for (size_t d = 0; d < A[n][l].size(); d++) {
        float diff = std::fabs(A[n][l][d] - B[n][l][d]);
        if (diff > tol) {
          return false;
        }
      }
    }
  }
  return true;
}


}  // namespace rsz
