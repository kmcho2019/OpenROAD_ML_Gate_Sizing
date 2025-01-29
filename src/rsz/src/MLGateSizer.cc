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

#include <algorithm>
#include <fstream>
#include <set>

#include <iomanip>
#include <sstream>
#include <numeric>  // For std::accumulate

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


void MLGateSizer::addToken(const std::vector<float>& pin_data,
                           const std::string& gate_type)
{
  pin_tokens_.push_back(pin_data);
  gate_types_.push_back(gate_type);
}


void MLGateSizer::getEndpointAndCriticalPaths(const std::string& output_base_path,
                                              const std::string& tech_embedding_file_path,
                                              const std::string& label_size_file_path,
                                              const std::string& model_weight_file_path)
{
  // Output file base directory
  //const std::string output_base_path = "/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/NV_NVDLA_partition_m";
  // Ensure that output directory exists
  std::filesystem::create_directories(output_base_path);


  // Input files
  // Tech embedding load file path
  //const std::string tech_embedding_file_path = "/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/embedding_generation/ASAP7_libcell_embeddings.bin";
  // Label .size file load file path
  //const std::string label_size_file_path = "/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/NV_NVDLA_partition_m.size";
  // Transsizer model weight load file path
  //const std::string model_weight_file_path = "/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/transformer_params.bin";

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
      10,//5 * endpoints->size(),//100, //10 * endpoints->size(), // group_count
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

    std::cout << "All Libcell Names: " << std::endl;
    for (const auto& cell : ordered_libcells_) {
      std::cout << cell << std::endl;
    }

    std::cout << std::endl;


    // Attempt to load libcell embeddings from a binary file
    size_t embedding_size = 768; // Embedding size of deberta-v3-base model

    // Load the embeddings based on the number of libcells
    // Currently supports ASAP7 and Nangate45 temporary fix needs to be more robust and general
    if (libcell_to_id_.size() == 216) { // ASAP7
      //loadEmbeddingsBinary("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/embedding_generation/ASAP7_libcell_embeddings.bin", embedding_size);
      loadEmbeddingsBinary(tech_embedding_file_path, embedding_size);
      updateLibcellTypeEmbeddings();
      //const std::string type_embedding_filename = "/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/ASAP7_libcell_type_embeddings.bin";
      //const std::string weight_filename = "/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/transformer_params.bin";
      std::filesystem::path type_embedding_export_path = std::filesystem::path(output_base_path) / "ASAP7_libcell_type_embeddings.bin";
      exportTypeEmbeddings(type_embedding_export_path.string());
      loadWeights(model_weight_file_path); // Load the transformer weights
    }
    else if (libcell_to_id_.size() == 135) {  // Nangate45
      //loadEmbeddingsBinary("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/embedding_generation/nangate45_libcell_embeddings.bin", embedding_size);
      loadEmbeddingsBinary(tech_embedding_file_path, embedding_size);
      updateLibcellTypeEmbeddings();
      //const std::string type_embedding_filename = "/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/nangate45_libcell_type_embeddings.bin";
      std::filesystem::path type_embedding_export_path = std::filesystem::path(output_base_path) / "nangate45_libcell_type_embeddings.bin";
      exportTypeEmbeddings(type_embedding_export_path.string());
      const std::string weight_filename = "";
      //loadWeights(model_weight_file_path); // Load the transformer weights
      std::cout << "Unsupported, unable to load model weights." << std::endl;
      std::cout << "Not implemented for Nangate45 (only ASAP7 so far)" << std::endl;
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

    std::cout << std::endl;
    checkDataConsistencyAndPrint();


    // Normalization constants for PinMetrics
    // Tempoarily commented out, will be used for normalization
    /*
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

    */




  


    

    
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
                                      libcell_to_id_,
                                      libcell_to_type_id_,
                                      libcell_id_to_embedding_);

    auto [data_array_temp, pin_ids_temp, cell_ids_temp, libcell_ids_temp, libcell_type_ids_temp] = builder.build();

    // Take ownership by moving the data (rvalue references) into new variables
    std::vector<std::vector<std::vector<float>>> data_array = std::move(data_array_temp);
    std::vector<std::vector<int>> pin_ids = std::move(pin_ids_temp);
    std::vector<std::vector<int>> cell_ids = std::move(cell_ids_temp);
    std::vector<std::vector<int>> libcell_ids = std::move(libcell_ids_temp);
    std::vector<std::vector<int>> libcell_type_ids = std::move(libcell_type_ids_temp);

    
    // Debugging print statements to check data_array and libcell_type_ids
    
    // Print shape of the data_array, pin_ids, cell_ids, and libcell_type_ids
    // N = number of sequences, L = max sequence length, D = number of features(token dimensions = pin_data + embedding)
    // data_array shape: (N, L, D)
    // pin_ids shape: (N, L): Pin IDs are used to lookup the pin name
    // cell_ids shape: (N, L): Cell IDs are used to lookup the cell name
    // libcell_type_ids shape: (N, L): Libcell IDs are used to lookup the libcell name
    // libcell_type_ids shape: (N, L): Libcell Type IDs are used to lookup the cell type which are used to input the available embeddings
    std::cout << "Data Array Shape: (" << data_array.size() << ", " << data_array[0].size() << ", " << data_array[0][0].size() << ")" << std::endl;
    std::cout << "Pin IDs Shape: (" << pin_ids.size() << ", " << pin_ids[0].size() << ")" << std::endl;
    std::cout << "Cell IDs Shape: (" << cell_ids.size() << ", " << cell_ids[0].size() << ")" << std::endl;
    std::cout << "Libcell Type IDs Shape: (" << libcell_type_ids.size() << ", " << libcell_type_ids[0].size() << ")" << std::endl;

    // Print example of the data_array (N, L, D), first 5 tokens(20 dim of token or less) of the first 2 sequence
    for (size_t o = 0; o < ((2 > data_array.size()) ? data_array.size() : 2); o++) {
      std::cout << "Sequence " << o << ": " << std::endl;
      for (size_t i = 0; i < ((5 > data_array[o].size()) ? data_array[o].size() : 5); i++) {
        std::cout << "Token " << i << ": ";
        for (size_t j = 0; j < ((20 > data_array[0][0].size()) ? data_array[0][0].size() : 20); j++) {
          std::cout << data_array[o][i][j] << " ";
        }
        std::cout << std::endl;
      }
    }

    // Print libcell_type_ids for the first 5 tokens of the first 2 sequences
    for (size_t o = 0; o < ((2 > data_array.size()) ? data_array.size() : 2); o++) {
      std::cout << "Sequence " << o << ": " << std::endl;
      for (size_t i = 0; i < ((5 > data_array[o].size()) ? data_array[o].size() : 5); i++) {
        std::cout << "Libcell Type ID " << i << ": " << libcell_type_ids[o][i] << std::endl;
      }
    }

    // Print Pin, Cell and Cell Type IDs for the first 10 sequences
    for (size_t o = 0; o < ((10 > pin_ids.size()) ? pin_ids.size() : 10); o++) {
      std::cout << "Sequence " << o << ": " << std::endl;
      std::cout << "Pin ID: " << std::endl;
      for (size_t i = 0; i < pin_ids[0].size(); i++) {
        std::cout << pin_ids[o][i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Cell ID: " << std::endl;
      for (size_t i = 0; i < cell_ids[0].size(); i++) {
        std::cout << cell_ids[o][i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Libcell ID: " << std::endl;
      for (size_t i = 0; i < libcell_ids[0].size(); i++) {
        std::cout << libcell_ids[o][i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Libcell Type ID: " << std::endl;
      for (size_t i = 0; i < libcell_type_ids[0].size(); i++) {
        std::cout << libcell_type_ids[o][i] << " ";
      }
      std::cout << std::endl;
    }

    // Check if Cell ID and Libcell Type ID are consistent and if they occur in pairs (e.g. 4 4 0 0 1 1 2 2 3 3)





    
    /*
    // Print out the slack of each path (for debugging), remove later
    for (size_t i = 0; i < path_slacks.size(); i++) {
      std::cout << "Path " << i << " Slack: " << path_slacks[i] << std::endl;
    }
    */

    int num_heads = 8;//4;
    size_t N = data_array.size(); // Number of sequences (100)
    size_t L = data_array[0].size();  // Max sequence length (78)
    // check if L is even or not, assert if not with error L must be even
    assert(L % 2 == 0 && "L must be even");
    // size_t L_2 = L/2;  // Half of the sequence length, corresponds to the 2nd encoder layer input sequence length (commented out as not used)
    size_t D_in = data_array[0][0].size();  // Input feature dimensions (17)
    size_t D_out = 50; // Output feature dimensions (50), maximum number of classes
    size_t D_emb = embedding_size_; // Embedding dimensions (768)
    size_t D_model = 128;//64; // Transformer model hidden dimensions (64)
    size_t FF_hidden_dim = 4 * D_model; // Feedforward hidden dimensions (256)
    int num_encoder_layers = 6;//2;
    int num_encoder_layers_2 = 2; // Number of encoder layers for the 2nd encoder layer

    // Calculate the number of parameters in the model
    // Initial Projection weights: (D_in x D_model)
    // Transformer Encoder layers: num_encoder_layers
    // - Q, K, V, O weights: (D_model x D_model) x 4
    // - FF weights: (D_model x FF_hidden_dim) x 2 + FF_hidden_dim + D_model
    // Projection for the 2nd encoder layer: (D_emb x D_model)
    // Final Output Projection weights: (D_model x D_out) (optional consider switching to classification head, in that case D_out becomes the number of classes)
    size_t projection_param_num = D_in * D_model + D_emb * D_model + D_model * D_out;
    size_t encoder_1_param_num = ((D_model * D_model) * 4 + (D_model * FF_hidden_dim) * 2 + FF_hidden_dim + D_model) * num_encoder_layers;
    size_t encoder_2_param_num = ((D_model * D_model) * 4 + (D_model * FF_hidden_dim) * 2 + FF_hidden_dim + D_model) * num_encoder_layers_2;
    size_t model_params = projection_param_num + encoder_1_param_num + encoder_2_param_num;

    // encoder_2_input input to 2nd encoder layer shape: (N, L/2, D_emb)
    // consists of the embeddings of the libcell type IDs
    // from every other libcell type ID in the libcell_type_ids use the corresponding embedding from libcell_type_id_to_embedding_
    std::vector<std::vector<int>> encoder_2_input_libcell_type_ids;
    encoder_2_input_libcell_type_ids = std::vector<std::vector<int>>(N, std::vector<int>(L/2, -1));

    // Fill in the encoder_2_input_libcell_type_ids with every other libcell type ID from libcell_type_ids
    // As the libcell_type_ids are in the shape (N, L) representing pin's libcell type, we need to extract every other libcell type ID to get the cell's libcell type ID
    // Due to every 2 pin corresponds to a cell.
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < L/2; j++) {
        encoder_2_input_libcell_type_ids[i][j] = libcell_type_ids[i][j*2];
      }
    }

    // encoder_2_output_avail_libcell_num shape: (N, L/2)
    // The number of available libcell types for each cell when resizing
    // The number of available libcell types can be found by taking the length libcell_type_id_to_libcell_ids_ vector for given libcell_type_id
    std::vector<std::vector<int>> encoder_2_output_avail_libcell_num;
    encoder_2_output_avail_libcell_num = std::vector<std::vector<int>>(N, std::vector<int>(L/2, 0));
    // Initialize the encoder_2_output_avail_libcell_num with the number of available libcell types for each cell
    // by searching through encoder_2_input_libcell_type_ids and finding the corresponding libcell type ID in libcell_type_id_to_libcell_ids_
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < L/2; j++) {
        int libcell_type_id = encoder_2_input_libcell_type_ids[i][j];
        if (libcell_type_id >= 0 && libcell_type_id < libcell_type_id_to_libcell_ids_.size()) {
          encoder_2_output_avail_libcell_num[i][j] = libcell_type_id_to_libcell_ids_[libcell_type_id].size();
        }
      }
    }


    std::vector<std::vector<std::vector<float>>> encoder_2_input;
    encoder_2_input = std::vector<std::vector<std::vector<float>>>(N, std::vector<std::vector<float>>(L/2, std::vector<float>(embedding_size_, 0.0)));
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < L/2; j++) {
        int libcell_type_id = libcell_type_ids[i][j*2];
        // only perform lookup if the libcell type id is valid
        if (libcell_type_id >= 0 && libcell_type_id < libcell_type_id_to_embedding_.size()) {
          encoder_2_input[i][j] = libcell_type_id_to_embedding_[libcell_type_id];
        }
      }
    }


    // Generate the labels shape (N, L) for the transformer model integer classification labels corresponding to correct libcell ID for each cell
    // Read .size file to get cell->libcell mapping
    std::unordered_map<std::string, std::string> cell_name_to_libcell_name = readSizeFile(label_size_file_path);
    std::unordered_map<int, int> cell_id_to_libcell_id;
    for (const auto& [cell_name, libcell_name] : cell_name_to_libcell_name) {
      // Lookup cell ID - skip if not found
      auto cell_it = cell_name_to_id.find(cell_name);
      if (cell_it == cell_name_to_id.end()) continue;

      // Lookup libcell ID - skip if not found
      auto libcell_it = libcell_to_id_.find(libcell_name);
      if (libcell_it == libcell_to_id_.end()) continue;

      // Only insert if both IDs exist
      cell_id_to_libcell_id[cell_it->second] = libcell_it->second;
    }
    std::vector<std::vector<int>> labels; // (N, L/2) for each cell, the corresponding libcell's order within same type
    // Example if libcells 1, 4, and 7 are the same type, then the labels for the cells are 0, 1, 2
    // So if libcell 1 is the libcell in .size file, then the label for the cell is 0 and so on
    // If the cell is not found in the libcell_to_id_ mapping, then the label is 0
    // This is because the model is trained to predict within each libcell type.
    labels = std::vector<std::vector<int>>(N, std::vector<int>(L/2, 0));
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < L/2; j++) {
        int cell_id = cell_ids[i][j*2];
        // only perform lookup if the cell id is valid
        if (cell_id_to_libcell_id.find(cell_id) != cell_id_to_libcell_id.end()) {
          int libcell_id = cell_id_to_libcell_id[cell_id];
          // Get the libcell type id
          if (libcell_id_to_libcell_type_id_.find(libcell_id) != libcell_id_to_libcell_type_id_.end()) {
            int libcell_type_id = libcell_id_to_libcell_type_id_[libcell_id];
            // Find the index of the libcell_id within the type's vector
            if (libcell_type_id_to_libcell_ids_.find(libcell_type_id) != libcell_type_id_to_libcell_ids_.end()) {
                const std::vector<int>& libcells_of_type = libcell_type_id_to_libcell_ids_[libcell_type_id];
                auto it = std::find(libcells_of_type.begin(), libcells_of_type.end(), libcell_id);
                if (it != libcells_of_type.end()) {
                  labels[i][j] = std::distance(libcells_of_type.begin(), it);
                }
            }
          }
        }
      }
    }

    // Save the data_array, encoder_2_input, pin_ids, cell_ids, libcell_ids, libcell_type_ids to a binary file
    // The data_array is saved as a 3D array, pin_ids, cell_ids, libcell_ids, libcell_type_ids are saved as 2D arrays
    // The data_array is saved as a float32 array, pin_ids, cell_ids, libcell_ids, libcell_type_ids are saved as int32 arrays
    // Also save the labels as a 2D array



    // Save the 3D arrays
    // Group files by type using tuples (filename_suffix, data, writer_function)
    const std::vector<std::tuple<const char*, 
                                std::vector<std::vector<std::vector<float>>>&, 
                                void (rsz::MLGateSizer::*)(const std::string&, const std::vector<std::vector<std::vector<float>>>&)>> float3d_files = {
        {"./data_array.bin", data_array, &rsz::MLGateSizer::writeBinaryFile3DFloat},
        {"./encoder_2_input.bin", encoder_2_input, &rsz::MLGateSizer::writeBinaryFile3DFloat}
    };

    //writeBinaryFile3DFloat("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/data_array.bin", data_array);
    //writeBinaryFile3DFloat("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/encoder_2_input.bin", encoder_2_input);
    // Save the 2D arrays
    // Group files by type using tuples (filename_suffix, data, writer_function)
    const std::vector<std::tuple<const char*, 
                                std::vector<std::vector<int>>&, 
                                void (rsz::MLGateSizer::*)(const std::string&, const std::vector<std::vector<int>>&)>> int2d_files = {
        {"./pin_ids.bin", pin_ids, &rsz::MLGateSizer::writeBinaryFile2DInt},
        {"./cell_ids.bin", cell_ids, &rsz::MLGateSizer::writeBinaryFile2DInt},
        {"./libcell_ids.bin", libcell_ids, &rsz::MLGateSizer::writeBinaryFile2DInt},
        {"./libcell_type_ids.bin", libcell_type_ids, &rsz::MLGateSizer::writeBinaryFile2DInt},
        {"./encoder_2_input_libcell_type_ids.bin", encoder_2_input_libcell_type_ids, &rsz::MLGateSizer::writeBinaryFile2DInt},
        {"./encoder_2_output_avail_libcell_num.bin", encoder_2_output_avail_libcell_num, &rsz::MLGateSizer::writeBinaryFile2DInt},
        {"./labels.bin", labels, &rsz::MLGateSizer::writeBinaryFile2DInt}
    };
    
    //writeBinaryFile2DInt("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/pin_ids.bin", pin_ids);
    //writeBinaryFile2DInt("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/cell_ids.bin", cell_ids);
    //writeBinaryFile2DInt("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/libcell_ids.bin", libcell_ids);
    //writeBinaryFile2DInt("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/libcell_type_ids.bin", libcell_type_ids);
    //writeBinaryFile2DInt("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/encoder_2_input_libcell_type_ids.bin", encoder_2_input_libcell_type_ids);
    //writeBinaryFile2DInt("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/encoder_2_output_avail_libcell_num.bin", encoder_2_output_avail_libcell_num);
    //writeBinaryFile2DInt("/home/kmcho/2_Project/ML_GateSizing_OpenROAD/dev_repo/test_scripts/pytorch_transsizer_training_code/labels.bin", labels);

    // Process output for all files
    this->batch_process_files(output_base_path, float3d_files);
    this->batch_process_files(output_base_path, int2d_files);


    //auto outputs = runTransformer(data_array, encoder_2_input, num_heads, N, L, D_in, embedding_size_, D_model, FF_hidden_dim, num_encoder_layers, num_encoder_layers_2);  // Run the transformer model
    //auto eigen_outputs = runTransformerEigen(data_array, encoder_2_input, num_heads, N, L, D_in, embedding_size_, D_model, FF_hidden_dim, num_encoder_layers, num_encoder_layers_2);  // Run the transformer model using Eigen

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

    // 3) Print the speed & correctness info
    // total tokens processed: N*L
    size_t total_tokens = data_array.size() * data_array[0].size();
    size_t total_tokens_encoder2 = encoder_2_input.size() * encoder_2_input[0].size();

    // Calculate total padded tokens for encoder1 and encoder2
    // From encoder1 count the number of -1s in libcell_type_ids
    // From encoder2 count the number of -1s in encoder_2_input_libcell_type_ids
    size_t total_padded_tokens_encoder1 = 0;
    size_t total_padded_tokens_encoder2 = 0;
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < L/2; j++) {
        if (libcell_type_ids[i][2*j] == -1) {
          total_padded_tokens_encoder1++;
        }
        if (libcell_type_ids[i][2*j+1] == -1) {
          total_padded_tokens_encoder1++;
        }
        if (encoder_2_input_libcell_type_ids[i][j] == -1) {
          total_padded_tokens_encoder2++;
        }
      }
    }



    //std::cout << "Projection parameters: " << D_in << " x " << D_model << " + " << D_emb << " x " << D_model << " + " << D_model << " x " << D_out << " = " << projection_param_num << std::endl;
    std::cout << "Projection parameters: " << projection_param_num << std::endl;
    std::cout << "Encoder 1 parameters: " << encoder_1_param_num << std::endl;
    std::cout << "Encoder 2 parameters: " << encoder_2_param_num << std::endl;
    std::cout << "Total model parameters: " << model_params << std::endl;

    std::cout << "Total tokens for Encoder 1: " << total_tokens << ": " << data_array.size() << "X" << data_array[0].size() << std::endl;
    std::cout << "Total tokens for Encoder 2: " << total_tokens_encoder2 << ": " << encoder_2_input.size() << "X" << encoder_2_input[0].size() << std::endl;
    std::cout << "Total non-padded tokens for Encoder 1: " << total_tokens - total_padded_tokens_encoder1 << std::endl;
    std::cout << "Total non-padded tokens for Encoder 2: " << total_tokens_encoder2 - total_padded_tokens_encoder2 << std::endl;
    std::cout << "Total padded tokens for Encoder 1: " << total_padded_tokens_encoder1 << std::endl;
    std::cout << "Total padded tokens for Encoder 2: " << total_padded_tokens_encoder2 << std::endl;
    std::cout << "Average valid length for Encoder 1: " << (total_tokens - total_padded_tokens_encoder1) / N << std::endl;
    std::cout << "Average valid length for Encoder 2: " << (total_tokens_encoder2 - total_padded_tokens_encoder2) / N << std::endl;

    std::cout << "1st Encoder Input shape: (" << data_array.size() << ", " << data_array[0].size() << ", " << data_array[0][0].size() << ")" << std::endl;
    std::cout << "2nd Encoder Input shape: (" << encoder_2_input.size() << ", " << encoder_2_input[0].size() << ", " << encoder_2_input[0][0].size() << ")" << std::endl;

    bool compare_naive_eigen = false;

    if (compare_naive_eigen) {
      // 3-1) Actually store the outputs to compare correctness
      auto start = std::chrono::steady_clock::now();
      auto out_naive  = runTransformer(data_array, encoder_2_input, num_heads, N, L, D_in, D_out, embedding_size_, D_model, FF_hidden_dim, num_encoder_layers, num_encoder_layers_2);
      auto end = std::chrono::steady_clock::now();
      auto naive_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
      start = std::chrono::steady_clock::now();
      auto out_eigen  = runTransformerEigen(data_array, encoder_2_input, num_heads, N, L, D_in, D_out, embedding_size_, D_model, FF_hidden_dim, num_encoder_layers, num_encoder_layers_2);
      end = std::chrono::steady_clock::now();
      auto eigen_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
      bool ok = compareOutputs(out_naive, out_eigen, 1e-3f);

      std::cout << "Output shape: (" << out_naive.size() << ", " << out_naive[0].size() << ", " << out_naive[0][0].size() << ")" << std::endl;

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


    // Test using loaded model if transformer_weights_ is not empty
    if (transformer_weights_.loaded) {

      // Check if transformer_weights_ are loaded properly and actually contain the weights
      // Iterate through the weights and print the shape of each weight
      // Print projection matrices dimensions
      std::cout << "Projection Matrices:" << std::endl;
      std::cout << "W_in1: " << transformer_weights_.W_in1.rows() << "x" << transformer_weights_.W_in1.cols() << std::endl;
      std::cout << "W_in2: " << transformer_weights_.W_in2.rows() << "x" << transformer_weights_.W_in2.cols() << std::endl;
      std::cout << "W_out: " << transformer_weights_.W_out.rows() << "x" << transformer_weights_.W_out.cols() << std::endl;

      // Print encoder1 weights dimensions (6 layers)
      for (int layer = 0; layer < num_encoder_layers; layer++) {
          std::cout << "\nEncoder1 Layer " << layer << ":" << std::endl;
          std::cout << "Wq: " << transformer_weights_.encoder1_0_Wq_1[layer].rows() << "x" 
                    << transformer_weights_.encoder1_0_Wq_1[layer].cols() << std::endl;
          std::cout << "Wk: " << transformer_weights_.encoder1_0_Wk_1[layer].rows() << "x"
                    << transformer_weights_.encoder1_0_Wk_1[layer].cols() << std::endl;
          std::cout << "Wv: " << transformer_weights_.encoder1_0_Wv_1[layer].rows() << "x"
                    << transformer_weights_.encoder1_0_Wv_1[layer].cols() << std::endl;
          std::cout << "FF_W1: " << transformer_weights_.encoder1_0_FF_W1_1[layer].rows() << "x"
                    << transformer_weights_.encoder1_0_FF_W1_1[layer].cols() << std::endl;
          std::cout << "FF_b1: " << transformer_weights_.encoder1_0_FF_b1_1[layer].size() << std::endl;
          std::cout << "FF_W2: " << transformer_weights_.encoder1_0_FF_W2_1[layer].rows() << "x"
                    << transformer_weights_.encoder1_0_FF_W2_1[layer].cols() << std::endl;
          std::cout << "FF_b2: " << transformer_weights_.encoder1_0_FF_b2_1[layer].size() << std::endl;
      }

      // Print encoder2 weights dimensions (2 layers)
      for (int layer = 0; layer < num_encoder_layers_2; layer++) {
          std::cout << "\nEncoder2 Layer " << layer << ":" << std::endl;
          std::cout << "Wq: " << transformer_weights_.encoder2_0_Wq_2[layer].rows() << "x"
                    << transformer_weights_.encoder2_0_Wq_2[layer].cols() << std::endl;
          std::cout << "Wk: " << transformer_weights_.encoder2_0_Wk_2[layer].rows() << "x"
                    << transformer_weights_.encoder2_0_Wk_2[layer].cols() << std::endl;
          std::cout << "Wv: " << transformer_weights_.encoder2_0_Wv_2[layer].rows() << "x"
                    << transformer_weights_.encoder2_0_Wv_2[layer].cols() << std::endl;
          std::cout << "FF_W1: " << transformer_weights_.encoder2_0_FF_W1_2[layer].rows() << "x"
                    << transformer_weights_.encoder2_0_FF_W1_2[layer].cols() << std::endl;
          std::cout << "FF_b1: " << transformer_weights_.encoder2_0_FF_b1_2[layer].size() << std::endl;
          std::cout << "FF_W2: " << transformer_weights_.encoder2_0_FF_W2_2[layer].rows() << "x"
                    << transformer_weights_.encoder2_0_FF_W2_2[layer].cols() << std::endl;
          std::cout << "FF_b2: " << transformer_weights_.encoder2_0_FF_b2_2[layer].size() << std::endl;
      }

      auto start_ = std::chrono::steady_clock::now();
      auto loaded_eigen_output = runTransformerEigen(data_array, encoder_2_input, num_heads, N, L, D_in, D_out, embedding_size_, D_model, FF_hidden_dim, num_encoder_layers, num_encoder_layers_2, transformer_weights_);
      auto end_ = std::chrono::steady_clock::now();

      auto loaded_eigen_us = std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count();

      std::cout << "Loaded Eigen  time: " << loaded_eigen_us << " us   => "
                << (1e6 * double(total_tokens) / double(loaded_eigen_us))
                << " tokens/sec\n";

      // Based on loaded_eigen_output, determine gate sizes, calculate accuracy (compared to label), and apply the gate sizes to the design
      // The loaded_eigen_output is in the shape (N, L/2, D_out) where D_out is the number of classes
      // The labels is in the shape (N, L/2) where the value is the index of the correct libcell among the same type
      // The gate sizes are determined by the index of the maximum value in the D_out dimension of the loaded_eigen_output
      // The accuracy is calculated by comparing the maximum value index to the label value
      // The gate sizes are applied to the design by updating the libcell name of the cell with the corresponding libcell name from the libcell_id_to_libcell_ mapping

      // Peform argmax on the D_out dimension of the loaded_eigen_output to get the predicted libcell index
      // Use encoder_2_output_avail_libcell_num to only consider the available libcell types for each cell
      // For example if encoder_2_output_avail_libcell_num[i][j] is 3, then only consider the first 3 elements of the D_out dimension
      std::vector<std::vector<int>> predicted_libcell_indices;
      predicted_libcell_indices = std::vector<std::vector<int>>(N, std::vector<int>(L/2, 0));
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < L/2; j++) {
          int max_index = 0;
          float max_value = loaded_eigen_output[i][j][0];
          int num_avail_libcells = encoder_2_output_avail_libcell_num[i][j]; // Number of available libcell types for the cell
          for (size_t k = 1; k < static_cast<size_t> (num_avail_libcells); k++) {
            if (loaded_eigen_output[i][j][k] > max_value) {
              max_index = k;
              max_value = loaded_eigen_output[i][j][k];
            }
          }
          predicted_libcell_indices[i][j] = max_index;
        }
      }

      // Calculate the accuracy by comparing the predicted libcell index to the label
      // Ignore the padding tokens by checking for -1 in the encoder_2_input_libcell_type_ids
      size_t num_correct = 0;
      size_t total_cells = 0;
      float accuracy = 0.0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < L/2; j++) {
          if (encoder_2_input_libcell_type_ids[i][j] != -1) {
            total_cells++;
            if (predicted_libcell_indices[i][j] == labels[i][j]) {
              num_correct++;
            }
          }
        }
      }
      if (total_cells > 0) {
        accuracy = (float)num_correct / (float)total_cells;
      }

      std::cout << "Number of non-padding token cells: " << total_cells << std::endl;
      std::cout << "Number of correct predictions: " << num_correct << std::endl;
      std::cout << "Accuracy: " << accuracy << std::endl;

      // Generate a map of cell ID to a vector of libcell IDs that the TransSizer model predicted
      std::unordered_map<int, std::vector<int>> cell_id_to_predicted_libcell_ids;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < L/2; j++) {
          if (encoder_2_input_libcell_type_ids[i][j] != -1) {
            int cell_id = cell_ids[i][j*2];
            int predicted_libcell_id = 0;
            // Check if predicted_libcell_indices[i][j] exceeds the bounds of the libcell_type_id_to_libcell_ids_ vector
            // This happens as the maximum number of classes is set to 50, but the actual number of classes is less
            if (predicted_libcell_indices[i][j] >= libcell_type_id_to_libcell_ids_[encoder_2_input_libcell_type_ids[i][j]].size()) {
              predicted_libcell_indices[i][j] = 0; // Set to 0 if out of bounds to prevent crash
            }
            else {
              predicted_libcell_id = libcell_type_id_to_libcell_ids_[encoder_2_input_libcell_type_ids[i][j]][predicted_libcell_indices[i][j]];
            }
            cell_id_to_predicted_libcell_ids[cell_id].push_back(predicted_libcell_id);
          }
        }
      }

      // If there are multiple predicted libcell IDs for a cell, the correct way to break ties is by selecting the libcell with the lower fo4 delay
      // Generate a map of cell ID to the libcell ID with ties broken by the lower fo4 delay
      // Need to have a way of looking up the fo4 delay of a libcell ID
      // Use dynamic programming to store the fo4 delay of each libcell ID and then use the map to lookup the fo4 delay of the predicted libcell IDs
      std::unordered_map<int, int> cell_id_to_predicted_libcell_id;
      

    }

  }



  // Next step is extracting the data from each pin
  // Then, transform the data into the format expected by the transformer model
  // Finally, apply the transformer model to classify and resize the gates

  // Use PinMetrics and PinDataSequence to store the extracted data

}


void MLGateSizer::writeBinaryFile3DFloat(const std::string& filename, 
                                 const std::vector<std::vector<std::vector<float>>>& data)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        logger_->error(utl::RSZ, 1001, "Cannot open file {} for writing (writeBinaryFile3DFloat)", filename);
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
        logger_->error(utl::RSZ, 1002, "Error writing to file {} (writeBinaryFile3DFloat)", filename);
    }
}

void MLGateSizer::writeBinaryFile2DInt(const std::string& filename, 
                                 const std::vector<std::vector<int>>& data)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        logger_->error(utl::RSZ, 1003, "Cannot open file {} for writing (writeBinaryFile2DInt)", filename);
        return;
    }

    // Write dimensions
    size_t N = data.size();
    size_t L = N > 0 ? data[0].size() : 0;

    out.write(reinterpret_cast<const char*>(&N), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&L), sizeof(size_t));

    // Write data in flattened format
    for (const auto& seq : data) {
        out.write(reinterpret_cast<const char*>(seq.data()), L * sizeof(int));
    }

    if (!out) {
        logger_->error(utl::RSZ, 1004, "Error writing to file {} (writeBinaryFile2DInt)", filename);
    }
}

// Read .size file to generate labels
// .size format: libcell_name, libcell_type_id for each line of the file
std::unordered_map<std::string, std::string> MLGateSizer::readSizeFile(const std::string& filename)
{
    std::unordered_map<std::string, std::string> libcell_to_type;
    std::ifstream in(filename);
    if (!in) {
        logger_->error(utl::RSZ, 1005, "Cannot open file {} for reading (readSizeFile)", filename);
        return libcell_to_type;
    }

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string libcell_name, libcell_type;
        if (!(iss >> libcell_name >> libcell_type)) {
            logger_->error(utl::RSZ, 1006, "Error reading line from file {} (readSizeFile)", filename);
            return libcell_to_type;
        }
        libcell_to_type[libcell_name] = libcell_type;
    }

    return libcell_to_type;
}


void MLGateSizer::exportTypeEmbeddings(const std::string& filename)
{
  const auto& embedding_map = libcell_type_id_to_embedding_;
  
  if (embedding_map.empty()) {
    logger_->error(utl::RSZ, 1007, "Embedding map is empty (exportTypeEmbeddings)");
    return;
  }

  // Collect and sort keys
  std::vector<int> keys;
  for (const auto& pair : embedding_map) {
    keys.push_back(pair.first);
  }
  std::sort(keys.begin(), keys.end());

  // Check if keys are contiguous starting from 0
  size_t N = keys.size();
  for (size_t i = 0; i < N; ++i) {
    if (keys[i] != static_cast<int>(i)) {
      logger_->error(utl::RSZ, 1008, "libcell_type_ids are not contiguous starting from 0 (exportTypeEmbeddings)");
      return;
    }
  }

  // Check all embeddings have the same dimension
  size_t D = embedding_map.at(keys[0]).size();
  for (const auto& key : keys) {
  if (embedding_map.at(key).size() != D) {
    logger_->error(utl::RSZ, 1009, "Inconsistent embedding dimensions (exportTypeEmbeddings)");
    return;
  }
  }

  // Open file
  std::ofstream out(filename, std::ios::binary);
  if (!out) {
    logger_->error(utl::RSZ, 1010, "Cannot open file {} for writing (exportTypeEmbeddings)", filename);
    return;
  }

  // Write dimensions (N and D)
  out.write(reinterpret_cast<const char*>(&N), sizeof(size_t));
  out.write(reinterpret_cast<const char*>(&D), sizeof(size_t));

  // Write embeddings in order of sorted keys
  for (const auto& key : keys) {
    const auto& embedding = embedding_map.at(key);
    out.write(reinterpret_cast<const char*>(embedding.data()), D * sizeof(float));
  }

  if (!out) {
    logger_->error(utl::RSZ, 1011, "Error writing to file {} (exportTypeEmbeddings)", filename);
  }
}

std::vector<std::vector<std::vector<float>>> MLGateSizer::readBinaryFile(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        logger_->error(utl::RSZ, 1012, "Cannot open file {} for reading (readBinaryFile)", filename);
        return {};
    }

    // Read dimensions
    size_t N, L, D;
    in.read(reinterpret_cast<char*>(&N), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&L), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&D), sizeof(size_t));

    if (!in) {
        logger_->error(utl::RSZ, 1013, "Error reading dimensions from file {} (readBinaryFile)", filename);
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
                logger_->error(utl::RSZ, 1014, "Error reading data from file {} (readBinaryFile)", filename);
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
void MLGateSizer::updateLibcellTypeMap(const std::unordered_map<std::string, int>& libcell_to_type_id_temp)
{
    libcell_to_type_id_.clear();
    libcell_type_id_to_libcell_ids_.clear();
    libcell_id_to_libcell_type_id_.clear();

	for (const auto& libcell : libcell_to_type_id_temp) {
		const std::string& libcell_name = libcell.first;
		int libcell_type_id = libcell.second;
		int libcell_id = libcell_to_id_[libcell_name];

    // Update the type maps
    libcell_to_type_id_[libcell_name] = libcell_type_id;
		libcell_type_id_to_libcell_ids_[libcell_type_id].push_back(libcell_id);
		libcell_id_to_libcell_type_id_[libcell_id] = libcell_type_id;
	}

}


void MLGateSizer::saveEmbeddingsBinary(const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        logger_->error(utl::RSZ, 1015, "Cannot open embeddings file {} (saveEmbeddingsBinary)", filename);
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
        logger_->error(utl::RSZ, 1016, "Cannot open embeddings file {} (loadEmbeddingsBinary)", filename);
        return;
    }

    // Read dimensions
    size_t num_items, file_embedding_size;
    in.read(reinterpret_cast<char*>(&num_items), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&file_embedding_size), sizeof(size_t));
    
    if (num_items != ordered_libcells_.size() || file_embedding_size != embedding_size) {
        logger_->error(utl::RSZ, 1017, 
            "Embedding file dimensions mismatch. Expected {} items of size {}, got {} items of size {} (loadEmbeddingsBinary)",
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
            logger_->error(utl::RSZ, 1018, "Error reading embedding {} from file (loadEmbeddingsBinary)", i);
            return;
        }
        libcell_id_to_embedding_[i] = embedding;
    }
}

void MLGateSizer::checkDataConsistencyAndPrint()
{
  bool consistent = true; // Flag to keep track of any inconsistency found

  // 1. Check that the size of ordered_libcells_ matches libcell_to_id_.size().
  if (ordered_libcells_.size() != libcell_to_id_.size()) {
    std::cerr << "[Error] Mismatch in number of libcells: "
              << "ordered_libcells_.size()=" << ordered_libcells_.size()
              << ", libcell_to_id_.size()=" << libcell_to_id_.size() << "\n";
    consistent = false;
  }

  // 2. Check each libcell in ordered_libcells_ to see if it maps back correctly in libcell_to_id_.
  for (int id = 0; id < static_cast<int>(ordered_libcells_.size()); ++id) {
    const std::string &name = ordered_libcells_[id];

    // (a) Ensure name exists in libcell_to_id_.
    auto it_id = libcell_to_id_.find(name);
    if (it_id == libcell_to_id_.end()) {
      std::cerr << "[Error] Libcell name '" << name 
                << "' not found in libcell_to_id_ map.\n";
      consistent = false;
    } else {
      // (b) Check that the ID matches the index in ordered_libcells_.
      if (it_id->second != id) {
        std::cerr << "[Error] Inconsistent ID for libcell name '" << name
                  << "': expected " << id << ", but found " << it_id->second 
                  << " in libcell_to_id_.\n";
        consistent = false;
      }
    }

    // (c) Ensure name exists in libcell_to_type_id_.
    auto it_type = libcell_to_type_id_.find(name);
    if (it_type == libcell_to_type_id_.end()) {
      std::cerr << "[Error] Libcell name '" << name
                << "' not found in libcell_to_type_id_ map.\n";
      consistent = false;
      continue; // Can't check type consistency below without this
    }
    int type_from_name = it_type->second;

    // (d) Check if we have an entry in libcell_id_to_libcell_type_id_ for this ID.
    auto it_id_to_type = libcell_id_to_libcell_type_id_.find(id);
    if (it_id_to_type == libcell_id_to_libcell_type_id_.end()) {
      std::cerr << "[Error] Libcell ID " << id
                << " not found in libcell_id_to_libcell_type_id_ map.\n";
      consistent = false;
      continue;
    }
    int type_from_id = it_id_to_type->second;

    // (e) Check that the type IDs match.
    if (type_from_name != type_from_id) {
      std::cerr << "[Error] Inconsistent type for libcell name '" << name
                << "' (type_from_name=" << type_from_name
                << ", type_from_id=" << type_from_id << ").\n";
      consistent = false;
    }

    // (f) Check that libcell_type_id_to_libcell_ids_[type] contains this ID.
    auto it_type_to_ids = libcell_type_id_to_libcell_ids_.find(type_from_name);
    if (it_type_to_ids == libcell_type_id_to_libcell_ids_.end()) {
      std::cerr << "[Error] Type ID " << type_from_name
                << " not found in libcell_type_id_to_libcell_ids_ map.\n";
      consistent = false;
    } else {
      const auto &ids_for_type = it_type_to_ids->second;
      if (std::find(ids_for_type.begin(), ids_for_type.end(), id) == ids_for_type.end()) {
        std::cerr << "[Error] Libcell ID " << id << " (name='" << name
                  << "') not found in libcell_type_id_to_libcell_ids_["
                  << type_from_name << "].\n";
        consistent = false;
      }
    }
  }

  // After all checks, if inconsistent, report and exit early or continue to print partial data.
  if (!consistent) {
    std::cerr << "[Warning] Inconsistencies found in the data structures.\n";
    // Can either return early or proceed to print partial data anyway.
    // return; 
  }

  // 3. Print out the organized data.
  //    a) Number of libcells
  std::cout << "Total number of libcells: " << ordered_libcells_.size() << "\n";
  
  //    b) Number of libcell types
  //       We take the distinct keys from libcell_type_id_to_libcell_ids_ as the count of types.
  std::cout << "Total number of libcell types: "
            << libcell_type_id_to_libcell_ids_.size() << "\n";

  //    c) Print each type in ascending order, then the libcells belonging to that type.
  //       First, gather the type IDs in a vector and sort them.
  std::vector<int> type_ids;
  type_ids.reserve(libcell_type_id_to_libcell_ids_.size());
  for (const auto &pair : libcell_type_id_to_libcell_ids_) {
    type_ids.push_back(pair.first);
  }
  std::sort(type_ids.begin(), type_ids.end());

  //    d) For each type ID in ascending order, print the libcells in that type.
  for (int type_id : type_ids) {
    std::cout << "\nType " << type_id << ":\n";
    const auto &ids_for_type = libcell_type_id_to_libcell_ids_.at(type_id);
    for (int libcell_id : ids_for_type) {
      if (libcell_id >= 0 && libcell_id < static_cast<int>(ordered_libcells_.size())) {
        std::cout << "  - ID = " << libcell_id
                  << ", Name = " << ordered_libcells_[libcell_id] << "\n";
      } else {
        std::cerr << "  [Error] Invalid libcell_id " << libcell_id
                  << " out of range [0, " << ordered_libcells_.size() - 1 << "].\n";
      }
    }
  }

  if (consistent) {
    std::cout << "\n[Info] Data structures are consistent.\n";
  } else {
    std::cout << "\n[Warning] Finished checking, but inconsistencies were found.\n";
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

// For a 2D param with shape [rows, cols], fill Eigen::MatrixXf:
//   - If transpose_if_needed is false, expect pd.shape == mat.shape
//   - If transpose_if_needed is true, expect pd.shape == [cols, rows] (Pytorch convention)
bool MLGateSizer::fillEigenMatrix(const ParamData& pd,
                     Eigen::MatrixXf& mat,
                     bool transpose_if_needed,
                     std::string& errMsg)
{
    // Must be exactly 2D
    if (pd.shape.size() != 2) {
        errMsg = "Expected 2D, but shape has " + std::to_string(pd.shape.size()) + " dims.";
        return false;
    }

    // Debugging
    //std::cout << "Param shape: [" << pd.shape[0] << "x" << pd.shape[1] << "]\n";
    //std::cout << "Eigen shape: [" << mat.rows() << "x" << mat.cols() << "]\n";
    //std::cout << "Data size: " << pd.values.size() << "\n";

    // The shape from the file
    size_t rows = pd.shape[0];
    size_t cols = pd.shape[1];

    if (!transpose_if_needed) {
        // Expect pd.shape == mat.shape
        if (rows != (size_t)mat.rows() || cols != (size_t)mat.cols()) {
            errMsg = "Dimension mismatch. Param shape is [" 
                      + std::to_string(rows) + "x" + std::to_string(cols) 
                      + "], but Eigen matrix is ["
                      + std::to_string(mat.rows()) + "x"
                      + std::to_string(mat.cols()) + "].";
            return false;
        }

        // Fill in row-major order
        size_t idx = 0;
        for (int r = 0; r < mat.rows(); r++) {
            for (int c = 0; c < mat.cols(); c++) {
                mat(r, c) = pd.values[idx++];
            }
        }
    }
    else {
        // If transposing, we expect param shape reversed:
        //    param is [rows, cols] but we want [cols, rows]
        if (cols != (size_t)mat.rows() || rows != (size_t)mat.cols()) {
            errMsg = 
                "Transpose mismatch. Param shape is [" 
                + std::to_string(rows) + "x" + std::to_string(cols) 
                + "], but after transpose we expected ["
                + std::to_string(mat.rows()) + "x" 
                + std::to_string(mat.cols()) + "].";
            return false;
        }

        // Fill: interpret the original data as if it’s transposed
        // mat(r, c) = param(c, r), which in flattened row-major is param[c * cols + r].
        for (int r = 0; r < mat.rows(); r++) {
            for (int c = 0; c < mat.cols(); c++) {
              // 'cols' here is pd.shape[1], i.e. the original param's "columns".
              mat(r, c) = pd.values[ static_cast<size_t>(c) * cols + r ];
            }
        }
    }

    return true;
}


// For a 1D param of shape [length], fill Eigen::VectorXf, etc.
bool MLGateSizer::fillEigenVector(const ParamData& pd,
                     Eigen::VectorXf& vec,
                     std::string& errMsg)
{
    // Must be exactly 1D
    if (pd.shape.size() != 1) {
        errMsg = "Expected 1D, but shape has " + std::to_string(pd.shape.size()) + " dims.";
        return false;
    }

    size_t length = pd.shape[0];
    if (length != (size_t)vec.size()) {
        errMsg = "Dimension mismatch. Param shape is ["
                  + std::to_string(length) + "], but Eigen vector is ["
                  + std::to_string(vec.size()) + "].";
        return false;
    }

    for (int i = 0; i < vec.size(); i++) {
        vec(i) = pd.values[i];
    }
    return true;
}



// Example function that loads all parameters into a name->ParamData map
bool MLGateSizer::loadModelWeightsRobust(const std::string& filename, 
                            std::unordered_map<std::string, ParamData>& out_map,
                            std::string& errMsg)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        errMsg = "Could not open file: " + filename;
        return false;
    }

    // 1) param_count
    size_t param_count = 0;
    in.read(reinterpret_cast<char*>(&param_count), sizeof(size_t));
    if (!in) {
        errMsg = "Error reading param_count";
        return false;
    }

    for (size_t p = 0; p < param_count; p++) {
        // a) name_length
        size_t name_len = 0;
        in.read(reinterpret_cast<char*>(&name_len), sizeof(size_t));
        if (!in) { errMsg = "Error reading name_len"; return false; }

        // b) name
        std::string name;
        name.resize(name_len);
        in.read(reinterpret_cast<char*>(&name[0]), name_len);
        if (!in) { errMsg = "Error reading param name"; return false; }

        // c) ndims
        size_t ndims = 0;
        in.read(reinterpret_cast<char*>(&ndims), sizeof(size_t));
        if (!in) { errMsg = "Error reading ndims"; return false; }

        // d) dims
        std::vector<size_t> dims(ndims);
        for(size_t i = 0; i < ndims; i++) {
            size_t dsize;
            in.read(reinterpret_cast<char*>(&dsize), sizeof(size_t));
            if(!in) { errMsg = "Error reading dimension size"; return false; }
            dims[i] = dsize;
        }

        // compute total size
        size_t total_size = 1;
        for (auto dsz : dims) {
            total_size *= dsz;
        }

        // e) data
        std::vector<float> values(total_size);
        in.read(reinterpret_cast<char*>(values.data()), total_size * sizeof(float));
        if (!in) { errMsg = "Error reading float data for param " + name; return false; }

        // store in map
        ParamData pd;
        pd.shape = dims;
        pd.values = std::move(values);
        out_map[name] = std::move(pd);
    }

    return true;
}

// Load weights from a file to update transformer_weights_
void MLGateSizer::loadWeights(const std::string& weight_file) {
    std::string errMsg;
    std::unordered_map<std::string, ParamData> param_map;
    bool ok = loadModelWeightsRobust(weight_file, param_map, errMsg);
    
    // Print param_map structure for debugging
    for (const auto& pair : param_map) {
        const std::string& name = pair.first;
        const ParamData& pd = pair.second;
        std::cout << "Param: " << name << " Shape: [";
        for (size_t i = 0; i < pd.shape.size(); i++) {
            std::cout << pd.shape[i] << (i < pd.shape.size() - 1 ? ", " : "");
        }
        std::cout << "]\n";
    }

    // Print Model Summary of the loaded weights based on param_map
    printModelSummary(param_map);

    
    if (!ok) {
        logger_->error(utl::RSZ, 1019, "Cannot load weights: {} (loadWeights)", errMsg);
        return;
    }
    
    // Load projection layers
    const std::vector<std::string> proj_layers = {
        "proj_in1.weight", "proj_in2.weight", "proj_out.weight"
    };
    
    for (const auto& name : proj_layers) {
      
      auto it = param_map.find(name);
      if (it == param_map.end()) {
          logger_->error(utl::RSZ, 1020, "Missing param {} (loadWeights)", name);
          return;
      }

      const ParamData& pd = it->second;
      bool transpose = true;

      // Determine target dimensions based on transpose flag
      int target_rows = transpose ? pd.shape[1] : pd.shape[0];
      int target_cols = transpose ? pd.shape[0] : pd.shape[1];
      
      // Debugging
      //std::cout << "Loading " << name << " with shape [" << target_rows << "x" << target_cols << "]\n";
      //std::cout << "Data size: " << pd.values.size() << "\n";

      if (name == "proj_in1.weight") {
        transformer_weights_.W_in1.resize(target_rows, target_cols);
        if (!fillEigenMatrix(pd, transformer_weights_.W_in1, transpose, errMsg)) {
          logger_->error(utl::RSZ, 1021, "Failed to fill {}: {} (loadWeights)", name, errMsg);
          return;
        }
      } else if (name == "proj_in2.weight") {
        transformer_weights_.W_in2.resize(target_rows, target_cols);
        if (!fillEigenMatrix(pd, transformer_weights_.W_in2, transpose, errMsg)) {
          logger_->error(utl::RSZ, 1022, "Failed to fill {}: {} (loadWeights)", name, errMsg);
          return;
        }
      } else if (name == "proj_out.weight") {
        transformer_weights_.W_out.resize(target_rows, target_cols);
        if (!fillEigenMatrix(pd, transformer_weights_.W_out, transpose, errMsg)) {
          logger_->error(utl::RSZ, 1023, "Failed to fill {}: {} (loadWeights)", name, errMsg);
          return;
        }          
      }
    }
    // Determine the max layer indices from the keys in param_map.
    int max_enc1_layer = -1;
    int max_enc2_layer = -1;

    for (const auto& kv : param_map) {
      const std::string& key = kv.first;
      
      // Check if the key starts with "encoders1."
      static const std::string enc1_prefix = "encoders1.";
      if (key.rfind(enc1_prefix, 0) == 0) {  // rfind(..., 0) == 0 means "starts with"
        // Find the dot after the layer index
        // "encoders1." is length 10, so the layer index starts at key[10]
        size_t dot_pos = key.find('.', enc1_prefix.size());
        // Extract substring representing the layer number
        std::string layer_str = key.substr(enc1_prefix.size(), dot_pos - enc1_prefix.size());
        int layer_idx = std::stoi(layer_str);
        max_enc1_layer = std::max(max_enc1_layer, layer_idx);
      }

      // Check if the key starts with "encoders2."
      static const std::string enc2_prefix = "encoders2.";
      if (key.rfind(enc2_prefix, 0) == 0) {
        size_t dot_pos = key.find('.', enc2_prefix.size());
        std::string layer_str = key.substr(enc2_prefix.size(), dot_pos - enc2_prefix.size());
        int layer_idx = std::stoi(layer_str);
        max_enc2_layer = std::max(max_enc2_layer, layer_idx);
      }
    }

    // Convert to a total layer count
    int num_encoder1_layers = max_enc1_layer + 1; // e.g. if max layer index was 5, that means 6 layers.
    int num_encoder2_layers = max_enc2_layer + 1; // e.g. if max layer index was 1, that means 2 layers.

    // Resize the relevant vectors in transformer_weights_ to the correct size
    transformer_weights_.encoder1_0_Wq_1.resize(num_encoder1_layers);
    transformer_weights_.encoder1_0_Wk_1.resize(num_encoder1_layers);
    transformer_weights_.encoder1_0_Wv_1.resize(num_encoder1_layers);
    transformer_weights_.encoder1_0_Wo_1.resize(num_encoder1_layers);

    transformer_weights_.encoder1_0_FF_W1_1.resize(num_encoder1_layers);
    transformer_weights_.encoder1_0_FF_b1_1.resize(num_encoder1_layers);
    transformer_weights_.encoder1_0_FF_W2_1.resize(num_encoder1_layers);
    transformer_weights_.encoder1_0_FF_b2_1.resize(num_encoder1_layers);

    transformer_weights_.encoder2_0_Wq_2.resize(num_encoder2_layers);
    transformer_weights_.encoder2_0_Wk_2.resize(num_encoder2_layers);
    transformer_weights_.encoder2_0_Wv_2.resize(num_encoder2_layers);
    transformer_weights_.encoder2_0_Wo_2.resize(num_encoder2_layers);

    transformer_weights_.encoder2_0_FF_W1_2.resize(num_encoder2_layers);
    transformer_weights_.encoder2_0_FF_b1_2.resize(num_encoder2_layers);
    transformer_weights_.encoder2_0_FF_W2_2.resize(num_encoder2_layers);
    transformer_weights_.encoder2_0_FF_b2_2.resize(num_encoder2_layers);



    
    // Load first encoder layers
    for (int layer = 0; layer < num_encoder1_layers; layer++) {
      
      std::string base = "encoders1." + std::to_string(layer);
      
      // Load attention weights
      std::vector<std::string> attn_weights = {".self_attn.Wq.weight", ".self_attn.Wk.weight", 
                                              ".self_attn.Wv.weight", ".self_attn.Wo.weight"};
      
      for (const auto& w : attn_weights) {
        
          auto it = param_map.find(base + w);
          if (it == param_map.end()) {
              logger_->error(utl::RSZ, 1024, "Missing encoder1 attention weight {} (loadWeights)", base + w);
              return;
          }
          
          const ParamData& pd = it->second;
          bool transpose = true;
          int target_rows = transpose ? pd.shape[1] : pd.shape[0];
          int target_cols = transpose ? pd.shape[0] : pd.shape[1];

          // Debugging
          //std::cout << "Loading " << base + w << " with shape [" << target_rows << "x" << target_cols << "]\n";
          //std::cout << "Data size: " << pd.values.size() << "\n";


          if (w.find("Wq") != std::string::npos) {
            transformer_weights_.encoder1_0_Wq_1[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder1_0_Wq_1[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1025, "Failed to fill encoder1 layer %d Wq: %s (loadWeights)", layer, errMsg);
              return;
            }
          } else if (w.find("Wk") != std::string::npos) {
            transformer_weights_.encoder1_0_Wk_1[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder1_0_Wk_1[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1026, "Failed to fill encoder1 layer %d Wk: %s (loadWeights)", layer, errMsg);
              return;
            }
          } else if (w.find("Wv") != std::string::npos) {
            transformer_weights_.encoder1_0_Wv_1[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder1_0_Wv_1[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1027, "Failed to fill encoder1 layer %d Wv: %s (loadWeights)", layer, errMsg);
              return;
            }
          } else if (w.find("Wo") != std::string::npos) {
            transformer_weights_.encoder1_0_Wo_1[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder1_0_Wo_1[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1028, "Failed to fill encoder1 layer %d Wo: %s (loadWeights)", layer, errMsg);
              return;
            }              
          }
        }
        
        // Load FF weights
        std::vector<std::string> ff_weights = {".ff.net.0.weight", ".ff.net.0.bias",
                                              ".ff.net.2.weight", ".ff.net.2.bias"};
        
        for (const auto& w : ff_weights) {
          
          auto it = param_map.find(base + w);
          if (it == param_map.end()) {
            logger_->error(utl::RSZ, 1029, "Missing encoder1 FF weight {} (loadWeights)", base + w);
            return;
          }

          const ParamData& pd = it->second;

          // Debugging
          //std::cout << "Loading " << base + w << " with shape [" << pd.shape[0] << "]\n";
          //std::cout << "Data size: " << pd.values.size() << "\n";
          
          if (w.find("net.0.weight") != std::string::npos) {
            bool transpose = true;
            int target_rows = transpose ? pd.shape[1] : pd.shape[0];
            int target_cols = transpose ? pd.shape[0] : pd.shape[1];
            transformer_weights_.encoder1_0_FF_W1_1[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder1_0_FF_W1_1[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1030, "Failed to fill encoder1 layer %d FF_W1: %s (loadWeights)", layer, errMsg);
              return;
            }
          } else if (w.find("net.0.bias") != std::string::npos) {
            transformer_weights_.encoder1_0_FF_b1_1[layer].resize(pd.shape[0]);
            if (!fillEigenVector(pd, transformer_weights_.encoder1_0_FF_b1_1[layer], errMsg)) {
              logger_->error(utl::RSZ, 1031, "Failed to fill encoder1 layer %d FF_b1: %s (loadWeights)", layer, errMsg);
              return;
            }
          } else if (w.find("net.2.weight") != std::string::npos) {
            bool transpose = true;
            int target_rows = transpose ? pd.shape[1] : pd.shape[0];
            int target_cols = transpose ? pd.shape[0] : pd.shape[1];
            transformer_weights_.encoder1_0_FF_W2_1[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder1_0_FF_W2_1[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1032, "Failed to fill encoder1 layer %d FF_W2: %s (loadWeights)", layer, errMsg);
              return;
            }
          } else if (w.find("net.2.bias") != std::string::npos) {
            transformer_weights_.encoder1_0_FF_b2_1[layer].resize(pd.shape[0]);
            if (!fillEigenVector(pd, transformer_weights_.encoder1_0_FF_b2_1[layer], errMsg)) {
              logger_->error(utl::RSZ, 1033, "Failed to fill encoder1 layer %d FF_b2: %s (loadWeights)", layer, errMsg);
              return;
            }
          }
        }
    }
    
    // Load second encoder layers
    for (int layer = 0; layer < num_encoder2_layers; layer++) {
      
      std::string base = "encoders2." + std::to_string(layer);
      
      // Load cross attention weights
      std::vector<std::string> attn_weights = {".cross_attn.Wq.weight", ".cross_attn.Wk.weight",
                                              ".cross_attn.Wv.weight", ".cross_attn.Wo.weight"};
      
      for (const auto& w : attn_weights) {
          auto it = param_map.find(base + w);
          if (it == param_map.end()) {
            logger_->error(utl::RSZ, 1034, "Missing encoder2 attention weight {} (loadWeights)", base + w);
            return;
          }

          const ParamData& pd = it->second;
          bool transpose = true;
          int target_rows = transpose ? pd.shape[1] : pd.shape[0];
          int target_cols = transpose ? pd.shape[0] : pd.shape[1];

          // Debugging
          //std::cout << "Loading " << base + w << " with shape [" << target_rows << "x" << target_cols << "]\n";
          //std::cout << "Data size: " << pd.values.size() << "\n";
          
          if (w.find("Wq") != std::string::npos) {
            transformer_weights_.encoder2_0_Wq_2[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder2_0_Wq_2[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1035, "Failed to fill encoder2 layer %d Wq: %s (loadWeights)", layer, errMsg);
              return;
            }
          } 
          else if (w.find("Wk") != std::string::npos) {
            transformer_weights_.encoder2_0_Wk_2[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder2_0_Wk_2[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1036, "Failed to fill encoder2 layer %d Wk: %s (loadWeights)", layer, errMsg);
              return;
            }
          } 
          else if (w.find("Wv") != std::string::npos) {
            transformer_weights_.encoder2_0_Wv_2[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder2_0_Wv_2[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1037, "Failed to fill encoder2 layer %d Wv: %s (loadWeights)", layer, errMsg);
              return;
            }
          } 
          else if (w.find("Wo") != std::string::npos) {
            transformer_weights_.encoder2_0_Wo_2[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder2_0_Wo_2[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1038, "Failed to fill encoder2 layer %d Wo: %s (loadWeights)", layer, errMsg);
              return;
            }
          }
        }
        
        // Load FF weights
        std::vector<std::string> ff_weights = {".ff.net.0.weight", ".ff.net.0.bias",
                                              ".ff.net.2.weight", ".ff.net.2.bias"};
        
        for (const auto& w : ff_weights) {
          auto it = param_map.find(base + w);
          if (it == param_map.end()) {
              logger_->error(utl::RSZ, 1039, "Missing encoder2 FF weight {} (loadWeights)", base + w);
              return;
          }

          const ParamData& pd = it->second;

          // Debugging
          //std::cout << "Loading " << base + w << " with shape [" << pd.shape[0] << "]\n";
          //std::cout << "Data size: " << pd.values.size() << "\n";
          
          if (w.find("net.0.weight") != std::string::npos) {
            bool transpose = true;
            int target_rows = transpose ? pd.shape[1] : pd.shape[0];
            int target_cols = transpose ? pd.shape[0] : pd.shape[1];
            transformer_weights_.encoder2_0_FF_W1_2[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder2_0_FF_W1_2[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1040, "Failed to fill encoder2 layer %d FF_W1: %s (loadWeights)", layer, errMsg);
              return;
            }
          } 
          else if (w.find("net.0.bias") != std::string::npos) {
            transformer_weights_.encoder2_0_FF_b1_2[layer].resize(pd.shape[0]);
            if (!fillEigenVector(pd, transformer_weights_.encoder2_0_FF_b1_2[layer], errMsg)) {
              logger_->error(utl::RSZ, 1041, "Failed to fill encoder2 layer %d FF_b1: %s (loadWeights)", layer, errMsg);
              return;
            }
          } 
          else if (w.find("net.2.weight") != std::string::npos) {
            bool transpose = true;
            int target_rows = transpose ? pd.shape[1] : pd.shape[0];
            int target_cols = transpose ? pd.shape[0] : pd.shape[1];
            transformer_weights_.encoder2_0_FF_W2_2[layer].resize(target_rows, target_cols);
            if (!fillEigenMatrix(pd, transformer_weights_.encoder2_0_FF_W2_2[layer], transpose, errMsg)) {
              logger_->error(utl::RSZ, 1042, "Failed to fill encoder2 layer %d FF_W2: %s (loadWeights)", layer, errMsg);
              return;
            }
          } 
          else if (w.find("net.2.bias") != std::string::npos) {
            transformer_weights_.encoder2_0_FF_b2_2[layer].resize(pd.shape[0]);
            if (!fillEigenVector(pd, transformer_weights_.encoder2_0_FF_b2_2[layer], errMsg)) {
              logger_->error(utl::RSZ, 1043, "Failed to fill encoder2 layer %d FF_b2: %s (loadWeights)", layer, errMsg);
              return;
            }
          }
        }
    }
    // Successfully loaded all weights
    transformer_weights_.loaded = true;
}

// Print Model Summary for Loaded Weights
void MLGateSizer::printModelSummary(const std::unordered_map<std::string, ParamData>& param_map) 
{
  // Data extraction variables
  int max_encoder1_layer = -1;
  int max_encoder2_layer = -1;
  size_t D_model = 0;
  size_t FF_hidden_dim = 0;
  size_t D_in = 0;
  size_t D_out = 0;
  size_t D_emb = 0;
  bool dims_initialized = false;

  // Temporary storage
  size_t total_params = 0;
  std::unordered_map<std::string, size_t> projection_params;
  size_t encoder1_params = 0;
  size_t encoder2_params = 0;

  // First pass: Extract dimensions
  for (const auto& [name, pd] : param_map) {
    // Extract encoder layers
    if (name.compare(0, 10, "encoders1.") == 0) {
      size_t dot_pos = name.find('.', 10);
      if (dot_pos != std::string::npos) {
        int layer = std::stoi(name.substr(10, dot_pos-10));
        max_encoder1_layer = std::max(max_encoder1_layer, layer);
      }
    }
    else if (name.compare(0, 10, "encoders2.") == 0) {
      size_t dot_pos = name.find('.', 10);
      if (dot_pos != std::string::npos) {
        int layer = std::stoi(name.substr(10, dot_pos-10));
        max_encoder2_layer = std::max(max_encoder2_layer, layer);
      }
    }

    // Extract model dimensions
    if (!dims_initialized) {
      if (name.find("Wq.weight") != std::string::npos && pd.shape.size() == 2) {
        D_model = pd.shape[0];
      }
      else if (name.find("net.0.weight") != std::string::npos && pd.shape.size() == 2) {
        FF_hidden_dim = pd.shape[0];
      }
      else if (name == "proj_in1.weight" && pd.shape.size() == 2) {
        D_in = pd.shape[1];
      }
      else if (name == "proj_in2.weight" && pd.shape.size() == 2) {
        D_emb = pd.shape[1];
      }
      else if (name == "proj_out.weight" && pd.shape.size() == 2) {
        D_out = pd.shape[0];
      }
      
      if (D_model && FF_hidden_dim && D_in && D_emb && D_out) {
        dims_initialized = true;
      }
    }
  }

  // Calculate layer counts
  const int num_encoder1 = max_encoder1_layer + 1;
  const int num_encoder2 = max_encoder2_layer + 1;

  // Second pass: count parameters
  for (const auto& [name, pd] : param_map) {
    const size_t params = std::accumulate(
      pd.shape.begin(), pd.shape.end(), 1, 
      [](size_t a, size_t b) { return a * b; }
    );
    
    if (name.find("proj") != std::string::npos) {
      projection_params[name] = params;
    }
    else if (name.find("encoders1") != std::string::npos) {
      encoder1_params += params;
    }
    else if (name.find("encoders2") != std::string::npos) {
      encoder2_params += params;
    }
    total_params += params;
  }

  // Calculate component params
  const size_t attn_per_layer = 4 * D_model * D_model;
  const size_t ffw_per_layer = D_model * FF_hidden_dim * 2;
  const size_t ffb_per_layer = FF_hidden_dim + D_model;

  // Print summary
  std::cout << "\nModel Summary                 \n";
  std::cout << "==============================\n";
  std::cout << "Dimensions:\n";
  std::cout << "  D_model:       " << D_model << "\n";
  std::cout << "  FF_hidden_dim: " << FF_hidden_dim << "\n";
  std::cout << "  D_in:          " << D_in << "\n";
  std::cout << "  D_out:         " << D_out << "\n";
  std::cout << "  D_emb:         " << D_emb << "\n";
  std::cout << "  Encoder1:      " << num_encoder1 << " layers\n";
  std::cout << "  Encoder2:      " << num_encoder2 << " layers\n";
  std::cout << "==========================\n\n";

  std::cout << "1. Projections\n";
  std::cout << "==========================\n";
  std::cout << std::setw(25) << std::left << "Name" 
            << std::setw(15) << "Shape" 
            << "Params\n";
  std::cout << "==========================\n";
  for (const auto& [name, params] : projection_params) {
    const auto& shape = param_map.at(name).shape;
    std::string shape_str;
    for (size_t i = 0; i < shape.size(); ++i) {
      shape_str += std::to_string(shape[i]);
      if (i != shape.size()-1) shape_str += "×";
    }
    std::cout << std::setw(25) << name 
              << std::setw(15) << ("[" + shape_str + "]")
              << params << "\n";
  }
  std::cout << "Total Projections: " 
            << std::accumulate(projection_params.begin(), projection_params.end(), 0, 
                             [](auto sum, const auto& p) { return sum + p.second; })
            << "\n";
  std::cout << "==========================\n\n";
  std::cout << "2. Encoder1\n";
  std::cout << "==========================\n";
  std::cout << "Component       Per Layer    Total\n";
  std::cout << "==========================\n";
  std::cout << "Attention       " 
            << std::setw(11) << attn_per_layer 
            << std::setw(11) << attn_per_layer * num_encoder1 << "\n";
  std::cout << "FF Weights      " 
            << std::setw(11) << ffw_per_layer 
            << std::setw(11) << ffw_per_layer * num_encoder1 << "\n";
  std::cout << "FF Biases       " 
            << std::setw(11) << ffb_per_layer 
            << std::setw(11) << ffb_per_layer * num_encoder1 << "\n";
  std::cout << "Actual Total:   " 
            << std::setw(22) << encoder1_params << "\n";
  std::cout << "==========================\n\n";

  std::cout << "3. Encoder2\n";
  std::cout << "==========================\n";
  std::cout << "Component       Per Layer    Total\n";
  std::cout << "==========================\n";
  std::cout << "Attention       " 
            << std::setw(11) << attn_per_layer 
            << std::setw(11) << attn_per_layer * num_encoder2 << "\n";
  std::cout << "FF Weights      " 
            << std::setw(11) << ffw_per_layer 
            << std::setw(11) << ffw_per_layer * num_encoder2 << "\n";
  std::cout << "FF Biases       " 
            << std::setw(11) << ffb_per_layer 
            << std::setw(11) << ffb_per_layer * num_encoder2 << "\n";
  std::cout << "Actual Total:   " 
            << std::setw(22) << encoder2_params << "\n";
  std::cout << "==========================\n\n";

  std::cout << "Grand Total\n";
  std::cout << "==========================\n";
  std::cout << "Projections:  " << std::setw(10) << projection_params.begin()->second << "\n";
  std::cout << "Encoder1:     " << std::setw(10) << encoder1_params << "\n";
  std::cout << "Encoder2:     " << std::setw(10) << encoder2_params << "\n";
  std::cout << "-----------------------------\n";
  std::cout << "Total:        " << std::setw(10) << total_params << "\n";
  std::cout << "==========================\n\n";
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

static void applyRoPE(std::vector<std::vector<float>>& Q_or_K, float base = 10000.0f) {
    // Q_or_K has shape [L x Dh]. We rotate pairs: (0,1), (2,3), ...
    const int L = Q_or_K.size();
    if (L == 0) return;
    
    const int Dh = Q_or_K[0].size();
    assert(Dh % 2 == 0 && "Dimension Dh must be even for RoPE");
    const int half = Dh / 2;

    for (int i = 0; i < L; i++) {
        const float pos = static_cast<float>(i);
        std::vector<float>& row = Q_or_K[i];
        
        for (int j = 0; j < half; j++) {
            const float theta = std::pow(base, -2.0f * j / static_cast<float>(Dh));
            const float angle = pos * theta;
            
            const float cos_theta = std::cos(angle);
            const float sin_theta = std::sin(angle);
            
            const int idx0 = 2 * j;
            const int idx1 = 2 * j + 1;
            
            const float x0 = row[idx0];
            const float x1 = row[idx1];
            
            // Apply rotation
            row[idx0] = x0 * cos_theta - x1 * sin_theta;
            row[idx1] = x1 * cos_theta + x0 * sin_theta;
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

  //std::vector<std::vector<float>> Wq = randomMatrix(D, D);
  //std::vector<std::vector<float>> Wk = randomMatrix(D, D);
  //std::vector<std::vector<float>> Wv = randomMatrix(D, D);

  // 1) Flatten seq from [L x D] to a form we can do matMul with:
  // For matMul, we treat it as [L x D].
  // Q = seq * Wq => shape [L x D]
  // K = seq * Wk => shape [L x D]
  // V = seq * Wv => shape [L x D]

  // We'll define a small helper to do matMul (LxD) * (DxD) => (LxD)
  auto Q = matMul(seq, Wq);
  auto K = matMul(seq, Wk);
  auto V = matMul(seq, Wv);

  // 2) Reshape Q,K,V into [L x H x Dh], do attention per head
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

  // Apply RoPE to each head's Q and K for positional encoding
  for (size_t h = 0; h <static_cast<size_t> (num_heads); h++) {
    applyRoPE(Q3[h]);
    applyRoPE(K3[h]);
  }

  // 3) Self-attention per head
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

  // 4) Concatenate heads back => shape [L x D]
  std::vector<std::vector<float>> outSeq(L, std::vector<float>(D, 0.0f));
  for (size_t l_ = 0; l_ < L; l_++) {
    for (size_t h_ = 0; h_ < (size_t)num_heads; h_++) {
      for (size_t d_ = 0; d_ < Dh; d_++) {
        outSeq[l_][h_ * Dh + d_] = outHeads[h_][l_][d_];
      }
    }
  }

  // 5) Output projection
  // We'll define Wo: [D x D] for simplicity
  //std::vector<std::vector<float>> Wo = randomMatrix(D, D);

  // outSeq * Wo: [L x D] * [D x D] => [L x D]
  auto finalOut = matMul(outSeq, Wo);

  return finalOut;
}

// ------------------------------------------------------------
// Helper function: naive multi-head cross-attention for one sequence
// We'll do a single pass of multi-head cross-attention for a shape [LQ x D].
// Q_seq: [LQ x D]
// M_seq: [LM x D] (sequence from 1st encoder used for K, V)
// Heads: H
// We split D into H heads each of size Dh = D/H.
// For simplicity, we do random Wq, Wk, Wv. Real code would load from TransSizer.
// Model loading has not been implemented yet.
// We do the standard formula: softmax(Q*K^T / sqrt(Dh)) * V
// Then we re-concatenate the heads
// ------------------------------------------------------------
static std::vector<std::vector<float>> multiHeadCrossAttention(
    const std::vector<std::vector<float>>& Q_seq,
    const std::vector<std::vector<float>>& M_seq,
    int num_heads,
    size_t D,
    size_t LQ,
    size_t LM,
    std::vector<std::vector<float>>& Wq,
    std::vector<std::vector<float>>& Wk,
    std::vector<std::vector<float>>& Wv,
    std::vector<std::vector<float>>& Wo)
{
    // Check dimensions
    assert(D % num_heads == 0 && "D must be divisible by num_heads");
    size_t Dh = D / num_heads;

    // Project Q, K, V
    auto Q = matMul(Q_seq, Wq); // [LQ x D]
    auto K = matMul(M_seq, Wk); // [LM x D]
    auto V = matMul(M_seq, Wv); // [LM x D]

    // Split into heads for Q, K, V
    auto to3D = [&](const std::vector<std::vector<float>>& X, size_t L) {
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

    auto Q3 = to3D(Q, LQ);
    auto K3 = to3D(K, LM);
    auto V3 = to3D(V, LM);

    // Compute attention per head
    auto softmax = [&](std::vector<float>& logits) {
        float max_val = *std::max_element(logits.begin(), logits.end());
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
        num_heads, std::vector<std::vector<float>>(LQ, std::vector<float>(Dh, 0.0f)));

    for (size_t h = 0; h < (size_t)num_heads; h++) {
        for (size_t lq = 0; lq < LQ; lq++) {
            std::vector<float> attn_logits(LM, 0.0f);
            for (size_t lm = 0; lm < LM; lm++) {
                float dot_val = 0.0f;
                for (size_t d = 0; d < Dh; d++) {
                    dot_val += Q3[h][lq][d] * K3[h][lm][d];
                }
                attn_logits[lm] = dot_val / std::sqrt((float)Dh);
            }
            softmax(attn_logits);

            for (size_t lm = 0; lm < LM; lm++) {
                for (size_t d = 0; d < Dh; d++) {
                    outHeads[h][lq][d] += attn_logits[lm] * V3[h][lm][d];
                }
            }
        }
    }

    // Concatenate heads
    std::vector<std::vector<float>> outSeq(LQ, std::vector<float>(D, 0.0f));
    for (size_t l_ = 0; l_ < LQ; l_++) {
        for (size_t h_ = 0; h_ < (size_t)num_heads; h_++) {
            for (size_t d_ = 0; d_ < Dh; d_++) {
                outSeq[l_][h_ * Dh + d_] = outHeads[h_][l_][d_];
            }
        }
    }

    // Output projection
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
    const std::vector<std::vector<std::vector<float>>>& data_array_1,
    const std::vector<std::vector<std::vector<float>>>& data_array_2,
    int num_heads,
    size_t N,
    size_t L,
    size_t D_in,
    size_t D_out,
    size_t D_emb,
    size_t D_model,
    size_t FF_hidden_dim,
    int num_encoder_layers,
    int num_encoder_layers_2)
{

  size_t L2 = L / 2;

  if (N == 0 || L == 0 || data_array_1[0].size() != L || data_array_2[0].size() != L2) {
      return {};
  }

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
  // W_out: [D_model x D_out]
  auto W_in  = randomMatrix(D_in, D_model);
  auto W_in2 = randomMatrix(D_emb, D_model);
  auto W_out = randomMatrix(D_model, D_out);

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
      std::vector<std::vector<float>>(L2, std::vector<float>(D_out, 0.f)));

  // For each sequence in the batch:
  for (size_t n = 0; n < N; n++) {
    // current seq: [L x D_in]
    std::vector<std::vector<float>> seq = data_array_1[n];

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

    // Process data_array_2 through cross-attention layers
    const std::vector<std::vector<float>>& seq2 = data_array_2[n];
    std::vector<std::vector<float>> seq_proj2 = matMul(seq2, W_in2); // shape [L2 x D_model]
    std::vector<std::vector<float>> x2 = seq_proj2; // shape [L2 x D_model]

    for (int layer_2 = 0; layer_2 < num_encoder_layers_2; layer_2++) {
        auto attn_out2 = multiHeadCrossAttention(x2, x, num_heads, D_model, L2, L, Wq, Wk, Wv, Wo);

        // Residual and LayerNorm
        for (size_t l = 0; l < L2; l++) {
            for (size_t d = 0; d < D_model; d++) {
                attn_out2[l][d] += x2[l][d];
            }
        }
        layerNorm(attn_out2);

        auto ff_out2 = feedForward(attn_out2, L2, D_model, FF_hidden_dim, FF_W1, FF_W2, FF_b1, FF_b2);

        for (size_t l = 0; l < L2; l++) {
            for (size_t d = 0; d < D_model; d++) {
                ff_out2[l][d] += attn_out2[l][d];
            }
        }
        layerNorm(ff_out2);
        x2 = ff_out2;
    }

    // x2 is now [L/2 x D_model]. Project back to [L/2 x D_out]
    auto final_out = matMul(x2, W_out);

    // Apply softmax to each row like a classification head same as runTransformerEigen
    for (size_t l = 0; l < L2; l++) {
        float max_val = *std::max_element(final_out[l].begin(), final_out[l].end());
        float sum = 0.0f;
        for (auto& v : final_out[l]) {
            v = std::exp(v - max_val);
            sum += v;
        }
        for (auto& v : final_out[l]) {
            v /= sum;
        }
    }

    // Store final_out in output[n]
    output[n] = final_out; // shape [L/2 x D_out]
  }

  // Return shape: [N x L/2 x D_out]
  return output;
}

static void applyRoPEEigen(Eigen::MatrixXf & Q_or_K, float base = 10000.0f)
{
    // Q_or_K has shape [L x Dh].  We rotate pairs: (0,1), (2,3), ...
    int L = Q_or_K.rows();
    int Dh = Q_or_K.cols();

    // For each position i in [0..L-1]:
    for (int i = 0; i < L; i++) {
        float pos = float(i);  // position index
        // For each pair j in [0..Dh/2-1]:
        int half = Dh / 2;  // must be even dimension for pairs
        for (int j = 0; j < half; j++) {
            // 2*j and 2*j+1 form a pair
            float x0 = Q_or_K(i, 2*j);
            float x1 = Q_or_K(i, 2*j+1);

            // Compute angle = pos * theta, where:
            //   theta = base^(-2*j / float(Dh))
            // This is a typical RoPE frequency progression.
            float theta = std::pow(base, -2.0f * j / float(Dh));
            float angle = pos * theta;

            float c = std::cos(angle);
            float s = std::sin(angle);

            // Apply 2D rotation:
            //   x0' =  x0 * cos(angle) - x1 * sin(angle)
            //   x1' =  x1 * cos(angle) + x0 * sin(angle)
            float rx0 = x0 * c - x1 * s;
            float rx1 = x1 * c + x0 * s;

            // Store back
            Q_or_K(i, 2*j)   = rx0;
            Q_or_K(i, 2*j+1) = rx1;
        }
    }
}


static Eigen::MatrixXf eigenSelfAttention(
  const Eigen::MatrixXf& seq, 
  int num_heads,
  size_t D,
  size_t L,
  const Eigen::MatrixXf& Wq,
  const Eigen::MatrixXf& Wk,
  const Eigen::MatrixXf& Wv,
  const Eigen::MatrixXf& Wo)
{
  // seq shape: [L x D]
  // Return shape: [L x D]
  //const int L = seq.rows();
  //const int D = seq.cols();

  assert(D % num_heads == 0 && "D must be divisible by num_heads");
  size_t Dh = D / num_heads;



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

    // Apply RoPE to Qsplit[h] and Ksplit[h] for position encoding
    applyRoPEEigen(Qsplit[h]); // shape [L x Dh] modified in place
    applyRoPEEigen(Ksplit[h]); // shape [L x Dh] modified in place
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

static Eigen::MatrixXf eigenCrossAttention(
    const Eigen::MatrixXf & Qseq, // shape [LQ x D]
    const Eigen::MatrixXf & Mseq, // shape [LM x D] used for K,V
    int num_heads,
    size_t D,
    size_t LQ,
    size_t LM,
    const Eigen::MatrixXf & Wq,
    const Eigen::MatrixXf & Wk,
    const Eigen::MatrixXf & Wv,
    const Eigen::MatrixXf & Wo)
{
    // 1) Project Q from Qseq, K and V from Mseq:
    //    Qseq: [LQ x D], Mseq: [LM x D]
    //    Wq, Wk, Wv: [D x D]
    Eigen::MatrixXf Q = Qseq * Wq; // [LQ x D]
    Eigen::MatrixXf K = Mseq * Wk; // [LM x D]
    Eigen::MatrixXf V = Mseq * Wv; // [LM x D]

    // 2) Split into heads.  Let Dh = D / num_heads.
    //    We gather them into Qsplit[h] shape [LQ x Dh], Ksplit[h] shape [LM x Dh], etc.
    size_t Dh = D / num_heads;
    std::vector<Eigen::MatrixXf> Qsplit(num_heads), Ksplit(num_heads), Vsplit(num_heads);
    for (int h = 0; h < num_heads; h++)
    {
        Qsplit[h].resize(LQ, Dh);
        Ksplit[h].resize(LM, Dh);
        Vsplit[h].resize(LM, Dh);

        for (int lq = 0; lq < (int)LQ; lq++) {
            for (int d_ = 0; d_ < (int)Dh; d_++) {
                Qsplit[h](lq, d_) = Q(lq, h*Dh + d_);
            }
        }
        for (int lm = 0; lm < (int)LM; lm++) {
            for (int d_ = 0; d_ < (int)Dh; d_++) {
                Ksplit[h](lm, d_) = K(lm, h*Dh + d_);
                Vsplit[h](lm, d_) = V(lm, h*Dh + d_);
            }
        }

        // (Optional) If we RoPE for cross‐attention as well,
        // you can apply it to Qsplit[h] and Ksplit[h] here.
        // Undecided if this is necessary for cross-attention.
        // applyRoPE(Qsplit[h]);
        // applyRoPE(Ksplit[h]);
    }

    // 3) Perform attention for each head
    auto softmax = [&](Eigen::VectorXf & logits) {
        float max_val = logits.maxCoeff();
        logits = (logits.array() - max_val).exp();
        float sum_ = logits.sum();
        logits /= sum_;
        return logits;
    };

    std::vector<Eigen::MatrixXf> outHeads(num_heads, Eigen::MatrixXf::Zero(LQ, Dh));
    for (int h = 0; h < num_heads; h++)
    {
        for (int lq = 0; lq < (int)LQ; lq++) {
            // Dot Q[lq] with each K[lm]
            Eigen::VectorXf attn_logits(LM);
            for (int lm = 0; lm < (int)LM; lm++) {
                float dot_val = Qsplit[h].row(lq).dot(Ksplit[h].row(lm));
                attn_logits(lm) = dot_val / std::sqrt((float)Dh);
            }
            // Softmax
            softmax(attn_logits);

            // Weighted sum of V
            for (int lm = 0; lm < (int)LM; lm++) {
                outHeads[h].row(lq).noalias() += attn_logits(lm) * Vsplit[h].row(lm);
            }
        }
    }

    // 4) Concatenate heads => [LQ x D]
    Eigen::MatrixXf out(LQ, D);
    for (int lq = 0; lq < (int)LQ; lq++) {
        for (int h = 0; h < num_heads; h++) {
            for (int d_ = 0; d_ < (int)Dh; d_++) {
                out(lq, h*Dh + d_) = outHeads[h](lq, d_);
            }
        }
    }

    // 5) Apply final projection Wo => [D x D]
    Eigen::MatrixXf finalOut = out * Wo; // shape [LQ x D]
    return finalOut;
}


static Eigen::MatrixXf eigenFF(
  const Eigen::MatrixXf& seq,
  size_t L,
  size_t D,
  size_t H,
  const Eigen::MatrixXf& W1,
  const Eigen::MatrixXf& W2,
  const Eigen::VectorXf& b1,
  const Eigen::VectorXf& b2)
{
  // seq: L x D, do 2D hidden layer
  //int L = seq.rows();
  //int D = seq.cols();
  //int H = 2*D;

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
// data_array_1: shape [N x L x D_in]
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
// RoPE is applied to Q/K in self-attention.
// 
// Also implement encoder-encoder architecture where the 2nd sequence consisting of libcell types are given.
// data_array_2: shape [N x L2 x D_emb], L2 is the length of the libcell sequence, L2 = L / 2.
// As L is the pin sequence length, L2 is the number of libcells which have a 2:1 ratio.
// The encoder output from data_array_1
// Uses Eigen for matrix operations.
//
// Currently uses output projection where Z : [L x D_model] is projected back to [L x D_in]
// Z : [L x D_model] * W_out [D_model x D_in]=> [L x D_in]
// For Gate sizing need to take encoder output and pass it through a classification head.
// Each class corresponding to the libcell.
// --------------------------------------------------------------------
std::vector<std::vector<std::vector<float>>> MLGateSizer::runTransformerEigen(  // Random weights
    const std::vector<std::vector<std::vector<float>>>& data_array_1, // shape [N x L   x D_in]
    const std::vector<std::vector<std::vector<float>>>& data_array_2, // shape [N x L/2 x D_emb]
    int num_heads,
    size_t N,
    size_t L,
    size_t D_in,
    size_t D_out,
    size_t D_emb,
    size_t D_model,
    size_t FF_hidden_dim,
    int num_encoder_layers,
    int num_encoder_layers_2)
{
  //size_t N = data_array_1.size();
  if (N == 0) return {};

  //size_t L = data_array_1[0].size();
  if (L == 0) return {};

  size_t L2 = L / 2;

  // Original input dimension
  //size_t D_in = data_array_1[0][0].size();

  // Decide how many heads we want
  //const int num_heads = 4;

  // Decide a model dimension that is divisible by num_heads
  //size_t D_model = 64;

  // Number of encoder layers
  //const int num_encoder_layers = 2;

  // Check if D_model is divisible by num_heads
  assert(D_model % num_heads == 0 && "D_model must be divisible by num_heads");


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
  // W_out: [D_model x D_out]
  auto W_in  = randomMatrix(D_in, D_model);
  auto W_in2 = randomMatrix(D_emb, D_model);
  auto W_out = randomMatrix(D_model, D_out);

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
      std::vector<std::vector<float>>(L2, std::vector<float>(D_out, 0.f)));
  
  
  // For each sequence in the batch:
  for (size_t n = 0; n < N; n++) {
    // current seq: [L x D_in]
    Eigen::MatrixXf seq(L, D_in);
    for (size_t l = 0; l < L; l++) {
      for (size_t d = 0; d < D_in; d++) {
        seq(l, d) = data_array_1[n][l][d];
      }
    }

    // First encoder: self-attention
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

    // Second encoder: cross-attention
    // 1) Project seq => [L/2 x D_model]
    Eigen::MatrixXf seq_proj2(L2, D_emb);
    for (size_t l = 0; l < L2; l++) {
      for (size_t d = 0; d < D_emb; d++) {
        seq_proj2(l, d) = data_array_2[n][l][d];
      }
    }

    // [L/2 x D_emb] x [D_emb x D_model] => [L/2 x D_model]
    Eigen::MatrixXf x2 = seq_proj2 * W_in2;
    for (int layer_2 = 0; layer_2 < num_encoder_layers_2; layer_2++) {
      // (a) Multi-head cross-attention
      auto attn_out2 = eigenCrossAttention(x2, x, num_heads, D_model, L2, L, Wq, Wk, Wv, Wo);
      // (b) Residual
      x2 += attn_out2;
      // (c) LayerNorm
      eigenLayerNorm(x2);

      // (d) FeedForward
      auto ff_out2 = eigenFF(x2, L2, D_model, FF_hidden_dim, FF_W1, FF_W2, FF_b1, FF_b2);
      // (e) Residual
      x2 += ff_out2;
      // (f) LayerNorm
      eigenLayerNorm(x2);
    }




    // x is now [L/2 x D_model]. Project back to [L/2 x D_out]
    Eigen::MatrixXf final_out = x2 * W_out;

    // Conduct softmax over D_out dimension to get class probabilities
    for (size_t l = 0; l < L2; l++) {
      float max_val = final_out.row(l).maxCoeff();
      final_out.row(l) = (final_out.row(l).array() - max_val).exp();
      float sum_ = final_out.row(l).sum();
      final_out.row(l) /= sum_;
    }

    // Store final_out in output[n]
    for (size_t l = 0; l < L2; l++) {
      for (size_t d = 0; d < D_out; d++) {
        output[n][l][d] = final_out(l, d);
      }
    }
  }


  // Return shape: [N x L/2 x D_out]
  return output;
}


std::vector<std::vector<std::vector<float>>> MLGateSizer::runTransformerEigen(  // Uses given weights
    const std::vector<std::vector<std::vector<float>>>& data_array_1, // shape [N x L   x D_in]
    const std::vector<std::vector<std::vector<float>>>& data_array_2, // shape [N x L/2 x D_emb]
    int num_heads,
    size_t N,
    size_t L,
    size_t D_in,
    size_t D_out,
    size_t D_emb,
    size_t D_model,
    size_t FF_hidden_dim,
    int num_encoder_layers,
    int num_encoder_layers_2,
    const TransformerWeights& weights)
{
  //size_t N = data_array_1.size();
  if (N == 0) return {};

  //size_t L = data_array_1[0].size();
  if (L == 0) return {};

  size_t L2 = L / 2;

  // Original input dimension
  //size_t D_in = data_array_1[0][0].size();

  // Decide how many heads we want
  //const int num_heads = 4;

  // Decide a model dimension that is divisible by num_heads
  //size_t D_model = 64;

  // Number of encoder layers
  //const int num_encoder_layers = 2;

  // Check if D_model is divisible by num_heads
  assert(D_model % num_heads == 0 && "D_model must be divisible by num_heads");

  //std::cout << "Debug 1" << std::endl;
  // Allocate output (N x L/2 x D_out) or (N x L/2 x num_classes)
  std::vector<std::vector<std::vector<float>>> output(N,
      std::vector<std::vector<float>>(L2, std::vector<float>(D_out, 0.f)));
  
  //std::cout << "Debug 2" << std::endl;
  // For each sequence in the batch:
  for (size_t n = 0; n < N; n++) {
    // current seq: [L x D_in]
    Eigen::MatrixXf seq(L, D_in);
    for (size_t l = 0; l < L; l++) {
      for (size_t d = 0; d < D_in; d++) {
        seq(l, d) = data_array_1[n][l][d];
      }
    }
    //std::cout << "Debug 3" << std::endl;
    // First encoder: self-attention
    // 1) Project seq => [L x D_model]
    Eigen::MatrixXf seq_proj = seq * weights.W_in1;
    //std::cout << "Debug 4" << std::endl;
    // 2) Run through encoder blocks
    Eigen::MatrixXf x = seq_proj; // shape [L x D_model]
    for (int layer = 0; layer < num_encoder_layers; layer++) {
      // (a) Multi-head self-attention
      auto attn_out = eigenSelfAttention(x, num_heads, D_model, L, weights.encoder1_0_Wq_1[layer], weights.encoder1_0_Wk_1[layer], weights.encoder1_0_Wv_1[layer], weights.encoder1_0_Wo_1[layer]);
      // (b) Residual
      x += attn_out;
      // (c) LayerNorm
      eigenLayerNorm(x);

      // (d) FeedForward
      auto ff_out = eigenFF(x, L, D_model, FF_hidden_dim, weights.encoder1_0_FF_W1_1[layer], weights.encoder1_0_FF_W2_1[layer], weights.encoder1_0_FF_b1_1[layer], weights.encoder1_0_FF_b2_1[layer]);
      // (e) Residual
      x += ff_out;
      // (f) LayerNorm
      eigenLayerNorm(x);
    }
    //std::cout << "Debug 5" << std::endl;
    // Second encoder: cross-attention
    // 1) Project seq => [L/2 x D_model]
    Eigen::MatrixXf seq_proj2(L2, D_emb);
    for (size_t l = 0; l < L2; l++) {
      for (size_t d = 0; d < D_emb; d++) {
        seq_proj2(l, d) = data_array_2[n][l][d];
      }
    }
    //std::cout << "Debug 6" << std::endl;
    // [L/2 x D_emb] x [D_emb x D_model] => [L/2 x D_model]
    Eigen::MatrixXf x2 = seq_proj2 * weights.W_in2;
    for (int layer_2 = 0; layer_2 < num_encoder_layers_2; layer_2++) {
      // (a) Multi-head cross-attention
      auto attn_out2 = eigenCrossAttention(x2, x, num_heads, D_model, L2, L, weights.encoder2_0_Wq_2[layer_2], weights.encoder2_0_Wk_2[layer_2], weights.encoder2_0_Wv_2[layer_2], weights.encoder2_0_Wo_2[layer_2]);
      // (b) Residual
      x2 += attn_out2;
      // (c) LayerNorm
      eigenLayerNorm(x2);

      // (d) FeedForward
      auto ff_out2 = eigenFF(x2, L2, D_model, FF_hidden_dim, weights.encoder2_0_FF_W1_2[layer_2], weights.encoder2_0_FF_W2_2[layer_2], weights.encoder2_0_FF_b1_2[layer_2], weights.encoder2_0_FF_b2_2[layer_2]);
      // (e) Residual
      x2 += ff_out2;
      // (f) LayerNorm
      eigenLayerNorm(x2);
    }



    //std::cout << "Debug 7" << std::endl;
    // x is now [L/2 x D_model]. Project back to [L/2 x D_out] or [L/2 x num_classes]
    Eigen::MatrixXf final_out = x2 * weights.W_out;

    // Conduct softmax over D_out dimension to get class probabilities
    for (size_t l = 0; l < L2; l++) {
      float max_val = final_out.row(l).maxCoeff();
      final_out.row(l) = (final_out.row(l).array() - max_val).exp();
      float sum_ = final_out.row(l).sum();
      final_out.row(l) /= sum_;
    }

    // Store final_out in output[n]
    for (size_t l = 0; l < L2; l++) {
      for (size_t d = 0; d < D_out; d++) {
        output[n][l][d] = final_out(l, d);
      }
    }
  }


  // Return shape: [N x L/2 x D_out]/ [N x L/2 x num_classes]
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
  float total_diff = 0.0f;
  if (A.size() != B.size()) return false;
  for (size_t n = 0; n < A.size(); n++) {
    if (A[n].size() != B[n].size()) return false;
    for (size_t l = 0; l < A[n].size(); l++) {
      if (A[n][l].size() != B[n][l].size()) return false;
      for (size_t d = 0; d < A[n][l].size(); d++) {
        float diff = std::fabs(A[n][l][d] - B[n][l][d]);
        total_diff += diff;
        if (diff > tol) {
          return false;
        }
      }
    }
  }

  std::cout << "Success! Total diff less than " << tol << std::endl;
  std::cout << "Total diff: " << total_diff << std::endl;
  std::cout << "Average diff: " << total_diff / (A.size() * A[0].size() * A[0][0].size()) << std::endl;
  return true;
}


}  // namespace rsz
