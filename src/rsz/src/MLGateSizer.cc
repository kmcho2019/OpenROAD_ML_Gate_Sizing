#include "MLGateSizer.hh"
#include <iostream>

namespace rsz {

MLGateSizer::MLGateSizer()
{
  // Initialize model, load weights, etc.
}

void MLGateSizer::loadWeights(const std::string& weight_file)
{
  // Load transformer model weights from PyTorch model and store in Eigen matrix
  // Example: transformer_weights_ = LoadFromPyTorch(weight_file);
}

void MLGateSizer::addToken(const std::vector<float>& pin_data, const std::string& gate_type)
{
  pin_tokens_.push_back(pin_data);
  gate_types_.push_back(gate_type);
}

void MLGateSizer::classifyAndResize()
{
  // Extract worst paths from the design
  extractWorstPaths();
  // Transform pin data into the format expected by the transformer model
  applyTransformerModel();
  
  // Once classification is done, resize gates based on predictions
  resizeGates();
}

void MLGateSizer::extractWorstPaths()
{
  // Example logic to extract worst timing paths using STA
  VertexSeq worst_paths;
  sta_->worstPaths(3, /*per-endpoint*/ true, worst_paths);

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

void MLGateSizer::transformPinDataToTokens(const Pin* pin, float slew, float cap, float slack, const std::string& gate_type)
{
  // Turn the extracted pin data into tokens that can be used by the transformer model
  std::vector<float> tokenized_data = {slew, cap, slack};
  
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
  // Implement feedforward network using Eigen operations and loaded transformer weights
  // This could include positional encoding, multihead attention, and the like
  Eigen::MatrixXf output = input * transformer_weights_; // Example of matrix multiplication
  // Apply activation functions, attention mechanisms, etc.
}

void MLGateSizer::resizewithML()
{
  // Based on the output from the transformer model, resize gates
  for (size_t i = 0; i < gate_types_.size(); ++i) {
    std::string gate_type = gate_types_[i];
    // Output the classification result (e.g., gate size)
    std::cout << "Resizing gate of type: " << gate_type << " to predicted size" << std::endl;
    // Resize the gate accordingly in the design
  }
}

void MLGateSizer::getEndpointAndCriticalPaths()
{
  // Retrieve endpoints and critical paths for debugging or further analysis

  Sta* sta = getSta();
  Network* network = sta->network();
  Graph* graph = sta->ensureGraph();

  // Retrieve endpoints
  PinSet endpoints = sta->endpoints();

  if (endpoints.empty()) {
    std::cout << "No endpoints found." << std::endl;

    return;
  } else {
    std::cout << "Found " << endpoints.size() << " endpoints." << std::endl;

    // Iterate through each endpoint and retrieve the critical path for each endpoint
    for (const Pin* endpoint : endpoints) {
      std::string endpoint_name = network->name(endpoint);
      std::cout << "Endpoint: " << endpoint_name << std::endl;

      PathEndSeq path_ends = sta->findPathEnds(
        nullptr, nullptr, new ExceptionTo(endpoint), false, nullptr,
        MinMaxAll::max(), 1, 1, false, -1e30, 1e30, true, nullptr,
        true, false, false, false, false, false);

      if (path_ends.empty()) {
        std::cout << "No critical path found for endpoint " << endpoint_name << std::endl;
      } else {
        for (const PathEnd* path_end : path_ends) {
          path_end->reportPath(std::cout, network, graph, 2);
        }
      }
    }

  }
}

} // namespace rsz

