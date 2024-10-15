#include "TransformerResizer.hh"
#include "MLGateSizer.hh"

namespace rsz {

TransformerResizer::TransformerResizer()
{
  // Initialize the gateSizer (load transformer weights, etc.)
  gateSizer.loadWeights("path/to/weights");
}

void TransformerResizer::optimizeTNS()
{
  extractWorstPaths();
  // Perform TNS optimization using ML model
  gateSizer.classifyAndResize();
}

void TransformerResizer::extractWorstPaths()
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

void TransformerResizer::transformPinDataToTokens(const Pin* pin, float slew, float cap, float slack, const std::string& gate_type)
{
  // Turn the extracted pin data into tokens that can be used by the transformer model
  std::vector<float> tokenized_data = {slew, cap, slack};
  
  // Add token to the gateSizer for later classification
  gateSizer.addToken(tokenized_data, gate_type);
}

} // namespace rsz

