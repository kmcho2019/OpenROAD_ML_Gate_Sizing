#ifndef RSZ_TRANSFORMER_RESIZER_HH
#define RSZ_TRANSFORMER_RESIZER_HH

#include "MLGateSizer.hh"
#include "Resizer.hh"
#include <vector>
#include <string>

namespace rsz {

class TransformerResizer : public Resizer {
 public:
  TransformerResizer();
  
  void optimizeTNS();
  void optimizeWNS();
  void optimizeLeakagePower();
  void optimizeSlew();
  void optimizeCapacitance();

 private:
  void extractWorstPaths();
  void transformPinDataToTokens();
  MLGateSizer gateSizer;
};

} // namespace rsz

#endif

