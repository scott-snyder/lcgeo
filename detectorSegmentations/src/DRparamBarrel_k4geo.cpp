#include "detectorSegmentations/DRparamBarrel_k4geo.h"

#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

  double DRparamBarrel_k4geo::currentInnerR (double innerX, double thetaOfCenter) const
  {
    return innerX / std::cos(thetaOfCenter);
  }


  int DRparamBarrel_k4geo::getBETowerNo (int signedTowerNo, int /*BEtrans*/) const
  {
    return unsignedTowerNo(signedTowerNo);
  }
} // namespace DDSegmentation
} // namespace dd4hep
