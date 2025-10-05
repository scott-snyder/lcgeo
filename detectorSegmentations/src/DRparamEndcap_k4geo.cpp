#include "detectorSegmentations/DRparamEndcap_k4geo.h"

#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

  double DRparamEndcap_k4geo::currentInnerR (double innerX, double thetaOfCenter) const
  {
    return innerX / std::sin(thetaOfCenter);
  }

  int DRparamEndcap_k4geo::getBETowerNo (int signedTowerNo, int BEtrans) const
  {
    return unsignedTowerNo(signedTowerNo) - unsignedTowerNo(BEtrans);
  }

} // namespace DDSegmentation
} // namespace dd4hep
