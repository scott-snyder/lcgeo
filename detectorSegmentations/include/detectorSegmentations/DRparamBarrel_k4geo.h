#ifndef DETSEGMENTATION_DRPARAMBARREL_H
#define DETSEGMENTATION_DRPARAMBARREL_H

#include "detectorSegmentations/DRparamBase_k4geo.h"

#include "DD4hep/DetFactoryHelper.h"

namespace dd4hep {
namespace DDSegmentation {
  class DRparamBarrel_k4geo : public DRparamBase_k4geo {
  public:
    DRparamBarrel_k4geo() = default;
    virtual ~DRparamBarrel_k4geo() = default;

  protected:
    virtual double currentInnerR (double innerX, double thetaOfCenter) const override;
    virtual int getBETowerNo (int signedTowerNo, int BEtrans) const override;
  };
} // namespace DDSegmentation
} // namespace dd4hep

#endif
