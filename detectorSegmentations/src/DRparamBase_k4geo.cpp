#include "detectorSegmentations/DRparamBase_k4geo.h"

#include "Math/GenVector/RotationZYX.h"

#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

  DRparamBase_k4geo::DRparamBase_k4geo() = default;
  DRparamBase_k4geo::~DRparamBase_k4geo() = default;

  dd4hep::RotationZYX DRparamBase_k4geo::GetRotationZYX(int numPhi) const {
    double xRot = m_isRHS ? -m_thetaOfCenter : m_thetaOfCenter;
    double zRot = m_isRHS ? -M_PI / 2. : M_PI / 2.;
    dd4hep::RotationZYX rot = dd4hep::RotationZYX(zRot, M_PI / 2. + xRot, 0.);
    ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(numPhi * m_phiZRot);
    return rotZ * rot;
  }

  dd4hep::Position DRparamBase_k4geo::GetTowerPos(int numPhi) const {
    double x = std::cos(numPhi * m_phiZRot) * m_currentCenter.X();
    double y = std::sin(numPhi * m_phiZRot) * m_currentCenter.X();
    double z = m_isRHS ? m_currentCenter.Z() : -m_currentCenter.Z();
    return dd4hep::Position(x, y, z);
  }

  dd4hep::Position DRparamBase_k4geo::GetAssemblePos(int numPhi) const {
    double alpha = (m_currentCenter.Mag() + m_sipmHeight / 2.) / m_currentCenter.Mag();
    return GetTowerPos(numPhi) * alpha;
  }

  dd4hep::Position DRparamBase_k4geo::GetSipmLayerPos(int numPhi) const {
    double alpha = (m_currentCenter.Mag() + m_towerH / 2. + m_sipmHeight / 2.) / m_currentCenter.Mag();
    return GetTowerPos(numPhi) * alpha;
  }

  dd4hep::Transform3D DRparamBase_k4geo::GetTransform3D(int numPhi) const {
    auto rot = GetRotationZYX(numPhi);
    auto pos = GetTowerPos(numPhi);
    return dd4hep::Transform3D(rot, pos);
  }

  dd4hep::Transform3D DRparamBase_k4geo::GetAssembleTransform3D(int numPhi) const {
    auto rot = GetRotationZYX(numPhi);
    auto pos = GetAssemblePos(numPhi);
    return dd4hep::Transform3D(rot, pos);
  }

  dd4hep::Transform3D DRparamBase_k4geo::GetSipmTransform3D(int numPhi) const {
    auto rot = GetRotationZYX(numPhi);
    auto pos = GetSipmLayerPos(numPhi);
    return dd4hep::Transform3D(rot, pos);
  }

  void DRparamBase_k4geo::SetFullLengthFibers(int rmin, int rmax, int cmin, int cmax) {
    if (m_filled)
      throw std::runtime_error(
          "DRparamBase_k4geo: An attempt to modify the geometry outside the detector construction is forbidden!");

    m_fullLengthFibers[m_currentTowerNum] = fullLengthFibers(rmin, rmax, cmin, cmax);
  }

  double DRparamBase_k4geo::shortFibers::retrieveFiberLength(const int row, const int col) const {
    auto it = m_fiberLengths.find(std::make_pair(row, col));
    if (it == m_fiberLengths.end())
      return m_towerH;

    return it->second;
  }

  void DRparamBase_k4geo::SetShortFibers(const shortFibers& input) {
    if (m_filled)
      throw std::runtime_error(
          "DRparamBase_k4geo: An attempt to modify the geometry outside the detector construction is forbidden!");

    m_shortFibers.emplace (m_currentTowerNum, input);
  }

  void DRparamBase_k4geo::init()
  {
    m_currentInnerR = currentInnerR (m_innerX, m_thetaOfCenter);

    const double r1 = m_currentInnerR;
    const double r2 = r1 + m_towerH;
    const double r3 = r2 + m_sipmHeight;

    const double trnsLength = m_towerH / 2. + m_currentInnerR;

    m_currentCenter = TVector3(std::cos(m_thetaOfCenter) * trnsLength, 0., std::sin(m_thetaOfCenter) * trnsLength);

    m_currentInnerHalf     = r1 * std::tan(m_deltaTheta / 2.);
    m_currentOuterHalf     = r2 * std::tan(m_deltaTheta / 2.);
    m_currentOuterHalfSipm = r3 * std::tan(m_deltaTheta / 2.);

    m_v1 = TVector3(std::cos(m_thetaOfCenter) * r1 +
                    std::sin(m_thetaOfCenter) * r1 * std::tan(m_deltaTheta / 2.),
                    0.,
                    std::sin(m_thetaOfCenter) * r1 -
                    std::cos(m_thetaOfCenter) * r1 * std::tan(m_deltaTheta / 2.));

    m_v2 = TVector3(std::cos(m_thetaOfCenter) * r2 +
                    std::sin(m_thetaOfCenter) * r2 * std::tan(m_deltaTheta / 2.),
                    0.,
                    std::sin(m_thetaOfCenter) * r2 -
                    std::cos(m_thetaOfCenter) * r2 * std::tan(m_deltaTheta / 2.));

    m_v3 = TVector3(std::cos(m_thetaOfCenter) * r1 -
                    std::sin(m_thetaOfCenter) * r1 * std::tan(m_deltaTheta / 2.),
                    0.,
                    std::sin(m_thetaOfCenter) * r1 +
                    std::cos(m_thetaOfCenter) * r1 * std::tan(m_deltaTheta / 2.));

    m_v4 = TVector3(std::cos(m_thetaOfCenter) * r2 -
                    std::sin(m_thetaOfCenter) * r2 * std::tan(m_deltaTheta / 2.),
                    0.,
                    std::sin(m_thetaOfCenter) * r2 +
                    std::cos(m_thetaOfCenter) * r2 * std::tan(m_deltaTheta / 2.));

    m_v2sipm =
        TVector3(std::cos(m_thetaOfCenter) * r3 +
                 std::sin(m_thetaOfCenter) * r3 * std::tan(m_deltaTheta / 2.),
                 0.,
                 std::sin(m_thetaOfCenter) * r3 -
                 std::cos(m_thetaOfCenter) * r3 * std::tan(m_deltaTheta / 2.));

    m_v4sipm =
        TVector3(std::cos(m_thetaOfCenter) * r3 -
                 std::sin(m_thetaOfCenter) * r3 * std::tan(m_deltaTheta / 2.),
                 0.,
                 std::sin(m_thetaOfCenter) * r3 +
                 std::cos(m_thetaOfCenter) * r3 * std::tan(m_deltaTheta / 2.));

    if (!m_filled) {
      m_deltaThetaVec.push_back(m_deltaTheta);
      m_thetaOfCenterVec.push_back(m_thetaOfCenter);
    }
  }

  void DRparamBase_k4geo::SetDeltaThetaByTowerNo(int signedTowerNo, int BEtrans) {
    if (!m_filled)
      throw std::runtime_error("Attempt to set by tower num while endcap parameter is not filled!");

    m_deltaTheta = m_deltaThetaVec.at(getBETowerNo(signedTowerNo, BEtrans));
  }

  void DRparamBase_k4geo::SetThetaOfCenterByTowerNo(int signedTowerNo, int BEtrans) {
    if (!m_filled)
      throw std::runtime_error("Attempt to set by tower num while endcap parameter is not filled!");

    m_thetaOfCenter = m_thetaOfCenterVec.at(getBETowerNo(signedTowerNo, BEtrans));
  }

} // namespace DDSegmentation
} // namespace dd4hep
