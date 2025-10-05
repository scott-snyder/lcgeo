#ifndef DETSEGMENTATION_DRPARAMBASE_H
#define DETSEGMENTATION_DRPARAMBASE_H

#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"

#include <cmath>
#include <vector>

namespace dd4hep {
namespace DDSegmentation {
  class DRparamBase_k4geo {
  public:
    DRparamBase_k4geo();
    virtual ~DRparamBase_k4geo();

    void SetIsRHS(bool isRHS) { m_isRHS = isRHS; }
    void SetInnerX(double innerX) { m_innerX = innerX; }
    void SetTowerH(double towerH) { m_towerH = towerH; }
    void SetNumZRot(int num) {
      m_numZRot = num;
      m_phiZRot = 2 * M_PI / (double)num;
    }
    void SetDeltaTheta(double theta) { m_deltaTheta = theta; }
    void SetThetaOfCenter(double theta) { m_thetaOfCenter = theta; }
    void SetSipmHeight(double SipmHeight) { m_sipmHeight = SipmHeight; }

    bool GetIsRHS() const { return m_isRHS; }
    int GetNumZRot() const { return m_numZRot; }
    double GetCurrentInnerR() const { return m_currentInnerR; }
    double GetTowerH() const { return m_towerH; }
    double GetSipmHeight() const { return m_sipmHeight; }
    double GetH1() const { return m_currentInnerHalf; }
    double GetBl1() const { return m_v3.X() * std::tan(m_phiZRot / 2.); }
    double GetTl1() const { return m_v1.X() * std::tan(m_phiZRot / 2.); }
    double GetH2() const { return m_currentOuterHalf; }
    double GetBl2() const { return m_v4.X() * std::tan(m_phiZRot / 2.); }
    double GetTl2() const { return m_v2.X() * std::tan(m_phiZRot / 2.); }

    double GetH2sipm() const { return m_currentOuterHalfSipm; }
    double GetBl2sipm() const { return m_v4sipm.X() * std::tan(m_phiZRot / 2.); }
    double GetTl2sipm() const { return m_v2sipm.X() * std::tan(m_phiZRot / 2.); }

    dd4hep::RotationZYX GetRotationZYX(int numPhi) const;
    dd4hep::Position GetTowerPos(int numPhi) const;
    dd4hep::Position GetAssemblePos(int numPhi) const;
    dd4hep::Position GetSipmLayerPos(int numPhi) const;

    dd4hep::Transform3D GetTransform3D(int numPhi) const;
    dd4hep::Transform3D GetAssembleTransform3D(int numPhi) const;
    dd4hep::Transform3D GetSipmTransform3D(int numPhi) const;

    int signedTowerNo(int unsignedTowerNo) const { return m_isRHS ? unsignedTowerNo : -unsignedTowerNo - 1; }
    int unsignedTowerNo(int signedTowerNo) const { return signedTowerNo >= 0 ? signedTowerNo : -signedTowerNo - 1; }

    void SetDeltaThetaByTowerNo(int, int);
    void SetThetaOfCenterByTowerNo(int, int);
    void SetIsRHSByTowerNo(int signedTowerNo) { m_isRHS = (signedTowerNo >= 0 ? true : false); }

    int GetTotTowerNum() const { return m_totNum; }
    void SetTotTowerNum(int totNum) { m_totNum = totNum; }

    int GetCurrentTowerNum() const { return m_currentTowerNum; }
    void SetCurrentTowerNum(int numEta) { m_currentTowerNum = numEta; }

    // the size of the neighborhood (radius) in units of fiber-fiber distance (i.e. grid size)
    void SetNeighborSize(double nsize) { m_neighborSize = nsize; }
    double GetNeighborSize() const { return m_neighborSize; }
    // the margin in units of fiber-fiber distance when incorporating satelite cells at the edge of the tower
    void SetMargin(int margin) { m_margin = margin; }
    int GetMargin() const { return m_margin; }

    void init();
    void filled() { m_filled = true; }
    void finalized() { m_finalized = true; }
    bool IsFinalized() const { return m_finalized; }

    // store information of which fibers have full length or not
    // full length fibers within rmin <= n_row <= rmax and cmin <= n_column <= cmax
    struct fullLengthFibers {
    public:
      fullLengthFibers(int rmin_, int rmax_, int cmin_, int cmax_)
          : rmin(rmin_), rmax(rmax_), cmin(cmin_), cmax(cmax_) {}

      fullLengthFibers() = default;

      int rmin = 0; // min n_row with full length fibers
      int rmax = 0; // max n_row
      int cmin = 0; // min n_column with full length fibers
      int cmax = 0; // max n_column
    };

    const fullLengthFibers& GetFullLengthFibers(int numEta) const { return m_fullLengthFibers.at(unsignedTowerNo(numEta)); }
    void SetFullLengthFibers(int rmin, int rmax, int cmin, int cmax);

    // length of fibers that don't have full length
    // it is a map containing (row, col) as a key and the length of the fiber as a value
    // and each tower in eta has one of it
    struct shortFibers {
    public:
      shortFibers(const double towerH) : m_towerH(towerH) {}

      void addShortFibers(const int row, const int col, const double len) {
        m_fiberLengths[std::make_pair(row, col)] = len;
      }
      double retrieveFiberLength(const int row, const int col) const;

    private:
      std::map<std::pair<int, int>, double> m_fiberLengths;
      const double m_towerH;
    };

    const shortFibers& GetShortFibers(int numEta) const { return m_shortFibers.at(unsignedTowerNo(numEta)); }
    void SetShortFibers(const shortFibers& input);

  protected:
    virtual double currentInnerR (double innerX, double thetaOfCenter) const = 0;
    virtual int getBETowerNo (int signedTowerNo, int BEtrans) const = 0;

  private:
    bool m_isRHS = 0;
    double m_phiZRot = 0;
    double m_innerX = 0;
    double m_towerH = 0;
    int m_numZRot = 0;
    double m_deltaTheta = 0;
    double m_thetaOfCenter = 0;
    double m_currentInnerR = 0;
    TVector3 m_currentCenter{};
    TVector3 m_v1{};
    TVector3 m_v2{};
    TVector3 m_v3{};
    TVector3 m_v4{};
    TVector3 m_v2sipm{};
    TVector3 m_v4sipm{};
    double m_sipmHeight = 0;

    double m_currentInnerHalf = 0;
    double m_currentOuterHalf = 0;
    double m_currentOuterHalfSipm = 0;

    // parameters for neighbour finding algorithm
    double m_neighborSize = 0;
    int m_margin = 0;

    int m_totNum = 0;
    int m_currentTowerNum = 0;
    std::vector<double> m_deltaThetaVec;
    std::vector<double> m_thetaOfCenterVec;
    std::map<int, fullLengthFibers> m_fullLengthFibers;
    std::map<int, shortFibers> m_shortFibers;
    bool m_filled = false;
    bool m_finalized = false;
  };
} // namespace DDSegmentation
} // namespace dd4hep

#endif
