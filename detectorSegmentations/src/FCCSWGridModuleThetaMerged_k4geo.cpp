#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

#include <iostream>
#include "DD4hep/Detector.h"
#include "DD4hep/VolumeManager.h"
#include "DDRec/DetectorData.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
FCCSWGridModuleThetaMerged_k4geo::FCCSWGridModuleThetaMerged_k4geo(const std::string& cellEncoding) : GridTheta_k4geo(cellEncoding) {
  // define type and description
  _type = "FCCSWGridModuleThetaMerged_k4geo";
  _description = "Module-theta segmentation with per-layer merging along theta and/or module";

  // register all necessary parameters (additional to those registered in GridTheta_k4geo)
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerIdentifier("identifier_module", "Cell ID identifier for readout module", m_moduleID, "module");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedModules", "Numbers of merged modules per layer", m_mergedModules, std::vector<int>());
  GetNModulesFromGeom();
  GetNLayersFromGeom();
}

FCCSWGridModuleThetaMerged_k4geo::FCCSWGridModuleThetaMerged_k4geo(const BitFieldCoder* decoder) : GridTheta_k4geo(decoder) {
  // define type and description
  _type = "FCCSWGridModuleThetaMerged_k4geo";
  _description = "Module-theta segmentation with per-layer merging along theta and/or module";

  // register all necessary parameters (additional to those registered in GridTheta_k4geo)
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerIdentifier("identifier_module", "Cell ID identifier for module", m_moduleID, "module");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedModules", "Numbers of merged modules per layer", m_mergedModules, std::vector<int>());
  GetNModulesFromGeom();
  GetNLayersFromGeom();
}

void FCCSWGridModuleThetaMerged_k4geo::GetNModulesFromGeom() {
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  try {
    m_nModules = dd4hepgeo->constant<int>("ECalBarrelNumPlanes");
  }
  catch(...) {
    std::cout << "Number of modules not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  std::cout << "Number of modules read from detector metadata and used in readout class: " << m_nModules << std::endl;
}

void FCCSWGridModuleThetaMerged_k4geo::GetNLayersFromGeom() {
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  try {
    m_nLayers = dd4hepgeo->constant<int>("ECalBarrelNumLayers");
  }
  catch(...) {
    std::cout << "Number of layers not found in detector metadata, exiting..." << std::endl;
    exit(1);
  }
  std::cout << "Number of layers read from detector metadata and used in readout class: " << m_nLayers << std::endl;
}

/// Return the per-layer information corresponding to a volume.
const FCCSWGridModuleThetaMerged_k4geo::LayerInfo&
FCCSWGridModuleThetaMerged_k4geo::getLayerInfo(VolumeID vID) const
{
  using LayerInfoData = dd4hep::rec::StructExtension<LayerInfo>;

  // Look up the DetElement for this volume and see if we've already
  // added the per-layer information.
  dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());
  VolumeManager vman = VolumeManager::getVolumeManager(*dd4hepgeo);
  DetElement de = vman.lookupDetElement(vID);
  const LayerInfo* li = de.extension<LayerInfoData>(false);

  if (!li) {
    // Not there.  We need to make it.
    // Look up this volume in the volume manager, and find its radius
    // by transforming the origin in the local coordinate system to global
    // coordinates.
    VolumeManagerContext* vc = vman.lookupContext(vID);
    Position wpos = vc->localToWorld({0,0,0});
    double rho = wpos.Rho();

    // If different modules are ganged together, we want to put hits
    // in the center of all the ganged modules.  This represents
    // a rotation in phi (in global coordinates).  Convert this rotation
    // to displacement in local coordinates.  This will be in the x-z plane,
    // and will be constant for all modules in a layer, but will be different
    // for different layers (even with identical ganging).
    double xloc = 0;
    double zloc = 0;
    double phioff = phi(vID);
    if (phioff > 0) {
      // We need to apply a phi offset.  Calculate it by rotating
      // the global position in phi and converting back to local
      // coordinates.
      //
      // For the Allegro calorimeter, this can also be calculated as:
      //
      //   d = rho * 2 * sin(phioff/2)
      //   beta = pi/2 - phioff/2 - asin(sin(InclinationAngle)*EMBarrel_rmin/rho)
      //   xloc = - d * sin(beta)
      //   yloc =   d * cos(beta)
      //
      // but we prefer to do the calculation via explicit rotations because
      // it's easier to see that that was is correct, and it also avoids
      // the explicit dependencies on the geometry parameters.
      Position wpos2 = RotateZ(wpos, phioff);
      Position lpos2 = vc->worldToLocal(wpos2);
      xloc = lpos2.X();
      zloc = lpos2.Z();
    }

    // Remember this by adding it to the DE as an extension.
    li = de.addExtension<LayerInfoData>(new LayerInfoData(LayerInfo(rho, xloc, zloc)));
  }

  return *li;
}


/// determine the local position based on the cell ID
Vector3D FCCSWGridModuleThetaMerged_k4geo::position(const CellID& cID) const {

  VolumeID vID = cID;
  _decoder->set(vID, m_thetaID, 0);
  const LayerInfo& li = this->getLayerInfo(vID);

  // debug
  // std::cout << "cellID: " << cID << std::endl;

  // Calculate the position in local coordinates.
  // The volume here has the cross-section of a cell in the x-z plane;
  // it extends the length of the calorimeter along the y-axis.
  // We set the y-coordinate based on the theta bin, and x and z
  // based on the phi offset required for this layer.
  return Vector3D(li.xloc,
                  -li.rho / tan(theta(cID)),
                  li.zloc);
}

/// determine the cell ID based on the global position
CellID FCCSWGridModuleThetaMerged_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;

  // retrieve layer (since merging depends on layer)
  int layer = this->layer(vID);

  // retrieve theta
  double lTheta = thetaFromXYZ(globalPosition);

  // calculate theta bin with original segmentation
  int thetaBin = positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta);

  // adjust theta bin if cells are merged along theta in this layer
  // assume that m_mergedCellsTheta[layer]>=1
  thetaBin -= (thetaBin % m_mergedCellsTheta[layer]);

  // set theta field of cellID
  _decoder->set(cID, m_thetaID, thetaBin);

  // retrieve module number
  int module = _decoder->get(vID, m_moduleID);

  // adjust module number if modules are merged in this layer
  // assume that m_mergedModules[layer]>=1
  module -= (module % m_mergedModules[layer]);

  // set module field of cellID
  _decoder->set(cID, m_moduleID, module);

  return cID;
}

/// determine the azimuth based on the cell ID
/// the value returned is the relative shift in phi
/// with respect to the first module in the group of
/// merged ones - which will be then added on top of
/// the phi of the volume containing the first cell
/// by the positioning tool
double FCCSWGridModuleThetaMerged_k4geo::phi(const CellID& cID) const {

  // retrieve layer
  int layer = this->layer(cID);

  // calculate phi offset due to merging
  // assume that m_mergedModules[layer]>=1
  double phi = (m_mergedModules[layer]-1) * M_PI / m_nModules ;

  // debug
  // std::cout << "layer: " << layer << std::endl;
  // std::cout << "merged modules: " << m_mergedModules[layer] << std::endl;
  // std::cout << "phi: " << phi << std::endl;

  return phi;
}

/// determine the polar angle based on the cell ID and the
/// number of merged theta cells
double FCCSWGridModuleThetaMerged_k4geo::theta(const CellID& cID) const {

  // retrieve layer
  int layer = this->layer(cID);

  // retrieve theta bin from cellID and determine theta position
  CellID thetaValue = _decoder->get(cID, m_thetaID);
  double _theta = binToPosition(thetaValue, m_gridSizeTheta, m_offsetTheta);

  // adjust return value if cells are merged along theta in this layer
  // shift by (N-1)*half theta grid size
  // assume that m_mergedCellsTheta[layer]>=1
  _theta += (m_mergedCellsTheta[layer]-1) * m_gridSizeTheta / 2.0 ;

  // debug
  // std::cout << "layer: " << layer << std::endl;
  // std::cout << "theta bin: " << thetaValue << std::endl;
  // std::cout << "gridSizeTheta, offsetTheta: " << m_gridSizeTheta << " , " << m_offsetTheta << std::endl;
  // std::cout << "merged cells: " << m_mergedCellsTheta[layer] << std::endl;
  // std::cout << "theta: " << _theta << std::endl;

  return _theta;
}

/// Extract the layer index fom a cell ID.
int FCCSWGridModuleThetaMerged_k4geo::layer(const CellID& cID) const {
  return _decoder->get(cID, m_layerID);
}

/// Determine the volume ID from the full cell ID by removing all local fields
VolumeID FCCSWGridModuleThetaMerged_k4geo::volumeID(const CellID& cID) const {
  VolumeID vID = cID;
  _decoder->set(vID, m_thetaID, 0);
  return vID;
}

}
}
