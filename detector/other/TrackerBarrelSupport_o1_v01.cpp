#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"
#include "TGeoTrd2.h"
#include "XML/Utilities.h"


using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)
{
  static double tolerance = 0e0;

  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  DetElement    sdet(det_name, det_id);

  Volume envelope = XML::createPlacedEnvelope(lcdd,  e , sdet) ;

  if (lcdd.buildType() == BUILD_ENVELOPE) return sdet ;

  Material air = lcdd.air();
  PlacedVolume pv;
  int n = 0;

  for (xml_coll_t i(x_det, _U(layer)); i; ++i, ++n)  {
    xml_comp_t x_layer = i;
    string  l_name = det_name + _toString(n, "_layer%d");
    double  z    = x_layer.outer_z();
    double  rmin = x_layer.inner_r();
    double  r    = rmin;
    DetElement layer(sdet, _toString(n, "layer%d"), x_layer.id());
    Tube    l_tub(rmin, 2 * rmin, z);
    Volume  l_vol(l_name, l_tub, air);
    int m = 0;

    for (xml_coll_t j(x_layer, _U(slice)); j; ++j, ++m)  {
      xml_comp_t x_slice = j;
      Material mat = lcdd.material(x_slice.materialStr());
      string s_name = l_name + _toString(m, "_slice%d");
      double thickness = x_slice.thickness();

      Tube   s_tub(r, r + thickness, z, 2 * M_PI);
      Volume s_vol(s_name, s_tub, mat);

      r += thickness;

      // Set Attributes
      s_vol.setAttributes(lcdd, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      pv = l_vol.placeVolume(s_vol);
      // Slices have no extra id. Take the ID of the layer!
      pv.addPhysVolID("slice", m);
    }
    l_tub.setDimensions(rmin, r, z);
    //cout << l_name << " " << rmin << " " << r << " " << z << endl;
    l_vol.setVisAttributes(lcdd, x_layer.visStr());

    pv = envelope.placeVolume(l_vol);
    pv.addPhysVolID("layer", n);
    layer.setPlacement(pv);
  }



  return sdet;

}

DECLARE_DETELEMENT(TrackerBarrelSupport_o1_v01, create_detector)
