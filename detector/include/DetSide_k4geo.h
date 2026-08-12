#ifndef K4GEO_DETSIDE_K4GEO_H
#define K4GEO_DETSIDE_K4GEO_H

namespace k4geo {

/** Values of the "side" field of the tracker cellID encoding
 *  "system:5,side:-2,layer:9,module:8,sensor:8".
 *
 *  Replaces lcio::ILDDetID::{bwd,barrel,fwd} (UTIL/ILDConf.h), which was one of
 *  the last two reasons k4geo linked against LCIO. DD4hep ships no equivalent.
 */
namespace DetSide {
  inline constexpr int bwd = -1;
  inline constexpr int barrel = 0;
  inline constexpr int fwd = 1;
} // namespace DetSide

} // namespace k4geo

#endif // K4GEO_DETSIDE_K4GEO_H
