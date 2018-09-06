#ifndef DPG_DEFINES_HEADER
#define DPG_DEFINES_HEADER

#define DPG_NAMESPACE_BEGIN namespace dpgrid {
#define DPG_NAMESPACE_END }

#define DPG_IMAGE_NAMESPACE_BEGIN namespace image {
#define DPG_IMAGE_NAMESPACE_END }

//#define DPG_GEOM_NAMESPACE_BEGIN namespace geom {
//#define DPG_GEOM_NAMESPACE_END }

#ifndef STD_NAMESPACE_BEGIN
#   define STD_NAMESPACE_BEGIN namespace std {
#   define STD_NAMESPACE_END }
#endif

/** dpgrid library. */
DPG_NAMESPACE_BEGIN
/** Image tools, loading and processing functions. */
DPG_IMAGE_NAMESPACE_BEGIN DPG_IMAGE_NAMESPACE_END
/** Geometric tools, loading and processing functions. */
//DPG_GEOM_NAMESPACE_BEGIN DPG_GEOM_NAMESPACE_END
DPG_NAMESPACE_END

#endif /* DPG_DEFINES_HEADER */

