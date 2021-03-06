#ifndef FEA_DEFINES_HEADER
#define FEA_DEFINES_HEADER

#define FEA_NAMESPACE_BEGIN namespace fea {
#define FEA_NAMESPACE_END }

#define FEA_MATCH_NAMESPACE_BEGIN namespace fea_mch {
#define FEA_MATCH_NAMESPACE_END }

#define DPG_MATCH_NAMESPACE_BEGIN namespace dpg_mch {
#define DPG_MATCH_NAMESPACE_END }

#ifndef STD_NAMESPACE_BEGIN
#   define STD_NAMESPACE_BEGIN namespace std {
#   define STD_NAMESPACE_END }
#endif

/** Structure-from-Motion library. */
FEA_NAMESPACE_BEGIN
FEA_MATCH_NAMESPACE_BEGIN 
DPG_MATCH_NAMESPACE_BEGIN DPG_MATCH_NAMESPACE_END
FEA_MATCH_NAMESPACE_END
FEA_NAMESPACE_END

#endif /* MVE_DEFINES_HEADER */
