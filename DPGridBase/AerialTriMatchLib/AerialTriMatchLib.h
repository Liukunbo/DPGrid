#pragma once

/** 支持跨平台的导出宏定义 */
#if _WIN32
# ifndef MAKE_DLL
#   define ArialTriMatchLib __declspec(dllimport)
# else
#	define ArialTriMatchLib __declspec(dllexport)
# endif
#else
#	define ArialTriMatchLib
#endif

#include "image\image.h"
#include "image\image_io.h"
#include "image\image_tools.h"
#include "feature\dpgrid_common.h"
#include "feature\dpgrid_matching.h"
#include "util\timer.h"
#include "math\algo.h"
