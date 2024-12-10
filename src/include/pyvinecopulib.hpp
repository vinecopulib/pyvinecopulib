#if defined(_MSC_VER) && _MSC_VER < 1500 // VC++ 8.0 and below
#define snprintf _snprintf
#endif

#include "bicop/class.hpp"
#include "bicop/family.hpp"
#include "bicop/fit_controls.hpp"
#include "misc/stats.hpp"
#include "vinecop/class.hpp"
#include "vinecop/fit_controls.hpp"
#include "vinecop/rvine_structure.hpp"