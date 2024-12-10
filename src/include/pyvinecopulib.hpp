#if (defined(_MSC_VER) && (_MSC_VER >= 1900)) || defined(__MINGW32__)
#define HAVE_SNPRINTF 1
#endif

#include "bicop/class.hpp"
#include "bicop/family.hpp"
#include "bicop/fit_controls.hpp"
#include "misc/stats.hpp"
#include "vinecop/class.hpp"
#include "vinecop/fit_controls.hpp"
#include "vinecop/rvine_structure.hpp"