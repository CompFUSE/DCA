#include "dca/config/haves_defines.hpp"
#include "dca/config/config_defines.hpp"
#include "dca/config/hpx_defines.hpp"
#include "dca/parallel/hpx/hpx.hpp"

#if defined(DCA_HAVE_HPX_LOGGING) && !defined(NDEBUG)
log_mutex_type log_mutex;
#endif

namespace dca {
namespace parallel {

//hpxthread::hpxthread() : threads_(0), data_(0) {}

}}
