
#ifndef DCA_MODEL_TRAITS_HPP
#define DCA_MODEL_TRAITS_HPP

#include <type_traits>

namespace dca {
namespace phys {
namespace models {

/** type trait indicating legacy return of pointer to static function variable behavior
 *  \todo remove this type trait entirely once this antipattern is purged from codebase
 */
template <class Lattice, typename = void>
struct AbusesStatic_t : public std::true_type{};

}
}  // namespace phys
}  // namespace dca

#endif
