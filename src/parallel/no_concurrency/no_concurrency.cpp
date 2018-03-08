#include "dca/parallel/no_concurrency/no_concurrency.hpp"

namespace dca {
namespace parallel {

constexpr char NoConcurrency::parallel_type_str_[];

std::ostream& operator<<(std::ostream& o, const NoConcurrency& c)
{
  o << '\n' << "parallel type:" << c.parallel_type_str_
    << '\n' << "number of processors:" << c.number_of_processors()
    << '\n' << "grouping first:" << c.first()
    << '\n' << "grouping last::" << c.last();
  return o;
}

} // parallel
} // dca
