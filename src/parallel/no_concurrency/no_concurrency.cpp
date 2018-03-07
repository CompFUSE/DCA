#include "dca/parallel/no_concurrency/no_concurrency.hpp"

namespace dca {
namespace parallel {

constexpr char NoConcurrency::concurrency_type_str_[];

std::ostream& operator<<(std::ostream& o, const NoConcurrency& c)
{
  o << '\n' << "concurrency type:" << c.concurrency_type_str_
    << '\n' << "number of processors:" << c.number_of_processors()
    << '\n' << "grouping first:" << c.first()
    << '\n' << "grouping last::" << c.last();
  return o;
}

} // parallel
} // dca
