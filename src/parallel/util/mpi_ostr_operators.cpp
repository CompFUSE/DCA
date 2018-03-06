#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"

namespace dca {
namespace parallel {

constexpr char MPIConcurrency::concurrency_type_str_[];

std::ostream& operator<<(std::ostream &o, const MPIConcurrency &c)
{
  o << '\n' << "concurrency type:" << c.concurrency_type_str_
    << '\n' << "number of processors:" << c.number_of_processors()
    << '\n' << "grouping first:" << c.first()
    << '\n' << "grouping last::" << c.last();
  return o;
}

}
}
