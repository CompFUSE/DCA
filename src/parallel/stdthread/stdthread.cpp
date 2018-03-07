#include "dca/parallel/stdthread/stdthread.hpp"

namespace dca {
namespace parallel {

constexpr char stdthread::concurrency_type_str_[];

std::ostream& operator<<(std::ostream& o, const stdthread& c)
{
  o << '\n' << "concurrency type:" << c.concurrency_type_str_
    << '\n' << "number of std::threads:" << c.threads_.size();
  return o;
}

} // parallel
} // dca
