#include "dca/parallel/no_threading/no_threading.hpp"

namespace dca {
namespace parallel {

constexpr char NoThreading::concurrency_type_str_[];

std::ostream& operator<<(std::ostream& o, const NoThreading& c) {
  o << '\n' << "concurrency type:" << c.concurrency_type_str_
    << '\n' << "number of threads:" << c.data_.num_threads;
  return o;
}

} // parallel
} // dca
