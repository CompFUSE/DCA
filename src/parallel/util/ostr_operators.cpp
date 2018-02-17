#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/parallel/pthreading/pthreading.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"

namespace dca {
namespace parallel {

constexpr char NoConcurrency::concurrency_type_str[];
constexpr char NoThreading::concurrency_type_str[];
constexpr char Pthreading::concurrency_type_str[];
constexpr char stdthread::concurrency_type_str[];


std::ostream& operator<<(std::ostream &o, const NoConcurrency &c)
{
  o << '\n' << "concurrency type:" << c.concurrency_type_str
    << '\n' << "number of processors:" << c.number_of_processors()
    << '\n' << "grouping first:" << c.first()
    << '\n' << "grouping last::" << c.last();
  return o;
}

std::ostream& operator<<(std::ostream &o, const NoThreading &c) {
  o << '\n' << "concurrency type:" << c.concurrency_type_str
    << '\n' << "number of threads:" << c.data_.num_threads;
  return o;
}

std::ostream& operator<<(std::ostream &o, const Pthreading &c)
{
  o << '\n' << "concurrency type:" << c.concurrency_type_str
    << '\n' << "number of posix threads:" << c.pthreads_.size();
  return o;
}

std::ostream& operator<<(std::ostream &o, stdthread const  &c)
{
  o << '\n' << "concurrency type:" << c.concurrency_type_str
    << '\n' << "number of posix threads:" << c.threads_.size();
  return o;
}

}
}
