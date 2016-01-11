#ifndef GIT_VERSION_HPP
#define GIT_VERSION_HPP

#include <string>

struct GitVersion {
  static const std::string GIT_LOG;
  static const std::string GIT_STATUS;

  static void        print();
  static std::string string();
};

#endif  // GIT_VERSION_HPP
