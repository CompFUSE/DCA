//-*-C++-*-

#ifdef ALLOW_GNUPLOT

#include "gnuplot_interface.hpp"

#else

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

class Gnuplot
{
public:

  Gnuplot(const std::string &style = "points") {}

  Gnuplot(const std::vector<double> &x,
          const std::string &title = "",
          const std::string &style = "points",
          const std::string &labelx = "x",
          const std::string &labely = "y") {}

  Gnuplot(const std::vector<double> &x,
          const std::vector<double> &y,
          const std::string &title = "",
          const std::string &style = "points",
          const std::string &labelx = "x",
          const std::string &labely = "y") {}

  Gnuplot(const std::vector<double> &x,
          const std::vector<double> &y,
          const std::vector<double> &z,
          const std::string &title = "",
          const std::string &style = "points",
          const std::string &labelx = "x",
          const std::string &labely = "y",
          const std::string &labelz = "z") {}

  ~Gnuplot() {}
  
  void plot_xy(std::vector<double> x,
	       std::vector<double> y) 
  {}

  void showonscreen() {}
};

#endif 

#endif 

#include "show_function.h"
