//-*-C++-*-

#ifndef CUPRATE_THREE_BAND_DATA_H
#define CUPRATE_THREE_BAND_DATA_H

struct cuprate_three_band_data
{
  static std::vector<std::vector<double> >& get_t_ij();

  static void initialize(cuprate_compounds_type compound_name);
};

std::vector<std::vector<double> >& cuprate_three_band_data::get_t_ij()
{
  static std::vector<std::vector<double> > t_ij(0, std::vector<double>(0));
  return t_ij;
}

void cuprate_three_band_data::initialize(cuprate_compounds_type compound_name)
{
  switch(compound_name)
    {
    case La2CuO4:
      IO::reader<IO::CSV>::execute("t_ij_La2CuO4.txt", get_t_ij());
      break;

    case HgBa2CuO4:
      IO::reader<IO::CSV>::execute("t_ij_HgBa2CuO4.txt", get_t_ij());
      break;

    case Ca2CuO2Cl2:
      IO::reader<IO::CSV>::execute("t_ij_HgBa2CuO4.txt", get_t_ij());
      break;

    case Tl2Ba2CuO6:
      IO::reader<IO::CSV>::execute("t_ij_Tl2Ba2CuO6.txt", get_t_ij());
      break;

    default:
      cout << "\n\n\t\t\t this compund is not known!!!" << endl;
      throw std::logic_error(__FUNCTION__);
    }
}

#endif
