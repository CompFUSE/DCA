//-*-C++-*-

/*
  ________  _________     _____
  \______ \ \_   ___ \   /  _  \      .__
  |    |  \/    \  \/  /  /_\  \   __|  |___
  |    `   \     \____/    |    \ /__    __/
  /_______  /\______  /\____|__  /    |__|
  \/        \/         \/

  __________        __                   _________ __
  \______   \ _____/  |_  ___________   /   _____//  |______  _____ _______
  |     ___// __ \   __\/ __ \_  __ \  \_____  \\   __\__  \ \__  \\_  __ \
  |    |   \  ___/|  | \  ___/|  | \/  /        \|  |  / __ \_/ __ \|  | \/
  |____|    \___  >__|  \___  >__|    /_______  /|__| (____  (____  /__|
  \/          \/                \/           \/     \/
  _________________________ ___   __________            .__       .__
  \_   _____/\__    ___/   |   \  \____    /__ _________|__| ____ |  |__
  |    __)_   |    | /    ~    \   /     /|  |  \_  __ \  |/ ___\|  |  \
  |        \  |    | \    Y    /  /     /_|  |  /|  | \/  \  \___|   Y  \
  /_______  /  |____|  \___|_  /  /_______ \____/ |__|  |__|\___  >___|  /
  \/                 \/           \/                    \/     \/
*/

#include "compiler_directives.h"

#include "include_files.h"

void print_device_info(){}

std::string get_version()
{
  string str(
	"commit bdff9dd0cb195220993ff3be43accdf1a2fb4583\n"
	"Author: Peter Staar <staarp@itp.phys.ethz.ch>\n"
	"Date:   Wed Oct 23 16:05:31 2013 +0200\n"
	"\n"
	"updated include_file and typedefinitions	2013-10-23	16:05\n");
  return str;
}

void print_version()
{
  string str = get_version();

  cout << "\n\n\n";
  cout << "*************************************************************************************\n";
  cout << "***                                  VERSION                                      ***\n";
  cout << "*************************************************************************************\n";
  cout << "\n\n\n";

  cout << str << endl;
}

;

int main(int argc,char *argv[])
{
  

  if(argc < 2)
    {
      std::cout << "Usage: "<< argv[0] << " inputFileName\n";
      return -1;
    }

  std::string file_name(argv[1]);

  static const COMP_LIB::PARALLELIZATION_LIBRARY_NAMES PARALLELIZATION_LIBRARY_NAME = COMP_LIB::SERIAL_LIBRARY;

  // ============================================================ Configure the calculation by selecting type definitions.


  typedef COMP_LIB::parallelization<PARALLELIZATION_LIBRARY_NAME>  concurrency_type;

  typedef DFT::VASP::parameters<concurrency_type> parameters_type;
  typedef DFT::VASP::data<parameters_type>        vasp_data_type;
  typedef DFT::VASP::reader<parameters_type>      vasp_reader_type;

  // ====================================================================== Create the algorithms and parameters object from the input file

  concurrency_type concurrency(argc, argv);

  //parameters_type::profiler_type::start();

  std::string stamp = get_version();

  /*
  if(concurrency.id() == concurrency.first())
    {
      cout << "DCA main: starting (MPI-world set up).\n\n";

      print_device_info();

      print_version();
    }
  */

  parameters_type parameters(stamp, concurrency);

  parameters.read_input(file_name);
  parameters.write_input("./check.json");
  //parameters.write_input(std::cout);

  parameters.update_domains();

  vasp_data_type   vasp_data(parameters);

  vasp_reader_type vasp_reader(parameters, vasp_data);

  vasp_reader.execute();

  vasp_data.execute();

//   {
//     vasp_data.construct_bloch_hamiltonians();

//     vasp_data.check_bloch_hamiltonians();    

//     vasp_data.downfold_bloch_hamiltonians();

//     vasp_data.construct_t_ij();
//   }

  //parameters_type::profiler_type::stop(concurrency, parameters.get_profiling_file_name());

  if(concurrency.id() == concurrency.last())
    cout << "\n\nmain: ending. \n\n";

  return 0;
}
