//-*-C++-*-

#ifndef CSV_WRITER_HEADER_H
#define CSV_WRITER_HEADER_H

namespace IO
{
  
  template<>
  class writer<IO::CSV>
  {
  public:

    template<typename scalartype>
    static void execute(std::string file_name,  std::vector<std::vector<scalartype> >& data);
  };
  
  template<typename scalartype>
  void writer<IO::CSV>::execute(std::string file_name,  std::vector<std::vector<scalartype> >& data)
  {
    ofstream myfile;
    myfile.open(file_name.c_str());

    for(size_t j=0; j<data.size(); ++j)
      {
	for(size_t i=0; i<data[j].size(); ++i)
	  {
	    myfile << data[i][j];
	    
	    if(i==data[j].size()-1)
	      myfile << "\n";	    
	    else
	      myfile << ", ";
	  }
      }

    myfile.close();
  }

}

#endif
