//-*-C++-*-

#ifndef KOSHEVNIKOV_PARSER_H_
#define KOSHEVNIKOV_PARSER_H_

/*!
 * parser.h
 *
 *  Created on: 21 December, 2009
 *      Author: peterstaar
 */
class Koshevnikov_parser 
{
public:

  static double read_Fermi_energy(std::string filename);
  static int    read_Nb_Wannier_functions(std::string filename);

  static void   read_grid_size(std::string filename, std::vector<int>& v);

  static void   read_k_basis(std::string filename, double* ptr);
  static void   read_r_basis(std::string filename, double* ptr);

  static void   read_DCA_k_basis(std::string filename, double* ptr);
  static void   read_DCA_r_basis(std::string filename, double* ptr);

  template<class domain>
  static void read_interaction_Hamiltonian(std::string filename,
					   function<double, domain>& H,
					   int BANDS);
  template<class domain>
  static void read_symmetries(std::string            filename,
			      function<int, domain>& H_symmetries,
			      int                    BANDS);

  template<class domain>
  static void   read_LDA_Hamiltonians(std::string filename,
				      function<std::complex<double>, domain>& H,
				      double Fermi_energy);


private:  
 
  template<class MultiOrbitalMultiSiteStructure>
  static void read_Hamiltonians(const char* filename, MultiOrbitalMultiSiteStructure& MOMS);

  static double              read_double(const char* p);
  static std::vector<double> read_vector(int Nb_doubles, const char* p);

};

double Koshevnikov_parser::read_Fermi_energy(std::string filename)
{
//   cout<<scientific;
//   cout.precision(6);
//   cout << endl << __FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);

  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;
      
      searched_string = "# fermi energy";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  file_op.getline(line,2000);
	  double value = read_double(line);
// 	  cout << "\t" << value << endl;

	  file_op.close();
	  return value;
	}
    }
  
  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("fermi energy not read !!!");
}

int Koshevnikov_parser::read_Nb_Wannier_functions(std::string filename)
{
//   cout << endl << __FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);

  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;
      
      searched_string = "# number of Wannier functions";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  file_op.getline(line,2000);
	  int value = read_double(line);
// 	  cout << "\t" << value << endl;
	  file_op.close();
	  return value;
	}
    }
  
  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("Nb_Wannier_functions not read !!!");
}

void Koshevnikov_parser::read_grid_size(std::string filename, std::vector<int>& v)
{
//   cout << endl << __FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);
  
  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# k-grid size";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  file_op.getline(line,2000);
	  std::vector<double> tmp = read_vector(3, line);
	  for(int i=0; i<3; i++){
	    v[i] = int(tmp[i]);
// 	    cout << "\t" << v[i] ;
	  }

	  file_op.close();
	  return;
	}
    }

  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("LDA grid-size not read !!!");
}

void Koshevnikov_parser::read_k_basis(std::string filename, double* ptr)
{
//   cout << endl << __FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);
  
  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# reciprocal lattice vectors (3 rows)";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  for(int i=0;i<3 ;i++)
	    {
	      file_op.getline(line,2000);
	      std::vector<double> k_vec  = read_vector(3 , &line[0]);
	      
	      for(int j=0;j<3 ;j++)
		{
// 		  cout << "\t" << k_vec[j];
		  ptr[j+3*i] = k_vec[j];
		}
// 	      cout << endl;
	    }

	  file_op.close();
	  return;
	}
    }

  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("k-basis not read !!!");
}

void Koshevnikov_parser::read_r_basis(std::string filename, double* ptr)
{
//   cout << endl <<__FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);
  
  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# lattice vectors (3 rows)";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  for(int i=0;i<3 ;i++)
	    {
	      file_op.getline(line,2000);
	      std::vector<double> r_vec  = read_vector(3 , &line[0]);
	      
	      for(int j=0;j<3 ;j++)
		{
		  ptr[j+3*i] = r_vec[j];
// 		  cout << "\t" << r_vec[j];
		}
// 	      cout << endl;
	    }

	  file_op.close();
	  return;
	}
    }

  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("r-basis not read !!!");
}

void Koshevnikov_parser::read_DCA_k_basis(std::string filename, double* ptr)
{
//   cout << endl << __FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);
  
  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# DCA reciprocal lattice vectors (3 rows)";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  for(int i=0;i<3 ;i++)
	    {
	      file_op.getline(line,2000);
	      std::vector<double> k_vec  = read_vector(3 , &line[0]);
	      
	      for(int j=0;j<3 ;j++)
		{
// 		  cout << "\t" << k_vec[j];
		  ptr[j+3*i] = k_vec[j];
		}
// 	      cout << endl;
	    }

	  file_op.close();
	  return;
	}
    }

  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("DCA k-basis not read !!!");
}

void Koshevnikov_parser::read_DCA_r_basis(std::string filename, double* ptr)
{
//   cout << endl <<__FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);
  
  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# DCA lattice vectors (3 rows)";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  for(int i=0;i<3 ;i++)
	    {
	      file_op.getline(line,2000);
	      std::vector<double> r_vec  = read_vector(3 , &line[0]);
	      
	      for(int j=0;j<3 ;j++)
		{
		  ptr[j+3*i] = r_vec[j];
// 		  cout << "\t" << r_vec[j];
		}
// 	      cout << endl;
	    }

	  file_op.close();
	  return;
	}
    }

  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("DCA r-basis not read !!!");
}


template<class domain>
void Koshevnikov_parser::read_interaction_Hamiltonian(std::string filename,
						      function<double, domain>& H_interaction,
						      int BANDS)
{
//   cout << endl <<__FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);
  
  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# interaction";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  for(int i=0;i<2*BANDS ;i++)
	    {
	      file_op.getline(line,2000);
	      std::vector<double> H_i_  = read_vector(2*BANDS , &line[0]);
	      
	      for(int j=0;j<2*BANDS ;j++)
		{
		  H_interaction(i,j,0) = H_i_[j];
// 		  cout << "\t" << H_i_[j];
		}
// 	      cout << endl;
	    }

	  file_op.close();
	  return;
	}
    }

  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("interaction not read !!!");
}

template<class domain>
void Koshevnikov_parser::read_symmetries(std::string            filename,
					 function<int, domain>& H_symmetries,
					 int                    BANDS)
{
//   cout << endl <<__FUNCTION__ << endl;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  fstream file_op(&filename[0],ios::in);
  
  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# Hamiltonian-symmetries";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  for(int i=0;i<2*BANDS ;i++)
	    {
	      file_op.getline(line,2000);
	      std::vector<double> H_symm_i  = read_vector(2*BANDS , &line[0]);

	      for(int j=0;j<2*BANDS ;j++)
		{
		  H_symmetries(i,j) = H_symm_i[j];
// 		  cout << "\t" << H_symmetries(i,j);
		}
// 	      cout << endl;

	    }

	  file_op.close();
	  return;
	}
    }

  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error("symmetries not read !!!");
}


template<class domain>
void Koshevnikov_parser::read_LDA_Hamiltonians(std::string filename,
					       function<std::complex<double>, domain>& H,
					       double Fermi_energy)
{
  cout.precision(5);
  cout << scientific;

  char line[2000];
  string l;
  string searched_string;
  size_t found;

  int Nb_Wannier_functions = read_Nb_Wannier_functions(filename);

  fstream file_op(&filename[0],ios::in);
  int K_ind=0;

  while(!file_op.eof()) 
    {
      file_op.getline(line,2000);
      l = line;

      searched_string = "# k-point :";
      found=l.find(searched_string);
      if (found!=string::npos)
	{
	  //cout << searched_string << "\t" << linind << endl;

	  file_op.getline(line,2000);
	  l = line;

	  searched_string = "# weight";
	  found=l.find(searched_string);
	  if (found!=string::npos)
	    {
	      file_op.getline(line,2000);
	      //MOMS.K_points_LDA_weight(K_ind) = read_double( &line[0]);
	    }


	  file_op.getline(line,2000);
	  l = line;
	  
	  searched_string = "# lattice coordinates";
	  found=l.find(searched_string);
	  if (found!=string::npos)
	    {
	      //cout << searched_string << endl;
	      file_op.getline(line,2000);
	      //std::vector<double> lattice_vec  = read_vector(3, &line[0]);
	      //MOMS_struct.K_points.set_element(K_ind,lattice_vec);
	    }

	  file_op.getline(line,2000);
	  l = line;

	  searched_string = "# Cartesian coordinates";
	  found=l.find(searched_string);
	  if (found!=string::npos)
	    {
	      file_op.getline(line,2000);
	      std::vector<double> lattice_vec  = read_vector(3, &line[0]);
// 	      for(int j=0;j<3 ;j++)
// 		cout << lattice_vec[j] << "\t";
	      
	      /*if( fabs(lattice_vec[0]-k(K_ind)[0]) > 1.e-6 || fabs(lattice_vec[1]-k(K_ind)[1]) > 1.e-6 || fabs(lattice_vec[2]-k(K_ind)[2]) > 1.e-6)
		{
		  cout << K_ind << endl;

		  cout << lattice_vec[0] << "\t" << lattice_vec[1] << "\t" << lattice_vec[2] << "\t <---> \t";
		  cout << k(K_ind)[0] << "\t" << k(K_ind)[1] << "\t" << k(K_ind)[2] << "\n";

		  cout << __FUNCTION__ << endl;

		  throw std::logic_error("LDA-k-grid does not match Koshevnikov-grid !!!");
		  }*/
	    }

	  file_op.getline(line,2000);
	  l = line;
	 
	  searched_string = "# real part of H";
	  found=l.find(searched_string);
	  if (found!=string::npos)
	    {
	      for(int i=0;i<Nb_Wannier_functions ;i++)
		{
		  file_op.getline(line,2000);
		  std::vector<double> hamil  = read_vector(Nb_Wannier_functions , &line[0]);
		  
		  for(int j=0;j<Nb_Wannier_functions ;j++)
		    {
		      real(H(i,0,j,0,K_ind)) = hamil[j];
		      real(H(i,1,j,1,K_ind)) = hamil[j];

		      if(i == j)
			{
			  real(H(i,0,j,0,K_ind)) -= Fermi_energy;
			  real(H(i,1,j,1,K_ind)) -= Fermi_energy;
			}
		    }

		}
	    }

	  file_op.getline(line,2000);
	  l = line;

	  searched_string = "# imaginary part of H";
	  found=l.find(searched_string);
	  if (found!=string::npos)
	    {
	      
	      for(int i=0;i<Nb_Wannier_functions ;i++)
		{
		  file_op.getline(line,2000);
		  std::vector<double> hamil  = read_vector(Nb_Wannier_functions , &line[0]);
		  
		  for(int j=0;j<Nb_Wannier_functions ;j++)
		    {
		      imag(H(i,0,j,0,K_ind)) = hamil[j];
		      imag(H(i,1,j,1,K_ind)) = hamil[j];
		    }
		}
	      
	    }
	  
// 	  if(K_ind < 3)
// 	    {	  
// 	      cout.precision(5);
// 	      cout << scientific;
// 	      for(int i=0;i<Nb_Wannier_functions*2 ;i++){
// 		for(int j=0;j<Nb_Wannier_functions*2 ;j++){
// 		  cout << H(i,j,K_ind) << "\t";
// 		}
// 		cout << endl;
// 	      }
// 	      cout << endl;
// 	    }

	  K_ind++;
	}
    }
  
  // close files
  file_op.close();

  // renormalization Hartree --> eV !!
  H *= 27.2113838600; 
}



double Koshevnikov_parser::read_double(const char* p)
{
  char* pEnd;
  double result = strtod(p,&pEnd);
  return result;
}



std::vector<double> Koshevnikov_parser::read_vector(int Nb_doubles, const char* p )
{
  std::vector<double> result(Nb_doubles);

  char* pEnd;

  for(int i=0;i<Nb_doubles;i++)
    {
      result[i] = strtod(p,&pEnd);
      p = pEnd;
    }

  return result;
}

#endif 
