//-*-C++-*-

#ifndef KOSHEVNIKOV_MODEL_H
#define KOSHEVNIKOV_MODEL_H

#include "Koshevnikov_parser.h"

/*!
 *  \author peter staar
 */
class Koshevnikov_model 
{

public:
  
  const static int DIMENSION = 3;

  const static cluster_shape_type DCA_cluster_shape = BETT_CLUSTER;
  const static cluster_shape_type LDA_cluster_shape = PARALLELEPIPED;

  typedef no_symmetry<DIMENSION> LDA_point_group;
  typedef no_symmetry<DIMENSION> DCA_point_group;

  static int    BANDS;
  static string filename;

  static double Fermi_energy;
  static double density;

  template<class parameters_type>
  static void initialize(parameters_type& parameters);
  
  template<class concurrency_type>
  static void read(concurrency_type& concurrency);

  static std::vector<int>& DCA_grid_size();
  static std::vector<int>& LDA_grid_size();

  static double* get_r_DCA_basis();
  static double* get_k_DCA_basis();

  static double* get_r_LDA_basis();
  static double* get_k_LDA_basis();

  static std::vector<int>                  get_flavors();
  static std::vector<std::vector<double> > get_a_vectors();

  static int* Hamiltonian_symmetries();

  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > get_orbital_permutations();

  template<class domain, class parameters_type>
  static void                 initialize_H_LDA(FUNC_LIB::function<std::complex<double> , domain >& H_LDA,
					       parameters_type&                          parameters);

  template<class domain, class parameters_type>
  static void                 initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interactions,
						       parameters_type&            parameters);

  template<class domain>
  static void                 initialize_H_symmetries(FUNC_LIB::function<int , domain >& H_symmetries);

  static void                 symmetrize_Hamiltonian(std::complex<double>* H_matrix);

  template<class concurrency_type>
  static int get_buffer_size(concurrency_type& concurrency);

  template<class concurrency_type>
  static void pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  static void unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);
};

int    Koshevnikov_model::BANDS = -1;
string Koshevnikov_model::filename = "";

double Koshevnikov_model::Fermi_energy = 0.;
double Koshevnikov_model::density = 0.;

template<class parameters_type>
void Koshevnikov_model::initialize(parameters_type& parameters)
{  
  Koshevnikov_model::filename = parameters.get_H_LDA_file();
  Koshevnikov_model::BANDS    = parameters.get_LDA_bands();
}

template<class concurrency_type>
void Koshevnikov_model::read(concurrency_type& concurrency)
{
  if(concurrency.id() == 0)
    {
      int bands = Koshevnikov_parser::read_Nb_Wannier_functions(filename);
      
      if( bands != BANDS)
	throw std::logic_error("LDA_BANDS is not correct !!");

      Koshevnikov_model::Fermi_energy = Koshevnikov_parser::read_Fermi_energy(filename);
      
      Koshevnikov_parser::read_grid_size(filename, Koshevnikov_model::LDA_grid_size());
      
      Koshevnikov_parser::read_DCA_k_basis(filename, Koshevnikov_model::get_k_DCA_basis());
      Koshevnikov_parser::read_DCA_r_basis(filename, Koshevnikov_model::get_r_DCA_basis());
      
      Koshevnikov_parser::read_k_basis(filename, Koshevnikov_model::get_k_LDA_basis());
      Koshevnikov_parser::read_r_basis(filename, Koshevnikov_model::get_r_LDA_basis());
    }

  Koshevnikov_model Koshevnikov;
  concurrency.broadcastObj(Koshevnikov);
}


 std::vector<int>& Koshevnikov_model::DCA_grid_size()
{
  static std::vector<int> grid(3,0);
  return grid;
}

 std::vector<int>& Koshevnikov_model::LDA_grid_size()
{
  static std::vector<int> grid(3,0);
  return grid;
}

 double* Koshevnikov_model::get_r_DCA_basis()
{
  static double* basis = new double[3*3];
  return basis;
}

 double* Koshevnikov_model::get_k_DCA_basis()
{
  static double* basis = new double[3*3];
  return basis;
}


 double* Koshevnikov_model::get_r_LDA_basis()
{
  static double* basis = new double[3*3];
  return basis;
}

double* Koshevnikov_model::get_k_LDA_basis()
{
  static double* basis = new double[3*3];
  return basis;
}

std::vector<int> Koshevnikov_model::get_flavors()
{
  static std::vector<int> flavors(BANDS, 0.);

  static bool tmp=true;
  if(tmp){
    //cout << "\n\n\n\t INITIALIZE this !!! " << __PRETTY_FUNCTION__<< "\n\n\n";
    
    for(int l=0; l<BANDS; ++l)
      flavors[l] = l;

    tmp = false;
  }

  return flavors;
}

std::vector<std::vector<double> > Koshevnikov_model::get_a_vectors()
{
  static std::vector<std::vector<double> > a_vecs(BANDS, std::vector<double>(DIMENSION,0.));

  static bool tmp=true;
  if(tmp){
    //cout << "\n\n\n\t INITIALIZE this !!! " << __PRETTY_FUNCTION__<< "\n\n\n";
    tmp = false;
  }

  return a_vecs;
}


int* Koshevnikov_model::Hamiltonian_symmetries()
{
  static bool tmp=true;
  if(tmp){
    cout << __FUNCTION__ << "\t" << Koshevnikov_model::BANDS << endl;
    tmp = false;
  }

  static int* symm = new int[(2*BANDS)*(2*BANDS)];
  return symm;
}

std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > Koshevnikov_model::get_orbital_permutations()
{
  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > permutations(0);
  return permutations;
}

template<class concurrency_type>
int Koshevnikov_model::get_buffer_size(concurrency_type& concurrency)
{
  int result = 0;

  result += concurrency.getBufferSize(BANDS);
  result += concurrency.getBufferSize(filename);
  result += concurrency.getBufferSize(Fermi_energy);

  result += concurrency.getBufferSize(DCA_grid_size()); //cout << result << endl;
  result += concurrency.getBufferSize(LDA_grid_size()); //cout << result << endl;

  result += concurrency.getBufferSize(std::vector<double>(9)); //cout << result << endl;
  result += concurrency.getBufferSize(std::vector<double>(9)); //cout << result << endl;

  result += concurrency.getBufferSize(std::vector<double>(9)); //cout << result << endl;
  result += concurrency.getBufferSize(std::vector<double>(9)); //cout << result << endl;

  return result;
}

template<class concurrency_type>
void Koshevnikov_model::pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, BANDS);
  concurrency.pack(buffer, buffer_size, position, filename);
  concurrency.pack(buffer, buffer_size, position, Fermi_energy);
  
  concurrency.pack(buffer, buffer_size, position, DCA_grid_size());
  concurrency.pack(buffer, buffer_size, position, LDA_grid_size());
  
  std::vector<double> tmp(9);
  memcpy(&tmp[0], get_r_DCA_basis(), sizeof(double)*9);
  concurrency.pack(buffer, buffer_size, position, tmp);
  
  memcpy(&tmp[0], get_k_DCA_basis(), sizeof(double)*9);
  concurrency.pack(buffer, buffer_size, position, tmp);
  
  memcpy(&tmp[0], get_r_LDA_basis(), sizeof(double)*9);
  concurrency.pack(buffer, buffer_size, position, tmp);
  
  memcpy(&tmp[0], get_k_LDA_basis(), sizeof(double)*9);
  concurrency.pack(buffer, buffer_size, position, tmp);
}

template<class concurrency_type>
void Koshevnikov_model::unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, BANDS);
  concurrency.unpack(buffer, buffer_size, position, filename);
  concurrency.unpack(buffer, buffer_size, position, Fermi_energy);

  concurrency.unpack(buffer, buffer_size, position, DCA_grid_size());
  concurrency.unpack(buffer, buffer_size, position, LDA_grid_size());

  std::vector<double> tmp(9);
  
  concurrency.unpack(buffer, buffer_size, position, tmp);
  memcpy(get_r_DCA_basis(), &tmp[0], sizeof(double)*9);
  
  concurrency.unpack(buffer, buffer_size, position, tmp);
  memcpy(get_k_DCA_basis(), &tmp[0], sizeof(double)*9);
  
  concurrency.unpack(buffer, buffer_size, position, tmp);
  memcpy(get_r_LDA_basis(), &tmp[0], sizeof(double)*9);
  
  concurrency.unpack(buffer, buffer_size, position, tmp);
  memcpy(get_k_LDA_basis(), &tmp[0], sizeof(double)*9);
}

template<class domain, class parameters_type>
void Koshevnikov_model::initialize_H_LDA(FUNC_LIB::function<std::complex<double> , domain >& H_LDA,
					 parameters_type&                          parameters)
{
  Koshevnikov_parser::read_LDA_Hamiltonians(filename, H_LDA, Fermi_energy);
}

template<class domain, class parameters_type>
void Koshevnikov_model::initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
						 parameters_type&            parameters)
{
  Koshevnikov_parser::read_interaction_Hamiltonian(filename, H_interaction, BANDS);
}

template<class domain>
void Koshevnikov_model::initialize_H_symmetries(FUNC_LIB::function<int , domain >& H_symmetries)
{
  Koshevnikov_parser::read_symmetries(filename, H_symmetries, BANDS);
}

void Koshevnikov_model::symmetrize_Hamiltonian(std::complex<double>* H_matrix)
{
  // spin -symmetrization

  for(int i=0; i<BANDS; i++){
    for(int j=0; j<BANDS; j++){
      std::complex<double> tmp = (H_matrix[i+j*2*BANDS] + H_matrix[(i+BANDS)+(j+BANDS)*2*BANDS])/2.;
      H_matrix[i+j*2*BANDS]                 = tmp;
      H_matrix[(i+BANDS)+(j+BANDS)*2*BANDS] = tmp; 
    }
  }
}

#endif
