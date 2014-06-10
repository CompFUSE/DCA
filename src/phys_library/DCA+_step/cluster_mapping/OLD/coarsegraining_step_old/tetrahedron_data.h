//-*-C++-*-

#ifndef DCA_TETRAHEDRON_INTEGRATION_DATA_H
#define DCA_TETRAHEDRON_INTEGRATION_DATA_H

namespace DCA
{

  template<typename scalar_type>
  class tetrahedron_integration_data
  {
  public:
    
    tetrahedron_integration_data(int N);
    ~tetrahedron_integration_data();

  public:

    std::complex<scalar_type>* G_inv_0;
    std::complex<scalar_type>* G_inv_1;
    std::complex<scalar_type>* G_inv_2;
    std::complex<scalar_type>* G_inv_3;

    std::complex<scalar_type>* VR_0;
    std::complex<scalar_type>* VR_1;
    std::complex<scalar_type>* VR_2;
    std::complex<scalar_type>* VR_3;

    std::complex<scalar_type>* VR_inv_0;
    std::complex<scalar_type>* VR_inv_1;
    std::complex<scalar_type>* VR_inv_2;
    std::complex<scalar_type>* VR_inv_3;

    std::complex<scalar_type>* W_0;
    std::complex<scalar_type>* W_1;
    std::complex<scalar_type>* W_2;
    std::complex<scalar_type>* W_3;
  };

  template<typename scalar_type>
  tetrahedron_integration_data<scalar_type>::tetrahedron_integration_data(int N):
    G_inv_0(NULL),
    G_inv_1(NULL),
    G_inv_2(NULL),
    G_inv_3(NULL),

    VR_0(NULL),
    VR_1(NULL),
    VR_2(NULL),
    VR_3(NULL),

    VR_inv_0(NULL),
    VR_inv_1(NULL),
    VR_inv_2(NULL),
    VR_inv_3(NULL),

    W_0(NULL),
    W_1(NULL),
    W_2(NULL),
    W_3(NULL)
  {
    G_inv_0 = new std::complex<scalar_type>[N*N];
    G_inv_1 = new std::complex<scalar_type>[N*N];
    G_inv_2 = new std::complex<scalar_type>[N*N];
    G_inv_3 = new std::complex<scalar_type>[N*N];

    VR_0 = new std::complex<scalar_type>[N*N];
    VR_1 = new std::complex<scalar_type>[N*N];
    VR_2 = new std::complex<scalar_type>[N*N];
    VR_3 = new std::complex<scalar_type>[N*N];

    VR_inv_0 = new std::complex<scalar_type>[N*N];
    VR_inv_1 = new std::complex<scalar_type>[N*N];
    VR_inv_2 = new std::complex<scalar_type>[N*N];
    VR_inv_3 = new std::complex<scalar_type>[N*N];

    W_0 = new std::complex<scalar_type>[N];
    W_1 = new std::complex<scalar_type>[N];
    W_2 = new std::complex<scalar_type>[N];
    W_3 = new std::complex<scalar_type>[N];
  }

  template<typename scalar_type>
  tetrahedron_integration_data<scalar_type>::~tetrahedron_integration_data()
  {
    delete [] G_inv_0;
    delete [] G_inv_1;
    delete [] G_inv_2;
    delete [] G_inv_3;

    delete [] VR_0;
    delete [] VR_1;
    delete [] VR_2;
    delete [] VR_3;

    delete [] VR_inv_0;
    delete [] VR_inv_1;
    delete [] VR_inv_2;
    delete [] VR_inv_3;

    delete [] W_0;
    delete [] W_1;
    delete [] W_2;
    delete [] W_3;
  }

}

#endif
