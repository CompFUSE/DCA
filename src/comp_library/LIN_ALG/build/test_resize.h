//-*-C++-*-

void test_resize()
{
    cout << __FUNCTION__ << endl;

    int Nx = 8;
    int Ny = 8;


    LIN_ALG::matrix<double, CPU> A(std::pair<int,int>(Nx, Ny));

    A.print();

    for(int i=0; i<Nx; ++i)
      for(int j=0; j<Ny; ++j)
	A(i,j) = i+j*Nx;

    A.print();
  
    LIN_ALG::matrix<double, GPU> B(A);
    B.print();

    A.resize(std::pair<int, int>(8,4));
    B.resize(std::pair<int, int>(8,4));
    
    assert(fabs(A.difference(B))<1.e-12);

    LIN_ALG::REMOVE<CPU>::row(A, 3);
    LIN_ALG::REMOVE<GPU>::row(B, 3);

    assert(fabs(A.difference(B))<1.e-12);

    LIN_ALG::REMOVE<CPU>::col(A, 2);
    LIN_ALG::REMOVE<GPU>::col(B, 2);

    assert(fabs(A.difference(B))<1.e-12);

    A.resize(std::pair<int, int>(10,11));
    B.resize(std::pair<int, int>(10,11));
    
    assert(fabs(A.difference(B))<1.e-12);

    A.print();
    B.print();
}
