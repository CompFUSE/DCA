//-*-C++-*-

void test_GEINV()
{
    int N = 8;

    LIN_ALG::matrix<double, CPU> A(N);
    
    for(int i=0; i<N; ++i)
	for(int j=0; j<N; ++j)
	  A(i,j) = double(rand())/double(RAND_MAX);

    A.print();

    LIN_ALG::matrix<double, CPU> A_cpy(A);
    LIN_ALG::matrix<double, CPU> U(N);

    LIN_ALG::GEINV<CPU>::execute(A);

    LIN_ALG::GEMM<CPU>::execute(1., A_cpy, A, 0., U);

    U.print();
}
