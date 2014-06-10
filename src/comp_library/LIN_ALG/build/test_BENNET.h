//-*-C++-*-

void do_GETRF(LIN_ALG::matrix<double, LIN_ALG::CPU>& A, LIN_ALG::matrix<double, LIN_ALG::CPU>& A_LU)
{
    LIN_ALG::matrix<double, LIN_ALG::CPU> L(A);
    LIN_ALG::matrix<double, LIN_ALG::CPU> U(A);

    int INFO=0;
    int IPIV[A_LU.get_current_size().first];

    int m = A_LU.get_current_size().first;
    int n = A_LU.get_current_size().second;
	
    LIN_ALG::GETRF<LIN_ALG::CPU>::execute(m, n, A_LU.get_ptr(), A_LU.get_global_size().first, IPIV, INFO);
	
    for(int l=0; l<A.get_current_size().first; l++)
	LIN_ALG::SWAP<LIN_ALG::CPU>::row(A, l, IPIV[l]-1);
    
    for(int i=0; i<A_LU.get_current_size().first; ++i){
	for(int j=0; j<A_LU.get_current_size().second; ++j){
	    if(i<j){
		U(i,j) = A_LU(i,j);
		L(i,j) = 0;
	    }
	    else{
		if(i==j){
		    L(i,i)=1;
		    U(i,i)=A_LU(i,i);
		}
		else{
		    U(i,j) = 0;
		    L(i,j) = A_LU(i,j);   
		}
	    }
	}
    }

    LIN_ALG::matrix<double, LIN_ALG::CPU> A_check(A);
    LIN_ALG::GEMM<LIN_ALG::CPU>::execute(1., L, U, 0., A_check);

    cout << A.difference(A_check) << endl;
}

void test_BENNET()
{
    cout << __FUNCTION__ << endl;

    int N=16;

    LIN_ALG::matrix<double, LIN_ALG::CPU> A(std::pair<int,int>(N, N));

    for(int i=0; i<N; ++i)
	for(int j=0; j<N; ++j)
	    A(i,j) = double(rand())/double(RAND_MAX);

    LIN_ALG::matrix<double, LIN_ALG::CPU> A_LU(A);
    LIN_ALG::matrix<double, LIN_ALG::CPU> C(std::pair<int,int>(N,1));
    LIN_ALG::matrix<double, LIN_ALG::CPU> R(std::pair<int,int>(1,N));

    do_GETRF(A, A_LU);
    
    for(int i=0; i<N; ++i){
	R(0,i) = double(rand())/double(RAND_MAX);
	C(i,0) = double(rand())/double(RAND_MAX);
    }

    LIN_ALG::matrix<double, LIN_ALG::CPU> B(A);

    for(int i=0; i<N; ++i)
	for(int j=0; j<N; ++j)
	    B(i,j) += C(i,0)*R(0,j);

    LIN_ALG::matrix<double, LIN_ALG::CPU> B_LU(B);

    do_GETRF(B, B_LU);

    //B_LU.print();

    LIN_ALG::matrix<double, LIN_ALG::GPU> A_GPU_LU(A_LU);
    LIN_ALG::matrix<double, LIN_ALG::GPU> C_GPU(C);
    LIN_ALG::matrix<double, LIN_ALG::GPU> R_GPU(R);

    LIN_ALG::BENNET<LIN_ALG::CPU>::execute(A_LU, C.get_ptr(), R.get_ptr());
    //A_LU.print();

    LIN_ALG::BENNET<LIN_ALG::GPU>::execute(A_GPU_LU, C_GPU.get_ptr(), R_GPU.get_ptr());
    //A_GPU_LU.print();

    cout << A_LU.difference(A_GPU_LU) << endl;
}
