//-*-C++-*-

void test_GEMM()
{
    cout << __FUNCTION__ <<endl;

    double time_CPU, time_GPU;

    int K  = 256;
    int Nx = 512;

    LIN_ALG::matrix<double, LIN_ALG::CPU> A1(std::pair<int,int>(Nx, K));    
    LIN_ALG::matrix<double, LIN_ALG::CPU> A2(std::pair<int,int>(K, Nx));
    LIN_ALG::matrix<double, LIN_ALG::CPU> A3(std::pair<int,int>(Nx, Nx));
    
    for(int i=0; i<Nx; ++i){
	for(int j=0; j<K; ++j){
	    A1(i,j) = rand();
	    A2(j,i) = j;
	}
    }
    
    LIN_ALG::matrix<double, LIN_ALG::GPU> B1(A1);
    LIN_ALG::matrix<double, LIN_ALG::GPU> B2(A2);
    LIN_ALG::matrix<double, LIN_ALG::GPU> B3(std::pair<int,int>(Nx, Nx));
    
    {
	clock_t start = clock();

	for(int l=0; l<100; l++)
	    LIN_ALG::GEMM<LIN_ALG::CPU>::execute(1., A1, A2, 1., A3);
	clock_t end = clock();
    
	time_CPU = double(end-start)/double(CLOCKS_PER_SEC)/100.;
    }

    {
	clock_t start = clock();

	for(int l=0; l<100; l++){
	    LIN_ALG::GEMM<LIN_ALG::GPU>::execute(1., B1, B2, 1., B3);
	    cuda_thread_synchronize();
	}

	clock_t end = clock();
    
	time_GPU = double(end-start)/double(CLOCKS_PER_SEC)/100.;
    }

    cout<<scientific;
    cout.precision(6);
    cout << "\t" << time_CPU << "\t" << time_GPU << "\t" << time_CPU/time_GPU << "\n";

    cout << A3.difference(B3) << endl;
}
