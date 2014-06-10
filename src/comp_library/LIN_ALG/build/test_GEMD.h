//-*-C++-*-

void test_GEMD()
{
    cout << __FUNCTION__ << endl;

    int m = 1000;
    int n = 1000;

    LIN_ALG::matrix<double, LIN_ALG::CPU> M(std::pair<int,int>(m, n));
    LIN_ALG::matrix<double, LIN_ALG::CPU> D(std::pair<int,int>(n, 1));
    LIN_ALG::matrix<double, LIN_ALG::CPU> N(std::pair<int,int>(m, n));

    for(int i=0; i<n; ++i){
	
	for(int j=0; j<m; ++j)
	    M(j,i) = 1;
	
	D(i,0) = i;
    }
    
    {
	clock_t start = clock();

	for(int l=0; l<10; l++)
	    LIN_ALG::GEMD<LIN_ALG::CPU>::execute(M, D.get_ptr(), N);
	clock_t end = clock();
    
	cout << double(end-start)/double(CLOCKS_PER_SEC)/10. << endl;
    }
    
    LIN_ALG::matrix<double, LIN_ALG::GPU> M_GPU(M);
    LIN_ALG::matrix<double, LIN_ALG::GPU> D_GPU(D);
    LIN_ALG::matrix<double, LIN_ALG::GPU> N_GPU(std::pair<int,int>(m, n));
    
    {
	clock_t start = clock();

	for(int l=0; l<100; l++){
	    LIN_ALG::GEMD<LIN_ALG::GPU>::execute(M_GPU, D_GPU.get_ptr(), N_GPU);
	    cuda_thread_synchronize();
	}
	clock_t end = clock();
    
	cout << double(end-start)/double(CLOCKS_PER_SEC)/100. << endl;
    }

    cout << N.difference(N_GPU) << endl;
}
