//-*-C++-*-

void test_operator()
{
    int N=1;

    LIN_ALG::matrix<double, GPU> A(N);

    clock_t start = clock();

    for(int l=0; l<10000; ++l)
	for(int i=0; i<N; ++i){
	    for(int j=0; j<N; ++j){
		A.set(A.get_ptr(i,j), double(j));
	    }
	}

    clock_t end = clock();
    
    cout << double(end-start)/double(CLOCKS_PER_SEC)/10000. << endl;
    //A.print();
}
