//-*-C++-*-

void test_reading()
{
    
    int N  = 16;
    int Nv = 6;

    std::vector<int> v_indices(Nv);
    std::vector<int> M_r_indices(Nv);
    std::vector<int> M_c_indices(Nv);

    LIN_ALG::matrix<double, LIN_ALG::CPU> M(std::pair<int,int>(N, N));
    LIN_ALG::matrix<double, LIN_ALG::CPU> v(std::pair<int,int>(1,Nv));

    for(int i=0; i<Nv; ++i){
	M_r_indices[i] = N*double(rand())/double(RAND_MAX); 
	M_c_indices[i] = N*double(rand())/double(RAND_MAX); 
    
	v_indices[i] = i;
    }

    for(int i=0; i<N; ++i)
	for(int j=0; j<N; ++j)
	    M(i,j) = double(rand())/double(RAND_MAX);


    LIN_ALG::COPY_FROM<LIN_ALG::CPU, LIN_ALG::CPU>::read_at_indices(Nv, 
								    &M_r_indices[0], &M_c_indices[0], M.get_ptr(), M.get_global_size().first, 
								    &v_indices[0]                   , v.get_ptr(), v.get_global_size().first);


    M.print();
    v.print();

    thrust::device_vector<int> v_indices_GPU = v_indices;

    thrust::device_vector<int> M_r_GPU = M_r_indices;
    thrust::device_vector<int> M_c_GPU = M_c_indices;

    for(int l=0; l<Nv; l++){
	cout << M_r_indices[l] << "\t" << M_r_GPU[l] << "\n";
    } 

    LIN_ALG::matrix<double, LIN_ALG::GPU> M_GPU(M);
    LIN_ALG::matrix<double, LIN_ALG::GPU> v_GPU(std::pair<int,int>(1,Nv));

    LIN_ALG::COPY_FROM<LIN_ALG::GPU, LIN_ALG::GPU>::read_at_indices(Nv, 
								    &M_r_GPU[0], &M_c_GPU[0], M_GPU.get_ptr(), M_GPU.get_global_size().first, 
								    &v_indices_GPU[0]       , v_GPU.get_ptr(), v_GPU.get_global_size().first);

    M_GPU.print();
    v_GPU.print();
}
