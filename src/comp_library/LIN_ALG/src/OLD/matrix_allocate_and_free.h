//-*-C++-*-

#ifndef QMC_MATRIX_MEMORY_MANAGEMENT_H
#define QMC_MATRIX_MEMORY_MANAGEMENT_H

template<device_type device_name>
struct MATRIX_MEMORY_MANAGEMENT
{
    template<typename scalartype>
    static scalartype& get(scalartype*& ptr, int index){
	return ptr[index];
    }

    template<typename scalartype>
    static void allocate(scalartype*& ptr, std::pair<int,int>& global_size){
	assert(ptr==NULL);
	
	cout << ptr << endl;

	ptr = (scalartype*) malloc(global_size.first*global_size.second*sizeof(scalartype));
    
	cout << ptr << endl;
    }

    template<typename scalartype>
    static void deallocate(scalartype*& ptr){
	//cout << ptr << endl;
	
	free(ptr);
    }

    template<typename scalartype>
    static void memcopy(scalartype*& target_ptr, scalartype*& source_ptr, int size){
	memcpy(target_ptr, source_ptr, sizeof(scalartype)*size);
    }

    template<typename scalartype>
    static void set_to_zero(scalartype*& ptr, int size){
	cout << ptr << endl;
	
	for(int l=0; l<size; l++)  
	    ptr[l] = scalartype(0);
    }

    template<typename scalartype>
    static void print(scalartype*& ptr, std::pair<int, int>& c_s, std::pair<int, int>& g_s){

	cout << ptr << endl;

	cout.precision(6);
	cout<<scientific;

	cout << "\n\n";
	cout << "\t current-size : " << c_s.first << "\t" << c_s.second << "\n";
	cout << "\t global -size : " << g_s.first << "\t" << g_s.second << "\n";
	cout << "\n\n";

	for(int i=0; i<c_s.second; i++){  

	    cout << "\n";
	    for(int j=0; j<c_s.first; j++)  
		cout << "\t" << ptr[i + g_s.first*j];
	}

	cout << "\n\n\n";
    }
};

template<>
struct MATRIX_MEMORY_MANAGEMENT<GPU>
{
    template<typename scalartype>
    static void allocate(scalartype* ptr, std::pair<int,int>& global_size){
        assert(ptr==NULL);
	cudaError_t ret = cudaMalloc((void**) ptr, global_size.first*global_size.second*sizeof(scalartype));
    
	if( ret != cudaSuccess)	    
	    std::cout << "NOT ENOUGH GPU MEMORY :: ret code: \n\n\n" << ret << std::endl;
    }

    template<typename scalartype>
    static void deallocate(scalartype* ptr){
        cudaFree(ptr);
    }

    template<typename scalartype>
    static void memcopy(scalartype* target_ptr, scalartype* source_ptr, int size){
        cudaMemcpy(target_ptr, source_ptr, sizeof(scalartype)*size);
    }

    template<typename scalartype>
    static void set_to_zero(scalartype*& ptr, int size){

	dim3 dim_G(128);
	dim3 dim_B(128);
	
	set_to_zero_kernel<<<dim_G,dim_B>>>(ptr, size);
        //for(int l=0; l<size; l++)
	//ptr[l] = scalartype(0);
    }

    template<typename scalartype>
    static void print(scalartype*& ptr, std::pair<int, int>& c_s, std::pair<int, int>& g_s){
/*
	cudaMemcpy(A_d, A_h, sizeof(float) * VAR, cudaMemcpyHostToDevice);

	cout.precision(6);
        cout<<scientific;

        cout << "\n\n";
        cout << "\t current-size : " << c_s.first << "\t" << c_s.second << "\n";
        cout << "\t global -size : " << g_s.first << "\t" << g_s.second << "\n";
        cout << "\n\n";

        for(int i=0; i<c_s.second; i++){

            cout << "\n";
            for(int j=0; j<c_s.first; j++)
                cout << "\t" << ptr[i + g_s.first*j];
        }

        cout << "\n\n\n";
*/  
    }
    
};


#endif
