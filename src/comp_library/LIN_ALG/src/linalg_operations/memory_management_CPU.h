//-*-C++-*-

#ifndef LIN_ALG_MEMORY_MANAGEMENT_CPU_H
#define LIN_ALG_MEMORY_MANAGEMENT_CPU_H

namespace LIN_ALG {

    namespace MEMORY_MANAGEMENT_ON_GPU {

      template<typename scalartype> void allocate_pinned_host_memory(scalartype*& ptr, int global_size);
      template<typename scalartype> void allocate_pinned_host_memory(scalartype*& ptr, std::pair<int,int> global_size);

      template<typename scalartype> void deallocate_pinned_host_memory(scalartype*& ptr);
    }

    template<>
    class MEMORY_MANAGEMENT<CPU>
    {

    public:

	template<typename scalartype>
	inline static scalartype& get(scalartype* ptr){
	    return ptr[0];
	}

	template<typename scalartype>
	inline static scalartype& get(scalartype* ptr, int index){
	    return ptr[index];
	}

	template<typename scalartype>
	inline static void set(scalartype* ptr, scalartype val){
	    ptr[0] = val;
	}

	template<typename scalartype>
	inline static void add(scalartype* ptr, scalartype val){
	    ptr[0] += val;
	}
	
	template<typename scalartype>
	inline static void allocate(scalartype*& ptr, int global_size){
	  assert(ptr==NULL);

#ifdef ENABLE_PINNED_MEMORY_ALLOCATION
	  MEMORY_MANAGEMENT_ON_GPU::allocate_pinned_host_memory(ptr, global_size);
#else
	  posix_memalign((void**) &ptr, 128, global_size*sizeof(scalartype));
#endif

	  assert(ptr!=NULL);
	}

	template<typename scalartype>
	inline static void allocate(scalartype*& ptr, std::pair<int,int>& global_size){
	  assert(ptr==NULL);

#ifdef ENABLE_PINNED_MEMORY_ALLOCATION
	  MEMORY_MANAGEMENT_ON_GPU::allocate_pinned_host_memory(ptr, global_size);
#else
	  posix_memalign((void**) &ptr, 128, global_size.first*global_size.second*sizeof(scalartype));
#endif

	  assert(ptr!=NULL);
	}
	
	template<typename scalartype>
	inline static void deallocate(scalartype*& ptr){
	  assert(ptr!=NULL);	  

#ifdef ENABLE_PINNED_MEMORY_ALLOCATION
	  MEMORY_MANAGEMENT_ON_GPU::deallocate_pinned_host_memory(ptr);
#else
	  free(ptr);
#endif

	  ptr = NULL;
	}
	
      /*
	template<typename scalartype>
	inline static void memcopy(scalartype* target_ptr, scalartype* source_ptr, int size){
	    memcpy(target_ptr, source_ptr, sizeof(scalartype)*size);
	}
      */

	template<typename scalartype>
	inline static void set_to_zero(scalartype* ptr, int size){
	    for(int l=0; l<size; ++l)  
		ptr[l] = scalartype(0);
	}

	template<typename scalartype>
	inline static void set_to_zero(scalartype* ptr, int LD, int size){
	    for(int l=0; l<size; ++l)  
		ptr[l*LD] = scalartype(0);
	}

	template<typename scalartype>
	static void print(scalartype* ptr, int c_s, int g_s){
	
	    cout.precision(6);
	    cout<<scientific;
	    
	    cout << "\n\n";
	    cout << "\t current-size : " << c_s << "\n";
	    cout << "\t global -size : " << g_s << "\n";
	    cout << "\n\n";
	    
	    for(int i=0; i<c_s; i++)  
		cout << "\t" << ptr[i];
	    cout << "\n";
	    	    
	    cout << "\n\n\n";
	}

	template<typename scalartype>
	static void print(scalartype*& ptr, std::pair<int, int>& c_s, std::pair<int, int>& g_s){
	    cout.precision(6);
	    cout<<scientific;
	    
	    cout << "\n\n";
	    cout << "\t current-size : " << c_s.first << "\t" << c_s.second << "\n";
	    cout << "\t global -size : " << g_s.first << "\t" << g_s.second << "\n";
	    cout << "\n\n";
	    
	    for(int i=0; i<c_s.first; i++){  
		for(int j=0; j<c_s.second; j++)  
		    cout << "\t" << ptr[i + g_s.first*j];
		cout << "\n";
	    }
	    
	    cout << "\n\n\n";
	}
    };
}
#endif
