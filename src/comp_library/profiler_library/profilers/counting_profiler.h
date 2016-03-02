//-*-C++-*-

#ifndef DCA_COUNTING_PROFILER_HEADER_H
#define DCA_COUNTING_PROFILER_HEADER_H

namespace PROFILER  
{
  /*! 
   *   \defgroup COUNTING-PROFILING
   *   \ingroup  PROFILING
   */

  /*! 
   *   \class    CountingProfiler
   *   \ingroup  COUNTING
   */
  template<typename time_event_type>
  class CountingProfiler 
  {

    const static int NB_COUNTERS = time_event_type::NB_COUNTERS;

    typedef typename time_event_type::scalar_type            scalar_type;

    typedef std::map<std::string,std::vector<scalar_type> >  profiling_table_type;
    typedef typename profiling_table_type::iterator          profiling_iterator_type;

    typedef std::map<std::string, std::string>  category_table_type;

    typedef typename category_table_type::iterator  category_table_iterator_type;

  public:

    static int MAX_THREADS();

    CountingProfiler(const char*        functionName_,
		     const char*        fileName_,
		     int                line);

    CountingProfiler(const char*        functionName_,
		     const char*        fileName_,
		     int                line,
		     int                thread_id);

    ~CountingProfiler();

    static void start();

    static void stop(std::string fileName);

    static void start_pthreading(int id);

    static void stop_pthreading(int id);


    template<typename concurrency_type>
    static void stop(concurrency_type& parallelProcessing, std::string fileName);

    static void to_TSV (std::string fileName); 
    static void to_JSON(std::string fileName);

  private:

    inline void finish();
    inline void constructionDone();
    
    static void merge_profiling_table();
    static void merge_categories();
    static void merge_categories_table();
    
    static inline profiling_table_type& get_profiling_table(int thread_id=0);

    static inline std::vector<std::string>& get_categories(int thread_id=0);

    static inline category_table_type& get_category_map(int thread_id=0);

    static inline std::vector<scalar_type>& start_running_counter(std::string name,
								  std::string category,
								  int         thread_id);

    static inline void print_counter(std::ostream& os, 
				     std::string   name_str,
				     std::vector<scalar_type>& counts);

    template<typename counter_value_type>
    static void get_counter_value(std::string         counter_name,
				  std::string         function_name,
				  counter_value_type& counter_value,
				  bool average=false);

  private:

    const std::string             functionName;

    int                           lineNumber;
      int                         thread_id;

    bool                          finished;
    
    time_event_type               startEvent;

    std::vector<std::string>      categories;

    bool                          using_GPU;
  };

  template<typename time_event_type>
  int CountingProfiler<time_event_type>::MAX_THREADS()
  {
    static int NB_THREADS = 32;
    return NB_THREADS;
  }

  template<typename time_event_type>
  CountingProfiler<time_event_type>::CountingProfiler(const char*        function_name,
						      const char*        category_name,
						      int                line):
    functionName      (function_name),
    lineNumber        (line),
    thread_id         (0),
    finished          (false),
    startEvent        (start_running_counter(function_name, category_name, thread_id), thread_id),
    categories(0)
  {
    assert(thread_id>-1 and thread_id<MAX_THREADS());
  }

  template<typename time_event_type>
  CountingProfiler<time_event_type>::CountingProfiler(const char*        function_name,
						      const char*        category_name,
						      int                line,
						      int                id):
    functionName      (function_name),
    lineNumber        (line),
    thread_id         (1+id),
    finished          (false),
    startEvent        (start_running_counter(function_name, category_name, thread_id), thread_id),
    categories(0)
  {
    assert(thread_id>-1 and thread_id<MAX_THREADS());
  }
  
  template<typename time_event_type>
  CountingProfiler<time_event_type>::~CountingProfiler() 
  {
    if (!finished) 
      finish();
  }
    
  template<typename time_event_type>
  void CountingProfiler<time_event_type>::finish() 
  {
    if(finished == true) 
      return; // already finished
    
    startEvent.end();

    finished = true;
  }

  template<typename time_event_type>
  typename CountingProfiler<time_event_type>::profiling_table_type& 
  CountingProfiler<time_event_type>::get_profiling_table(int thread_id) 
  {
    static std::vector<profiling_table_type> profiling_table(MAX_THREADS());
    return profiling_table[thread_id];
  }

  template<typename time_event_type>
  std::vector<std::string>& CountingProfiler<time_event_type>::get_categories(int thread_id)
  {
    static std::vector<std::vector<std::string> > categories(MAX_THREADS());
    return categories[thread_id];
  }

  template<typename time_event_type>
  std::map<std::string, std::string>& CountingProfiler<time_event_type>::get_category_map(int thread_id)
  {
    static std::vector<category_table_type> category_table(MAX_THREADS());
    return category_table[thread_id];
  }
  
  template<typename time_event_type>
  std::vector<typename CountingProfiler<time_event_type>::scalar_type>& 
  CountingProfiler<time_event_type>::start_running_counter(std::string name,
							   std::string category,
							   int         thread_id) 
  {
    size_t number_of_events = NB_COUNTERS;
    
    std::vector<scalar_type>& counts = get_profiling_table(thread_id)[name];

    if(counts.size() < number_of_events)
      {
	get_categories  (thread_id).push_back(category);
	get_category_map(thread_id)[name] = category;
	
	counts.resize(number_of_events, 0);
      }
         
    return counts;
  }

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::start_pthreading(int id)
  {
    time_event_type::start_pthreading(id+1);
  }

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::stop_pthreading(int id)
  {
    time_event_type::stop_pthreading(id+1);
  }
      
  template<typename time_event_type>
  void CountingProfiler<time_event_type>::start()
  {
    time_event_type::start();
  }

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::stop(std::string file_name) 
  {
    time_event_type::stop();
        
    merge_profiling_table();
    merge_categories();
    merge_categories_table();

    {
      to_JSON(file_name);
      
      to_TSV(file_name);
    }
  }

  template<typename time_event_type>
  template<typename concurrency_type>
  void CountingProfiler<time_event_type>::stop(concurrency_type& concurrency,
					       std::string       file_name) 
  {
    time_event_type::stop();
        
    merge_profiling_table();
    merge_categories();
    merge_categories_table();

    profiling_table_type& profiling_table = get_profiling_table(0);

    concurrency.sum(profiling_table);

    if(concurrency.id() == concurrency.last()) 
      {
	to_JSON(file_name);
	
	to_TSV(file_name);
    }
  }

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::merge_profiling_table()
  {
    size_t                number_of_events  = NB_COUNTERS;//time_event_type::names().size();
    profiling_table_type& profiling_table_0 = get_profiling_table(0);
    
    for(int id=1; id<MAX_THREADS(); ++id)
      {
	profiling_table_type& profiling_table_i = get_profiling_table(id);
	
	if(profiling_table_i.size()>0)
	  {
	    for(profiling_iterator_type itr = profiling_table_i.begin(); itr != profiling_table_i.end(); ++itr) 
	      {		
		std::string name = ((*itr).first);
		
		std::vector<scalar_type>& counts = profiling_table_0[name];
		
		if(counts.size() < number_of_events)
		  profiling_table_0[name] = ((*itr).second);
		else
		  for(size_t i=0; i<number_of_events; ++i)
		    profiling_table_0[name][i] += ((*itr).second)[i];
	      }
	    
// 	    cout << "\n\t merging : " << id << "\n";
// 	    for(profiling_iterator_type itr = profiling_table_0.begin(); itr != profiling_table_0.end(); ++itr) 
// 	      cout << "\t\t" << ((*itr).first) << "\t" << ((*itr).second)[0] <<"\n";
// 	    cout << "\n";
	  }
      }
  }    

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::merge_categories()
  {
    std::vector<std::string>& categories_0 = get_categories(0); 
    
    for(int id=1; id<MAX_THREADS(); ++id)
      for(size_t l=0; l<get_categories(id).size(); ++l)
	categories_0.push_back(get_categories(id)[l]); 
  }

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::merge_categories_table()
  {
    category_table_type& category_table_0 = get_category_map(0);
    
    for(int id=1; id<MAX_THREADS(); ++id)
      {
	category_table_type& category_table_i = get_category_map(id);
	
	if(category_table_i.size()>0)
	  for(category_table_iterator_type itr = category_table_i.begin(); itr != category_table_i.end(); ++itr) 
	    {		
	      std::string name = ((*itr).first);
	      category_table_0[name] = ((*itr).second);
	    }
	  
      }
  }    

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::to_JSON(std::string fileName) 
  {
    profiling_table_type& profiling_table = get_profiling_table();

    std::ofstream fileStream(fileName.c_str());

    int size  = profiling_table.size();

    fileStream << "{";
    
    int index = 0;
    for(profiling_iterator_type itr = profiling_table.begin(); itr != profiling_table.end(); ++itr) 
      {
	const std::string&              function_name = ((*itr).first);
	const std::vector<scalar_type>& counters      = ((*itr).second);
	
	std::vector<scalar_type> counter_normalized = counters;
	
	time_event_type::normalize(counter_normalized);
	
	fileStream << "\"" << index << "\" : ";
	print_counter(fileStream, function_name, counter_normalized);
	
	if(index==size-1)
	  fileStream << " \n";
	else
	  fileStream << ",\n";

	index += 1;
      }
    
    fileStream << "}";
    
    fileStream.close();
  }

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::to_TSV(std::string fileName) 
  {
    profiling_table_type& profiling_table = get_profiling_table();

    std::string file = fileName;
    file = file + std::string(".txt");

    std::ofstream of(file.c_str());

    of<<std::fixed;
    of.precision(6);

    {
      std::sort(get_categories().begin(), get_categories().end());    
      get_categories().resize(std::unique(get_categories().begin(), get_categories().end()) - get_categories().begin());
    }

    for(size_t l=0; l<get_categories().size(); ++l)
      {
	std::string category_name = get_categories()[l];

	std::vector<double>      number_of_counts (0);
	std::vector<double>      number_of_threads(0);
	std::vector<double>      function_times(0);
	std::vector<std::string> function_names(0);

	size_t max_size = 0;
	for(profiling_iterator_type itr = profiling_table.begin(); itr != profiling_table.end(); ++itr) 
	  {
	    const std::string&              function_name = ((*itr).first);
	    const std::vector<scalar_type>& counters      = ((*itr).second);

	    if(get_category_map()[function_name] == category_name)
	      {
		if(max_size<function_name.size())
		  max_size = function_name.size();

		number_of_counts .push_back(counters[6]);
		number_of_threads.push_back(counters[7]);

		function_names.push_back(function_name);
		function_times.push_back(double(counters[0])+1.e-6*double(counters[1]));
	      }
	  }

	double total_time=1.e-6;
	for(size_t l=0; l<function_times.size(); ++l){
	  function_names[l].resize(int(max_size*1.1), ' ');

	  total_time += function_times[l];
	}

	of << "\n\n\t" + category_name + " :\t" << total_time << "\n";
	for(size_t l=0; l<function_times.size(); ++l)
	  of << "\t" << function_names[l] 
	     << " | \t" << function_times[l]/number_of_threads[l] 
	     << " | \t" << function_times[l]/total_time 
	     << "\n";
      }
	
    of.close();
  }

  template<typename time_event_type>
  void CountingProfiler<time_event_type>::print_counter(std::ostream&             os, 
							std::string               name,
							std::vector<scalar_type>& counts) 
  {
    const std::vector<std::string> names = time_event_type::names();
    
    os << "{\n";
    
    os << "\"name\"     : \"" <<                    name  << "\",\n";
    os << "\"category\" : \"" << get_category_map()[name] << "\",\n";
    
    for (size_t i=0; i<names.size(); i++){ 
      os << "\"" << names[i] << "\" : " << counts[i];

      if(i==names.size()-1)
	os << "\n";
      else
	os << ",\n";
    }

    os << "}";
  }
  
}//namespace 

#endif



