## Random number library

We provide three types of random number generators

* a wrapper class for the random number library of C++11
* a wrapper class for the SPRNG library
* Ranq2 from Numerical Recipes (not to be released)

More details on them can be found below.  
However, any other type of random number generator can also be used as long as it uses the interface shared by these three.

### Interface

##### Contructor
	RandomNumberGenerator(const int proc_id, const int num_procs, const uint64_t seed = 0);

##### operator()
	double operator()();
	
Generates the next random number in the interval [0,1).

##### Non-copyable, but movable
Peferably, the random number generators are non-copyable. However, to store them in a STL-container like `std::vector` they have to be move-constructible.  
A `std::vector` containing `num_rngs` random number generators can be most efficiently created by using `emplace_back`:

	std::vector<RandomNumberGenerator> rngs;
	for (int i = 0; i < num_rngs; ++i)
    	rngs.emplace_back(proc_id, num_procs, seed);



### Random number library of C++11
Defined in header `std_random_wrapper.hpp`. See also <http://en.cppreference.com/w/cpp/numeric/random>.  
The wrapper class for the random number library of C++11 is templated on the type of the random number engine. For example, a random number generator that uses the *64-bit Mersenne Twister* algorithm is constructed by

	StdRandomWrapper<std::mt19937_64> rng(proc_id, num_procs, seed);


### SPRNG library
Defined in header `sprng_wrapper.hpp`. See also <http://www.sprng.org>.  
The wrapper class for the SPRNG library is templated on the type of the generator (= engine in C++11 language). Out of the six generators that SPRNG provides, we only support the *Modified Lagged Fibonacci Generator* (LFG) and the *Multiplicative Lagged Fibonacci Generator* (MLFG). For instance, a random number generator that uses the *Modified Lagged Fibonacci* algorithm is created by

	SprngWrapper<LFG> rng(proc_id, num_procs, seed);


### Ranq2
Defined in header `ranq2.hpp`.  
By using bitshift and xor operations Ranq2 is a very fast random number algorithm. A random number generator of type *Ranq2* is constructed by

	Ranq2 rng(proc_id, num_procs, seed);