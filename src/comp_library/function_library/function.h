//-*-C++-*-

#ifndef FUNCTION_H
#define FUNCTION_H

#include "scalar_cast_methods.h"
#include "function_operations/set_to_zero.h"
#include "function_operations/copy_from.h"
#include "print_type.h"
#include "function_operations/function_operators_collection.h"
#include "function_operations/subind_2_linind_collection.h"


namespace FUNC_LIB
{

/*!
 *  \defgroup FUNCTION
 */

/*!
 *  \class   function
 *  \ingroup FUNCTION
 *
 *  \brief  This class connects the function-values to the domains
 *  \author Peter Staar
 *
 *  \version 1.0
 *  \date    2009-2014
 */
    template<typename scalartype, class domain>
    class function
    {
    public:

        typedef scalartype this_scalar_type;
        typedef domain     this_domain_type;

    public:

        function();
        function(std::string name);
        function(const function<scalartype, domain>& other_one);

        ~function();

        void reset();

/**
   @name domain-parameters
   @{
*/
        domain&      get_domain();

        std::string&      get_name();
        int          signature();
        int          size();
        int          operator[](int index);
/**@}*/

        scalartype*  values();
        scalartype*  values() const;

/**
   @name linear index versus coordinate index
   @{
*/
        void linind_2_subind(int i, int*& subind);
        void linind_2_subind(int linind, std::vector<int>& subind);
        void subind_2_linind(int* subind, int& i);

        int  subind_2_linind(int i);

#include "function_operations/subind_2_linind.h"
/**@}*/

/**
   @name opertor(...) implementations
   @{
*/
        scalartype& operator()(int i);
        scalartype& operator()(int* subind);

#include "function_operations/function_operators.h"
/**@}*/

        function<scalartype, domain>& operator =(function<scalartype, domain>& f_other);
        void                          operator+=(function<scalartype, domain>& f_other);
        void                          operator-=(function<scalartype, domain>& f_other);
        void                          operator*=(function<scalartype, domain>& f_other);
        void                          operator/=(function<scalartype, domain>& f_other);

/**
   @name atomic operations
   @{
*/
        template<typename new_scalartype>
        void operator=(new_scalartype c);

        template<typename new_scalartype>
        void operator+=(new_scalartype c);

        template<typename new_scalartype>
        void operator-=(new_scalartype c);

        template<typename new_scalartype>
        void operator*=(new_scalartype c);

        template<typename new_scalartype>
        void operator/=(new_scalartype c);
/**@}*/

/**
   @name data-management
   @{
*/
        template<typename new_scalartype>
        void slice     (int sbdm_index, int* subind, new_scalartype* fnc_vals);
        template<typename new_scalartype>
        void distribute(int sbdm_index, int* subind, new_scalartype* fnc_vals);

        template<typename new_scalartype>
        void slice     (int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals);
        template<typename new_scalartype>
        void distribute(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals);
/**@}*/

/**
   @name I/O-interface
   @{
*/
        void print_fingerprint();
        void print_2_file(const char* file_name);

        void to_JSON(std::ofstream& ss);
        void to_JSON(std::stringstream& ss);

        void from_JSON(std::fstream& input_file);

        template<typename concurrency_t>
        void from_JSON(std::fstream& input_file, concurrency_t& concurrency);
/**@}*/

/**
   @name MPI-interface
   @{
*/
        template<typename concurrency_t>
        int get_buffer_size(concurrency_t& concurrency);

        template<class concurrency_t>
        void pack(concurrency_t& concurrency,
                  int* buffer,
                  int buffer_size,
                  int& position);

        template<class concurrency_t>
        void unpack(concurrency_t& concurrency,
                    int* buffer,
                    int buffer_size,
                    int& position);
/**@}*/

    private:

        void read_value_from_JSON(std::string& line, int index);

    private:

        std::string            name;
        std::string            function_type;

        domain            dmn;
        int               Nb_elements;

        int               Nb_sbdms;
        std::vector<int>& size_sbdm;
        std::vector<int>  step_sbdm;

        scalartype*      fnc_values;
    };

//++++++++++++++++++++++++++++++++++++++++//
//+++  CONSTRUCTOR & DESTRUCTOR        +++//
//++++++++++++++++++++++++++++++++++++++++//

    template<typename scalartype, class domain>
    function<scalartype, domain>::function():
        name("no name"),
        function_type(__PRETTY_FUNCTION__),
        dmn(),
        Nb_elements(dmn.get_size()),
        Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
        size_sbdm(dmn.get_leaf_domain_sizes()),
        step_sbdm(Nb_sbdms,1)
    {
        //std::cout << __PRETTY_FUNCTION__ << std::endl;

        for(int i=0; i<Nb_sbdms; i++)
            for(int j=0; j<i; j++)
                step_sbdm[i] *= dmn.get_subdomain_size(j);

        fnc_values = new scalartype[Nb_elements];

        for(int linind=0; linind<Nb_elements; linind++)
            set_to_zero::execute(fnc_values[linind]);
    }

    template<typename scalartype, class domain>
    function<scalartype, domain>::function(std::string fnc_name):
        name(fnc_name),
        function_type(__PRETTY_FUNCTION__),
        dmn(),
        Nb_elements(dmn.get_size()),
        Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
        size_sbdm(dmn.get_leaf_domain_sizes()),
        step_sbdm(Nb_sbdms,1)
    {

        for(int i=0; i<Nb_sbdms; i++)
            for(int j=0; j<i; j++)
                step_sbdm[i] *= dmn.get_subdomain_size(j);

        fnc_values = new scalartype[Nb_elements];

        for(int linind=0; linind<Nb_elements; linind++)
            set_to_zero::execute(fnc_values[linind]);
    }


    template<typename scalartype, class domain>
    function<scalartype, domain>::function(const function<scalartype, domain>& other_one):
        name("no_name"),
        function_type(__PRETTY_FUNCTION__),
        dmn(),
        Nb_elements(dmn.get_size()),
        Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
        size_sbdm(dmn.get_leaf_domain_sizes()),
        step_sbdm(Nb_sbdms,1)
    {
        //std::cout << __PRETTY_FUNCTION__ << std::endl;

        for(int i=0; i<Nb_sbdms; i++)
            for(int j=0; j<i; j++)
                step_sbdm[i] *= dmn.get_subdomain_size(j);

        fnc_values = new scalartype[Nb_elements];

        copy_from<scalartype>::execute(Nb_elements, fnc_values, other_one.values());
    }


    template<typename scalartype, class domain>
    function<scalartype, domain>::~function()
    {
        delete [] fnc_values;
    }

    template<typename scalartype, class domain>
    void function<scalartype, domain>::reset()
    {
        dmn.reset();

        for(int i=0; i<Nb_sbdms; i++)
            size_sbdm[i] = dmn.get_subdomain_size(i);

        for(int i=0; i<Nb_sbdms; i++){
            step_sbdm[i] = 1;
            for(int j=0; j<i; j++)
                step_sbdm[i] *= dmn.get_subdomain_size(j);
        }

        Nb_sbdms    = dmn.get_leaf_domain_sizes().size();
        Nb_elements = dmn.get_size();

        delete [] fnc_values;
        fnc_values = new scalartype[Nb_elements];

        for(int linind=0; linind<Nb_elements; linind++)
            set_to_zero::execute(fnc_values[linind]);
    }

//++++++++++++++++++++++++++++++++++++++++//
//+++  GET                             +++//
//++++++++++++++++++++++++++++++++++++++++//

    template<typename scalartype, class domain>
    domain& function<scalartype, domain>::get_domain()
    {
        return dmn;
    }

    template<typename scalartype, class domain>
    std::string& function<scalartype, domain>::get_name()
    {
        return name;
    }

    template<typename scalartype, class domain>
    int function<scalartype, domain>::signature()
    {
        return Nb_sbdms;
    }


    template<typename scalartype, class domain>
    int function<scalartype, domain>::size()
    {
        return Nb_elements;
    }

    template<typename scalartype, class domain>
    int function<scalartype, domain>::operator[](int index)
    {
        return size_sbdm[index];
    }

    template<typename scalartype, class domain>
    scalartype*  function<scalartype, domain>::values()
    {
        return fnc_values;
    }


    template<typename scalartype, class domain>
    scalartype* function<scalartype, domain>::values() const
    {
        return fnc_values;
    }

//++++++++++++++++++++++++++++++++++++++++//
//+++  LINEAR INDEX <=> SUB INDEX      +++//
//++++++++++++++++++++++++++++++++++++++++//


    template<typename scalartype, class domain>
    void function<scalartype, domain>::linind_2_subind(int linind, int*& subind)
    {
        int tmp = linind;
        for(int i=0; i<int(size_sbdm.size()); i++)
        {
            subind[i] = tmp % size_sbdm[i];
            tmp = (tmp - subind[i])/size_sbdm[i];
        }
    }

    template<typename scalartype, class domain>
    void function<scalartype, domain>::linind_2_subind(int linind, std::vector<int>& subind)
    {
        assert(int(subind.size()) == signature());

        int tmp = linind;
        for(int i=0; i<int(size_sbdm.size()); i++)
        {
            subind[i] = tmp % size_sbdm[i];
            tmp = (tmp - subind[i])/size_sbdm[i];
        }
    }

    template<typename scalartype, class domain>
    void function<scalartype, domain>::subind_2_linind(int* subind, int& linind)
    {
        linind = 0;
        for(int i=0; i<int(step_sbdm.size()); i++)
            linind += subind[i]*step_sbdm[i];
    }

    template<typename scalartype, class domain>
    int function<scalartype, domain>::subind_2_linind(int i)
    {
        assert(i<Nb_elements);
        return i;
    }

//++++++++++++++++++++++++++++++++++++++++//
//+++  FUNCTION HANDLING               +++//
//++++++++++++++++++++++++++++++++++++++++//

    template<typename scalartype, class domain>
    scalartype& function<scalartype, domain>::operator()(int* subind)
    {
        int linind;
        subind_2_linind(subind, linind);

        assert(linind >=0 && linind<Nb_elements);
        return fnc_values[linind];
    }

    template<typename scalartype, class domain>
    scalartype& function<scalartype, domain>::operator()(int i)
    {
        assert(i >=0 && i<Nb_elements);
        return fnc_values[i];
    }

    template<typename scalartype, class domain>
    function<scalartype, domain>& function<scalartype, domain>::operator=(function<scalartype, domain>& f_other)
    {
        domain& dmn_other = f_other.get_domain();

        if(dmn.get_size() != dmn_other.get_size()) // domains were not initialized when function was created !
        {
            dmn.get_size()                = dmn_other.get_size();
            dmn.get_branch_domain_sizes() = dmn_other.get_branch_domain_sizes();
            dmn.get_leaf_domain_sizes()   = dmn_other.get_leaf_domain_sizes();

            for(int i=0; i<Nb_sbdms; i++)
                size_sbdm[i] = dmn.get_subdomain_size(i);

            for(int i=0; i<Nb_sbdms; i++)
                for(int j=0; j<i; j++)
                    step_sbdm[i] *= dmn.get_subdomain_size(j);

            Nb_sbdms    = dmn.get_leaf_domain_sizes().size();
            Nb_elements = dmn.get_size();

            delete [] fnc_values;
            fnc_values = new scalartype[Nb_elements];
        }

        memcpy(fnc_values, f_other.values(), Nb_elements*sizeof(scalartype));

        return *this;
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::operator=(new_scalartype c)
    {
        scalartype c_new(c);

        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] = c_new;
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::operator+=(new_scalartype c)
    {
        scalartype c_new(c);

        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] += c_new;
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::operator-=(new_scalartype c)
    {
        scalartype c_new(c);

        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] -= c_new;
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::operator*=(new_scalartype c)
    {
        scalartype c_new(c);

        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] *= c_new;
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::operator/=(new_scalartype c)
    {
        scalartype c_new(c);

        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] /= c_new;
    }

    template<typename scalartype, class domain>
    void function<scalartype, domain>::operator+=(function<scalartype, domain>& f_other)
    {
        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] += f_other(linind);
    }

    template<typename scalartype, class domain>
    void function<scalartype, domain>::operator-=(function<scalartype, domain>& f_other)
    {
        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] -= f_other(linind);
    }

    template<typename scalartype, class domain>
    void function<scalartype, domain>::operator*=(function<scalartype, domain>& f_other)
    {
        for(int linind=0; linind<Nb_elements; linind++)
            fnc_values[linind] *= f_other(linind);
    }

    template<typename scalartype, class domain>
    void function<scalartype, domain>::operator/=(function<scalartype, domain>& f_other)
    {
        for(int linind=0; linind<Nb_elements; linind++){
            assert(ASSERT_NON_ZERO(f_other(linind)));
            fnc_values[linind] /= f_other(linind);
        }
    }


//++++++++++++++++++++++++++++++++++++++++//
//+++  SLICING && DISTRIBUTING         +++//
//++++++++++++++++++++++++++++++++++++++++//


    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::slice(int sbdm_index, int* subind, new_scalartype* fnc_vals)
    {
        assert(sbdm_index >= 0);
        assert(sbdm_index < Nb_sbdms);

        int linind=0;
        subind[sbdm_index] = 0;
        subind_2_linind(subind, linind);

        for(int i=0; i<size_sbdm[sbdm_index]; i++)
            fnc_vals[i] = do_cast<new_scalartype>::execute(fnc_values[linind + i*step_sbdm[sbdm_index] ]);
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::distribute(int sbdm_index, int* subind, new_scalartype* fnc_vals)
    {
        assert(sbdm_index >= 0);
        assert(sbdm_index < Nb_sbdms);

        int linind=0;
        subind[sbdm_index] = 0;
        subind_2_linind(subind, linind);

        for(int i=0; i<size_sbdm[sbdm_index]; i++)
            fnc_values[linind + i*step_sbdm[sbdm_index] ] = do_cast<scalartype>::execute(fnc_vals[i]);
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::slice(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals)
    {
        assert(sbdm_index_1 >= 0);
        assert(sbdm_index_2 >= 0);
        assert(sbdm_index_1 < Nb_sbdms);
        assert(sbdm_index_2 < Nb_sbdms);

        int linind=0;
        subind[sbdm_index_1] = 0;
        subind[sbdm_index_2] = 0;
        subind_2_linind(subind, linind);

        int size_sbdm_1 = size_sbdm[sbdm_index_1];
        int size_sbdm_2 = size_sbdm[sbdm_index_2];

        int step_sbdm_1 = step_sbdm[sbdm_index_1];
        int step_sbdm_2 = step_sbdm[sbdm_index_2];

        new_scalartype* fnc_ptr_left  = NULL;
        new_scalartype* fnc_ptr_right = NULL;

        for(int j=0; j<size_sbdm_2; j++){

            fnc_ptr_left  = &fnc_vals  [0      + j*size_sbdm_1];
            fnc_ptr_right = &fnc_values[linind + j*step_sbdm_2];

            for(int i=0; i<size_sbdm_1; i++)
                fnc_ptr_left[i] = fnc_ptr_right[i*step_sbdm_1];
//       fnc_vals[i+j*size_sbdm[sbdm_index_1]] = fnc_values[linind + i*step_sbdm[sbdm_index_1] + j*step_sbdm[sbdm_index_2]];
        }
    }

    template<typename scalartype, class domain>
    template<typename new_scalartype>
    void function<scalartype, domain>::distribute(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals)
    {
        assert(sbdm_index_1 >= 0);
        assert(sbdm_index_2 >= 0);
        assert(sbdm_index_1 < Nb_sbdms);
        assert(sbdm_index_2 < Nb_sbdms);

        int linind=0;
        subind[sbdm_index_1] = 0;
        subind[sbdm_index_2] = 0;
        subind_2_linind(subind, linind);

        for(int i=0; i<size_sbdm[sbdm_index_1]; i++)
            for(int j=0; j<size_sbdm[sbdm_index_2]; j++)
                fnc_values[linind + i*step_sbdm[sbdm_index_1] + j*step_sbdm[sbdm_index_2]] = fnc_vals[i+j*size_sbdm[sbdm_index_1]];
    }


//++++++++++++++++++++++++++++++++++++++++//
//+++  PRINT                           +++//
//++++++++++++++++++++++++++++++++++++++++//



    template<typename scalartype, class domain>
    void function<scalartype, domain>::print_fingerprint()
    {
        std::cout << std::endl << std::endl << "function : " << name <<std::endl;

        std::cout <<"*********************************"<<std::endl;

        std::cout << "# subdomains        : " << Nb_sbdms << std::endl;

        printTL<domain>::print();

        std::cout << "size of subdomains  : " << std::endl;
        for(int i=0; i<Nb_sbdms; i++)
            std::cout << size_sbdm[i] << "\t";
        std::cout << std::endl;

        std::cout << "memory step         : " << std::endl;
        for(int i=0; i<Nb_sbdms; i++)
            std::cout << step_sbdm[i] << "\t";
        std::cout << std::endl;

        std::cout << "# elements          : " << Nb_elements << std::endl;
        std::cout << "# size              : " << Nb_elements*sizeof(scalartype)*(1.e-6) << " (mega-bytes)" << std::endl;
        std::cout <<"*********************************"<<std::endl;
    }

    template<typename scalartype, class domain>
    template<typename concurrency_t>
    int function<scalartype, domain>::get_buffer_size(concurrency_t& concurrency)
    {
        int result = 0;
        result += concurrency.get_buffer_size(*this);
        return result;
    }

    template<typename scalartype, class domain>
    template<class concurrency_t>
    void function<scalartype, domain>::pack(concurrency_t& concurrency,
                                            int* buffer,
                                            int buffer_size,
                                            int& position)
    {
        concurrency.pack(buffer, buffer_size, position, *this);
    }

    template<typename scalartype, class domain>
    template<class concurrency_t>
    void function<scalartype, domain>::unpack(concurrency_t& concurrency,
                                              int* buffer,
                                              int buffer_size,
                                              int& position)
    {
        concurrency.unpack(buffer, buffer_size, position, *this);
    }


} // namespace FUNC_LIB

#endif
