#ifndef __DCA_UTILS_TYPE_UTILS_HPP_
#define __DCA_UTILS_TYPE_UTILS_HPP_

#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#   include <cxxabi.h>
#endif
#include <memory>
#include <string>

#include <dca/util/type_list.hpp>

// forward declare these templates so we can use them in print functions
// @TODO, move the domain print type functions to the domain classes
template<typename parameters>
class dmn_0;
template<typename ... domain_list>
class dmn_variadic;

namespace dca {
    namespace util {

        //----------------------------------------------------------------------------
        /// an assertion is thrown if the two types do not match
        /// extends std::is_same<> by forcing the compiler to print the types in
        /// the error message which helps with debugging
        //----------------------------------------------------------------------------
        template<typename T1, typename T2>
        struct assert_same
        {
            assert_same() {
                static_assert(std::is_same<T1, T2>::value, "Types must be equal");
            }
            static_assert(std::is_same<T1, T2>::value, "Types must be equal");
        };

        //----------------------------------------------------------------------------
        /// print a type cleanly if possible
        /// StackOverflow http://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
        //----------------------------------------------------------------------------
        template<class T>
        std::string type_name()
        {
            typedef typename std::remove_reference<T>::type TR;
            std::unique_ptr<char, void (*)(void*)> own
            (
                #ifndef _MSC_VER
                abi::__cxa_demangle(typeid(TR).name(), nullptr,
                    nullptr, nullptr),
                #else
                nullptr,
#endif
                std::free
                );
            std::string r = own != nullptr ? own.get() : typeid(TR).name();
            if (std::is_const<TR>::value)
                r += " const";
            if (std::is_volatile<TR>::value)
                r += " volatile";
            if (std::is_lvalue_reference<T>::value)
                r += "&";
            else if (std::is_rvalue_reference<T>::value)
                r += "&&";
            return r;
        }

        //----------------------------------------------------------------------------
        // print_type
        //----------------------------------------------------------------------------
        template<class T>
        struct print_type_name {
            static void print() {
                print(std::cout);
            }

            static void print(std::ostream& ss) {
                ss << "\t" << T::get_name() << std::endl;
            }

            template<class stream_type>
            static void to_JSON(stream_type& ss) {
                ss << "\t" << T::get_name() << std::endl;
            }
        };

        //----------------------------------------------------------------------------
        //
        //----------------------------------------------------------------------------
        // basic print displays the type of the template instantiation
        template<typename D>
        struct print_type {
            static void print() {
                print(std::cout);
            }
            //
            static void print(std::ostream &stream) {
                stream << "\t" << type_name<D>().c_str() << "\n";
            }
            static void to_JSON(std::ostream &stream) {
                stream << "\t" << type_name<D>().c_str();
            }
        };

        // dmn_0 override is actually the same as basic, but provided
        // for future customization
        template<typename Domain>
        struct print_type<dmn_0<Domain>> {
            static void print()
            {
                print(std::cout);
            }
            static void print(std::ostream &stream)
                {
                stream << "\t" << type_name<dmn_0<Domain>>().c_str() << "\n";
            }
        };

        // dmn_variadic prints out all subdomains recursively via pack expansion
        template<typename ... Domains>
        struct print_type<dmn_variadic<Domains...>> {
            static void print() {
                print(std::cout);
            }
            // we can't expand a pack out without passing it as a parameter
            // so expand the pack as a parameter list, and drop dummy return values
            // use func(),0 because func() returns void
            static void print(std::ostream &s) {
                ignore_returnvalues((print_type<Domains>::print(s),0)...);
            }
        };

        // dmn_variadic prints out all subdomains recursively via pack expansion
        template<typename Domain, typename ... Domains>
        struct print_type<dca::util::Typelist<Domain, Domains...>> {
            static void print() {
                print(std::cout);
            }
            // we can't expand a pack out without passing it as a parameter
            // so expand the pack as a parameter list, and drop dummy return values
            // use func(),0 because func() returns void
            static void print(std::ostream &s) {
                ignore_returnvalues((print_type<Domains>::print(s),0)...);
            }

            static void to_JSON(std::ostream &s) {
                s << "\"";
                print_type<Domain>::to_JSON(s);
                if (sizeof...(Domains)==0) {
                    s << "\"\n";
                }
                else {
                    s << "\",\n";
                    print_type<dca::util::Typelist<Domains...>>::to_JSON(s);
                }
            }
        };

    } // namespace util
} //namespace dca
#endif
