#ifndef __DCA_UTILS_TYPE_UTILS_HPP_
#define __DCA_UTILS_TYPE_UTILS_HPP_

#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#   include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>

namespace dca { namespace util {

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

} // namespace util
} //namespace dca
#endif