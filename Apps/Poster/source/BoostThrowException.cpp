#include <exception>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

namespace boost {

[[noreturn]] void throw_exception(const std::exception& e) {
    std::fprintf(stderr, "boost::throw_exception called: %s (type=%s)\n",
                 e.what(), typeid(e).name());
    std::abort();
    }
}
