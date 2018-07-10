#ifndef OPTIONS_H
#define OPTIONS_H

#include <memory>
#include <vector>

#include "candock/candockexport.hpp"

namespace candock {

namespace help {
        class CANDOCK_EXPORT Options {
        private:
                static std::unique_ptr<Options> __current_options;
        public:

                virtual ~Options();
                static const Options* get_options();
                static void  set_options(Options* opts);

                virtual const std::string& get_string_option (const std::string& option) const = 0;
                virtual bool        get_bool_option   (const std::string& option) const = 0;
                virtual int         get_int_option    (const std::string& option) const = 0;
                virtual double      get_double_option (const std::string& option) const = 0;

                virtual const std::vector<std::string>& get_string_vector (const std::string& option) const = 0;

                virtual std::string program_name() const = 0;
                virtual int ncpu() const = 0;

                virtual std::string configuration_file() const = 0;
        };
}

#define cmdl (*candock::help::Options::get_options())

}

#endif