#ifndef OPTIONS_H
#define OPTIONS_H

#include <memory>
#include <vector>

namespace help {
        class Options {
        private:
                static std::unique_ptr<Options> __current_options;
        public:

                static const Options* get_options();
                static void  set_options(Options* opts);
                
                virtual const std::string& get_string_option (const std::string& option) const = 0;
                virtual bool        get_bool_option   (const std::string& option) const = 0;
                virtual int         get_int_option    (const std::string& option) const = 0;
                virtual double      get_double_option (const std::string& option) const = 0;

                virtual const std::vector<std::string>& get_string_vector (const std::string& option) const = 0;
                
                virtual int ncpu() const = 0;
        };
}

#define cmdl (*help::Options::get_options())

#endif