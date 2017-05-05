#include "options.hpp"

namespace help {
        std::unique_ptr<Options> Options::__current_options(nullptr);

        Options::~Options() {
                // Prevent a double free
                __current_options.release();
        }

        const help::Options * Options::get_options() {
                return __current_options.get();
        }

        void Options::set_options(help::Options* opts) {
                __current_options.reset(opts);
        }

}
