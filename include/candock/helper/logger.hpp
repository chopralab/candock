#ifndef HELP_LOGGER_H
#define HELP_LOGGER_H

#include <memory>
#include <iostream>

#include "candock/candockexport.hpp"

namespace Inout {

        enum Severity {
                CD_NONE      = 0,
                CD_NOTE      = 1 << 1,
                CD_STEP      = 1 << 2,
                CD_BENCHMARK = 1 << 3,
                CD_WARNING   = 1 << 4,
                CD_ERROR     = 1 << 5
        };

        class CANDOCK_EXPORT Logger {

                static int __application_setting;
                static bool __all_to_stderr;

                template <Severity s>
                class LoggerImpl {
                        std::ostream& __log_out;
                public:

                        explicit LoggerImpl(std::ostream& out) : __log_out(out) {}

                        template<typename T>
                        const LoggerImpl &operator<< (const T &v) const {
                                if (__application_setting & s) {
                                        __log_out << v;
                                }

                                return *this;
                        }

                        LoggerImpl const &operator<< (std::ostream& (*F) (std::ostream &)) const {
                                if (__application_setting & s) {
                                        F(__log_out);
                                }
                                return *this;
                        }
                };

        public:

                static void flip_mode ( int s1 ) {
                        __application_setting = __application_setting^s1;
                }

                static void set_all_stderr(bool mode) {
                        __all_to_stderr = mode;
                }

                static bool all_stderr() {
                        return __all_to_stderr;
                }

                template<Severity log_type>
                static LoggerImpl<log_type> log() {
                        return log_type >= CD_WARNING || all_stderr()?
                            LoggerImpl<log_type>(std::cerr) : LoggerImpl<log_type>(std::cout);
                }

        };
}

#define log_note      Inout::Logger::log<Inout::Severity::CD_NOTE>()
#define log_step      Inout::Logger::log<Inout::Severity::CD_STEP>()
#define log_benchmark Inout::Logger::log<Inout::Severity::CD_BENCHMARK>()
#define log_warning   Inout::Logger::log<Inout::Severity::CD_WARNING>()
#define log_error     Inout::Logger::log<Inout::Severity::CD_ERROR>()

#endif // HELP_LOGGER_H
