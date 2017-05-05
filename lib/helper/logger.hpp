#ifndef HELP_LOGGER_H
#define HELP_LOGGER_H

#include <memory>
#include <iostream>

namespace Inout {

        enum Severity {
                NONE      = 0,
                NOTE      = 1 << 1,
                STEP      = 1 << 2,
                BENCHMARK = 1 << 3,
                WARNING   = 1 << 4,
                ERROR     = 1 << 5,
        };

        class Logger {

                static int __application_setting;

                template <Severity s>
                class LoggerImpl {
                        std::ostream& __out;
                public:

                        explicit LoggerImpl(std::ostream& out) : __out(out) {}

                        template<typename T>
                        const LoggerImpl &operator<< (const T &v) const {
                                if (__application_setting & s) {
                                        __out << v;
                                }

                                return *this;
                        }

                        LoggerImpl const &operator<< (std::ostream& (*F) (std::ostream &)) const {
                                if (__application_setting & s) {
                                        F(__out);
                                }
                                return *this;
                        }
                };

        public:

                static void flip_mode ( int s1 ) {
                        __application_setting = __application_setting^s1;
                }

                template<Severity log_type>
                static LoggerImpl<log_type> log() {
                        return log_type >= WARNING? LoggerImpl<log_type>(std::cerr) : LoggerImpl<log_type>(std::cout);
                }

        };
}

#define log_note Inout::Logger::log<Inout::NOTE>()
#define log_step Inout::Logger::log<Inout::STEP>()
#define log_benchmark Inout::Logger::log<Inout::BENCHMARK>()
#define log_warning Inout::Logger::log<Inout::WARNING>()
#define log_error Inout::Logger::log<Inout::ERROR>()

#endif // HELP_LOGGER_H
