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

                static Severity __application_setting;

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
                                F (__out);
                                return *this;
                        }
                };

        public:

                static void flip_mode ( Severity s1 ) {
                        __application_setting = static_cast<Severity>(__application_setting^s1);
                }

                template<Severity log_type>
                static LoggerImpl<log_type> log() {
                        return log_type >= WARNING? LoggerImpl<log_type>(std::cerr) : LoggerImpl<log_type>(std::cout);
                }

        };
}

#define log_note Inout::Logger::log<NOTE>()
#define log_step Inout::Logger::log<STEP>()
#define log_benchmark Inout::Logger::log<BENCHMARK>()
#define log_warning Inout::Logger::log<WARNING>()
#define log_error Inout::Logger::log<ERROR>()

#endif // HELP_LOGGER_H
