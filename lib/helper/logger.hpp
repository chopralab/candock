/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef HELP_LOGGER_H
#define HELP_LOGGER_H

#include <memory>
#include <iostream>

namespace Inout {

        enum Severity {
                CD_NONE      = 0,
                CD_NOTE      = 1 << 1,
                CD_STEP      = 1 << 2,
                CD_BENCHMARK = 1 << 3,
                CD_WARNING   = 1 << 4,
                CD_ERROR     = 1 << 5
        };

        class Logger {

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
