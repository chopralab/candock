#include "candock/helper/inout.hpp"

#include <regex>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <sys/stat.h> /* for S_* constants */
#include <string.h> /* for strerror(3) prototype */
#include <stdio.h> /* for fprintf(3),printf(3),stderr protype */
#include <errno.h> /* for errno prototype */

#ifndef _WINDOWS

#include <sys/file.h> /* for flock(2) */
#include <unistd.h> /* for close(2) prototypes */

#else

#include <Windows.h>

#endif

#include "candock/helper/error.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/logger.hpp"

namespace candock{
namespace Inout {

        int Logger::__application_setting = Severity::CD_ERROR;
        bool Logger::__all_to_stderr = false;

        void __mkdir(const string &dir_path) {
                // makes a path "janez/aska/mia from e.g. "janez/aska/mia/test.txt"
                boost::filesystem::path dir(dir_path);
                dir.remove_filename();
                if(!dir.string().empty()) {
                        if(!boost::filesystem::exists(dir)) {
                                if(!boost::filesystem::create_directories(dir)) {
                                        throw Error("die : Cannot create directory " + dir.string());
                                }
                        }
                }
        }

#ifndef _WINDOWS

        int  __lock(const string &name) {
                int fd = open(name.c_str(),O_RDWR|O_CREAT,S_IRUSR|S_IWUSR);
                if (flock(fd,LOCK_EX) == -1) {
                        dbgmsg("Could not lock file " + name + " because " + strerror(errno));
                } 
                return fd;
        }

        void __unlock(int fd) {
                if (flock(fd,LOCK_UN)==-1) {
                        dbgmsg("Could not unlock handle " + std::to_string(fd));
                }
                if (close(fd)==-1) {
                        dbgmsg("Could not close handle " + std::to_string(fd));
                }
        }
#endif

        size_t file_size(const string &name) {
                if ( boost::filesystem::exists(name) &&
                        boost::filesystem::is_regular_file(name) )
                {
                        return boost::filesystem::file_size(name);
                } else {
                        return 0;
                }
        }

        void read_file(const string &name, string &s, f_not_found w, 
                const int num_occur, const string &pattern) {
                streampos pos_in_file = 0;
                read_file(name, s, pos_in_file, w, num_occur, pattern);
        }
        void read_file(const string &name, vector<string> &s,
                f_not_found w, const int num_occur, const string &pattern) {
                streampos pos_in_file = 0;
                read_file(name, s, pos_in_file, w, num_occur, pattern);
        }

        void read_file(const string &name, string &s, streampos &pos_in_file, 
                f_not_found w, const int num_occur, const string &pattern) {
                vector<string> vs;
                read_file(name, vs, pos_in_file, w, num_occur, pattern);
                s = boost::algorithm::join(vs, "\n");
        }

        void read_file(const string &name, vector<string> &s, streampos &pos_in_file, 
                f_not_found w, const int num_occur, const string &pattern) {
                dbgmsg(pos_in_file);
#ifndef _WINDOWS
                int fd = __lock(name);
                ifstream in(name, ios::in);
#else
                ifstream in(name, ios::in, _SH_DENYRW);
#endif
                if (!in.is_open() && w == panic) {
#ifndef _WINDOWS
                        __unlock(fd);
#endif
                        throw Error("Cannot read " + name);
                }
                read_stream(in, s, pos_in_file, num_occur, pattern);
                in.close();
#ifndef _WINDOWS
                __unlock(fd);
#endif
        }

        void read_stream(std::istream &in, vector<string> &s, streampos &pos_in_file, 
                const int num_occur, const string &pattern) {
                in.seekg(pos_in_file);
                string line;
                int i = 0;
                streampos pos = in.tellg();
                while (getline(in, line)) {
                        if (num_occur != -1 && (line.find(pattern) != string::npos && i++ == num_occur)) {
                                dbgmsg("breaking on line = " << line);
                                break;
                        } else
                                pos = in.tellg();
                        s.push_back(line);
                }
                in.seekg(pos);
                pos_in_file = pos;   
        }

        void read_stream(std::istream &in, vector<string> &s, const int num_occur, const string &pattern) {
                string line;
                int i = 0;
                streampos pos;
                while (getline(in, line)) {
                        if (num_occur != -1 && (line.find(pattern) != string::npos && i++ == num_occur)) {
                                dbgmsg("breaking on line = " << line);
                                break;
                        } else
                                pos = in.tellg();
                        s.push_back(line);
                }
                in.seekg(pos);
        }

        void file_open_put_contents(const string &name, 
                const vector<string> &v, ios_base::openmode mode) {
                stringstream ss;
                for (auto &s : v) ss << s << endl;
                file_open_put_stream(name, ss, mode);
        }

        void file_open_put_stream(const string &name, 
                const stringstream &ss, ios_base::openmode mode) {
                __mkdir(name);
#ifndef _WINDOWS
                int fd = __lock(name);
                ofstream output_file(name, mode);
#else
                ofstream output_file(name, mode, _SH_DENYRW);
#endif
                if (!output_file.is_open()) {
#ifndef _WINDOWS
                        __unlock(fd);
#endif
                        throw Error("Cannot open output file: " + name);  // how to emulate $! ?
                }
                output_file << ss.str();
                output_file.close();
#ifndef _WINDOWS
                __unlock(fd);
#endif
        }

        vector<string> files_matching_pattern(const string &file_or_dir, 
                const string &pattern) {

                vector<string> fn;
                namespace fs = boost::filesystem;
                fs::path p(file_or_dir);
                if (p.extension().string() == ".lst") {
                        read_file(file_or_dir, fn);
                } 
                else if (fs::is_directory(p)) {
                        vector<fs::path> files;
                        copy(fs::directory_iterator(p), 
                                fs::directory_iterator(), back_inserter(files));  // get all files in directoy...
                        for (auto &p : files) { 
                                fn.push_back(p.string()); 
                        }
                } else {
                        throw Error("die : Cannot open directory \"" + p.string() + '"');
                }
                vector<string> filenames;
                for(string &s : fn) {
                        if (std::regex_search(s, std::regex(pattern))) {
                                filenames.push_back(s);
                        }
                }
                if (filenames.empty())
                        throw Error("die : should provide either a .lst file with each line having one pdb filename, or a directory in which pdb files are...");
                return filenames;
	}
};
}
