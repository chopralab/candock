#ifndef INOUT_H
#define INOUT_H
#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <sys/file.h> /* for flock(2) */
#include <sys/stat.h> /* for S_* constants */
#include <string.h> /* for strerror(3) prototype */
#include <stdio.h> /* for fprintf(3),printf(3),stderr protype */
#include <errno.h> /* for errno prototype */
#include <unistd.h> /* for close(2) prototypes */
using namespace std;

namespace inout {
	class Inout {
		static int __lock(const string &name);
		static void __unlock(int fd);
	public:
		enum f_not_found {panic=0, no_panic=1};
		static streampos read_file(const string &name, vector<string>&, 
			f_not_found=panic, const streampos &pos=0, const int num_occur=-1, 
			const string &pattern=""); // throws Error
		static streampos read_file(const string &name, string&, 
			f_not_found=panic, const streampos &pos=0, const int num_occur=-1, 
			const string &pattern=""); // throws Error
		static void file_open_put_contents(const string &name, 
			const vector<string> &v, ios_base::openmode=ios_base::out);
		static void file_open_put_stream(const string &name, 
			const stringstream &ss, ios_base::openmode=ios_base::out);
		static vector<string> files_matching_pattern(const string&, const string&);
	};	
	template<class T>
	void output_file(T &anything, const string &filename, ios_base::openmode mode=ios_base::out) {
		stringstream ss;
		ss << anything;
		Inout::file_open_put_stream(filename, ss, mode);
	}
	template<class T>
	void output_file(const T &anything, const string &filename, ios_base::openmode mode=ios_base::out) {
		stringstream ss;
		ss << anything;
		Inout::file_open_put_stream(filename, ss, mode);
	}
};
#endif
