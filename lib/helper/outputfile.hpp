#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

namespace inout {

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