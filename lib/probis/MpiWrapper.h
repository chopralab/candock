/* MIT License
*
* Copyright (c) 2017 Janez Konc
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#ifndef MPIWRAPPER_H_INCLUDED
#define MPIWRAPPER_H_INCLUDED

/// mpi header does not like being included after iostream, so special care should be taken using
/// this header. Either include mpi.h earlier, or include this header before iostream


#include <mpi.h>
#include <exception>
#include <vector>
#include "Array.h"
#include <cassert>


class BinaryStream {
public:
	std::vector<char> stream;
	mutable size_t outPos;
	
public:
	BinaryStream() : outPos(0) {}
	
	
	void* voidStream() {
		return &(stream[0]);
	}
	
	void* voidStream() const {
		return const_cast<void*>((void*)&(stream[0]));
	}
	
	size_t size() const {
		return stream.size();
	}
	
	void resize(size_t s) {
		stream.resize(s);
	}
	
	void clear() {
		stream.clear();
		outPos = 0;
	}
};

	
template<class Pod>
BinaryStream& operator<< (BinaryStream& bs, const Pod& pod) {
	size_t position = bs.stream.size();
	bs.stream.resize(bs.stream.size() + sizeof(Pod));
	*((Pod*)&bs.stream[position]) = pod;
	return bs;
}


template<class Pod>
BinaryStream& operator<< (BinaryStream& bs, const std::vector<Pod>& vec) {
	bs << vec.size();
	size_t position = bs.stream.size();
	bs.stream.resize(bs.stream.size() + sizeof(Pod) * vec.size());
	std::copy(vec.begin(), vec.end(), (Pod*)&bs.stream[position]);
	return bs;
}


template<class Pod>
BinaryStream& operator<< (BinaryStream& bs, const std::basic_string<Pod>& s) {
	bs << s.size();
	size_t position = bs.stream.size();
	bs.stream.resize(bs.stream.size() + sizeof(Pod) * s.size());
	std::copy(s.data(), s.data()+s.size(), (Pod*)&bs.stream[position]);
	return bs;
}


template<class Pod>
BinaryStream& operator>> (BinaryStream& bs, Pod& pod) {
	assert(bs.outPos+sizeof(Pod) <= bs.stream.size());
	pod = *reinterpret_cast<const Pod*>(&bs.stream[bs.outPos]);
	bs.outPos += sizeof(Pod);
	return bs;
}


template<class Pod>
BinaryStream& operator>> (BinaryStream& bs, std::vector<Pod>& vec) {
	size_t s;
	bs >> s;
	vec.resize(s);
	assert(bs.outPos+sizeof(Pod)*vec.size() <= bs.stream.size());
	std::copy((Pod*)&bs.stream[bs.outPos], ((Pod*)&bs.stream[bs.outPos]) 
		+ vec.size(), vec.begin());
	bs.outPos += sizeof(Pod)*vec.size();
	return bs;
}


template<class Pod>
BinaryStream& operator>> (BinaryStream& bs, std::basic_string<Pod>& st) {
	size_t s;
	bs >> s;
	st.resize(s);
	assert(bs.outPos+sizeof(Pod)*st.size() <= bs.stream.size());
	std::copy((Pod*)&bs.stream[bs.outPos], ((Pod*)&bs.stream[bs.outPos]) + st.size(), st.begin());
	bs.outPos += sizeof(Pod)*st.size();
	return bs;
}


namespace Mpi {
	
	class Communicator {
		MPI_Comm comm;
		
	public:
		Communicator(MPI_Comm c = MPI_COMM_WORLD) : comm(c) {}
		
		const Communicator& operator= (MPI_Comm c) {
			comm = c;
			return *this;
		}
		
		int getSize() const {
			int size;
			MPI_Comm_size(comm, &size);
			return size;
		}
		
		int getRank() const {
			int rank;
			MPI_Comm_rank(comm, &rank);
			return rank;
		}
		
		operator MPI_Comm() const {
			return comm;
		}
	};
	
	
	struct Status {
		MPI_Status status;
		bool messageWaiting;
		
		int source() const {return status.MPI_SOURCE;}
		int tag() const {return status.MPI_TAG;}
		int error() const {return status.MPI_ERROR;}
		int count() const {
			int c;
			MPI_Get_count(const_cast<MPI_Status*>(&status), MPI_CHAR, &c);
			return c;
		}
		
		Status() : messageWaiting(false) {
		}
		
		bool probeSpecific(int sourcep, int tagp, const Communicator& comm = MPI_COMM_WORLD) {
			int flag;
			MPI_Iprobe(sourcep, tagp, comm, &flag, &status);
			messageWaiting = (flag != 0);
			return messageWaiting;
		}
		
		bool probe(const Communicator& comm = MPI_COMM_WORLD) {
			return probeSpecific(MPI_ANY_SOURCE, MPI_ANY_TAG, comm);
		}
		
		/// blocking probe (wait for message)
		void waitSpecific(int sourcep, int tagp, const Communicator& comm = MPI_COMM_WORLD) {
			MPI_Probe(sourcep, tagp, comm, &status);
		}
		
		void wait(const Communicator& comm = MPI_COMM_WORLD) {
			waitSpecific(MPI_ANY_SOURCE, MPI_ANY_TAG, comm);
		}
	};
	
	
	struct Request {
		MPI_Request request;
		
		void wait() {
			Status status;
			MPI_Wait(&request, &status.status);
		}
	};
	
	
	template<class T>
	struct Streamify {
		T* data;
		
		Streamify(const T& t) : data(const_cast<T*>(&t)) {}
		void* ptr() {return (void*)(data);}
		size_t size() {return sizeof(T);}
		void getSize(int source, int tag, const Communicator& comm) {}
		void setSize(size_t ssize) {assert(ssize == sizeof(T));}
	};
	
	template<class T>
	struct Streamify<std::vector<T> > {
		std::vector<T>& vec;
		
		Streamify(const std::vector<T>& t) : vec(const_cast<std::vector<T>&>(t)) {}
		void* ptr() {return (void*)(&vec[0]);}
		size_t size() const {return sizeof(T)*vec.size();}
		void getSize(int source, int tag, const Communicator& comm) {
			Status status;
			status.probeSpecific(source, tag, comm);
			vec.resize(status.count() / sizeof(T));
		}
		void setSize(size_t ssize) {vec.resize(ssize / sizeof(T));}
	};

	template<class T>
	struct Streamify<std::basic_string<T> > {
		std::basic_string<T>& s;
						
		Streamify(const std::basic_string<T>& t) : s(const_cast<std::basic_string<T>&>(t)) {}
		void* ptr() {return (void*)(s.data());}
		size_t size() const {return sizeof(T)*s.size();}
		void getSize(int source, int tag, const Communicator& comm) {
			Status status;
			status.probeSpecific(source, tag, comm);
			s.resize(status.count() / sizeof(T));
		}
			
		void setSize(size_t ssize) {s.resize(ssize / sizeof(T));}
	};
	
	template<>
	struct Streamify<BinaryStream> {
		BinaryStream& data;
		
		Streamify(const BinaryStream& t) : data(const_cast<BinaryStream&>(t)) {}
		void* ptr() {return data.voidStream();}
		size_t size() const {return data.size();}
		void getSize(int source, int tag, const Communicator& comm) {
			Status status;
			status.probeSpecific(source, tag, comm);
			data.resize(status.count());
		}
		void setSize(size_t ssize) {data.resize(ssize);}
	};
	
	template<class T, size_t N>
	struct Streamify<Array<T, N> > {
		Array<T, N>& data;
		
		Streamify(Array<T, N>& t) : data(t) {}
		void* ptr() {return (void*)(&data[0]);}
		size_t size() const {return N;}
		void getSize(int source, int tag, const Communicator& comm) {
			Status status;
			status.probeSpecific(source, tag, comm);
			assert(status.count() / sizeof(T) == N);
		}
		void setSize(size_t ssize) {assert(ssize == N);}
	};


	template<class T>
	void send(const T& t, int dest, int tag, const Communicator& comm = MPI_COMM_WORLD) {
		Streamify<T> s(t);
		MPI_Send(s.ptr(), s.size(), MPI_CHAR, dest, tag, comm);
	}
	
	template<class T>
	void send(const T& t, int dest, int tag, Request& req, const Communicator& comm = MPI_COMM_WORLD) {
		Streamify<T> s(t);
		MPI_Ibsend(s.ptr(), s.size(), MPI_CHAR, dest, tag, comm, &req.request);
	}
	
	
	template<class T>
	void receive(T& t, int source, int tag, const Communicator& comm = MPI_COMM_WORLD) {
		MPI_Status status;
		Streamify<T> s(t);
		s.getSize(source, tag, comm);
		MPI_Recv(s.ptr(), s.size(), MPI_CHAR, source, tag, comm, &status);
	}
	
	template<class T>
	void receive(T& t, Status& status, const Communicator& comm = MPI_COMM_WORLD) {
		Streamify<T> s(t);
		s.setSize(status.count());
		MPI_Recv(s.ptr(), s.size(), MPI_CHAR, status.source(), status.tag(), comm, &status.status);
	}
	
	
	class Environment {
		Communicator comm;
		
	public:
		Environment(int &argc, char**& argv) {
			if (MPI_Init(&argc, &argv) != MPI_SUCCESS) throw std::exception();
		} 
		
		~Environment() {
			MPI_Finalize();
		}
		
		const Communicator& getCommunicator() const {
			return comm;
		}
		
		operator const Communicator& () const {
			return comm;
		}
	};
	
	
	// class that wraps process buffer allocation and deallocation
	class Buffer {
		char* buf;
		
	public:
		Buffer(size_t size) {
			size += MPI_BSEND_OVERHEAD;
			buf = new char[size];
			MPI_Buffer_attach((void*)buf, size);
		}
		
		~Buffer() {
			void* dummyAddr;
			int dummySize;
			MPI_Buffer_detach(&dummyAddr, &dummySize);
			delete[] buf;
		}
	};
	
};


#endif // MPIWRAPPER_H_INCLUDED
