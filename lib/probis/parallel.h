#ifndef _PARALLEL_H
#define _PARALLEL_H

#include "states.h"
#include "ParallelFramework.h"
#include <vector>
#include <string>
#include <sstream>

typedef std::string JobInput;
typedef NoSql::WriteStruct JobOutput;

BinaryStream& operator<< (BinaryStream& bs, const JobOutput& s);
BinaryStream& operator>> (BinaryStream& bs, JobOutput& s);

struct ParallelJob {
	Args* args;
	ParallelJob(Args* a) : args(a) {}
	JobOutput operator() (const JobInput& input);
};

struct Controller {
	std::vector<std::string> surfVector;
	size_t sentSurfs;
	std::vector<double> timerVals;
	
	void init();
	int moreJobs(int desiredNumber);
	JobInput generateJob(int id);
	void newResult(JobOutput& output, int id, int rank, double calcTime, double lifeTime);
//	void finalize();
};

#endif // _PARALLEL_H
