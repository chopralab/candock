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
