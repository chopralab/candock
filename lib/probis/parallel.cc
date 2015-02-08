#include "parallel.h"
#include <map>
#include <utility>
#include <sstream>
#include <algorithm>
#include "geo.h"


BinaryStream& operator<< (BinaryStream& bs, const JobOutput& s) {
	bs << s.suffix << s.mol2_pdb_id << s.mol2_chain_id << s.mol1_pdb_id << s.mol1_chain_id;
	return bs;
}

BinaryStream& operator>> (BinaryStream& bs, JobOutput& s) {
	bs >> s.suffix >> s.mol2_pdb_id >> s.mol2_chain_id >> s.mol1_pdb_id >> s.mol1_chain_id;
	return bs;
}

JobOutput ParallelJob::operator() (const JobInput& input) { 
	std::stringstream bufferStr(input);
	bufferStr >> PROTEIN2 >> CHAIN2;
	std::cerr << "comparing to " << PROTEIN2 << " " << CHAIN2 << "\n";
	JobOutput strings;
	state0(args, strings);
	return strings;
}

// init Control object form file (read proteins for comparison)
void Controller::init() {
	sentSurfs = 0;

	ifstream surf((add_end_slash(INDIR) + string(SURF_FILE)).c_str());
	if (!surf.is_open()) {
		cout << "Error (STATES) : Cannot open SURF_FILE " << SURF_FILE << "!" << endl;
		exit (1);
	}

	vector<pair<int, string> > surfComplVector;
	surfComplVector.reserve(30000);
	do {  
		string buffer;
		getline(surf, buffer);
		if (buffer.size() > 0) {
			surfComplVector.push_back(make_pair(1000000000, buffer));
		}
	}
	while (!surf.eof());

	// optionally (if file exists) read complexities of surfaces
	ifstream surfCompl((add_end_slash(INDIR) + string(SURF_FILE)+".complexity").c_str());
//	std::cout << "requested file " << (string(SURF_FILE)+".complexity").c_str() << "\n";
//	if (surfCompl) std::cout << "file found\n";
	map<string, int> surfComplMap;
	while (surfCompl) {
		string name;
		int complexity;
		surfCompl >> name >> complexity;
		if (surfCompl) {
			surfComplMap[name] = complexity;
		}
	}

	// assign complexities to surfaces that will be used in comparison	
	for (size_t i = 0; i < surfComplVector.size(); ++i) {
		istringstream stream(surfComplVector[i].second);
		string name;
		stream >> name;
		map<string, int>::const_iterator it = surfComplMap.find(name);
		if (it != surfComplMap.end())
			surfComplVector[i].first = it->second;
	}

	// sort by complexities, then write the result to surfVector)
	surfVector.resize(surfComplVector.size());
	sort(surfComplVector.begin(), surfComplVector.end());
	for (size_t i = 0; i < surfComplVector.size(); ++i)
		surfVector[surfVector.size()-i-1] = surfComplVector[i].second;

// //  debug sorting
//	for (size_t i = 0; i < surfComplVector.size(); ++i)
//		cout << surfComplVector[i].second << " " << surfComplVector[i].first << "\n";
//	exit(0);

	timerVals.reserve(surfVector.size());
}

// are there more jobs? If yes - return how many (return may be temporarily 0), if no - return -1
int Controller::moreJobs(int desiredNumber) {
	int moreToSend = surfVector.size() - sentSurfs;
	if (moreToSend == 0)
		moreToSend = -1;
	return std::min(moreToSend, desiredNumber);
}

// generate a new job (input for the job) with id (id may be safely ignored)
JobInput Controller::generateJob(int id) {
	//std::cerr << "generating job " << id << ": " << sentSurfs << "/" <<  surfVector.size() << " :" << surfVector[sentSurfs] << "\n";
	return surfVector[sentSurfs++];
}

// gather the result from finished job
void Controller::newResult(JobOutput& output, int id, int rank, double calcTime, double lifeTime) {
	timerVals.push_back(calcTime);
	if (output.suffix != "") {
		NoSql nosql(output.mol1_pdb_id, output.mol1_chain_id);
		nosql.write(output);
	}
}

//// finalize parallel processing
//void Controller::finalize() {
//	std::ofstream timerOut("timer.log");
//	for (size_t i = 0; i < timerVals.size(); ++i)
//		timerOut << timerVals[i] << "\n";
//}

