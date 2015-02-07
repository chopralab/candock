#include "optics.hpp"

int main(int argc, char* argv[]) {
	try {
		vector<string> db {"a", "b", "c", "d"};
		cluster::PairwiseDistances<string> pairwise_distances {
			{&db[0], {{&db[1], 10.6},{&db[2], 10.2},{&db[3], 10.3}}},
			{&db[1], {{&db[2], 0.1},{&db[3], 0.2}}},
			{&db[2], {{&db[3], 0.5}}}
		};
		//~ cluster::PairwiseDistances<string> pairwise_distances {
			//~ {&db[0], {{&db[1], 0.6},{&db[2], 0.2},{&db[3], 0.3}}},
			//~ {&db[1], {{&db[2], 0.1},{&db[3], 0.2},{&db[0], 0.6}}},
			//~ {&db[2], {{&db[3], 0.5},{&db[1], 0.1},{&db[0], 0.2}}},
			//~ {&db[3], {{&db[0], 0.3},{&db[1], 0.2},{&db[2], 0.5}}}
		//~ };
		cluster::MapD<string> scores {{&db[0], 1},{&db[1], 2},{&db[2], 3},{&db[3], 4}};
		//~ cluster::Optics<string, std::less<double>, std::greater<double>> optics(pairwise_distances, scores, 1.0, 2);
		cluster::Optics<string> optics(pairwise_distances, scores, 1.0, 2);
		//~ cluster::Optics<string> optics(pairwise_distances, scores, 1.0, 2);
		optics.extract_dbscan(0.6);
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}
