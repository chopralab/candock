#ifndef CONFORMATIONS_H
#define CONFORMATIONS_H

namespace Docker {

	typedef vector<pair<Molib::Atom*, Docker::Gpoint*>> OneConformation;
	typedef map<int, map<int, map<int, vector<int>>>> ConfMap;
	
	class Conformations {
		vector<OneConformation> __conf;
		ConfMap __confmap;
		void __init_conformations(const Molib::Molecule &molecule, Gpoints &gpoints, 
			const double &grid_spacing);
	public:
		Conformations(const Molib::Molecule &molecule, Gpoints &gpoints, const double &grid_spacing);
		vector<OneConformation> &get_conformations() { return __conf; }
		ConfMap &get_confmap() { return __confmap; }
		friend ostream& operator<<(ostream& os, const Conformations &conformations);
	};
};

ostream& operator<<(ostream& os, const Docker::OneConformation &conf);

#endif
