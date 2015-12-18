#include "help.hpp"
#include <string>
#include <iostream>
#include <stdio.h>

namespace help {

	ostream& operator<<(ostream& os, const smiles& edges) {
		for (auto &e : edges)
			os << "{" << e.atom_property1 
				<< "," << e.atom_property2
				<< "," << e.bond_property << "}";
		return os;
	}

	ostream& operator<<(ostream& os, const rename_rule& rule) {
		os << "PATTERN = " << rule.pattern << " RULE = {";
		for (auto &txt : rule.rule) os << txt << ",";
		os << "}";
		//~ os << endl;
		return os;
	}

	vector<string> gnuplot(const string &x1, const string &x2, const string &datapoints) {

		// fit 1/x^12 repulsion term onto datapoints formatted x y (one-per-line)
		const string cmd = "gnuplot << EOF 2>&1\n"
			"f(x) = a/x**12 + b\n"
			"a=1000000\n"
			"b=1\n"
			"fit [" + x1 + ":" + x2 + "] f(x) \"-\" u 1:2 via a,b\n"
			+ datapoints + "\n"
			"e\n"
			"print 'a = ', a, ' b = ', b\n"
			"EOF\n";
			
	    FILE* pipe = popen(cmd.c_str(), "r");
	    if (!pipe) throw Error("die : install gnuplot!");
	    char buffer[128];
	    vector<string> result;
	    while(!feof(pipe)) {
	    	if(fgets(buffer, 128, pipe) != NULL)
	    		result.push_back(buffer);
	    }
	    pclose(pipe);
	    return result;
		
	}

	string memusage(const string &msg) {
		const string cmd = "ps ax -o rss,command | sort -nr | head -n 10|grep test_link|cut -f1 -d' '";
	    FILE* pipe = popen(cmd.c_str(), "r");
	    if (!pipe) return "ERROR";
	    char buffer[128];
	    string result = "";
	    while(!feof(pipe)) {
	    	if(fgets(buffer, 128, pipe) != NULL)
	    		result += buffer;
	    }
	    pclose(pipe);
	    cerr << "Memusage:"  << msg << ":" << result << endl;
	    return result;
	}

	vector<vector<string>> get_replacement(const vector<string> &initial) {
		vector<vector<string>> result;
		if (initial.size() == 2) {
			const string &iniclass1 = initial[0];
			const string &iniclass2 = initial[1];
			for (auto &repclass1 : gaff_replacement.at(iniclass1)) {
				for (auto &repclass2 : gaff_replacement.at(iniclass2)) {
					if (!(repclass1 == iniclass1 && repclass2 == iniclass2)) {
						dbgmsg("replacement bond type [" << repclass1 << " " 
							<< repclass2 << "]");
						result.push_back({repclass1, repclass2});
					}
				}
			}
		} else if (initial.size() == 3) {
			const string &iniclass1 = initial[0];
			const string &iniclass2 = initial[1];
			const string &iniclass3 = initial[2];
			for (auto &repclass1 : gaff_replacement.at(iniclass1)) {
				for (auto &repclass2 : gaff_replacement.at(iniclass2)) {
					for (auto &repclass3 : gaff_replacement.at(iniclass3)) {
						if (!(repclass1 == iniclass1 && repclass2 == iniclass2 && repclass3 == iniclass3)) {
							dbgmsg("replacement angle type [" << repclass1 << " " 
								<< repclass2 << " " << repclass3 << "]");
							result.push_back({repclass1, repclass2, repclass3});
						}
					}
				}
			}
		} else if (initial.size() == 4) {
			const string &iniclass1 = initial[0];
			const string &iniclass2 = initial[1];
			const string &iniclass3 = initial[2];
			const string &iniclass4 = initial[3];
			for (auto &repclass1 : gaff_replacement.at(iniclass1)) {
				for (auto &repclass2 : gaff_replacement.at(iniclass2)) {
					for (auto &repclass3 : gaff_replacement.at(iniclass3)) {
						for (auto &repclass4 : gaff_replacement.at(iniclass4)) {
							if (!(repclass1 == iniclass1 && repclass2 == iniclass2 
								&& repclass3 == iniclass3 && repclass4 == iniclass4)) {
								dbgmsg("replacement dihedral type [" << repclass1 << " " 
									<< repclass2 << " " << repclass3 << " " << repclass4 << "]");
								result.push_back({repclass1, repclass2, repclass3, repclass4});
							}
						}
					}
				}
			}
		}	
		return result;
	}

};
