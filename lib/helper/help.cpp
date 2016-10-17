#include "help.hpp"
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>

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

	const IdatmEntry &get_info_map(const string &name) { 
		try { 
			return help::infoMap.at(name); 
		} catch(const std::out_of_range&) { 
			throw Error("cannot find idatm type " + name + " in IdatmInfoMap"); 
		}
	}
	

	std::tuple<double, double, double> gnuplot(const double &x1, const double &x2, const string &datapoints) {

		double coeffA = 0, coeffB = 0, WSSR = HUGE_VAL;
		
		// try a range of coefficients to get the best fit
		for (double a = 1e+5; a < 1e+10; a *= 2) { 
			
			// fit 1/x^12 repulsion term onto datapoints formatted x y (one-per-line)
			const string cmd = "gnuplot << EOF 2>&1\n"
				"f(x) = a/x**12 + b\n"
				"a=" + std::to_string(a) + "\n"
				"b=-1\n"
				"fit [" + std::to_string(x1) + ":" + std::to_string(x2) + "] f(x) \"-\" u 1:2 via a,b\n"
				+ datapoints + "\n"
				"e\n"
				"print 'JANEZ_FITTED ', a, ' ' , b, ' ', FIT_WSSR\n"
				"EOF\n";
				
		    FILE* pipe = popen(cmd.c_str(), "r");
		    if (!pipe) throw Error("die : install gnuplot!");
		    char buffer[128];
		    vector<string> output;
		    while(!feof(pipe)) {
				if(fgets(buffer, 128, pipe) != NULL)
					output.push_back(buffer);
			}
			pclose(pipe);

			// convert output to potential
			for (auto &line : output) {
				stringstream ss(line);
				string str1, str2, str3, str4;
				ss >> str1 >> str2 >> str3 >> str4;
				if (str1 == "JANEZ_FITTED") {
					dbgmsg("str1 = " << str1 << " str2 = " << str2 << " str3 = " << str3 << " str4 = " << str4);
					try { 
						double w = stod(str4); // this throws exceptions
						
						if (WSSR > w) {
							coeffA = stod(str2);
							coeffB = stod(str3);
							WSSR = w;
						}
						
					} catch (...) {}
				}
			}

		}
	    return std::make_tuple(coeffA, coeffB, WSSR);
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
