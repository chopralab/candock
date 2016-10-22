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

	const map<const string, const map<string,string>> standard_residues {
		{"A",{{"P","Pac"},{"OP1","O3-"},{"OP2","O3-"},{"O5'","O3"},{"C5'","C3"},{"C4'","C3"},{"O4'","O3"},{"C3'","C3"},{"O3'","O3"},{"C2'","C3"},{"O2'","O3"},{"C1'","C3"},{"N1","Npl"},{"N3","Npl"},{"N6","Npl"},{"N7","Npl"},{"N9","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"C8","Car"}}},
		{"C",{{"P","Pac"},{"OP1","O3-"},{"OP2","O3-"},{"O5'","O3"},{"C5'","C3"},{"C4'","C3"},{"O4'","O3"},{"C3'","C3"},{"O3'","O3"},{"C2'","C3"},{"O2'","O3"},{"C1'","C3"},{"N1","Npl"},{"N3","Npl"},{"N4","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"O2","O2"}}},
		{"G",{{"P","Pac"},{"OP1","O3-"},{"OP2","O3-"},{"O5'","O3"},{"C5'","C3"},{"C4'","C3"},{"O4'","O3"},{"C3'","C3"},{"O3'","O3"},{"C2'","C3"},{"O2'","O3"},{"C1'","C3"},{"N1","Npl"},{"N2","Npl"},{"N3","Npl"},{"N7","Npl"},{"N9","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"C8","Car"},{"O6","O2"}}},
		{"U",{{"P","Pac"},{"OP1","O3-"},{"OP2","O3-"},{"O5'","O3"},{"C5'","C3"},{"C4'","C3"},{"O4'","O3"},{"C3'","C3"},{"O3'","O3"},{"C2'","C3"},{"O2'","O3"},{"C1'","C3"},{"N1","Npl"},{"N3","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"O2","O2"},{"O4","O2"}}},
		{"ALA",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB1","HC"},{"HB2","HC"},{"HB3","HC"}}},
		{"ARG",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C3"},{"CG","C3"},{"CZ","C2"},{"NE","Ng+"},{"NH1","Ng+"},{"NH2","Ng+"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HD2","HC"},{"HD3","HC"},{"HE","H"},{"HH11","H"},{"HH12","H"},{"HH21","H"},{"HH22","H"}}},
		{"ASN",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG","C2"},{"ND2","Npl"},{"N","Npl"},{"OD1","O2"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD21","H"},{"HD22","H"}}},
		{"ASP",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG","Cac"},{"N","Npl"},{"OD1","O2-"},{"OD2","O2-"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"}}},
		{"CYS",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"O","O2"},{"SG","S3"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG","H"}}},
		{"CYX",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"O","O2"},{"SG","S3"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"}}},
		{"GLN",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C2"},{"CG","C3"},{"NE2","Npl"},{"N","Npl"},{"OE1","O2"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HE21","H"},{"HE22","H"}}},
		{"GLU",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","Cac"},{"CG","C3"},{"N","Npl"},{"OE1","O2-"},{"OE2","O2-"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"}}},
		{"GLY",{{"OXT","O2-"},{"CA","C3"},{"C","C2"}, {"N","Npl"},{"O","O2"},{"H","H"},{"HA2","HC"},{"HA3","HC"}}},
		{"HID",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","Car"},{"CE1","Car"},{"CG","Car"},{"ND1","Npl"},{"NE2","N2"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","H"},{"HE1","HC"},{"HD2","HC"}}},
		{"HIE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","Car"},{"CE1","Car"},{"CG","Car"},{"ND1","N2"},{"NE2","Npl"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HE1","HC"},{"HE2","H"},{"HD2","HC"}}},
		{"HIP",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","C2"},{"CE1","C2"},{"CG","C2"},{"ND1","Npl"},{"NE2","Npl"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","H"},{"HE1","HC"},{"HE2","H"},{"HD2","HC"}}},
		{"HIS",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","Car"},{"CE1","Car"},{"CG","Car"},{"ND1","Npl"},{"NE2","N2"},{"N","Npl"},{"O","O2"},{"H","H"},{"HN","H"},{"HA","HC"},{"HB1","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","H"},{"HD2","HC"},{"HE1","HC"},{"HE2","H"}}},
		{"ILE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","C3"},{"CG1","C3"},{"CG2","C3"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB","HC"},{"HG21","HC"},{"HG22","HC"},{"HG23","HC"},{"HG12","HC"},{"HG13","HC"},{"HD11","HC"},{"HD12","HC"},{"HD13","HC"}}},
		{"LEU",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","C3"},{"CD2","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG","HC"},{"HD11","HC"},{"HD12","HC"},{"HD13","HC"},{"HD21","HC"},{"HD22","HC"},{"HD23","HC"}}},
		{"LYS",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C3"},{"CE","C3"},{"CG","C3"},{"N","Npl"},{"NZ","N3+"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HD2","HC"},{"HD3","HC"},{"HE2","HC"},{"HE3","HC"},{"HZ1","H"},{"HZ2","H"},{"HZ3","H"}}},
		{"MET",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CE","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"SD","S3"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HE1","HC"},{"HE2","HC"},{"HE3","HC"}}},
		{"MSE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CE","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"SE","Se"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HE1","HC"},{"HE2","HC"},{"HE3","HC"}}},
		{"PHE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","Car"},{"CD2","Car"},{"CE1","Car"},{"CE2","Car"},{"CG","Car"},{"CZ","Car"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","HC"},{"HE1","HC"},{"HZ","HC"},{"HE2","HC"},{"HD2","HC"}}},
		{"PRO",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"HD2","HC"},{"HD3","HC"},{"HG2","HC"},{"HG3","HC"},{"HB2","HC"},{"HB3","HC"},{"HA","HC"}}},
		{"PTR",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","Car"},{"CD2","Car"},{"CE1","Car"},{"CE2","Car"},{"CG","Car"},{"CZ","Car"},{"N","Npl"},{"OH","O3"},{"O","O2"},{"O1P","O3-"},{"O2P","O3-"},{"O3P","O3-"},{"P","Pac"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","HC"},{"HE1","HC"},{"HH","H"},{"HE2","HC"},{"HD2","HC"}}},
		{"SER",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"OG","O3"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG","H"}}},
		{"SEP",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"OG","O3"},{"O","O2"},{"O1P","O3-"},{"O2P","O3-"},{"O3P","O3-"},{"P","Pac"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG","H"}}},
		{"THR",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG2","C3"},{"N","Npl"},{"OG1","O3"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB","HC"},{"HG21","HC"},{"HG22","HC"},{"HG23","HC"},{"HG1","H"}}},
		{"TPO",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG2","C3"},{"N","Npl"},{"OG1","O3"},{"O","O2"},{"O1P","O3-"},{"O2P","O3-"},{"O3P","O3-"},{"P","Pac"},{"H","H"},{"HA","HC"},{"HB","HC"},{"HG21","HC"},{"HG22","HC"},{"HG23","HC"},{"HG1","H"}}},
		{"TRP",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","Car"},{"CD2","Car"},{"CE2","Car"},{"CE3","Car"},{"CG","Car"},{"CH2","Car"},{"CZ2","Car"},{"CZ3","Car"},{"NE1","Npl"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","HC"},{"HE1","H"},{"HZ2","HC"},{"HH2","HC"},{"HZ3","HC"},{"HE3","HC"}}},
		{"TYR",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","Car"},{"CD2","Car"},{"CE1","Car"},{"CE2","Car"},{"CG","Car"},{"CZ","Car"},{"N","Npl"},{"OH","O3"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","HC"},{"HE1","HC"},{"HH","H"},{"HE2","HC"},{"HD2","HC"}}},
		{"VAL",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG1","C3"},{"CG2","C3"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB","HC"},{"HG11","HC"},{"HG12","HC"},{"HG13","HC"},{"HG21","HC"},{"HG22","HC"},{"HG23","HC"}}},

		{"HOH",{{"O","O3"},{"H1","H"},{"H2","H"}}},
	};

};
