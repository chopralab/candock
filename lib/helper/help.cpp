#include "help.hpp"
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "renamerules.hpp"

namespace help {

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

        const map<const string, const map<string,string>> cofactor_residues {
                {"FAD", {{"PA" ,"Pac"},
                         {"O1A","O3-"},
                         {"O2A","O3-"},
                         {"O5B","O3"},
                         {"C5B","C3"},
                         {"C4B","C3"},
                         {"O4B","O3"},
                         {"C3B","C3"},
                         {"O3B","O3"},
                         {"C2B","C3"},
                         {"O2B","O3"},
                         {"C1B","C3"},
                         {"N9A","Npl"},
                         {"C8A","Car"},
                         {"N7A","N2"},
                         {"C5A","Car"},
                         {"C6A","Car"},
                         {"N6A","Npl"},
                         {"N1A","N2"},
                         {"C2A","Car"},
                         {"N3A","N2"},
                         {"C4A","Car"},
                         {"N1" ,"N2"},
                         {"C2" ,"Car"},
                         {"O2" ,"O2"},
                         {"N3" ,"Npl"},
                         {"C4" ,"Car"},
                         {"O4" ,"O2"},
                         {"C4X","Car"},
                         {"N5" ,"N2"},
                         {"C5X","Car"},
                         {"C6" ,"Car"},
                         {"C7" ,"Car"},
                         {"C7M","C3"},
                         {"C8" ,"Car"},
                         {"C8M","C3"},
                         {"C9" ,"Car"},
                         {"C9A","Car"},
                         {"N10","Npl"},
                         {"C10","Car"},
                         {"C1'","C3"},
                         {"C2'","C3"},
                         {"O2'","O3"},
                         {"C3'","C3"},
                         {"O3'","O3"},
                         {"C4'","C3"},
                         {"O4'","O3"},
                         {"C5'","C3"},
                         {"O5'","O3"},
                         {"P"  , "Pac"},
                         {"O1P","O3-"},
                         {"O2P","O3-"},
                         {"O3P","O3"}}
                },
                {"NAD", {{"PA" ,"Pac"},
                         {"O1A","O3-"},
                         {"O2A","O3-"},
                         {"O5B","O3"},
                         {"C5B","C3"},
                         {"C4B","C3"},
                         {"O4B","O3"},
                         {"C3B","C3"},
                         {"O3B","O3"},
                         {"C2B","C3"},
                         {"O2B","O3"},
                         {"C1B","C3"},
                         {"N9A","Npl"},
                         {"C8A","Car"},
                         {"N7A","N2"},
                         {"C5A","Car"},
                         {"C6A","Car"},
                         {"N6A","Npl"},
                         {"N1A","N2"},
                         {"C2A","Car"},
                         {"N3A","N2"},
                         {"C4A","Car"},
                         {"O3" ,"O3"},
                         {"PN" ,"Pac"},
                         {"O1N","O3-"},
                         {"O2N","O3-"},
                         {"O5D","O3"},
                         {"C5D","C3"},
                         {"C4D","C3"},
                         {"O4D","O3"},
                         {"C3D","C3"},
                         {"O3D","O3"},
                         {"C2D","C3"},
                         {"O2D","O3"},
                         {"C1D","C3"},
                         {"N1N","Ng+"},
                         {"C2N","Car"},
                         {"C3N","Car"},
                         {"C7N","C2"},
                         {"O7N","O2"},
                         {"N7N","Npl"},
                         {"C4N","Car"},
                         {"C5N","Car"},
                         {"C6N","Car"}}
                },
                {"HEM", {{"FE" , "Fe"},
                         {"CHA","C2"},
                         {"CHB","C2"},
                         {"CHC","C2"},
                         {"CHD","C2"},
                         {"NA" ,"Npl"},
                         {"C1A","Car"},
                         {"C2A","Car"},
                         {"C3A","Car"},
                         {"C4A","Car"},
                         {"CMA","C3"},
                         {"CAA","C3"},
                         {"CBA","C3"},
                         {"CGA","Cac"},
                         {"O1A","O2-"},
                         {"O2A","O2-"},
                         {"NB" ,"N2"},
                         {"C1B","C2"},
                         {"C2B","C2"},
                         {"C3B","C2"},
                         {"C4B","C2"},
                         {"CMB","C3"},
                         {"CAB","C2"},
                         {"CBB","C2"},
                         {"NC" ,"Npl"},
                         {"C1C","C2"},
                         {"C2C","C2"},
                         {"C3C","C2"},
                         {"C4C","C2"},
                         {"CMC","C3"},
                         {"CAC","C2"},
                         {"CBC","C2"},
                         {"ND" ,"N2"},
                         {"C1D","C2"},
                         {"C2D","C2"},
                         {"C3D","C2"},
                         {"C4D","C2"},
                         {"CMD","C3"},
                         {"CAD","C3"},
                         {"CBD","C3"},
                         {"CGD","Cac"},
                         {"O1D","O2-"},
                         {"O2D","O2-"}}
                }
        };
        
	const map<const string, const pair<string,string>> non_specific_binders {
		{"BE7", {"(4-CARBOXYPHENYL)(CHLORO)MERCURY", "C7 H5 O2 CL1 HG1"}},
		{"MRD", {"(4R)-2-METHYLPENTANE-2,4-DIOL", "C6 H14 O2"}},
		{"MHA", {"(CARBAMOYLMETHYL-CARBOXYMETHYL-AMINO)-ACETIC ACID", "C6 H10 N2 O5"}},
		{"BU3", {"(R,R)-2,3-BUTANEDIOL", "C4 H10 O2"}},
		{"EDO", {"1,2-ETHANEDIOL", "C2 H6 O2"}},
		{"PGO", {"1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"BU2", {"1,3-BUTANEDIOL", "C4 H10 O2"}},
		{"PDO", {"1,3-PROPANEDIOL", "C3 H8 O2"}},
		{"BU1", {"1,4-BUTANEDIOL", "C4 H10 O2"}},
		{"PG6", {"1-(2-METHOXY-ETHOXY)-2-{2-[2-(2-METHOXY-ETHOXY]-ETHOXY}-ETHANE", "C12 H26 O6"}},
		{"1BO", {"1-BUTANOL", "C4 H10 O1"}},
		{"PE7", {"1-DEOXY-1-THIO-HEPTAETHYLENE GLYCOL", "C14 H30 O7 S1"}},
		{"PG5", {"1-METHOXY-2-[2-(2-METHOXY-ETHOXY]-ETHANE", "C8 H18 O4"}},
		{"TFP", {"10-[3-(4-METHYL-PIPERAZIN-1-YL)-PROPYL]-2-TRIFLUOROMETHYL-10H-PHENOTHIAZINE", "C21 H24 N3 F3 S1"}},
		{"DHD", {"2,4-DIHYDROXY-5-OXO-HEXA-2,3-DIENOIC ACID", "C5 H4 O6"}},
		{"PEU", {"2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80-HEPTACOSAOXADOOCTACONTAN-82-OL", "C55 H112 O28"}},
		{"TRS", {"2-AMINO-2-HYDROXYMETHYL-PROPANE-1,3-DIOL", "C4 H12 N1 O3"}},
		{"TAU", {"2-AMINOETHANESULFONIC ACID", "C2 H7 N1 O3 S1"}},
		{"SBT", {"2-BUTANOL", "C4 H10 O1"}},
		{"SAL", {"2-HYDROXYBENZOIC ACID", "C7 H6 O3"}},
		{"MPD", {"2-METHYL-2,4-PENTANEDIOL", "C6 H14 O2"}},
		{"IOH", {"2-PROPANOL, ISOPROPANOL", "C3 H8 O1"}},
		{"IPA", {"2-PROPANOL, ISOPROPANOL", "C3 H8 O1"}},
		{"PGE", {"2-[2-(2-HYDROXY-ETHOXY)-ETHOXY]-ETHANOL", "C6 H14 O4"}},
		{"PIG", {"2-[2-(2-HYDROXY-ETHOXY)-ETHOXY]-ETHANOL", "C6 H14 O4"}},
		{"B3P", {"2-[3-(2-HYDROXY-1,1-DIHYDROXYMETHYL-ETHYLAMINO)-PROPYLAMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL", "C11 H26 N2 O6"}},
		{"BTB", {"2-[BIS-(2-HYDROXY-ETHYL)-AMINO]-2-HYDROXYMETHYL- PROPANE-1,3-DIOL", "C8 H19 N1 O5"}},
		{"NHE", {"2-[N-CYCLOHEXYLAMINO]ETHANE SULFONIC ACID", "C8 H17 N1 O3 S1"}},
		{"C8E", {"2-{2-[2-(2-OCTYLOXYETHOXY)-ETHOXYL]-ETHOXY}ETHANOL", "C16 H34 O5"}},
		{"OTE", {"2-{2-[2-(2-OCTYLOXYETHOXY)-ETHOXYL]-ETHOXY}ETHANOL", "C16 H34 O5"}},
		{"PE4", {"2-{2-[2-(2-{2-[2-(2-ETHOXY-ETHOXY)-ETHOXY]-ETHOXY}-ETHOXY)-ETHOXY]-ETHOXY}-ETHANOL", "C16 H34 O8"}},
		{"XPE", {"3,6,9,12,15,18,21,24,27-NONAOXANONACOSANE-1,29-DIOL", "C20 H42 O11"}},
		{"PE8", {"3,6,9,12,15,18,21-HEPTAOXATRICOSANE-1,23-DIOL", "C16 H34 O9"}},
		{"P33", {"3,6,9,12,15,18-HEXAOXAICOSANE-1,20-DIOL", "C14 H30 O8"}},
		{"N8E", {"3,6,9,12,15-PENTAOXATRICOSAN-1-OL", "C18 H38 O6"}},
		{"2OS", {"3-N-OCTANOYLSUCROSE", "C20 H36 O12"}},
		{"1PS", {"3-PYRIDINIUM-1-YLPROPANE-1-SULFONATE", "C8 H11 N1 O3 S1"}},
		{"CPS", {"3-[(3-CHOLAMIDOPROPYL)DIMETHYLAMMONIO]-1-PROPANESULFONATE", "C32 H58 N2 O7 S1"}},
		{"DMX", {"3-[BENZYL(DIMETHYL)AMMONIO]PROPANE-1-SULFONATE", "C12 H19 N1 O3 S1"}},
		{"MPO", {"3[N-MORPHOLINE]PROPANE SULFONIC ACID", "C7 H15 N1 O4 S1"}},
		{"GCD", {"4,5-DEHYDRO-D-GLUCURONIC ACID", "C6 H8 O7"}},
		{"IDT", {"4,5-DEHYDRO-L-IDURONIC ACID", "C6 H8 O7"}},
		{"DXG", {"4-DEOXYGLUCARATE", "C6 H8 O7"}},
		{"CM5", {"5-CYCLOHEXYL-1-PENTYL-BETA-D-MALTOSIDE", "C23 H42 O11"}},
		{"ACA", {"6-AMINOHEXANOIC ACID", "C6 H13 N1 O2"}},
		{"ACT", {"ACETATE ION", "C2 H3 O2"}},
		{"ACN", {"ACETONE", "C3 H6 O1"}},
		{"CCN", {"ACETONITRILE", "C2 H3 N1"}},
		{"AGC", {"ALPHA-D-GLUCOSE", "C6 H12 O6"}},
		{"GLC", {"ALPHA-D-GLUCOSE", "C6 H12 O6"}},
		{"MAN", {"ALPHA-D-MANNOSE", "C6 H12 O6"}},
		{"DR6", {"ALPHA-[4-(1,1,3,3 - TETRAMETHYLBUTYL)PHENYL]- OMEGA-HYDROXY-POLY(OXY-1,2-ETHANEDIYL", "C74 H142 O31"}},
		{"AL", {"ALUMINUM ION", "AL1"}},
		{"NH4", {"AMMONIUM ION", "H4 N1"}},
		{"AZI", {"AZIDE ION", "N3"}},
		{"BNG", {"B-NONYLGLUCOSIDE", "C15H30O6"}},
		{"BOG", {"B-OCTYLGLUCOSIDE", "C14 H28 O6"}},
		{"BA", {"BARIUM ION", "BA1"}},
		{"LAK", {"BETA-D-GALACTOPYRANOSYL-1-6-BETA-D-GLUCOPYRANOSE", "C12 H22 O11"}},
		{"BGC", {"BETA-D-GLUCOSE", "C6 H12 O6"}},
		{"BMA", {"BETA-D-MANNOSE", "C6 H12 O6"}},
		{"BCN", {"BICINE", "C6 H13 N1 O4"}},
		{"BR", {"BROMIDE ION", "BR1"}},
		{"BR", {"BROMO GROUP", "BR1"}},
		{"BRO", {"BROMO GROUP", "BR1"}},
		{"CAC", {"CACODYLATE ION", "C2 H6 O2 AS1 1-"}},
		{"CD", {"CADMIUM ION", "CD1"}},
		{"CA", {"CALCIUM ION", "CA1"}},
		{"CBX", {"CARBOXY GROUP", "C1 H1 O2"}},
		{"FMT", {"CARBOXY GROUP", "C1 H1 O2"}},
		{"ACY", {"CARBOXYMETHYL GROUP", "C2 H3 O2"}},
		{"CBM", {"CARBOXYMETHYL GROUP", "C2 H3 O2"}},
		{"CM", {"CARBOXYMETHYL GROUP", "C2 H3 O2"}},
		{"CS", {"CESIUM ION", "CS1"}},
		{"CL", {"CHLORIDE ION", "CL1"}},
		{"CL", {"CHLORO GROUP", "CL1"}},
		{"CLO", {"CHLORO GROUP", "CL1"}},
		{"FCL", {"CITRATE ANION", "C6 H5 O7 3-"}},
		{"CIT", {"CITRIC ACID", "C6 H8 O7"}},
		{"CO", {"COBALT (II) ION", "CO1"}},
		{"3CO", {"COBALT (III) ION", "CO1"}},
		{"NCO", {"COBALT HEXAMMINE ION", "H18 N6 CO1"}},
		{"CU1", {"COPPER (I) ION", "CU1"}},
		{"CU", {"COPPER (II) ION", "CU1"}},
		{"CN", {"CYANIDE GROUP", "C1 N1"}},
		{"CYN", {"CYANIDE GROUP", "C1 N1"}},
		{"CYN", {"CYANIDE ION", "C1 N1"}},
		{"MA4", {"CYCLOHEXYL-HEXYL-BETA-D-MALTOSIDE", "C24 H44 O11"}},
		{"BTC", {"CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"CYS", {"CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"TAR", {"D(-)-TARTARIC ACID", "C4 H6 O6"}},
		{"GLO", {"D-GLUCOSE IN LINEAR FORM", "C6 H12 O6"}},
		{"MTL", {"D-MANNITOL", "C6 H14 O6"}},
		{"DPR", {"D-PROLINE", "C5 H9 N1 O2"}},
		{"SOR", {"D-SORBITOL", "C6 H14 O6"}},
		{"SYL", {"D-XYLITOL", "C7 H16 O5"}},
		{"DMU", {"DECYL-BETA-D-MALTOPYRANOSIDE", "C22 H42 O11"}},
		{"DDQ", {"DECYLAMINE-N,N-DIMETHYL-N-OXIDE", "C12 H27 N1 O1"}},
		{"DMS", {"DIMETHYL SULFOXIDE", "C2 H6 O1 S1"}},
		{"DMF", {"DIMETHYLFORMAMIDE", "C3 H7 N1 O1"}},
		{"DIO", {"DIOXANE", "C4 H8 O2"}},
		{"DOX", {"DIOXANE", "C4 H8 O2"}},
		{"12P", {"DODECAETHYLENE GLYCOL", "C24 H50 O13"}},
		{"SDS", {"DODECYL SULFATE", "C12 H26 O4 S1"}},
		{"LMT", {"DODECYL-BETA-D-MALTOSIDE", "C24 H46 O11"}},
		{"EOH", {"ETHANOL", "C2 H6 O1"}},
		{"EEE", {"ETHYL ACETATE", "C4 H8 O2"}},
		{"EDO", {"ETHYLENE GLYCOL", "C2 H6 O2"}},
		{"EGL", {"ETHYLENE GLYCOL", "C2 H6 O2"}},
		{"FE2", {"FE (II) ION", "FE1"}},
		{"FE", {"FE (III) ION", "FE1"}},
		{"F", {"FLUORIDE ION", "F1"}},
		{"F", {"FLUORO GROUP", "F1"}},
		{"FLO", {"FLUORO GROUP", "F1"}},
		{"TRT", {"FRAGMENT OF TRITON X-100", "C21 H36 O4"}},
		{"CYS", {"FREE CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"FCY", {"FREE CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"FRU", {"FRUCTOSE", " C6 H12 O6"}},
		{"GBL", {"GAMMA-BUTYROLACTONE", "C4 H6 O2"}},
		{"GLC", {"GLUCOSE", "C6 H12 O6"}},
		{"GOL", {"GLYCEROL", "C3 H8 O3"}},
		{"GLY", {"GLYCINE", "C2 H5 N1 O2"}},
		{"GPX", {"GUANOSINE 5'-DIPHOSPHATE 2':3'-CYCLIC MONOPHOSPHATE", "C10 H14 N5 O13 P3"}},
		{"HTO", {"HEPTANE-1,2,3-TRIOL", "C7 H16 O3"}},
		{"HTG", {"HEPTYL 1-THIOHEXOPYRANOSIDE", "C13 H26 O5 S1"}},
		{"B7G", {"HEPTYL-BETA-D-GLUCOPYRANOSIDE", "C13 H26 O6"}},
		{"C10", {"HEXAETHYLENE GLYCOL MONODECYL ETHER", "C22 H46 O7"}},
		{"16D", {"HEXANE-1,6-DIAMINE", "C6 H16 N2"}},
		{"HEZ", {"HEXANE-1,6-DIOL", "C6 H14 O2"}},
		{"IOD", {"IODIDE ION", "I1"}},
		{"IDO", {"IODO GROUP", "I1"}},
		{"IOD", {"IODO GROUP", "I1"}},
		{"ICI", {"ISOCITRIC ACID", "C6 H8 O7"}},
		{"ICT", {"ISOCITRIC ACID", "C6 H8 O7"}},
		{"IPA", {"ISOPROPYL ALCOHOL", "C3 H8 01"}},
		{"TLA", {"L(+)-TARTARIC ACID", "C4 H6 O6"}},
		{"LAT", {"LACTOSE", "C12 H22 O11"}},
		{"LBT", {"LACTOSE", "C12 H22 O11"}},
		{"LDA", {"LAURYL DIMETHYLAMINE-N-OXIDE", "C14 H31 N1 O1"}},
		{"PB", {"LEAD (II) ION", "PB1"}},
		{"LI", {"LITHIUM ION", "LI1"}},
		{"MG", {"MAGNESIUM ION", "MG1"}},
		{"MN", {"MANGANESE (II) ION", "MN1"}},
		{"MN3", {"MANGANESE (III) ION", "MN1"}},
		{"HG", {"MERCURY (II) ION", "HG1"}},
		{"MRY", {"MESO-ERYTRHITOL", "C4 H10 O4"}},
		{"MOH", {"METHANOL", "C1 H4 O1"}},
		{"BEQ", {"N-(CARBOXYMETHYL)-N,N-DIMETHYL-3-[(1-OXODODECYL)AMINO]-1-PROPANAMINIUM INNER SALT", "C19 H38 N2 O3"}},
		{"C15", {"N-DODECYL-N,N-DIMETHYL-3-AMMONIO-1-PROPANESULFONATE", "C17 H38 N1 O3 S1"}}, 
		{"MG8", {"N-OCTANOYL-N-METHYLGLUCAMINE", "C15 H31 N1 O6"}},
		{"POL", {"N-PROPANOL", "C3 H8 O1"}},
		{"NI", {"NICKEL (II) ION", "NI1"}},
		{"3NI", {"NICKEL (III) ION", "NI1"}},
		{"NO3", {"NITRATE ION", "N1 O3"}},
		{"JEF", {"O-(O-(2-AMINOPROPYL)-O'-(2-METHOXYETHYL)POLYPROPYLENEGLYCOL 500)", "C31 H65 N1 O10"}},
		{"P4C", {"O-ACETALDEHYDYL-HEXAETHYLENE GLYCOL", "C14 H28 O8"}},
		{"CE1", {"O-DODECANYL OCTAETHYLENE GLYCOL", "C28 H58 O9"}},
		{"DIA", {"OCTANE 1,8-DIAMINE", "C8 H20 N2"}},
		{"CXE", {"PENTAETHYLENE GLYCOL MONODECYL ETHER", "C20 H42 O6"}},
		{"IPH", {"PHENOL", "C6 H6 O1"}},
		{"PIN", {"PIPERAZINE-N,N'-BIS(2-ETHANESULFONIC ACID)", "C8 H18 N2 O6 S2"}},
		{"15P", {"POLYETHYLENE GLYCOL (N=34)", "C69 H140 O35"}},
		{"K", {"POTASSIUM ION", "K1"}},
		{"CRY", {"PROPANE-1,2,3-TRIOL", "C3 H8 O3"}},
		{"GOL", {"PROPANE-1,2,3-TRIOL", "C3 H8 O3"}},
		{"PGR", {"R-1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"RB", {"RUBIDIUM ION", "RB1"}},
		{"PGO", {"S-1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"PGQ", {"S-1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"AG", {"SILVER ION", "AG1 1+"}},
		{"NA", {"SODIUM ION", "NA1"}},
		{"SPD", {"SPERMIDINE", "C7 H19 N3"}},
		{"SPK", {"SPERMINE (FULLY PROTONATED FORM)", "C10 H30 N4 4+"}},
		{"SPM", {"SPERMINE", "C10 H26 N4"}},
		{"SR", {"STRONTIUM ION", "SR1"}},
		{"SUC", {"SUCROSE", "C12 H22 O11"}},
		{"SO4", {"SULFATE ANION", "O4 S1"}},
		{"SUL", {"SULFATE ANION", "O4 S1"}},
		{"SO4", {"SULFATE ION", "O4 S1"}},
		{"TBU", {"TERTIARY-BUTYL ALCOHOL", "C4 H10 O1"}},
		{"TMA", {"TETRAMETHYLAMMONIUM ION", "C4 H12 N1 1+"}},
		{"TEP", {"THEOPHYLLINE", "C7 H8 N4 O2"}},
		{"SCN", {"THIOCYANATE ION", "C1 N1 S1"}},
		{"TRE", {"TREHALOSE", "C12 H22 O11"}},
		{"PGE", {"TRIETHYLENE GLYCOL", "C6 H14 O4"}},
		{"ETF", {"TRIFLUOROETHANOL", "C2 H3 O1 F3"}},
		{"144", {"TRIS-HYDROXYMETHYL-METHYL-AMMONIUM", "C4 H12 N1 O3"}},
		{"UMQ", {"UNDECYL-MALTOSIDE", "C23 H44 O11"}},
		{"URE", {"UREA", "C1 H4 N2 O1"}},
		{"YT3", {"YTTRIUM (III) ION", "Y1"}},
		{"Y1", {"YTTRIUM ION", "Y1"}},
		{"ZN2", {"ZINC ION ON 3-FOLD CRYSTAL AXIS", "ZN1"}},
		{"ZN", {"ZINC ION ON 3-FOLD CRYSTAL AXIS", "ZN1"}},
		{"ZN", {"ZINC ION", "ZN1"}},
		{"PEG", {"DI(HYDROXYETHYL)ETHER", ""}},
		{"BEN", {"BENZAMIDINE", ""}},
		{"2PE", {"NONAETHYLENE GLYCOL", ""}},
		{"P6G", {"HEXAETHYLENE GLYCOL", ""}},
		{"1PE", {"PENTAETHYLENE GLYCOL", ""}},
		{"DOD", {"DEUTERATED WATER", ""}},
		{"ACE", {"ACETYL GROUP", ""}},
		{"BEZ", {"BENZOIC ACID", ""}},
		{"XYP", {"BETA-D-XYLOPYRANOSE", ""}},
		{"NAG", {"N-ACETYL-D-GLUCOSAMINE", ""}},
		{"NDG", {"2-(ACETYLAMINO)-2-DEOXY-A-D-GLUCOPYRANOSE", ""}},
		{"MSE", {"SELENOMETHIONINE", ""}},
		{"PO4", {"PHOSPHATE ION", ""}},
		{"VO4", {"VANADATE ION", ""}},
		{"EPE", {"4-(2-HYDROXYETHYL)-1-PIPERAZINE ETHANESULFONIC ACID", ""}},
		{"MES", {"2-(N-MORPHOLINO)-ETHANESULFONIC ACID", ""}},
		{"MXE", {"2-METHOXYETHANOL", ""}},
		{"XE", {"XENON", ""}},
	};
        
        string replace_str_char(string str, const string& replace, char ch) {
                size_t found = str.find_first_of(replace);
                while (found != string::npos) { // While our position in the sting is in range.
                        str[found] = ch; // Change the character at position.
                        found = str.find_first_of(replace, found+1); // Relocate again.
                }
                return str; // return our new string.
        }
        vector<string> ssplit(const string &source, const char*delimiter, bool keepEmpty) {
                vector<string> results;
                size_t prev = 0;
                size_t next = 0;
                while ((next = source.find_first_of(delimiter, prev)) != string::npos){
                        if (keepEmpty || (next - prev != 0)) {
                                results.push_back(source.substr(prev, next - prev));
                        }
                        prev = next + 1;
                }
                if (prev < source.size()) {
                        results.push_back(source.substr(prev));
                }
                return results;
        }

        char get_one_letter(string s) { 
                auto it = one_letter.find(s); 
                if (it == one_letter.end()) { 
                        dbgmsg("warn : residues name " << s << 
                                " not found in to_one_letter..."); 
                        return ' '; 
                } 
                return it->second; 
        }

};
