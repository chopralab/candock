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

#include "args.h"
#include "geo.h"
#include "motif.h"


void Args::read_constants(string par_file) {
  char buffer[XLARGE], ctrl_word[SMALL], const_str[MEDIUM];
  char const_chr;
  float const_float;
  ifstream par((add_end_slash(INDIR) + par_file).c_str());

  if (!par.is_open()) { 
    throw Err("Error (ARGS) : Cannot open the parameter file " + par_file, 14); 
  }
  cout << "Reading parameter file..." << endl;

  while (!par.eof()) {
    par.getline(buffer, XLARGE);
    strcpy(ctrl_word, "");
    strcpy(const_str, "");
    const_chr = ' ';
    sscanf(buffer, "%s %f", ctrl_word, &const_float);
    sscanf(buffer, "%s %s", ctrl_word, const_str);
    sscanf(buffer, "%s %c", ctrl_word, &const_chr);
    cout << "PARAM> " << buffer << endl;
    if      (strcmp(ctrl_word, "SRF_FILE") == 0)
      SRF_FILE = const_str;
    else if (strcmp(ctrl_word, "PROTEIN1") == 0)
      PROTEIN1 = const_str;
    else if (strcmp(ctrl_word, "PROTEIN2") == 0)
      PROTEIN2 = const_str;
    else if (strcmp(ctrl_word, "CHAIN1") == 0)
      CHAIN1 = const_str;
    else if (strcmp(ctrl_word, "CHAIN2") == 0)
      CHAIN2 = const_str;
    else if (strcmp(ctrl_word, "INTER_CHAIN_DIST") == 0)
      INTER_CHAIN_DIST = const_float;
    else if (strcmp(ctrl_word, "SURF_FILE") == 0)
      //~ strcpy(SURF_FILE, const_str);
      SURF_FILE = const_str;
    else if (strcmp(ctrl_word, "INDIR") == 0)
      INDIR = const_str;
    else if (strcmp(ctrl_word, "OUTDIR") == 0)
      OUTDIR = const_str;
    else if (strcmp(ctrl_word, "SPACING") == 0)
      SPACING = const_float;
#ifdef CILE
    else if (strcmp(ctrl_word, "CILER") == 0)
      CILER = const_float;
#endif
    else if (strcmp(ctrl_word, "PROBE") == 0)
      PROBE = const_float;
    else if (strcmp(ctrl_word, "MAXR") == 0)
      MAXR = const_float;
    else if (strcmp(ctrl_word, "CATOM") == 0)
      CATOM = const_float;
    else if (strcmp(ctrl_word, "NATOM") == 0)
      NATOM = const_float;
    else if (strcmp(ctrl_word, "OATOM") == 0)
      OATOM = const_float;
    else if (strcmp(ctrl_word, "SATOM") == 0)
      SATOM = const_float;
    else if (strcmp(ctrl_word, "HATOM") == 0)
      HATOM = const_float;
    else if (strcmp(ctrl_word, "EPS") == 0)
      EPS = const_float;
    else if (strcmp(ctrl_word, "SURF") == 0)
      SURF = const_float;
    else if (strcmp(ctrl_word, "DESC_PROBE_DIST") == 0)
      DESC_PROBE_DIST = const_float;
    else if (strcmp(ctrl_word, "BFACTOR") == 0)
      BFACTOR = const_float;
    else if (strcmp(ctrl_word, "POP_MEAN") == 0)
      POP_MEAN = const_float;
    else if (strcmp(ctrl_word, "POP_SD") == 0)
      POP_SD = const_float;
    else if (strcmp(ctrl_word, "NFP") == 0)
      NFP = const_float;
    else if (strcmp(ctrl_word, "RESOLUTION") == 0)
      RESOLUTION = const_float;
    else if (strcmp(ctrl_word, "OGA") == 0)
      OGA = const_float;
    else if (strcmp(ctrl_word, "OGB") == 0)
      OGB = const_float;
    else if (strcmp(ctrl_word, "OGC") == 0)
      OGC = const_float;
    else if (strcmp(ctrl_word, "OGV") == 0)
      OGV = const_float;
    else if (strcmp(ctrl_word, "THRMSD") == 0)
      THRMSD = const_float;
    else if (strcmp(ctrl_word, "CUTOFF_FIRST") == 0)
      CUTOFF_FIRST = const_float;
    else if (strcmp(ctrl_word, "CUTOFF_SECOND") == 0)
      CUTOFF_SECOND = const_float;
//    else if (strcmp(ctrl_word, "PCUT") == 0)
//      PCUT = const_float;
    else if (strcmp(ctrl_word, "NORM_EPS") == 0)
      NORM_EPS = const_float;
    else if (strcmp(ctrl_word, "MATRIX_CUTOFF") == 0)
      MATRIX_CUTOFF = const_float;
    else if (strcmp(ctrl_word, "VECTOR_CUTOFF") == 0)
      VECTOR_CUTOFF = const_float;
    else if (strcmp(ctrl_word, "RMSD_INCR") == 0)
      RMSD_INCR = const_float;
    else if (strcmp(ctrl_word, "RMSD_INCR_ADD") == 0)
      RMSD_INCR_ADD = const_float;
    else if (strcmp(ctrl_word, "CLIQUE_TLIMIT") == 0)
      CLIQUE_TLIMIT = const_float;
    else if (strcmp(ctrl_word, "CLUS_SPEC") == 0)
      CLUS_SPEC = (int) const_float;
    else if (strcmp(ctrl_word, "PAIRWISE") == 0)
      _pairwise = true;
    else if (strcmp(ctrl_word, "LOCAL") == 0)
      _local = true;
    else if (strcmp(ctrl_word, "NOFP") == 0)
      _nofp = true;
    else if (strcmp(ctrl_word, "CLUS_TO_OUTPUT") == 0)
      CLUS_TO_OUTPUT = (int) const_float;
    else if (strcmp(ctrl_word, "SCONS") == 0)
      SCONS = (int) const_float;
    else if (strcmp(ctrl_word, "NCPU") == 0)
      NCPU = (int) const_float;
    else if (strcmp(ctrl_word, "ALIGNMENT_NO") == 0)
      ALIGNMENT_NO = (int) const_float;
    else if (strcmp(ctrl_word, "BLOSUM") == 0)
      strcpy(BLOSUM, const_str);
    else if (strcmp(ctrl_word, "MOTIF1") == 0) {
      MOTIF1 = Motif::debracket(buffer);
      _motif1 = true;
    }
    else if (strcmp(ctrl_word, "MOTIF2") == 0) {
      MOTIF2 = Motif::debracket(buffer);
      _motif2 = true;
    }
    else if (strcmp(ctrl_word, "BSITE1") == 0) {
      BSITE1 = const_str;
      _bsite1 = true;
    }
    else if (strcmp(ctrl_word, "BSITE2") == 0) {
      BSITE2 = const_str;
      _bsite2 = true;
    }
    else if (strcmp(ctrl_word, "SURF_VECTOR_ANGLE") == 0)
      SURF_VECTOR_ANGLE = const_float;
    else if (strcmp(ctrl_word, "BLOSUM_SCORE") == 0)
      BLOSUM_SCORE = const_float;
    else if (strcmp(ctrl_word, "CALPHA_RMSD") == 0)
      CALPHA_RMSD = const_float;
    else if (strcmp(ctrl_word, "K") == 0)
      K = const_float;
    else if (strcmp(ctrl_word, "LAMBDA") == 0)
      LAMBDA = const_float;
    else if (strcmp(ctrl_word, "Z_SCORE") == 0)
      Z_SCORE = const_float;
    else if (strcmp(ctrl_word, "Z_SCORE_CONS") == 0)
      Z_SCORE_CONS = const_float;
    else if (strcmp(ctrl_word, "Z_SCORE_FP") == 0)
      Z_SCORE_FP = const_float;
    else if (strcmp(ctrl_word, "GLOBAL_DIST") == 0)
      GLOBAL_DIST = const_float;
    else if (strcmp(ctrl_word, "GLOBAL_SVAL") == 0)
      GLOBAL_SVAL = const_float;
    else if (strcmp(ctrl_word, "GLOBAL_ANGLE") == 0)
      GLOBAL_ANGLE = const_float;
    else if (strcmp(ctrl_word, "LIG_NEIGHB") == 0)
      LIG_NEIGHB = const_float;
    else if (strcmp(ctrl_word, "PP_DIST") == 0)
      PP_DIST = const_float;
    else if (strcmp(ctrl_word, "SL_DIST") == 0)
      SL_DIST = const_float;
    else if (strcmp(ctrl_word, "NU_DIST") == 0)
      NU_DIST = const_float;
    else if (strcmp(ctrl_word, "IO_DIST") == 0)
      IO_DIST = const_float;
    else if (strcmp(ctrl_word, "NCLTRIES") == 0)
      NCLTRIES = (int) const_float;
    else if (strcmp(ctrl_word, "NMOL") == 0)
      NMOL = (int) const_float;
    else if (strcmp(ctrl_word, "LIGDIR") == 0) {
      LIGDIR = const_str;
      _lig = true;
	}
    else if (strcmp(ctrl_word, "NOSQL_FILE") == 0)
      NOSQL_FILE = const_str;
    else if (strcmp(ctrl_word, "JSON_FILE") == 0)
      JSON_FILE = const_str;
  }
  par.close();
}

void Args::print_modifiers() {
  cout << "MODIFIERS USED: ";
  if (_noclus) cout << "noclus ";
  if (_noprune) cout << "noprune ";
  if (_srf) cout << "srf ";
  if (_verbose) cout << "verbose ";
  if (_local) cout << "local ";
  if (_nobb) cout << "nobb ";
  if (_nomarkbb) cout << "nomarkbb ";
  if (_super) cout << "super ";
  if (_caonly) cout << "caonly ";
  if (_motif1) cout << "motif1 ";
  if (_motif2) cout << "motif2 ";
  if (_bsite1) cout << "bsite1 ";
  if (_bsite2) cout << "bsite2 ";
  if (_longnames) cout << "longnames ";
  if (_database) cout << "database ";
  if (_pairwise) cout << "pairwise ";
  if (_nofp) cout << "nofp ";
  cout << endl;
}

void Args::print_constants() {
  cout << "PARAMETERS USED - To modify, copy rows between the" << endl;
  cout << "dashed lines to a separate input parameter file -"  << endl; 
  cout << "----------------------------------------------------------" << endl;
  cout << "NCPU		                " << NCPU               << endl;
  cout << "SRF_FILE                " << SRF_FILE               << endl;
  cout << "PROTEIN1                " << PROTEIN1               << endl;
  cout << "PROTEIN2                " << PROTEIN2               << endl;
  cout << "CHAIN1                  " << CHAIN1                 << endl;
  cout << "CHAIN2                  " << CHAIN2                 << endl;
  cout << "INTER_CHAIN_DIST        " << INTER_CHAIN_DIST       << endl;
  cout << "SURF_FILE               " << SURF_FILE              << endl;
  cout << "INDIR                   " << INDIR                  << endl;
  cout << "OUTDIR                  " << OUTDIR                 << endl;
  cout << "SPACING                 " << SPACING                << endl;
#ifdef CILE
  cout << "CILER                   " << CILER                  << endl;
#endif
  cout << "PROBE                   " << PROBE                  << endl;       
  cout << "MAXR                    " << MAXR                   << endl;        
  cout << "CATOM                   " << CATOM                  << endl;
  cout << "NATOM                   " << NATOM                  << endl;
  cout << "OATOM                   " << OATOM                  << endl;
  cout << "SATOM                   " << SATOM                  << endl;
  cout << "HATOM                   " << HATOM                  << endl;
  cout << "EPS                     " << EPS                    << endl;
  cout << "SURF                    " << SURF                   << endl;
  cout << "DESC_PROBE_DIST         " << DESC_PROBE_DIST        << endl;
  cout << "BFACTOR                 " << BFACTOR                << endl;
  cout << "POP_MEAN                " << POP_MEAN               << endl;
  cout << "POP_SD                  " << POP_SD                 << endl;
  cout << "NFP                     " << NFP                    << endl;
  cout << "SITE_RMSD               " << SITE_RMSD              << " *" << endl;
  cout << "RESOLUTION              " << RESOLUTION             << endl;
  cout << "OGA                     " << OGA                    << endl;
  cout << "OGB                     " << OGB                    << endl;
  cout << "OGC                     " << OGC                    << endl;
  cout << "OGV                     " << OGV                    << endl;
  cout << "THRMSD                  " << THRMSD                 << endl;
  cout << "CUTOFF_FIRST            " << CUTOFF_FIRST           << endl;
  cout << "CUTOFF_SECOND           " << CUTOFF_SECOND          << endl;
//  cout << "PCUT                    " << PCUT                   << endl;
  cout << "NORM_EPS                " << NORM_EPS               << endl;
  cout << "MAX_VERTICES            " << MAX_VERTICES           << " *" << endl;
  cout << "MATRIX_CUTOFF           " << MATRIX_CUTOFF          << endl;
  cout << "VECTOR_CUTOFF           " << VECTOR_CUTOFF          << endl;
  cout << "RMSD_INCR               " << RMSD_INCR              << endl;
  cout << "RMSD_INCR_ADD           " << RMSD_INCR_ADD          << endl;
  cout << "WORD                    " << WORD                   << " *" << endl;
  cout << "CLIQUE_TLIMIT           " << CLIQUE_TLIMIT          << endl;
  cout << "CLUS_SPEC               " << CLUS_SPEC              << endl;           
  cout << "CLUS_TO_OUTPUT          " << CLUS_TO_OUTPUT         << endl;           
  cout << "SCONS                   " << SCONS                  << endl;           
  cout << "ALIGNMENT_NO            " << ALIGNMENT_NO           << endl;           
  cout << "BLOSUM                  " << BLOSUM                 << endl;           
  cout << "MOTIF1                  " << MOTIF1                 << endl;           
  cout << "MOTIF2                  " << MOTIF2                 << endl;           
  cout << "BSITE1                  " << BSITE1                 << endl;           
  cout << "BSITE2                  " << BSITE2                 << endl;           
  cout << "SURF_VECTOR_ANGLE       " << SURF_VECTOR_ANGLE      << endl;           
  cout << "BLOSUM_SCORE            " << BLOSUM_SCORE           << endl;           
  cout << "CALPHA_RMSD             " << CALPHA_RMSD            << endl;           
  cout << "K                       " << K                      << endl;           
  cout << "LAMBDA                  " << LAMBDA                 << endl;           
  cout << "Z_SCORE                 " << Z_SCORE                << endl;           
  cout << "Z_SCORE_CONS            " << Z_SCORE_CONS           << endl;           
  cout << "Z_SCORE_FP              " << Z_SCORE_FP             << endl;           
  cout << "GLOBAL_DIST             " << GLOBAL_DIST            << endl;           
  cout << "GLOBAL_SVAL             " << GLOBAL_SVAL            << endl;           
  cout << "GLOBAL_ANGLE            " << GLOBAL_ANGLE           << endl;           
  cout  << "LIG_NEIGHB             " << LIG_NEIGHB             << endl;
  cout  << "PP_DIST                " << PP_DIST                << endl;
  cout  << "SL_DIST                " << SL_DIST                << endl;
  cout  << "NU_DIST                " << NU_DIST                << endl;
  cout  << "IO_DIST                " << IO_DIST                << endl;
  cout  << "NCLTRIES               " << NCLTRIES               << endl;
  cout  << "NMOL                   " << NMOL                   << endl;
  cout  << "LIGDIR                 " << LIGDIR                 << endl;
  cout  << "NOSQL_FILE             " << NOSQL_FILE             << endl;
  cout  << "JSON_FILE              " << JSON_FILE              << endl;
  cout << "----------------------------------------------------------" << endl;
  cout << "(*) These parameters cannot be changed." << endl;
}

Args::Args() {
  opts.insert(make_pair("-align","Read a rotational matrix of an alignment from an .nosql file and superimpose the two given proteins accordingly (first run -compare or -surfdb). Output the superimposed proteins' coordinates in a .pdb file. You need to provide both .pdb files that you want to superimpose (see -f1, -c1, -f2, c2 modifiers) and an alignment number (see -alno modifier)."));
  opts.insert(make_pair("-compare","Compare two protein surfaces (.pdb or .srf files). (If you use .pdb files, surfaces will be computed first.) Output their local structural alignments in an .nosql file. Each alignment consists of a rotational matrix, alignment scores, and aligned residues of the compared proteins."));
  opts.insert(make_pair("-extract","Calculate surface of a protein. Redirect the output (which is the surface) to a surface (.srf) file. Surface files can be used instead of .pdb files together with -compare or -surfdb options, which improves performance when doing repetitive comparisons (because surface does not need to be recalculated for each comparison). Option -surfdb works with .srf files exclusively!"));
  opts.insert(make_pair("-results","Read alignments from an .nosql file, filter them according to their scores, and calculate fingerprint residues (which can also be used as filter). Output results in Json format. In addition, replace B-factors in the query protein’s PDB file with the degrees of structural conservation. If used with -ligdir modifier, output ligands in Json format as well."));
  opts.insert(make_pair("-surfdb","Compare the query protein surface (.srf) with other protein surfaces listed in the SURF_FILE (see -sfile modifier). This does the same calculation as the -compare option, but faster, because protein surfaces are precalculated (see -extract option). This options also supports parallel computation on multiple CPUs. Output is the same as with -compare."));
  opts.insert(make_pair("-h","Show list of all parameters and their current values. You can copy/paste the parameters into a separate file and change their values (see -param modifier)."));
  mods.insert(make_pair("-alno ALIGNMENT_NO","Each comparison of a pair of proteins may result in many different local structural alignments.The alignment number can be 0 to 4 (see default CLUS_TO_OUTPUT parameter)."));
  mods.insert(make_pair("-ncpu NCPU","Use this number of concurrent threads."));
  mods.insert(make_pair("-bsite BSITE, -bsite1 BSITE, or -bsite2 BSITE","This selects protein residues in a certain radius (set by -dist) around the given ligand, and takes these residues as input (use with -extract or -compare options). If used with -compare, it only works with .pdb files (not .srf).For example: ‘-bsite ATP.305.A’ - ATP (residue name), 305 (residue number), A (chain id) or ‘-bsite *.*.B’ - chain B is the ligand (so you can select protein-protein binding sites as input)"));
  mods.insert(make_pair("-c1, c2","Chain identifiers of the compared proteins. You may give multiple chains, e.g., '-c1 ABC'."));
  mods.insert(make_pair("-database","Used by the web server (with -surfdb and -compare options). It will output an .nosql file as usual, and additional .nosql file with inversed rotational matrices, whose lines are marked with asterisks."));
  mods.insert(make_pair("-dist INTER_CHAIN_DIST","The distance between protein chains or between ligand and protein. Use with -bsite modifier or -mark and -results options. "));
  mods.insert(make_pair("-srffile","Surface file (.srf) for input or output."));
  mods.insert(make_pair("-f1, f2","Compared proteins' files (.pdb or .srf)."));
  mods.insert(make_pair("-in INDIR","Directory where input files are."));
  mods.insert(make_pair("-local","Use this to perform local alignments search (with -compare or -surfdb options). By default, after the local alignment is found (with maximum clique algorithm), an attempt is made to extend this alignment along the backbones of the compared proteins. In this way, parts of proteins that adopt different conformations (e.g., loops) can be aligned (these aligned residues are marked with 'flx' in alignments.json)."));
  mods.insert(make_pair("-longnames","Use this to allow long file names. By default the protein names are trimmed down to 4 letters. "));
  mods.insert(make_pair("-motif MOTIF, -motif1 MOTIF, or -motif2 MOTIF","This selects residues to be used as a query instead of the whole protein structure (use with -extract or -compare options). This will generate a .srf file with only the selected residues. To select some residues on chains A and B of the input protein, use -motif \"[:A and (14,57,69-71) or :B and (33,34,50)]\". Note that chain Iids are case sensitive. Square brackets are mandatory!"));
  mods.insert(make_pair("-nobb","Do not include descriptors originating from backbone atoms."));
  mods.insert(make_pair("-nomarkbb","Turns off the default action, which is to mark (not delete) backbone descriptors. In the first step of filtering, only non-backbone descriptors are used, while in the maximum clique step, all descriptors are used."));
  mods.insert(make_pair("-noclus","Local structural alignments found (maximum cliques) are not clustered."));
  mods.insert(make_pair("-nofp","Do not calculate fingerprint residues. Do not filter by fingerprint residues (use with -results option)."));
  mods.insert(make_pair("-noprune","Alignments are not pruned. By default bad scoring cliques are deleted (see SURF_VECTOR_ANGLE, BLOSUM_SCORE and CALPHA_RMSD parameters)."));
  mods.insert(make_pair("-out OUTDIR","Directory to write output files to."));
  mods.insert(make_pair("-param PAR_FILE","Read parameters from the specified parameter file."));
  mods.insert(make_pair("-sfile SURF_FILE","Specify file that contains names of .srf files to be compared with the query protein (see -surfdb option). Each line must contain one .srf file-name.Example:{newline}protein1.srf A{newline}protein2.srf B{newline}protein3.srf A{newline}etc."));
  mods.insert(make_pair("-super","Find local structural alignments between two proteins (use with -compare option) and superimpose the two proteins according all alignments found.  For each alignment, output the '.rota.pdb' file with the proteins superimposed according to this alignment."));
  mods.insert(make_pair("-verbose","Output debugging information. Use when testing the program."));
  mods.insert(make_pair("-z_score Z_SCORE","The cutoff value for z_score. Low z_score (<2) means that more insignificant alignments will be outputted (these can also occur by chance), higher z_score (>2) means only significant alignments will be outputted (use with -results option)."));
  mods.insert(make_pair("-ligdir DIRNAMES","Directory where the names of proteins as files: each file is PdbIdChainId, e.g., 1allA."));
  mods.insert(make_pair("-nosql NOSQL_FILE","Output / read nosql alignments to / from this file."));
  mods.insert(make_pair("-json JSON_FILE","Output json alignments to this file."));
}

ostream& operator<< (ostream &out, Args &arg) {
  //
  // Kopiraj iz google docs v opts.txt in potem: cat opts.txt|sed -r 's/^-(.*)/"));\nopts.insert(make_pair("-\1",/'|sed 's/^$/"/'
  // Izbrisi vse \n , nadomesti z {newline} !
  //
  cout << "Protein Binding Sites by Local Structural Alignments --- Ver. 2.4.4 --- Sep 11 2014 ---" << endl;
  cout << "Usage: probis {OPTION} [MODIFIERS] -f1 PDB_FILE -c1 CHAIN_ID [-f2 PDB_FILE -c2 CHAIN_ID]" << endl;
	unsigned int n = thread::hardware_concurrency();
	cout << endl << "Detected support for " << n << " concurrent threads." << endl;

  out << endl << endl << "--------------------------------------- OPTIONS --------------------------------------------------------------------------------" << endl<<endl;
  for (map<string,string>::iterator it = arg.opts.begin(); it != arg.opts.end(); it++) {
    vector<string> s;
    size_t i = 0, i_new = 0;
    while ((i_new = it->second.find_first_of(" \n", i + 60)) != string::npos) {
      out <<setw(47)<< (i==0 ? it->first: "") << setw(10) <<"" << it->second.substr(i, i_new - i) << endl;
      i = i_new;
    }
    out <<setw(47)<< (i==0 ? it->first: "") << setw(10) << "" << it->second.substr(i, it->second.size()) << endl << endl;
 }
  out << endl << "---------------------------------------- MODIFIERS --------------------------------------------------------------------------------" << endl<<endl;   
  for (map<string,string>::iterator it = arg.mods.begin(); it != arg.mods.end(); it++) {
    vector<string> s;
    size_t i = 0, i_new = 0;
    while ((i_new = it->second.find_first_of(" \n", i + 60)) != string::npos) {
      out <<setw(47)<< (i==0 ? it->first: "") << setw(10) <<"" << it->second.substr(i, i_new - i) << endl;
      i = i_new;
    }
    out <<setw(47)<< (i==0 ? it->first: "") << setw(10) << "" << it->second.substr(i, it->second.size()) << endl<< endl;
  }
  out << endl << endl;
  return out;
}

void Args::fill_args(int argc, char *argv[]) {
  char str[XLARGE];
  string par_file;
  char *pch;
  for (int i = 1; i < argc; i++) {
    strcat(str, argv[i]);
    strcat(str, " ");
  }
  if (argc < 2) {
    cout << *this;
    throw Err("Error (ARGS) : Wrong number of arguments.", 8); 
  }
  /*read options */
  if (strstr(str, "-compare")) {
    _state = TOTAL;
  }
  else if (strstr(str, "-extract")) {
    _state = WRITE;
  }
  else if (strstr(str, "-surfdb")) {
    _state = SURFDB;
  }
  else if (strstr(str, "-mark")) {
    _state = MARK;
  }
  else if (strstr(str, "-lig ")) { // <- presledek mora biti za "lig" drugace je zmesnjava z "ligdir"
    _state = LIGAND;
  }
  else if (strstr(str, "-header")) {
    _state = HEAD;
  }
  else if (strstr(str, "-allbio")) {
    _state = ALLBIO;
  }
#ifdef CILE
  else if (strstr(str, "-patch")) {
    _state = PATCH;
  }
#endif
  else if (strstr(str, "-results")) {
    _state = RESULTS;
  }
  else if (strstr(str, "-align")) {
    _state = ALIGN;
  }
  else if (strstr(str, "-h")) {
    cout << *this;
//    error_message();
    print_constants();
    throw Err("", 9);
  }

  /*read modifiers */
  if ((pch = strstr(str, "-param"))) {
    char pf[MEDIUM];
    sscanf(pch, "%*s %s", pf);
    read_constants(pf);
  }
  if ((pch = strstr(str, "-in"))) {
    char indir[MEDIUM];
    sscanf(pch, "%*s %s", indir);
    INDIR = indir;
  }
  if ((pch = strstr(str, "-out"))) {
    char outdir[MEDIUM];
    sscanf(pch, "%*s %s", outdir);
    OUTDIR = outdir;
  }
  if ((pch = strstr(str, "-ligdir"))) {
    char ligdir[MEDIUM];
    sscanf(pch, "%*s %s", ligdir);
    LIGDIR = ligdir;
    _lig = true;
  }
  if ((pch = strstr(str, "-nosql"))) {
    char nosql_file[MEDIUM];
    sscanf(pch, "%*s %s", nosql_file);
    NOSQL_FILE = nosql_file;
  }
  if ((pch = strstr(str, "-json"))) {
    char json_file[MEDIUM];
    sscanf(pch, "%*s %s", json_file);
    JSON_FILE = json_file;
  }
  if ((pch = strstr(str, "-motif ")) || (pch = strstr(str, "-motif1"))) {
    MOTIF1 = Motif::debracket(pch);
    _motif1 = true;
  }
  if ((pch = strstr(str, "-motif2"))) {
    MOTIF2 = Motif::debracket(pch);
    _motif2 = true;
  }
  if ((pch = strstr(str, "-bsite ")) || (pch = strstr(str, "-bsite1"))) {
    char bsite[MEDIUM];
    sscanf(pch, "%*s %s", bsite);
    BSITE1 = bsite;
    _bsite1 = true;
  }
  if ((pch = strstr(str, "-bsite2"))) {
    char bsite[MEDIUM];
    sscanf(pch, "%*s %s", bsite);
    BSITE2 = bsite;
    _bsite2 = true;
  }
  if ((pch = strstr(str, "-longnames"))) {
    _longnames = true;
  }
  if ((pch = strstr(str, "-nofp"))) {
    _nofp = true;
  }
  if ((pch = strstr(str, "-dist"))) {
    sscanf(pch, "%*s %f", &INTER_CHAIN_DIST);
  }
  if ((pch = strstr(str, "-z_score"))) {
    sscanf(pch, "%*s %f", &Z_SCORE);
  }
  if ((pch = strstr(str, "-bfactor"))) {
    sscanf(pch, "%*s %f", &BFACTOR);
  }
  if ((pch = strstr(str, "-ncpu"))) {
    sscanf(pch, "%*s %d", &NCPU);
  }
  if ((pch = strstr(str, "-srffile"))) {
    char tmp_srf_file[MEDIUM];
    sscanf(pch, "%*s %s", tmp_srf_file);
    SRF_FILE = tmp_srf_file;
  }

  if ((pch = strstr(str, "-alno"))) {
    sscanf(pch, "%*s %d", &ALIGNMENT_NO);
  }
  if ((pch = strstr(str, "-sfile"))) {
    //~ sscanf(pch, "%*s %s", SURF_FILE);
    char surf_file[MEDIUM];
    sscanf(pch, "%*s %s", surf_file);
    SURF_FILE = surf_file;
  }
  if ((pch = strstr(str, "-noclus"))) {
    _noclus = true;
  }
  if ((pch = strstr(str, "-noprune"))) {
    _noprune = true;
  }
  if ((pch = strstr(str, "-verbose"))) {
    _verbose = true;
  }
  if ((pch = strstr(str, "-local"))) {
    _local = true;
  }
  if ((pch = strstr(str, "-nobb"))) {
    _nobb = true;
  }
  if ((pch = strstr(str, "-nomarkbb"))) {
    _nomarkbb = true;
  }
  if ((pch = strstr(str, "-super"))) {
    _super = true;
  }
  if ((pch = strstr(str, "-caonly"))) {
    _caonly = true;
  }
  if ((pch = strstr(str, "-database"))) {
    _database = true;
  }

  /*read pdb file(s) or srf files */
  if ((pch = strstr(str, "-f1"))) {
    char tmp_PROTEIN1[MEDIUM];
    sscanf(pch, "%*s %s", tmp_PROTEIN1);
    PROTEIN1 = tmp_PROTEIN1;

    if ((pch = strstr(str, "-c1"))) {
      char tmp_CHAIN1[SMALL];
      sscanf(pch, "%*s %s", tmp_CHAIN1);
      CHAIN1 = tmp_CHAIN1;
    }
    else {
      throw Err("Error (ARGS) : Please provide chain identifier!", 15);
    }
  }

  if ((pch = strstr(str, "-f2"))) {
    char tmp_PROTEIN2[MEDIUM];
    sscanf(pch, "%*s %s", tmp_PROTEIN2);
    PROTEIN2 = tmp_PROTEIN2;

    if ((pch = strstr(str, "-c2"))) {
      char tmp_CHAIN2[SMALL];
      sscanf(pch, "%*s %s", tmp_CHAIN2);
      CHAIN2 = tmp_CHAIN2;

    }
    else {
      throw Err("Error (ARGS) : Please provide chain identifier! ", 13);
    }

  }


  bool srf1 = false, srf2 = false;
  if (PROTEIN1.find("srf")!=string::npos||PROTEIN1.find("SRF")!=string::npos)
    srf1 = true;
  if (PROTEIN2.find("srf")!=string::npos||PROTEIN2.find("SRF")!=string::npos)
    srf2 = true;
  if (!PROTEIN2.empty() && srf1 != srf2)
    throw Err("Error (ARGS) : Both proteins should be in the same format (.srf or .pdb) !", 11);
  else if (srf1)
    _srf = true;
  
}

void Args::print_state() {
  switch(_state) {
  case TOTAL  : cout << "--- Complete proteins' search ---" << endl; break;
//  case PPI    : cout << "--- Protein-protein interaction search ---" << endl; break;
//  case BIND   : cout << "--- Binding cavity search ---" << endl; break;
  case WRITE  : cout << "--- Writing out surface of a protein ---" << endl; break;
  case SURFDB : cout << "--- Database search for a complete protein ---" << endl; break;  // obsolete
  case MARK   : cout << "--- Find interface residues ---" << endl; break;
  case LIGAND : cout << "--- Find existing binding sites & ligands ---" << endl; break;
  case RESULTS  : cout << "--- Reading results from nosql and outputing alignments as Json ---" << endl; break;
//  case SEQ    : cout << "--- Reading sequence ---" << endl; break;
  case ALIGN  : cout << "--- Aligning the protein structure ---" << endl; break;
//  case BIOU   : cout << "--- Generating biological units ---" << endl; break;
  case HEAD   : cout << "--- Ligand names are output ---" << endl; break;
  case ALLBIO : cout << "--- Output all biological units for a given PDB file ---" << endl; break;
#ifdef CILE
  case PATCH   : cout << "--- Generate surface patches ---" << endl; break;
#endif
    //  case SUPER  : cout << "--- Superimposing two protein structures ---" << endl; break;
  }
}
