#include "states.h"
#include "args.h"
#include "molecule.h"
#include "grid.h"
#include "desc.h"
#include "score.h"
#include "probe.h"
#include "clique.h"
#include "product.h"
#include "subgraph.h"
#include "cluster.h"
#include "output.h"
#include "ligands.h"
#include "item.h"
#include "motif.h"
#include "nosql.h"
#include "debug.hpp"



void get_motif(Molecule *m, string motif_sele) {
	Motif motif(motif_sele);
	// oznacimo atome, ki niso v motivu za visited, tako jih all_triples ne bo uposteval pri izracunu povrsine
	m->mark_motif_atoms(&motif);
	}

void get_bsite(Molecule *m, string bsite_sele, string protein, string chain_id) {
	try {
		/* preberi kompleks v PROTEIN in izracunaj dejanski interface */
		Molecule m2(protein, chain_id);
		m2.read_plain_pdb(_firstModel);
		// preberemo celoten pdb (s hetero in nucleic) v novo molecule
		m2.read_PDB(_firstModel, _nucleic, _allChains, _saveMemLig);
		m2.append_coor_delete_comp();
		Grid g;
		g.init( m2.minmax_coor() );
		m2.init_grid_molecule(&g, 1.1 * INTER_CHAIN_DIST); // take a bit bigger boxes 
		// dolocimo aminokisline v okolici liganda - to je motif
		m2.interface_models(INTER_CHAIN_DIST);
		// oznacimo atome, ki niso v motif za visited
		m->mark_binding_site(&m2, bsite_sele);
	}
	catch (Err e) {
		throw e;
	}
}


void state0(string protein2, string chain2) {
	NoSql::WriteStruct strings;
	//~ cout << "before state0 protein2= " << protein2 << " chain2=" << chain2 << endl;
	state0(strings, protein2, chain2);
	//~ cout << "out of state0 protein2= " << protein2 << " chain2=" << chain2 << endl;
	if (strings.suffix != "") {
		//~ NoSql nosql(strings.mol1_pdb_id, strings.mol1_chain_id);
		NoSql nosql(NOSQL_FILE);
		nosql.write(strings);
	}
}

void state0(NoSql::WriteStruct& strings, string protein2, string chain2) {

	//~ cout << "in state0 protein2= " << protein2 << " chain2=" << chain2 << endl;
	Probe *probe1 = nullptr, *probe2 = nullptr;
	Descriptor *desc1 = nullptr, *desc2 = nullptr;
	try {
		//~ cout << "in state0 protein2= " << protein2 << " chain2=" << chain2 << endl;
		//~ cout << "desc2 = " << desc2 << endl;
		Grid grid1, grid2;
		//FIRST MOLECULE
		Molecule mol1(PROTEIN1, CHAIN1);
		if ( _srf ) { 
			/* reading from .srf file */
			//~ cout << "restoring " << PROTEIN1 << " " << CHAIN1 << " _srf = " << _srf << endl;
			//~ mol1.parse_srf_file(SRF_FILE);
			mol1.restore_atoms();
			mol1.restore_descriptors(desc1, _nobb);
			mol1.restore_probes(probe1);
			grid1.init( mol1.minmax_coor() );
		}
		else {	    
			/* reading from .pdb file */
			mol1.read_plain_pdb(_firstModel);
			mol1.read_PDB(_firstModel, _noHetero, _selectedChains);
			mol1.append_coor(); // ALL CHAINS
			
			/* prvic inicializiramo grid, zato moramo dolociti min in max koordinato */
			grid1.init( mol1.minmax_coor() );
			mol1.init_grid_molecule(&grid1, 2*PROBE);
			
			// ce je --motif modifier
			if (_motif1) get_motif(&mol1, MOTIF1);
			// ce je --bsite modifier, oznacimo atome, ki niso del tega binding site-a za visited
			if (_bsite1) get_bsite(&mol1, BSITE1, PROTEIN1, CHAIN1);
			
			mol1.all_triples(probe1);
			mol1.enumerate_clefts(probe1);
			mol1.surface_atoms(&grid1, probe1);
			mol1.set_descriptors(desc1, _nobb);
			
			/* izbrisemo atome iz grida */
			grid1.deallocate_content();
		}
		dbgmsg(protein2 << " " << chain2);

		/* ce smo brali srf file, je grid ze inicializiran (glej zgoraj) */
		desc1->init_grid_descriptors(&grid1, CUTOFF_FIRST, CUTOFF_FIRST);
		grid1.deallocate_content();
		
		//~ cout << "after grid deallocate" << endl;
		//SECOND MOLECULE
		Molecule mol2(protein2, chain2);
		if ( _srf ) {     
			/* reading from .srf file */
			mol2.restore_atoms();
			mol2.restore_descriptors(desc2, _nobb);
			//~ cout << "before score" << endl;
			mol2.restore_probes(probe2);
			grid2.init( mol2.minmax_coor() );
		}
		else {	    
			/* reading from .pdb file */
			mol2.read_plain_pdb(_firstModel);
			mol2.read_PDB(_firstModel, _noHetero, _selectedChains);
			mol2.append_coor(); // ALL CHAINS
			
			grid2.init( mol2.minmax_coor() );
			mol2.init_grid_molecule(&grid2, 2*PROBE);
			
			// ce je --motif modifier
			if (_motif2) get_motif(&mol2, MOTIF2);
			// ce je --bsite modifier, oznacimo atome, ki niso del tega binding site-a za visited
			if (_bsite2) get_bsite(&mol2, BSITE2, protein2, chain2);
			
			mol2.all_triples(probe2);
			mol2.enumerate_clefts(probe2);
			mol2.surface_atoms(&grid2, probe2);
			
			mol2.set_descriptors(desc2, _nobb);
			
			grid2.deallocate_content();
		}
	  
		desc2->init_grid_descriptors(&grid2, CUTOFF_FIRST, CUTOFF_FIRST);
		grid2.deallocate_content();
		
		Score score;
		score.get_internal_blosum();
		
		/* ce je podana opcija --markbackbone, oznacimo deskriptorje, ki pripadajo backbone atomom */
		if (!_nomarkbb) {
			desc1->mark_backbone();
			desc2->mark_backbone();
		}
		 
		/* naredimo produktni graf za oba seta deskriptorjev */
		Product product;
		product.tetrahedron(&score, desc1, desc2, TOTAL);
		
		desc1->init_grid_descriptors(&grid1, CUTOFF_SECOND, CUTOFF_SECOND);
		grid1.deallocate_content();
		desc2->init_grid_descriptors(&grid2, CUTOFF_SECOND, CUTOFF_SECOND);
		grid2.deallocate_content();
		
		product.init_product(&score);
		product.count_subgraphs();
		
		Subgraph subgraph(&product);

		unique_ptr<Clique> clique (new Clique());
		clique->max_clique(&subgraph);
		
		score.score_descriptor_ratio(&subgraph, desc1, desc2); 
		
		//PROBE CENTERS NEIGHBORING DESCRIPTORS ASSIGNED TO DESCRIPTORS
		desc1->init_grid_probe(&grid1, mol1.get_lcolor(), probe1, DESC_PROBE_DIST, DESC_PROBE_DIST);
		desc2->init_grid_probe(&grid2, mol2.get_lcolor(), probe2, DESC_PROBE_DIST, DESC_PROBE_DIST);
		
		score.score_surface_vector(&subgraph, probe1, probe2);
		score.score_probe(&subgraph, probe1);
		
		mol1.get_seq();
		mol2.get_seq();
		
		score.E_score(&mol1, &mol2, subgraph.qmax);
		score.cluster_score(subgraph.qmax);
		  
		if (!_noprune) subgraph.delete_bad_scoring();
		
		subgraph.sort_qmax(by_cluster_score);
		
		Output out;

		if (_verbose) subgraph.output_qmax(&out, &score);
		  
#ifndef NDEBUG
		for (unsigned int i = 0; i < subgraph.qmax.size(); i++)
			dbgmsg(i << "-th qmax size = " << subgraph.qmax[i]->vert.size());
#endif

		if (!_noclus) {
			/* cluster cliques */
			Cluster cluster;
			cluster.cluster(&subgraph);
			if (!_local) {
				if (subgraph.clus.size())
					subgraph.extend_all(&score, &mol1, &mol2, &product, desc1, desc2);
			}
			score.E_score(&mol1, &mol2, subgraph.clus);
			score.cluster_score(subgraph.clus);
			
			subgraph.sort_clus(by_cluster_score);
			if (_verbose) subgraph.output_clus(&out, &score);
		}

		if (!_noclus && !_super) {
			//~ NoSql nosql(&mol1);
			NoSql nosql(&mol1, NOSQL_FILE);
			nosql.writeStruct(&mol2, &subgraph, &score, strings);
		}

		/* prilegamo protein2 na protein1 kakor narekuje rota-trans matrika prvega klastra */
		if (!_srf && _super) {
		
			if (!_noclus) {
				for  (size_t i=0; i < subgraph.clus.size(); i++) {
					mol2.remark_pdb(subgraph.clus[i]); 
					mol2.output_rotated_asym(&mol1, subgraph.clus[i], i);         // input je lahko samo asymmetric unit!! (ne smes inputat vec kot en MODEL !!!!!!!!!!!)
				}
			}
			else {
				/* ce ni klastrov, pa prve klike */
				for  (size_t i=0; i < subgraph.qmax.size(); i++) {
					mol2.remark_pdb(subgraph.qmax[i]); 
					mol2.output_rotated_asym(&mol1, subgraph.qmax[i], i); 
				}
			}
		}
		desc1->free(); desc2->free();
		probe1->free(); probe2->free();
	}  // .. konec try bloka
	// poskrbimo za izjeme ..
	catch (Err e) {
		desc1->free(); desc2->free();
		probe1->free(); probe2->free();
		cout << e.what() << endl;
	}
}

void state3(Args *args) {
  /* 
     Naredimo .srf zapis iz .pdb . V .srf-ju je dovoljen le en chain in samo prvi model ! 
  */
	Probe *probe1 = NULL;
	Descriptor *desc1 = NULL;

  try {
    //FIRST MOLECULE
	Molecule mol1(PROTEIN1, CHAIN1);
    //~ mol1 = new Molecule(PROTEIN1, CHAIN1);
    mol1.read_plain_pdb(_firstModel);
    mol1.read_PDB(_firstModel, _noHetero, _selectedChains);
    mol1.append_coor(); // ALL CHAINS
	Grid grid1;
	grid1.init( mol1.minmax_coor() );
    mol1.init_grid_molecule(&grid1, 2*PROBE);

    // ce delamo srf samo za izbrane residueje (--motif modifier)
    if (_motif1) get_motif(&mol1, MOTIF1);
    // ce je --bsite modifier, oznacimo atome, ki niso del tega binding site-a za visited
    if (_bsite1) get_bsite(&mol1, BSITE1, PROTEIN1, CHAIN1);
    
    mol1.all_triples(probe1);
    mol1.enumerate_clefts(probe1);
    mol1.surface_atoms(&grid1, probe1);
    mol1.set_descriptors(desc1, _nobb);
    
    /* izpisemo .srf file */
    Output output;
    output.srf_file(SRF_FILE, &mol1, desc1, probe1, _motif1 || _bsite1);
    //  mol1.output(_motif);
    //~ mol1.output();
    //~ desc1->output();
    //~ probe1->output(mol1.get_lcolor(), _motif1 || _bsite1);
    
	  desc1->free(); 
	  probe1->free();
  }  // .. konec try bloka
  
  // poskrbimo za izjeme ..
  catch (Err e) {
		desc1->free(); 
		probe1->free();	    
		cout << e.what() << endl;
  }
}


void state4(Args *args) {
	//~ ifstream surf((add_end_slash(INDIR) + string(SURF_FILE)).c_str());
	ifstream surf((add_end_slash(INDIR) + SURF_FILE).c_str());
	if (!surf.is_open()) {
		throw Err("Error (STATES) : Cannot open SURF_FILE " + SURF_FILE + "!", 23);
	}
	
	/* This piece of code reads lines from SURF_FILE, on each line is name of the surface file */
	// contents of the file are read into a vector of strings
	
	vector<vector<string> > surface_files;
	surface_files.resize(NCPU, vector<string>());
	int i = 0;
	do {
		std::string buffer;
		std::getline(surf, buffer);
		if (buffer.size() > 0) {
			//~ cout << buffer << endl;
			surface_files[i % NCPU].push_back(buffer);
			//~ cout << i % NCPU << " " << NCPU << " " <<surface_files[i % NCPU].back() << endl;
			i++;
		}
	}
	while (!surf.eof());
	
	// vector of strings is traversed; individual strings represent jobs 
	vector<thread> threads;
	for(int j = 0; j < NCPU; ++j) {
		auto &sfiles = surface_files[j];
		cout  << "THREAD> Initializing thread " << j << endl;
		threads.push_back(
			thread([&sfiles] () {
				for (auto &str : sfiles) {
				stringstream ss(str);
				string protein2, chain2;
				ss >> protein2 >> chain2;
				//~ cout << "in loop = " <<  protein2 << " " << chain2 << endl;
				state0(protein2, chain2);
			}
		}));
	}
	for(auto& thread : threads) {
	thread.join();
	}
}

void state5(Args *args) {
  /*
    Izracunamo specificity in sensitivity napovedi. Input je proteinski kompleks, 
    z v b-faktorju oznacenimi napovedanimi binding sitei.
  */

  /* preberi kompleks v PROTEIN1 in izracunaj dejanski interface */
	Molecule mol1(PROTEIN1, CHAIN1);
  //~ mol1 = new Molecule(PROTEIN1, CHAIN1);
  

  mol1.read_plain_pdb(_firstModel);
  mol1.read_PDB(_firstModel, _hetero, _allChains);
//  mol1.read_hetatm_coor();


  mol1.append_coor();
	Grid grid1;
  //~ grid1 = new Grid( mol1.minmax_coor() );
  grid1.init( mol1.minmax_coor() );
  mol1.init_grid_molecule(&grid1, 1.1 * INTER_CHAIN_DIST); // take a bit bigger boxes 
  //~ delete grid1;


  mol1.interface(INTER_CHAIN_DIST);
  mol1.fn_important_residues();


//  /* moramo ponastaviti chain_id, ker sedaj racunamo samo za query verigo */
  mol1.goodness_of_prediction();
  mol1.output_iface_pdb();
  //~ delete mol1;

}


void state6(Args *args) {
  /*
    Prilegali smo vse v klastru (>95% sequence similarity) na predstavnika. Sedaj preberemo file z modeli, 
    prvi model je predstavnik, vsi ostali pa so prilegani proteini. Za vsako query aminokislino potem dobimo
    vse ligande (proteine, ione, ...), ki so dosti blizu (<3 A). Izpisemo query protein in za vsako aminokislino
    njene ligande.
  */


  try {
    /* preberi kompleks v PROTEIN1 in izracunaj dejanski interface */
	Molecule mol1(PROTEIN1, CHAIN1);
    //~ mol1 = new Molecule(PROTEIN1, CHAIN1);
    mol1.read_plain_pdb(_allModels);
    mol1.read_PDB(_allModels, _nucleic, _allChains, _saveMemLig);
    mol1.append_coor_delete_comp();
    
    //~ grid1 = new Grid( mol1.minmax_coor() );
    Grid grid1;
    grid1.init( mol1.minmax_coor() );
    mol1.init_grid_molecule(&grid1, 1.1 * INTER_CHAIN_DIST); // take a bit bigger boxes 
//    delete grid1;
    
    mol1.interface_models(INTER_CHAIN_DIST);
    mol1.output_interface_models();
//    delete mol1;
    //~ close_resources(); 
  }
  // poskrbimo za izjeme ..
  catch (Err e) {
    cout << e.what() << endl;
    //~ close_resources(); 
  }

}


void state7(Args *args) {
  /*
    Preberemo iz nosql baze rezultate prileganj in jih izpisemo.
  */

  try {
    //~ subgraph = new Subgraph();
    Subgraph subgraph;
    //~ clique = new Clique();
    //~ Clique clique;
    //~ Clique clique = new Clique();
    
    /* v PROTEIN1 je naslov kompleksa, v CHAIN1 pa identifikator obravnavane verige v tem kompleksu */
	Molecule mol1(PROTEIN1, CHAIN1);
    //~ mol1 = new Molecule(PROTEIN1, CHAIN1);
    mol1.read_plain_pdb(_firstModel);
    /* preberi REMARK 350, v katerih so rotacijske matrike za bioloske unite, kasneje rabis koordinate atomov pri gridu */  
    mol1.read_PDB(_firstModel, _nucleic, _allChains);
    // inicializiramo atom 1) zaradi get_seq in 2) zaradi minmax f-je 
    mol1.append_coor();
    // preberemo aminokislinsko zaporedje query proteina, samo query verig (v CHAIN1)
    mol1.get_seq();
    
    // beremo iz denormalizirane baze
    //~ nosql = new NoSql(&mol1);
    //~ NoSql nosql(&mol1);
    NoSql nosql(&mol1, NOSQL_FILE);
    nosql.read();

/* NOT FOR CANDOCK
    // dolocis kolikokrat je vsak query residue ohranjen (upostevas le dobro prilegane sekvence)
    if (! (_pairwise || _nofp)) { // ce je pairwise ali nofp (no fingerprint), ne filtriramo po fingerprint residue-jih
      mol1.conservation(Z_SCORE_FP); // z_score > 3
      mol1.binarize_cons_seq();
      //dolocis fingerprint residueje
      mol1.fingerprint();
    }
*/    

/* NOT FOR CANDOCK
    // se enkrat dolocis ohranjenost tokrat upostevas vse prilegane sekvence (kar jih je ostalo)
    mol1.conservation(Z_SCORE_CONS); // z_score > 2
    // izpisi aminokislinsko zaporedje proteina in za vsako aminokislino, kolikokrat je ohranjena
    mol1.binarize_cons_seq();
    // izpisemo asymetric unit
    mol1.output_cons_pdb();
    // generiramo biounite za query protein in jih izpisemo
    mol1.gen_cons_biounit();
*/    
    
    /* kadar server poganja probis, izpise spletno stran z rezultati */
//    if (_html) {
    mol1.out_alignments(&subgraph, JSON_FILE);
/* NOT FOR CANDOCK
    mol1.out_info();
    mol1.out_query();
*/
//    }
    
    /* samo, ce je _html flag in ni _pairwise in je v parameters.inp nastavljen LIGDIR (_lig = true), dodamo informacijo o ligandih v html results datoteko */
//    if ( _html && !_pairwise && _lig ) {
//    if (!_pairwise && _lig ) {
//      lig = new Ligands();
//      lig->read(mol1);
//      lig->generate(mol1);
//      
//      /* poiscemo sredisce vsakega liganda */
//      lig->center(mol1);
//      
//      /* najdemo sosednje ligande */
//      grid1 = new Grid( mol1.minmax_coor() );
//      lig->init_grid_ligands(grid1, LIG_NEIGHB);
//    
//      /* pripravimo in izpisemo binding site in pripadajoce ligande */
//      lig->initialize_data();
//      lig->zip_points();
//      lig->cluster();
//      lig->output();
//      /* rotiramo .bu.pdb datoteke prileganih proteinov kot narekujejo matrike in vektorji za vsak klaster ligandov */
//      lig->write_rotate_ligands(mol1, _caonly);
//    }
    //~ close_resources();
  }  // .. konec try bloka

  // poskrbimo za izjeme ..
  catch (Err e) {
    cout << e.what() << endl;
    //~ close_resources(); 
  }
}

void state9(Args *args) {
  /*
    Align two proteins according to the rotation and translation a cluster imposes.
  */

  Item *one_clq = new Item();

  try {
    /* preberi PROTEIN1, to je nas query protein, ki ostane na mestu */
	Molecule mol1(PROTEIN1, CHAIN1);
    //~ mol1 = new Molecule(PROTEIN1, CHAIN1);
    mol1.read_plain_pdb(_firstModel);
     /* vrtel bos PROTEIN2 */
	Molecule mol2(PROTEIN2, CHAIN2);
    //~ mol2 = new Molecule(PROTEIN2, CHAIN2);
    mol2.read_plain_pdb(_firstModel);
    //~ nosql = new NoSql(&mol1);
    //~ NoSql nosql(&mol1);
    NoSql nosql(&mol1, NOSQL_FILE);
    nosql.get_rota_trans_resi(&mol2, one_clq);
    // output both PROTEIN1 and rotated PROTEIN2 in one file
    mol2.output_rotate_pdb(&mol1, one_clq); 
    //~ close_resources(); 
    delete one_clq;
  } 
  catch (Err e) {
    cout << e.what() << endl;
    //~ close_resources(); 
    delete one_clq;
  }
}

void state11(Args *args) {
  /*
    Preberemo in izpisemo kljucne podatke iz HEADER sekcije dane PDB datoteke.
  */

  try {
    /* preberi PROTEIN1, to je nas query protein, ki ostane na mestu */
    //~ lig = new Ligands();
    Ligands lig;
    // v PROTEIN1 je PDB, katerega header bomo prebrali, OUTDIR je direktorij, kamor zapisemo imena ligandov
    //~ lig->pdb_header(PROTEIN1);  
    lig.pdb_header(PROTEIN1);  
    //~ close_resources();
  } 
  catch (Err e) {
    cout << e.what() << endl;
    //~ close_resources(); 
  }
//
//  delete lig;
}


void state12(Args *args) {
  /*
    Preberemo navaden PDB (asymmetric unit) in naredimo biological assemblije, ki jih izpisemo vsakega
    v svojo datoteko.
  */
  try {
	Molecule mol1(PROTEIN1, CHAIN1);
    //~ mol1 = new Molecule(PROTEIN1, CHAIN1);
    mol1.read_plain_pdb(_firstModel);
    /* preberi REMARK 350, v katerih so rotacijske matrike za bioloske unite, kasneje rabis koordinate atomov pri gridu */  
    mol1.read_PDB(_firstModel, _nucleic, _allChains);
    /* generiramo biounite za query protein in jih izpisemo */
    mol1.gen_cons_biounit();
    //~ close_resources();
  }
  catch (Err e) {
    cout << e.what() << endl;
    //~ close_resources(); 
  }
//  delete mol1;
}

#ifdef CILE
void stateX(Args *args) {
  /*
    Preberemo 1 protein in izpisemo patche residuejev na njegovi površini. Vsaka površinska aminokislina predstavlja 
    en patch. Radij (velikost patcha) je variabilen, sredisce patcha je CA atom, prav tako se lahko spreminja 
    SURF parameter, ki doloca, kako globoko pod povrsino vzamemo aminokisline kot 'povrsinske'. 
  */

	Probe *probe1 = NULL;
  /* preberi PROTEIN1 */
	Molecule mol1(PROTEIN1, CHAIN1);
  //~ mol1 = new Molecule(PROTEIN1, CHAIN1);
  mol1.read_plain_pdb(_firstModel);
  mol1.read_PDB(_firstModel, _noHetero, _selectedChains);
  mol1.append_coor();
  /* prvic inicializiramo grid, zato moramo dolociti min in max koordinato */
	Grid grid1;
  grid1.init( mol1.minmax_coor() );
  
  //  mol1.init_grid_molecule(grid1, 2*PROBE);
  mol1.init_grid_molecule(&grid1, 1.1 * INTER_CHAIN_DIST); // take a bit bigger boxes
  mol1.all_triples(probe1);
  //      mol1.enumerate_clefts(probe1);
  mol1.surface_atoms(&grid1, probe1);

  /* dolocimo patche aminokislin oziroma atomov na povrsini */
  mol1.patches_cile(&grid1, probe1);

//  /* izbrisemo atome iz grida */
//  grid1->deallocate_content();
//
//  /* probe-om dolocimo sosednje atome v mol1 (ki so v kontaktu s probe-om) */
//  probe1->init_grid(grid1, mol1, EPS);
//
  
  
	probe1->free();
  //~ delete grid1;
  //~ delete mol1;
  //~ delete probe1;
}
#endif
