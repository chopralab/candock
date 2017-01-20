#ifndef LINKFRAGMENTS_H
#define LINKFRAGMENTS_H

#include "programstep.hpp"
#include "pdbreader/molecule.hpp"
#include "pdbreader/molecules.hpp"
#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "pdbreader/pdbreader.hpp"

#include <mutex>

namespace Program {

	class LinkFragments : public ProgramStep
	{
		Molib::Molecules __all_top_poses;
		
		const Molib::Molecule& __receptor;
		const Molib::Score& __score;
		const OMMIface::ForceField& __ffield;
		const Molib::Atom::Grid& __gridrec;

		std::mutex __concurrent_numbering;
		int __ligand_cnt = 0;

		void __link_ligand ( Molib::Molecule& ligand, const CmdLnOpts& cmdl, const OMMIface::ForceField& ffield);
	protected:
		virtual bool __can_read_from_files ( const CmdLnOpts& cmdl );
		virtual void __read_from_files ( const CmdLnOpts& cmdl );
		virtual void __continue_from_prev ( const CmdLnOpts& cmdl );
	public:
		LinkFragments ( const Molib::Molecule& receptor,
						const Molib::Score& score,
						const OMMIface::ForceField& ffield,
						const Molib::Atom::Grid& gridrec ) : 
						__receptor(receptor), __score(score),
						__ffield(ffield), __gridrec(gridrec) {};

		void link_ligands (const Molib::Molecules& ligands, const CmdLnOpts &cmdl);
		const Molib::Molecules& top_poses() const { return __all_top_poses; }
		void clear_top_poses() { __all_top_poses.clear(); };
	};
}

#endif // LINKFRAGMENTSSTEP_H
