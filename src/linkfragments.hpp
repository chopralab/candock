#ifndef LINKFRAGMENTS_H
#define LINKFRAGMENTS_H

#include "programstep.hpp"
#include "pdbreader/molecule.hpp"
#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "pdbreader/pdbreader.hpp"

#include <mutex>

namespace Program {

	class LinkFragments : public ProgramStep
	{
		const Molib::Molecule& __receptor;
		const Molib::Score& __score;
		const OMMIface::ForceField& __ffield;
		const Molib::Atom::Grid& __gridrec;

		std::mutex __concurrent_numbering;
		int __ligand_cnt = 0;

		void __link_ligand_from_fragment( Molib::PDBreader& lpdb2, const CmdLnOpts& cmdl );
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
		
	};
}

#endif // LINKFRAGMENTSSTEP_H
