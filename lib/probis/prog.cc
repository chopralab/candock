#include "const.h"
#include "args.h"
#include "states.h"


int main_probis(int argc, char *argv[]) {
  try {
    Args args;
    args.fill_args(argc, argv);
    args.print_state();
    args.print_modifiers();
    args.print_constants();
    switch (_state) {
    case TOTAL  : state0(PROTEIN2, CHAIN2);       break;
//    case PPI    : state1(args);       break;
//    case BIND   : state2(args);       break;
    case WRITE  : state3(&args);       break;
    case SURFDB : state4(&args);       break;
    case MARK   : state5(&args);       break;
    case LIGAND : state6(&args);       break;
    case RESULTS: state7(&args);       break;
      //  case SEQ    : state8(args);       break;
    case ALIGN  : state9(&args);       break;
      //  case BIOU   : state10(args);       break;
    case HEAD   : state11(&args);       break;
    case ALLBIO : state12(&args);       break;
#ifdef CILE
    case PATCH  : stateX(&args);       break;
#endif
    }
  }  // .. konec try bloka
  catch (Err e) {
    cout << e.what() << endl;
  }

  return 0;
}
