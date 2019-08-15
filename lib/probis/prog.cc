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
