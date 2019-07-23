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

#include "candock/probis/probis.hpp"
#include "args.h"
#include "statchem/helper/benchmark.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/error.hpp"
#include "const.h"
#include "states.h"

using namespace std;

namespace probis {
void compare_against_bslib(const string& receptor_file,
                           const string& surface_file,
                           const string& receptor_chain_id,
                           const string& bslib_file, const int ncpu,
                           const string& nosql_file, const string& json_file) {
    try {
        statchem::Benchmark bench;
        cout << "Starting ProBiS for binding sites prediction" << endl;
        Args args;
        _longnames = true;
        _local = true;
        //~ args.fill_args(argc, argv);
        NCPU = ncpu;
        SRF_FILE = surface_file;
        PROTEIN1 = receptor_file;
        CHAIN1 = receptor_chain_id;
        //~ strcpy(SURF_FILE, bslib_file.c_str());
        SURF_FILE = bslib_file.c_str();
        NOSQL_FILE = nosql_file;
        JSON_FILE = json_file;
        args.print_state();
        args.print_modifiers();
        args.print_constants();
        state3(&args);
        _srf = true;
        PROTEIN1 = SRF_FILE;
        state4(&args);
        PROTEIN1 = receptor_file;
        state7(&args);
        cout << "Binding sites prediction took " << bench.seconds_from_start()
             << " wallclock seconds\n";

    }
    catch (Err e) {
        cout << e.what() << endl;
        throw statchem::Error("die : something went wrong in probis ...");
    }
}
};
