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

#include "motif.h"

bool Motif::is_in_motif(ChResi c) {
    /*
      Vrnemo true, ce se kombinacija chain id-ja in resi-ja pojavlja v motivu -
      "sele" setu.
    */

    return (sele.find(c) == sele.end() ? false : true);
}

string Motif::debracket(char str[]) {
    /*
      Vrnemo vsebino, zaobjeto v oglatih oklepajih, npr [:A or :B] --> :A or :B
    */
    string m = str;
    size_t found = m.find_first_of("[");
    return m.substr(found + 1, m.find_first_of("]", found + 1) - found - 1);
}

string Motif::delete_all_whitespaces(string str) {
    /*
       Zbrisemo blank znake kjerkoli v stringu
    */
    string whitespaces(" \t\f\v\n\r");
    size_t found = str.find_first_of(whitespaces);
    size_t first = 0;
    string new_str;
    // od zadaj
    while (found != string::npos) {
        new_str += str.substr(first, found - first);
        first = found + 1;
        found = str.find_first_of(whitespaces, found + 1);
    }

    new_str += str.substr(first, str.length() - first);
    return new_str;
}

pair<size_t, size_t> Motif::ffo(string str, set<string> n, size_t pos) {
    /*
      Poiscemo v stringu str tisti string v n, ki je najbolj na zacetku str-ja,
      ampak za pos.
      Vrnemo kazalec na najdeni string n v str-ju in pa dolzino stringa n.
      Na primer str=":Aand(32,33)" n="and|(|)|," or  ----> return value = <2,3>
      (nasel je separator
      "and" na poziciji 2, "and" pa je dolzine 3).
      Ce ne najde separotorja, potem vrne <string::npos,0>.
    */
    size_t pos1 = str.length();
    size_t lsep = 0;
    for (set<string>::iterator it = n.begin(); it != n.end(); it++) {
        size_t found = str.find(*it, pos);
        if (found != string::npos) {
            if (found < pos1) {
                pos1 = found;
                lsep = it->length();
            }
        }
    }
    return make_pair((pos1 == str.length() ? string::npos : pos1), lsep);
}

set<string> Motif::divide(string str, string needle) {
    /*
      Razdelimo str na chunk-e. V chunk vrnemo prvi del str-ja, ki je za indexom
      pos in se konca z enim izmed separatorjev
      navedenim v needle stringu. Separator stringi so loceni z |.
    */
    // zapisemo vse separatorje v set
    set<string> n;
    size_t ff = 0;
    size_t lf = needle.find_first_of("|");
    //  cout << lf << endl;
    n.insert(needle.substr(0, (lf == string::npos ? needle.length() : lf)));
    //  cout << needle.substr(0, (lf == string::npos ? needle.length() : lf)) <<
    //  endl;
    while (lf != string::npos) {
        ff = lf;
        lf = needle.find_first_of("|", ff + 1);
        n.insert(needle.substr(
            ff + 1, (lf == string::npos ? needle.length() : lf - ff - 1)));
        //    cout << needle.substr(ff + 1, (lf == string::npos ?
        //    needle.length() : lf - ff - 1)) << endl;
    }

    set<string> chunks;
    string chunk = "";
    // zapisemo vse chunke v set chunks
    size_t cf = 0;
    pair<size_t, size_t> r = ffo(str, n);
    while (r.first != string::npos) {
        chunk = str.substr(cf, r.first - cf);

        if (chunk != "") {
            //      cout << "chunk=[" << chunk << "]" << endl;
            chunks.insert(chunk);
        }
        cf = r.first + r.second;
        r = ffo(str, n, cf);
    }

    chunk = str.substr(cf, str.length() - cf);
    if (chunk != "") {
        //    cout << "chunk=[" << chunk << "]" << endl;
        chunks.insert(chunk);
    }

    return chunks;
}

Motif::Motif(string motif) {
    /*
      Interpretiramo zapis oblike :A and (14,57,69-71) or :B and (33,34,50) in
      zapisemo
      chain id-je in resi-je izbranih aminokislin v set.
    */
    sele.clear();
    string str = delete_all_whitespaces(motif);
    //  cout << "no spaces=" << str << endl;

    // najprej razdelimo string na odseke med or-i ali vejicami
    set<string> chunk1 = divide(str, "or");
    for (set<string>::iterator it1 = chunk1.begin(); it1 != chunk1.end();
         it1++) {
        //    cout << "chunk=" << *it1 << endl;
        // nato beremo odseke med and-i
        set<string> chunk2 = divide(*it1, "and|&");
        char chain_id = ' ';  // default
        set<int> resi;
        resi.clear();
        for (set<string>::iterator it2 = chunk2.begin(); it2 != chunk2.end();
             it2++) {
            size_t cp;
            // imamo chain_id
            if ((cp = it2->find(":")) != string::npos) {
                chain_id = it2->at(cp + 1);
                //        cout << "ChainID=" << chain_id << endl;
            } else {
                // ali pa residue id-je, locene z vejicami ali or-i ali oklepaji
                set<string> chunk3 = divide(*it2, "or|,|(|)");
                for (set<string>::iterator it3 = chunk3.begin();
                     it3 != chunk3.end(); it3++) {
                    // ce imamo array resi-jev, npr. 34-37
                    int fi = 0, li = 0, i = 0;

                    set<string> tmp_chunk4 = divide(
                        *it3, "-");  // stringe moramo sortirati kot int !!!
                    set<string, cmp> chunk4(tmp_chunk4.begin(),
                                            tmp_chunk4.end());
                    for (set<string, cmp>::iterator it4 = chunk4.begin();
                         it4 != chunk4.end(); it4++) {
                        istringstream ins;
                        ins.str(*it4);
                        if (i == 0)
                            ins >> fi;
                        else if (i == 1)
                            ins >> li;
                        else {
                            throw Err(
                                "Error (MOTIF) : Error in motif selection "
                                "string (wrong residue range)!",
                                6);
                        }
                        i++;
                    }
                    if (i == 1) li = fi;
                    //          cout << "li="<<li<<" fi="<<fi<<endl;
                    // napolnimo set chain id-jev in resi id-jev
                    for (int j = fi; j <= li; j++) {
                        resi.insert(j);
                    }
                }
            }
        }
        // zapisemo selekcijo v "sele"
        for (set<int>::iterator rit = resi.begin(); rit != resi.end(); rit++) {
            sele.insert(ChResi(chain_id, *rit));
        }
    }
}
