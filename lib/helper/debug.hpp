/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include <iostream>
#include <iomanip>

#ifdef NDEBUG

#define dbgprint(x)
#define dbgmsg(message)

#else

#define dbgprint(x) cerr << #x << ": " << x << " in "  << __FILE__ << ":" << __LINE__ << endl
#define dbgmsg(message) cerr << message << " in "  << __FILE__ << ":" << __LINE__ <<endl
#endif
