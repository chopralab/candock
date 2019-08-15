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

#ifndef VERTEX_H
#define VERTEX_H
#include "molib/it.hpp"
#include "helper/error.hpp"
#include "helper/debug.hpp"
#include <memory>
#include <queue>
#include <algorithm>
#include <functional>
#include <queue>
#include <set>
#include <boost/regex.hpp>

namespace Glib {
	
	template<typename T>
	class Vertex : public template_vector_container<Vertex<T>*, Vertex<T>> {
		T &__vertex;
	public:
		Vertex(const T &vertex) : __vertex(vertex) {}

		T& get_vertex() const { return __vertex; }
		friend ostream& operator<< (ostream& stream, const Vertex &v);

		bool compatible(const Vertex &other) const { 
			return __vertex.compatible(other.get_vertex()); 
		}
		const string& get_label() const { return __vertex.get_label(); }
		const string& print() const { return __vertex.get_label(); }
		const int weight() const { return __vertex.weight; }
	};
};

#endif