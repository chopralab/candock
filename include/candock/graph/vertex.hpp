#ifndef VERTEX_H
#define VERTEX_H
#include "candock/molib/it.hpp"
#include "candock/helper/error.hpp"
#include "candock/helper/debug.hpp"
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