#ifndef VERTEX_H
#define VERTEX_H
#include <algorithm>
#include <functional>
#include <memory>
#include <queue>
#include <queue>
#include <set>
#include "candock/helper/debug.hpp"
#include "candock/helper/error.hpp"
#include "candock/molib/it.hpp"

namespace candock {

namespace graph {

template <typename T>
class Vertex : public template_vector_container<Vertex<T>*, Vertex<T>> {
    T& __vertex;

   public:
    Vertex(const T& vertex) : __vertex(vertex) {}

    T& get_vertex() const { return __vertex; }
    friend ostream& operator<<(ostream& stream, const Vertex& v);

    bool compatible(const Vertex& other) const {
        return __vertex.compatible(other.get_vertex());
    }
    const string& get_label() const { return __vertex.get_label(); }
    const string& print() const { return __vertex.get_label(); }
    const int weight() const { return __vertex.weight; }
};
};
}

#endif
