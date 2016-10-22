#ifndef GRAPH_H
#define GRAPH_H
#include "mnts.hpp"
#include "ullsubstate.hpp"
#include "pdbreader/it.hpp"
#include "helper/error.hpp"
#include "helper/smiles.hpp"
#include "helper/debug.hpp"
#include "helper/benchmark.hpp"
#include <memory>
#include <queue>
#include <algorithm>
#include <functional>
#include <queue>
#include <set>
#include <boost/regex.hpp>

namespace Glib {
	typedef vector<vector<bool>> AdjacencyMatrix;
	
	template<class Vertex>
	class Graph : public template_vector_container<Vertex*, Vertex> {
	public:
		typedef set<Vertex*> VertexSet;
		typedef vector<Vertex*> Path;
		typedef vector<Path> Cliques;
		typedef set<VertexSet> Cycles;
		typedef pair<vector<node_id>, vector<node_id>> MatchedVertices;
		typedef vector<MatchedVertices> Matches;
	private:
		vector<unique_ptr<Vertex>> __vertices; // if graph owns the vertices, they (the unique_ptrs) are stored here
		AdjacencyMatrix __conn;
		vector<int> __num_edges;
		void __expand(Vertex&, Path&, Cycles&, VertexSet&);
		template<class Vertex2>
		bool __match(vector<node_id>&, vector<node_id>&, Matches&, 
			UllSubState<Graph<Vertex>, Graph<Vertex2>>*);
		template<typename T>
		void __init(const T&, const bool);
	public:
		Graph() {}
		template<typename T>
		Graph(const T &vertices, const bool ict=false) { __init(vertices, ict); }
		Graph(vector<unique_ptr<Vertex>> vertices, const bool ict, 
				const bool) : __vertices(std::move(vertices)) { // here graph owns the vertices
			__init(__vertices, ict);
		}
		void init_conn() { 
			__conn.resize(this->size()); 
			for (auto &row : __conn) 
				row.resize(this->size(), false); 
		}
		void set_conn(node_id i, node_id j) { __conn[i][j] = true; __conn[j][i] = true; }
		bool get_conn(node_id i, node_id j) const { return __conn[i][j]; }
		int get_num_edges(const node_id i) const { return __num_edges[i]; } // real number of edges of a vertex
		string get_smiles() const;
		Cycles find_cycles_connected_graph();
		Cycles find_fused_rings();
		Cycles find_rings();
		Cliques max_weight_clique(const int);
		template<class Vertex2>
		Matches match(Graph<Vertex2>&);
		bool isomorphic(Graph &g) { Matches m = match(g); if (m.size() > 0 
			&& m[0].first.size() == g.size() && this->size() == g.size()) 
			return true; return false; }
		template<class P>
		friend ostream& operator<< (ostream& stream, const Graph<P>& g);
	};

	template<typename T>
	T intersection(T p1, T p2) {
		T v;
		auto it = set_intersection(p1.begin(), p1.end(), 
			p2.begin(), p2.end(), 
			inserter(v, v.begin()));
		return v;
	}


	template<class Vertex> 
	template<class T> // vertices is any container of unique_ptr<Vertex> or Vertex*
	void Graph<Vertex>::__init(const T &vertices, const bool ict) {
		map<const Vertex*, node_id> idx;
		for (auto &v : vertices) {
			idx[&*v] = this->size();
			this->add(&*v);
#ifndef NDEBUG
			for (auto &adj_v : *v) {
				//~ dbgmsg("pv1 = " << &*v << " (" << v << ") " << " pv2 = " << &adj_v);
				dbgmsg("pv1 = " << &*v << " pv2 = " << &adj_v);
				//~ dbgmsg("v1 = " << v->atom_number() << " v2 = " << adj_v.atom_number());
			}
#endif
		}
		if (ict) { // init connection table if requested
			init_conn();
			__num_edges.resize(this->size());
			for (auto &v : *this) {
				const node_id i1 = idx[&v];
				for (auto &adj_v : v) {
					if (idx.count(&adj_v)) {
						const node_id i2 = idx.at(&adj_v);
						dbgmsg("new2");
						dbgmsg("__init : set_conn " << i1 << " " << i2
							<< " v1 " << v << " v2 " << adj_v);
						set_conn(i1, i2);
						__num_edges[i1]++; // this counts the real number of edges
											// if this is subgraph
					}
				}
			}
		}
		dbgmsg("exiting __init");
	}

	template<class Vertex>
	typename Graph<Vertex>::Path reconstruct_path(Vertex &goal, const map<Vertex*, Vertex*> &came_from) {
		typename Graph<Vertex>::Path path;
		path.push_back(&goal);
		dbgmsg("path = ");
		dbgmsg(*path.back());
		while (came_from.find(path.back()) != came_from.end()) {
			path.push_back(came_from.at(path.back()));
			dbgmsg(*path.back());
		}
		dbgmsg("----------------------");
		return path;
	}

	template<class Vertex>
	typename Graph<Vertex>::Path find_path(Vertex &start, Vertex &goal) {
		queue<Vertex*> openset;
		openset.push(&start);
		typename Graph<Vertex>::VertexSet closedset;
		map<Vertex*, Vertex*> came_from;
		while (!openset.empty()) {
			Vertex &curr = *openset.front();
			openset.pop();
			dbgmsg(curr);
			closedset.insert(&curr);
			if (&curr == &goal) {
				dbgmsg("finished successfully " << curr << " == " << goal);
				return reconstruct_path(curr, came_from);
			}
			for (auto &adj_v : curr) {
				if (closedset.find(&adj_v) == closedset.end()) {
					came_from[&adj_v] = &curr;
					openset.push(&adj_v);
				}
			}
		}
		return typename Graph<Vertex>::Path(); // failure - return empty path			
	}

	template<class Vertex>
	template<class Vertex2>
	bool Graph<Vertex>::__match(vector<node_id> &c1, vector<node_id> &c2, 
		Matches &m, UllSubState<Graph<Vertex>, Graph<Vertex2>> *s) {
		if (s->IsGoal()) {
			int n=s->CoreLen();
			s->GetCoreSet(c1, c2);
#ifndef NDEBUG
			for (auto &v : c1) dbgmsg(v);
#endif
			m.push_back(pair<vector<node_id>, vector<node_id>>(
				vector<node_id>(c1.begin(), c1.begin() + n), 
				vector<node_id>(c2.begin(), c2.begin() + n)));
			//~ MatchedVertices mv;
			//~ for (int i = 0; i < n; ++i) {
				//~ node_id nid1 = c1[i];
				//~ node_id nid2 = c2[i];
				//~ mv[c1[i]] = c2[i];
			//~ }
			//~ m.push_back(mv);
			return false;
		}
		if (s->IsDead())
			return false;
		node_id n1=NULL_NODE, n2=NULL_NODE;
		while (s->NextPair(n1, n2, n1, n2)) {
			dbgmsg("__match : n1 = " << n1 << " n2 = " << n2);
			if (s->IsFeasiblePair(n1, n2)) {
				UllSubState<Graph<Vertex>, Graph<Vertex2>> *s1=s->Clone();
				s1->AddPair(n1, n2);
				dbgmsg("__match : n1 = " << n1 << " n2 = " << n2 
					<< " s1->IsGoal() = " << boolalpha << s1->IsGoal()
					<< " c1.size() = " << c1.size()
					<< " c2.size() = " << c2.size());
				if (__match(c1, c2, m, s1)) {
					s1->BackTrack();
					delete s1;
					return true;
				}
				else {
					s1->BackTrack();
					delete s1;
				}
			}
		}
		return false;
	}
	
	template<class Vertex>
	template<class Vertex2>
	typename Graph<Vertex>::Matches Graph<Vertex>::match(Graph<Vertex2> &other) {
		UllSubState<Graph<Vertex>, Graph<Vertex2>> s0(*this, other);
		Graph<Vertex> &g1=s0.GetGraph1();
		Graph<Vertex2> &g2=s0.GetGraph();
		/* Choose a conservative dimension for the arrays */
		int n;
		if (g1.size()<g2.size())
			n=g2.size();
		else
			n=g1.size();
		dbgmsg(n);
		dbgmsg(g1);
		dbgmsg(g2);
		vector<node_id> c1(n), c2(n);
		Matches m;
		__match(c1, c2, m, &s0);
		return m;
	}
	
	template<class Vertex>
	typename Graph<Vertex>::Cliques Graph<Vertex>::max_weight_clique(const int iter) {
		Benchmark::reset();
		unique_ptr<int[]> weight(new int[this->size()]);
		for (int i = 0; i < this->size(); ++i) weight[i] = this->element(i).weight();
		vector<vector<int>> qmax;
		MNTS m(qmax, __conn, weight.get(), iter);
		Cliques clique;
		for (auto &rows : qmax) {
			dbgmsg("found max weight clique of " << to_string(rows.size()) 
				<< " vertices");
			clique.push_back(Path());
			for (auto &vnum : rows) {
				clique.back().push_back(&this->operator[](vnum));
				dbgmsg("clique push vertex = " << vnum);
			}
		}
		cout << "time to find max.weight clique " << Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		return clique;
	}

	//~ template<class Vertex>
	//~ typename Graph<Vertex>::Cliques Graph<Vertex>::max_clique() {
		//~ Benchmark::reset();
		//~ vector<vector<int>> qmax;
		//~ Maxclique m(__conn, this->size());
		//~ int qmax_sz;
		//~ m.mcq(qmax, qmax_sz);
		//~ Cliques clique;
		//~ for (auto &rows : qmax) {
			//~ dbgmsg("found max weight clique of " << to_string(rows.size()) 
				//~ << " vertices");
			//~ clique.push_back(Path());
			//~ for (auto &vnum : rows) {
				//~ clique.back().push_back(&this->operator[](vnum));
				//~ dbgmsg("clique push vertex = " << vnum);
			//~ }
		//~ }
		//~ cout << "time to find max.weight clique " << Benchmark::seconds_from_start() << " wallclock seconds" << endl;
		//~ return clique;
	//~ }

	template<class Vertex>
	typename Graph<Vertex>::Cycles Graph<Vertex>::find_fused_rings() {
		Cycles fused;
		Cycles cycles = find_cycles_connected_graph();
		bool mergeable = true;
		// merge cycles until no more can be merged
		while (mergeable) {
			mergeable = false;
			for (typename Cycles::iterator it = cycles.begin(); it != cycles.end(); ++it) {
				const VertexSet &first = *it;
				VertexSet current(first.begin(), first.end());
				typename Cycles::iterator it2 = it;
				for (++it2; it2 != cycles.end();) {
					const VertexSet &second = *it2;
					VertexSet inter;
					set_intersection(first.begin(), first.end(), 
						second.begin(), second.end(), inserter(inter, inter.begin()));
					// merge first with second if > 1 vertices in common
					if (inter.size() > 1) {
						mergeable = true;
						current.insert(second.begin(), second.end());
						it2 = cycles.erase(it2);
					}
					else {
						++it2;
					}
				}
				fused.insert(current);
			}
			if (mergeable) {
				cycles = fused;
				fused.clear();
			}
		}
#ifndef NDEBUG
		dbgmsg("+++++++++++++FUSED RINGS++++++++++++++++++");
		for (auto &cycle : fused) {
			dbgmsg("new fused ring : ");
			for (auto &pv : cycle)
				//~ dbgmsg(pv->print());
				//~ dbgmsg(*pv);
				dbgmsg(pv->get_label());
		}
#endif
		return fused;
	}

	template<class Vertex>
	typename Graph<Vertex>::Cycles Graph<Vertex>::find_rings() {
		Cycles rings;
		Cycles cycles = find_cycles_connected_graph();
		vector<VertexSet> v(cycles.begin(), cycles.end());
		sort(v.begin(), v.end(), [](const VertexSet &i, const VertexSet &j) { 
			return i.size() < j.size(); });
#ifndef NDEBUG
		for (typename vector<VertexSet>::iterator it = v.begin(); it != v.end(); it++) {
			dbgmsg("size = " << it->size());
			for (auto it2 = it->begin(); it2 != it->end(); it2++)
				//~ dbgmsg((*it2)->print() << " ");
				//~ dbgmsg(**it2);
				dbgmsg((*it2)->get_label() << " ");
		}
#endif
		VertexSet mapped;
		// go over cycles sorted by increasing size
		for (typename vector<VertexSet>::iterator it = v.begin(); it != v.end(); it++) {
			VertexSet &cycle = *it;
			VertexSet inter;
			// find cycles that are not made out of smaller cycles
			set_intersection(cycle.begin(), cycle.end(), mapped.begin(), 
						mapped.end(), inserter(inter, inter.begin()));
			if (inter.size() == cycle.size())
				continue;
			mapped.insert(cycle.begin(), cycle.end());
			rings.insert(cycle);
		}
		return rings;
	}
	
	template<class Vertex>
	void Graph<Vertex>::__expand(Vertex &v, Path &path, Cycles &cycles, VertexSet &visited) {
		path.push_back(&v);
		visited.insert(&v);
		dbgmsg("new vertex");
		for (auto &adj_v : v) {
			if (!visited.count(&adj_v)) {
				__expand(adj_v, path, cycles, visited);
			}
			else {
				typename Path::iterator it = find(path.begin(), path.end(), &adj_v);
				VertexSet cycle(it, path.end());
				if (cycle.size() > 2) {
					cycles.insert(cycle);
				}
			}
		}
		visited.erase(&v);
		path.pop_back();
	}
	
	template<class Vertex>
	typename Graph<Vertex>::Cycles Graph<Vertex>::find_cycles_connected_graph() {
		Path path;
		Cycles cycles;
		VertexSet visited;
		dbgmsg("expanding the vertices");
		__expand(this->first(), path, cycles, visited);
		dbgmsg("after expanding the vertices");
		dbgmsg("CYCLES FOUND : ");
		for (auto &cycle : cycles) {
			dbgmsg("CYCLE : ");
			for (auto &pvertex : cycle) {
				dbgmsg(*pvertex);
			}
		}

		//~ dbgmsg("CYCLES FOUND : " << cycles);
		return cycles;
	}
	
	template<class Vertex>
	string Graph<Vertex>::get_smiles() const {
		stringstream ss;
		map<Vertex*, int> idx;
		for (auto &v : *this)
			ss << v.get_label() << " ";
		return ss.str();
	}

	template<class P>
	ostream& operator<< (ostream& stream, const Graph<P>& g) {
		map<const P*, int> idx;
		int i = 0;
		for (auto &vertex : g) idx[&vertex] = i++;
		for (auto &vertex : g) {
			stream << "vertex label = " << vertex.get_label() << " index = " << idx.at(&vertex) << " edges = [";
			for (auto &adj_v : vertex) {
				if (idx.count(&adj_v)) {
					stream << "{label = " << adj_v.get_label() << " index = " << idx.at(&adj_v) << "} ";
				}
			}
			stream << "]" << endl;
		}
		if (!g.__conn.empty()) {
			int edge_size = 0;
			for (auto &v1 : g) {
				for (auto &v2 : g) {
					//~ if (g.__conn->size() > 0 && (*g.__conn)[idx[&v1]][idx[&v2]] == true) edge_size++;
					if (g.__conn[idx[&v1]][idx[&v2]] == true) edge_size++;
				}
			}
			stream << "p " << g.size() << " " << edge_size << endl;
			for (int i = 0; i < g.size(); i++) {
				for (int j = i + 1; j <g.size(); j++) {
					//~ if (g.__conn->size() > 0 && (*g.__conn)[i][j] == true) stream << "e " << i + 1 << " " << j + 1 << endl;
					if (g.__conn[i][j] == true) stream << "e " << i + 1 << " " << j + 1 << endl;
				}
			}
		}
		for (auto &vertex : g) {
			stream << "w " << idx[&vertex] << " " << vertex.weight() << endl;
		}
		return stream;
	}
};
#endif
