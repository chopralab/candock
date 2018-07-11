#ifndef ULL_SUB_STATE_H
#define ULL_SUB_STATE_H

#include <assert.h>
#include <algorithm>
#include "candock/helper/debug.hpp"

namespace candock {

namespace graph {
typedef unsigned short node_id;
static const node_id NULL_NODE = 0xFFFF;

/*----------------------------------------------------------
 * class UllSubState
 * A representation of the current search state
 * of the Ullmann's algorithm for graph-subgraph isomorphism
 ---------------------------------------------------------*/
template <class Graph1, class Graph2>
class UllSubState {
    int core_len;
    std::vector<node_id> core_1;
    std::vector<node_id> core_2;
    const Graph1& g1;
    const Graph2& g2;
    const int n1, n2;
    //~ typedef unsigned char byte;
    typedef unsigned int byte;
    std::unique_ptr<std::unique_ptr<byte[]>[]>
        M;  // Matrix encoding the compatibility of the nodes
    void __refine();

   public:
    UllSubState(const Graph1& g1, const Graph2& g2);
    UllSubState(const UllSubState& state);
    const Graph1& GetGraph1() { return g1; }
    const Graph2& GetGraph() { return g2; }
    bool NextPair(node_id& pn1, node_id& pn2, node_id prev_n1 = NULL_NODE,
                  node_id prev_n2 = NULL_NODE);
    bool IsFeasiblePair(node_id n1, node_id n2);
    void AddPair(node_id n1, node_id n2);
    bool IsGoal() { return core_len == n1; };
    bool IsDead() {
        if (n1 > n2) return true;

        for (int i = core_len; i < n1; i++) {
            bool dead = true;

            for (int j = 0; j < n2; j++)
                if (M[i][j] != 0) dead = false;

            if (dead) return true;
        }

        return false;
    };
    void BackTrack() {}
    int CoreLen() { return core_len; }
    void GetCoreSet(std::vector<node_id>& c1, std::vector<node_id>& c2);
    UllSubState* Clone();
    template <class P, class R>
    friend std::ostream& operator<<(std::ostream& stream, UllSubState<P, R>& s);
};

/*----------------------------------------------------------
 * UllSubState::UllSubState(g1, g2)
 * Constructor. Makes an empty state.
 ---------------------------------------------------------*/
template <class Graph1, class Graph2>
UllSubState<Graph1, Graph2>::UllSubState(const Graph1& ag1, const Graph2& ag2)
    : core_len(0),
      core_1(ag1.size(), NULL_NODE),
      core_2(ag2.size(), NULL_NODE),
      g1(ag1),
      g2(ag2),
      n1(ag1.size()),
      n2(ag2.size()) {
    dbgmsg(core_1.size());
    dbgmsg(core_2.size());
#ifndef NDEBUG

    for (auto& i : core_1) dbgmsg(i);

#endif
    //~ core_len=0;
    M = std::unique_ptr<std::unique_ptr<byte[]>[]>(
        new std::unique_ptr<byte[]>[n1]);

    for (int i = 0; i < n1; i++) M[i] = std::unique_ptr<byte[]>(new byte[n2]);

    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++) {
            // need to use function "get_num_edges" to get the real number of
            // edges, because g1 (or g2)
            // might be a subgraph - we want only the number of edges in the
            // subgraph (not edges
            // that extend out of the subgraph)
            M[i][j] = (g1.get_num_edges(i) <= g2.get_num_edges(j) &&
                               g1[i].compatible(g2[j])
                           ? 1
                           : 0);
            dbgmsg("M[" << i << "][" << j << "] = " << M[i][j] << " g1[" << i
                        << "].size() = " << g1[i].size() << " g2[" << j
                        << "].size() = " << g2[j].size() << " g1.get_num_edges("
                        << i << ") = " << g1.get_num_edges(i)
                        << " g2.get_num_edges(" << j << ") = "
                        << g2.get_num_edges(j) << " compatible = " << boolalpha
                        << g1[i].compatible(g2[j]));
        }
}
/*----------------------------------------------------------
 * UllSubState::UllSubState(state)
 * Copy constructor.
 ---------------------------------------------------------*/
template <class Graph1, class Graph2>
UllSubState<Graph1, Graph2>::UllSubState(
    const UllSubState<Graph1, Graph2>& state)
    : core_len(state.core_len),
      core_1(state.core_1.size()),
      core_2(state.core_2.size()),
      g1(state.g1),
      g2(state.g2),
      n1(state.n1),
      n2(state.n2) {
    copy(state.core_1.begin(), state.core_1.end(), core_1.begin());
    copy(state.core_2.begin(), state.core_2.end(), core_2.begin());
    M = std::unique_ptr<std::unique_ptr<byte[]>[]>(
        new std::unique_ptr<byte[]>[n1]);

    for (int i = 0; i < core_len; i++) M[i] = std::unique_ptr<byte[]>(nullptr);

    for (int i = core_len; i < n1; i++)
        M[i] = std::unique_ptr<byte[]>(new byte[n2]);

    for (int i = core_len; i < n1; i++)
        for (int j = 0; j < n2; j++) M[i][j] = state.M[i][j];
}
/*--------------------------------------------------------------------------
 * bool UllSubState::NextPair(pn1, pn2, prev_n1, prev_n2)
 * Puts in *pn1, *pn2 the next pair of nodes to be tried.
 * prev_n1 and prev_n2 must be the last nodes, or NULL_NODE (default)
 * to start from the first pair.
 * Returns false if no more pairs are available.
 -------------------------------------------------------------------------*/
template <class Graph1, class Graph2>
bool UllSubState<Graph1, Graph2>::NextPair(node_id& pn1, node_id& pn2,
                                           node_id prev_n1, node_id prev_n2) {
    if (prev_n1 == NULL_NODE) {
        prev_n1 = core_len;
        prev_n2 = 0;
    } else if (prev_n2 == NULL_NODE)
        prev_n2 = 0;
    else
        prev_n2++;

    if (prev_n2 >= n2) {
        prev_n1++;
        prev_n2 = 0;
    }

    if (prev_n1 != core_len) return false;

    while (prev_n2 < n2 && M[prev_n1][prev_n2] == 0) prev_n2++;

    if (prev_n2 < n2) {
        pn1 = prev_n1;
        pn2 = prev_n2;
        return true;
    } else
        return false;
}
/*---------------------------------------------------------------
 * bool UllSubState::IsFeasiblePair(node1, node2)
 * Returns true if (node1, node2) can be added to the state
 --------------------------------------------------------------*/
template <class Graph1, class Graph2>
bool UllSubState<Graph1, Graph2>::IsFeasiblePair(node_id node1, node_id node2) {
    assert(node1 < n1);
    assert(node2 < n2);
    return M[node1][node2] != 0;
}
/*--------------------------------------------------------------
 * void UllSubState::AddPair(node1, node2)
 * Adds a pair to the Core set of the state.
 * Precondition: the pair must be feasible
 -------------------------------------------------------------*/
template <class Graph1, class Graph2>
void UllSubState<Graph1, Graph2>::AddPair(node_id node1, node_id node2) {
    assert(node1 < n1);
    assert(node2 < n2);
    assert(core_len < n1);
    assert(core_len < n2);
    core_1[node1] = node2;
    core_2[node2] = node1;
    core_len++;
    int k;

    for (k = core_len; k < n1; k++) M[k][node2] = 0;

    dbgmsg("node1 " << node1 << " node2 " << node2 << " core_1.size() "
                    << core_1.size() << " core_2.size() " << core_2.size());
    __refine();
}
/*--------------------------------------------------------------
 * void UllSubState::GetCoreSet(c1, c2)
 * Reads the core set of the state into the arrays c1 and c2.
 * The i-th pair of the mapping is (c1[i], c2[i])
 --------------------------------------------------------------*/
template <class Graph1, class Graph2>
void UllSubState<Graph1, Graph2>::GetCoreSet(std::vector<node_id>& c1,
                                             std::vector<node_id>& c2) {
    for (int i = 0, j = 0; i < n1; i++)
        if (core_1[i] != NULL_NODE) {
            c1[j] = i;
            c2[j] = core_1[i];
            j++;
        }
}
/*------------------------------------------------------------
 * void UllSubState::refine()                             PRIVATE
 * Removes from the matrix M all the pairs which are not
 * compatible with the isomorphism condition
 -----------------------------------------------------------*/
template <class Graph1, class Graph2>
void UllSubState<Graph1, Graph2>::__refine() {
    dbgmsg("n1 " << n1 << " n2 " << n2 << " core_1.size() " << core_1.size()
                 << " core_2.size() " << core_2.size() << " core_len "
                 << core_len);

    for (int i = core_len; i < n1; i++)
        for (int j = 0; j < n2; j++) {
            //~ dbgmsg("i = " << i << " j = " << j);
            if (M[i][j]) {
                bool edge_ik, edge_ki, edge_jl, edge_lj;
                // The following (commented-out) for wasn't necessary
                // for(k=0; k<core_len; k++)
                int l;

                for (int k = core_len - 1; k < core_len; k++) {
                    l = core_1[k];
                    assert(l != NULL_NODE);
                    edge_ik = g1.get_conn(i, k);
                    edge_ki = g1.get_conn(k, i);
                    edge_jl = g2.get_conn(j, l);
                    edge_lj = g2.get_conn(l, j);

                    if (edge_ik != edge_jl || edge_ki != edge_lj) {
                        M[i][j] = 0;
                        break;
                    }
                    //~ else if (edge_ik  &&
                    //!g1->CompatibleEdge(g1->GetEdgeAttr(i,k),
                    //g2->GetEdgeAttr(j,l))) {
                    else if (edge_ik &&
                             !(g1[i].compatible(g2[j]) &&
                               g1[k].compatible(g2[l]))) {
                        M[i][j] = 0;
                        break;
                    }
                    //~ else if (edge_ki &&
                    //!g1->CompatibleEdge(g1->GetEdgeAttr(k,i),
                    //g2->GetEdgeAttr(l,j))) {
                    else if (edge_ik &&
                             !(g1[k].compatible(g2[l]) &&
                               g1[i].compatible(g2[j]))) {
                        M[i][j] = 0;
                        break;
                    }
                }
            }
        }
}
/*-----------------------------------------
 * Clones a state, allocating with new
 ----------------------------------------*/
template <class Graph1, class Graph2>
UllSubState<Graph1, Graph2>* UllSubState<Graph1, Graph2>::Clone() {
    return new UllSubState(*this);
}
template <class P, class R>
std::ostream& operator<<(std::ostream& stream, UllSubState<P, R>& s) {
    stream << "M = " << std::endl;

    for (int i = 0; i < s.n1; i++) {
        for (int j = 0; j < s.n2; j++) stream << s.M[i][j];

        stream << std::endl;
    }

    return stream;
}
}
}

#endif
