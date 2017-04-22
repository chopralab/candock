#ifndef OPTICS_H
#define OPTICS_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
//~ #include <unordered_map>
#include <iomanip>
#include <queue>
#include <cmath>
#include <algorithm>
#include <helper/debug.hpp>
#include <helper/benchmark.hpp>
#include <assert.h>
using namespace std;

/* OPTICS clustering algorithm C++ implementation
 * 
 * requires two parameters: epsilon (the maximum radius to consider) 
 * and min_pts (the number of points require to form a cluster)
 */

namespace cluster {
	//~ template <typename T> using PairwiseDistances = unordered_map<T*, unordered_map<T*, double>>;
	template <typename T> using PairwiseDistances = map<T*, map<T*, double>>;
	template <typename T> using Neighbors = map<T*, vector<T*>>;
	template <typename T> using MapD = map<T*, double>;
	template <typename T> using MapB = map<T*, bool>;
	template <typename T> using MapI = map<T*, int>;
	template <typename T> using Clusters = multimap<int, T*>;
	template <typename T> using Vector = vector<T*>;
		
	template<typename T, typename CompareDistanceLess = std::less<double>, typename CompareScoreLess = std::less<double>>
	class Optics {

		MapD<T> __reachability_distance;
		MapD<T> __core_distance;
		PairwiseDistances<T> &__distance; // distance from center object
		MapB<T> __processed;
		MapI<T> __cluster_id;
		Neighbors<T> __neighbors;
		Vector<T> __ordered_file;
		const MapD<T> &__scores;
		//~ static const MapD<T> &__scores;
		//~ MapD<T> &__scores;
		double __eps;
		size_t __min_pts;

		double __JANEZ_HV;

		void __set_core_distance(T &p, size_t min_pts) {
			if (__neighbors[&p].size() < min_pts) {
				__core_distance[&p] = __JANEZ_HV;
			} else {
				__core_distance[&p] = __distance[&p][__neighbors[&p][min_pts-1]];
				//~ dbgmsg("core_distance=" << __core_distance[&p]);
			}
		}

		void __expand_cluster_order(T &p) {
#ifndef NDEBUG
			assert(__processed[&p] == false);
#endif
			__processed[&p] = true;
			__reachability_distance[&p] = __JANEZ_HV; // reachability distance is already __JANEZ_HV
			__set_core_distance(p, __min_pts);
			__ordered_file.push_back(&p);
			Vector<T> order_seeds; // priority_queue sorted by reachability_distance from the core object (last element has the smallest reachability_distance)
			if (__core_distance[&p] != __JANEZ_HV) { // is p a core object?
				__update(order_seeds, p);
//~ #ifndef NDEBUG
				//~ dbgmsg("order_seeds = ");
				//~ for (auto &pp : order_seeds) 
					//~ dbgmsg("\t" << *pp << "\tprocessed=" << __processed[pp] 
						//~ << "\treachability_distance="<< __reachability_distance[pp]);
//~ #endif
				while (!order_seeds.empty()) {
					//~ pop_heap(order_seeds.begin(), order_seeds.end(), by_r_dist); // new
					pop_heap(order_seeds.begin(), order_seeds.end(), 
						[this] (T *i, T *j) { 
							//~ return __reachability_distance[i] > __reachability_distance[j]; 
							//~ dbgmsg("i = " << boolalpha << (i == nullptr) << " " << *i);
							//~ dbgmsg("j = " << boolalpha << (j == nullptr) << " " << *j);

							//~ dbgmsg("reachability distance = " << __reachability_distance.at(i) 
								//~ << " " << __reachability_distance.at(j));

							return !CompareDistanceLess()(__reachability_distance[i], __reachability_distance[j]); 
						});
					//~ const T &q = *(order_seeds.back()); // top of the heap
					T &q = *(order_seeds.back()); // top of the heap
//~ #ifndef NDEBUG
					//~ dbgmsg("order_seeds (before pop_heap)= ");
					//~ for (auto &pp : order_seeds) 
						//~ dbgmsg("\t" << *pp << "\tprocessed=" << __processed[pp] 
							//~ << "\treachability_distance="<< __reachability_distance[pp]);
//~ #endif
					order_seeds.pop_back(); // pop from heap
//~ #ifndef NDEBUG
					//~ dbgmsg("order_seeds (after pop_heap)= ");
					//~ for (auto &pp : order_seeds) 
						//~ dbgmsg("\t" << *pp << "\tprocessed=" << __processed[pp] 
							//~ << "\treachability_distance="<< __reachability_distance[pp]);
//~ #endif
//~ #ifndef NDEBUG
					//~ assert(__processed[&q] == false);
//~ #endif
					__processed[&q] = true;
					__set_core_distance(q, __min_pts);
					__ordered_file.push_back(&q);
					if (__core_distance[&q] != __JANEZ_HV) { // is q a core object?
						__update(order_seeds, q); // if yes, insert further candidates into the order_seeds
					}
				}			
			}
		}
		void __update(Vector<T> &order_seeds, T &p) { // p is center object
			double c_dist = __core_distance[&p];
			for (auto &po : __neighbors[&p]) {
				T &o = *po;
				if (!__processed[&o]) {
					//~ double new_r_dist = max(c_dist, __distance[&p][&o]); // o.distance() is the distance between o and p
					double new_r_dist = max(c_dist, __distance[&p][&o], 
						CompareDistanceLess()); // o.distance() is the distance between o and p
//~ #ifndef NDEBUG
					//~ dbgmsg("c_dist =" << c_dist << "\to.distance() =" << __distance[&p][&o]);
//~ #endif
					if (__reachability_distance[&o] == __JANEZ_HV) {
						__reachability_distance[&o] = new_r_dist;
						order_seeds.push_back(&o);
	
						//~ push_heap(order_seeds.begin(), order_seeds.end(), by_r_dist); // new
						push_heap(order_seeds.begin(), order_seeds.end(),
							[this] (T *i, T *j) { 
								//~ return __reachability_distance[i] > __reachability_distance[j]; 
								//~ dbgmsg("i = " << boolalpha << (i == nullptr) << " " << *i);
								//~ dbgmsg("j = " << boolalpha << (j == nullptr) << " " << *j);

								//~ dbgmsg("reachability distance = " << __reachability_distance.at(i) 
									//~ << " " << __reachability_distance.at(j));

								return !CompareDistanceLess()(__reachability_distance[i], __reachability_distance[j]); 
							});

#ifndef NDEBUG
						assert(__processed[&o] == false);
#endif
					}
					else { // object o already in order_seeds
						//~ if (new_r_dist < __reachability_distance[&o]) {
						if (CompareDistanceLess()(new_r_dist, __reachability_distance[&o])) {
//~ #ifndef NDEBUG
							//~ dbgmsg("order_seeds (before resorting)= ");
							//~ for (auto &pp : order_seeds) 
								//~ dbgmsg("\t" << *pp << "\tprocessed=" << __processed[pp] 
									//~ << "\treachability_distance="<< __reachability_distance[pp]);
//~ #endif
							__reachability_distance[&o] = new_r_dist; // o is pointed to from the order_seeds
//~ #ifndef NDEBUG
							//~ dbgmsg("order_seeds (before resorting2)= ");
							//~ for (auto &pp : order_seeds) 
								//~ dbgmsg("\t" << *pp << "\tprocessed=" << __processed[pp] 
									//~ << "\treachability_distance="<< __reachability_distance[pp]);
//~ #endif
							//~ make_heap(order_seeds.begin(), order_seeds.end(), by_r_dist); // new
							make_heap(order_seeds.begin(), order_seeds.end(),
								[this] (T *i, T *j) { 
									//~ return __reachability_distance[i] > __reachability_distance[j]; 

									//~ dbgmsg("reachability distance = " << __reachability_distance.at(i) 
										//~ << " " << __reachability_distance.at(j));

									return !CompareDistanceLess()(__reachability_distance[i], __reachability_distance[j]); 
								});

							//~ sort(order_seeds.begin(), order_seeds.end(), by_r_dist);  // old
//~ #ifndef NDEBUG
							//~ dbgmsg("order_seeds (after resorting)= ");
							//~ for (auto &pp : order_seeds) 
								//~ dbgmsg("\t" << *pp << "\tprocessed=" << __processed[pp]
									//~ << "\treachability_distance="<< __reachability_distance[pp]);
							//~ assert(__processed[&o] == false);
//~ #endif
						}
					}
				}
			}
		}
	public:
		Optics(PairwiseDistances<T> &pairwise_distances, const MapD<T> &scores, 
			double eps, size_t min_pts) : __distance(pairwise_distances), 
			__scores(scores), __eps(eps), __min_pts(min_pts) {

			Benchmark bench;
			cout << "starting clustering ..." << endl;

			__JANEZ_HV = (CompareDistanceLess()(-1.0, 1.0) ? HUGE_VAL : -HUGE_VAL);
			// make inverse pairs of pairwise distances if missing
			//~ for (auto &kv1 : pairwise_distances) {
			for (auto &kv1 : __distance) {
				auto &p1 = kv1.first;
				//~ for (auto &kv2 : pairwise_distances[p1]) {
				for (auto &kv2 : __distance[p1]) {
					auto &p2 = kv2.first;
					auto &distance = kv2.second;
					//~ pairwise_distances[p2][p1] = distance;
					__distance[p2][p1] = distance;
					//~ dbgmsg(*p1 << " " << *p2 << " " << pairwise_distances[p2][p1]);
					//~ dbgmsg(*p1 << " " << *p2);
				}
			}
			//~ throw Error("the end");
			//~ __distance = pairwise_distances;
			for (auto &kv1 : __distance) {
				auto &p1 = kv1.first;
				__reachability_distance[p1] = __core_distance[p1] = __JANEZ_HV;
				__cluster_id[p1] = 0;
				__processed[p1] = false;
				for (auto &kv2 : __distance[p1]) {
					auto &p2 = kv2.first;
					auto &distance = kv2.second;
					//~ if (distance < __eps) {
					if (CompareDistanceLess()(distance, __eps)) {
						__neighbors[p1].push_back(p2);
					}
				}
				sort(__neighbors[p1].begin(), __neighbors[p1].end(), 
					[&p1, this] (T *i, T *j) 
					//~ { return __distance[p1][i] < __distance[p1][j]; });
					{ return CompareDistanceLess()(__distance[p1][i], __distance[p1][j]); });
//~ #ifndef NDEBUG
				//~ dbgmsg("neighbors of " << *p1 << " are ");
				//~ for (auto &pp : __neighbors[p1])
					//~ dbgmsg("\t" << *pp);
//~ #endif
			}
			for (auto &kv : __distance) {
				auto &p = kv.first;
				if(!__processed[p])
					__expand_cluster_order(*p);
			}
			cout << "total time required for clustering was " 
				<< bench.seconds_from_start() << " wallclock seconds\n";
		}

		pair<Clusters<T>, Clusters<T>> extract_dbscan(const double eps_cur,
			const size_t max_num_clus=9999999, const bool do_sort=false) { // returns the original objects together with their cluster_ids

			assert(CompareDistanceLess()(eps_cur, __eps));  // eps_cur must be < eps
			Clusters<T> clustered_points, representative_points;
			int cluster_id = 0;
			
			for (auto &pp : __ordered_file) {
				if (!CompareDistanceLess()(__reachability_distance[pp], eps_cur)) {
					if (CompareDistanceLess()(__core_distance[pp], eps_cur)) {
						__cluster_id[pp] = ++cluster_id;
					}
					else {
						__cluster_id[pp] = -1; // noise
					}
				}
				else { 
					__cluster_id[pp] = cluster_id;
				}
			}
			
			// assign all yet unclustered points to clusters
			for (auto &pp : __ordered_file) {
				if (__cluster_id[pp] == -1) {
					// pass 1 : try to assign a noise point to a nearest cluster
					Vector<T> clustered_neighbors;
					for (auto &ppp : __neighbors[pp])
						if (__cluster_id[ppp] != -1)
							clustered_neighbors.push_back(ppp);
					if (!clustered_neighbors.empty()) {
						T *nearest = *min_element(clustered_neighbors.begin(), clustered_neighbors.end(),
							[&pp, this] (T *i, T *j) { return 
							CompareDistanceLess()(__distance[pp][i], __distance[pp][j]); });
						if (CompareDistanceLess()(__distance[pp][nearest], eps_cur)) {
							__cluster_id[pp] = __cluster_id[nearest];
							dbgmsg("pass 1 : assigned point " << *pp << " to cluster "
								<< __cluster_id[pp] << " because it was close (distance=" << 
								__distance[pp][nearest] << ") to point " << *nearest);
							continue;
						}
					}
					// pass 2 : find points that are close and assign them to a new cluster
					Vector<T> close_neighbors;
					close_neighbors.push_back(pp);
					++cluster_id;
					for (auto &ppp : __neighbors[pp])
						if (CompareDistanceLess()(__distance[pp][ppp], eps_cur))
							close_neighbors.push_back(ppp);
					for (auto &n : close_neighbors) {
						__cluster_id[n] = cluster_id;
						dbgmsg("pass 2 : assigned point " << *n << " to a new cluster "
							<< __cluster_id[n]);
					}
				}
			}
			if (do_sort) {
#ifndef NDEBUG
				for (auto &pp : __ordered_file) {
					dbgmsg("before sorting cluster_id = " << __cluster_id[pp] 
						<< " score = " << __scores.at(pp));
					dbgmsg("pp = " << *pp);
				}
#endif
				dbgmsg("sorting3 was requested");
				// sort according to scores
				for (size_t i = 0; i <__ordered_file.size(); ++i) {
					for (size_t j = i + 1; j <__ordered_file.size(); ++j) {
						//~ dbgmsg("before2 sorting : " << *__ordered_file[i] << " " << *__ordered_file[j]);
						if (CompareScoreLess()(__scores.at(__ordered_file[i]), __scores.at(__ordered_file[j]))) {
							T *tmp = __ordered_file[i];
							__ordered_file[i] = __ordered_file[j];
							__ordered_file[j] = tmp;
						}
					}
				}
				dbgmsg("immediately after sorting");
#ifndef NDEBUG
				for (auto &pp : __ordered_file) {
					dbgmsg("after sorting cluster_id = " << __cluster_id[pp] 
						<< " score = " << __scores.at(pp));
				}
#endif
			}
			// renumber clusters so that the best cluster is first
			map<int, int> added_clus;
			int i = 0;
			for (auto &pp : __ordered_file) {
				const int cluster_id = __cluster_id[pp];
				if (added_clus.size() < max_num_clus || added_clus.count(cluster_id)) {
					if (!added_clus.count(cluster_id))
						added_clus.insert({cluster_id, ++i});
					clustered_points.insert({added_clus[cluster_id], pp});
				}
			}
			
			// find cluster representatives
			dbgmsg("representatives : ");
			for(typename Clusters<T>::iterator it = clustered_points.begin(), end = clustered_points.end();
				it != end; it = clustered_points.upper_bound(it->first)) {
				
				auto range = clustered_points.equal_range(it->first);
				T *rep = max_element(range.first, range.second, [this] (
					const pair<const int, T*> &i, const pair<const int, T*> &j) { 
					return CompareScoreLess()(__scores.at(i.second), __scores.at(j.second)); })->second;
				representative_points.insert({it->first, rep});
				dbgmsg("\t[" << it->first << "]\t" << *rep << "\tscore="
					<< __scores.at(rep));
			}
#ifndef NDEBUG
			dbgmsg("clusters : ");
			for (auto pp : __ordered_file) {
				dbgmsg("\t[" << __cluster_id[pp] << "]\t" << *pp << "\tcore_distance=" 
					<< __core_distance[pp] << "\treachability_distance=" << __reachability_distance[pp]);
			}
#endif
			return {clustered_points, representative_points};
		}
	};
};
#endif
