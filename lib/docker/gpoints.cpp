#include "helper/inout.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "geom3d/geom3d.hpp"
#include "centro/centroids.hpp"
#include "score/score.hpp"
#include "gpoints.hpp"
#include <iostream>
#include <exception>
using namespace std;

namespace Docker {

	ostream& operator<<(ostream& os, const Docker::Gpoints::GpointVec &points)	{
		for (auto &point : points) {
			// to add : output of energies
			os << "ATOM      1   U  DIK     1    " << point.crd().pdb() << endl; 
			//~ os << "ATOM      1   U  DIK     1    " << point.crd().pdb() 
				//~ setw(6) << 1.0 << setw(6) << setprecision(2) << fixed
				//~ << point.energy(<< endl;
		}
		return os;
	}	


	ostream& operator<<(ostream& os, const Gpoints &gpoints)	{
		for (auto &kv : gpoints.__gridpoints) {
			const int bsite_id = kv.first;
			os << "MODEL" << setw(9) << right << bsite_id << endl;
			os << kv.second;
			os << "ENDMDL" << endl;
		}
		return os;
	}	
	
	Gpoints::Gpoints(const double &grid_spacing, const double &radial_check) 
		: __score(nullptr), __ligand_idatm_types(nullptr) {
		try {
			__identify_gridpoints(grid_spacing, radial_check);
		} catch(...) {
			dbgmsg("FAILURE: constructor of Gpoints failed ... cleaning up resources...");
			throw;
		}
	}
	Gpoints::Gpoints(const Molib::Score &score, const set<int> &ligand_idatm_types, 
		const Centro::Centroids &centroids, Molib::Atom::Grid &grid, const double &grid_spacing, 
		const int &dist_cutoff, const double &excluded_radius, const double &max_interatomic_distance)
		: __score(&score), __ligand_idatm_types(&ligand_idatm_types) {
		try {
			__identify_gridpoints(centroids, grid, grid_spacing, dist_cutoff, 
				excluded_radius, max_interatomic_distance);
		} catch(...) {
			dbgmsg("FAILURE: constructor of Gpoints failed ... cleaning up resources...");
			throw;
		}
	}

	Gpoints::Gpoints(const Centro::Centroids &centroids, Molib::Atom::Grid &grid, const double &grid_spacing, 
		const int &dist_cutoff, const double &excluded_radius, const double &max_interatomic_distance)
		: __score(nullptr), __ligand_idatm_types(nullptr) {
		try {
			__identify_gridpoints(centroids, grid, grid_spacing, dist_cutoff, 
				excluded_radius, max_interatomic_distance);
		} catch(...) {
			dbgmsg("FAILURE: constructor of Gpoints failed ... cleaning up resources...");
			throw;
		}
	}

	const Gpoints::Gpoint& Gpoints::get_center_point() const {
		Geom3D::Point center(0,0,0);
		double min_d = HUGE_VAL;
		const Gpoints::Gpoint *center_point = nullptr;
		for (auto &p : get_gridpoints0()) {
			const double d = p.crd().distance(center);
			if (d < min_d) {
				min_d = d;
				center_point = &p;
			}
		}
		if (!center_point)
			throw Error("die : central point cannot be found");
		return *center_point;
	}

	void Gpoints::__identify_gridpoints(const Centro::Centroids &centroids, Molib::Atom::Grid &grid, 
		const double &grid_spacing, const int &dist_cutoff, const double &excluded_radius, 
		const double &max_interatomic_distance) {
	
		Benchmark::reset();
	
		if (centroids.empty()) 
			throw Error("die : there are no centroids");
	
		// find the absolute minimium and maximum coordinates of all centroids
		auto &first = centroids.begin()->second;
		Geom3D::Coordinate min = first[0].get_centroid() - ceil(first[0].get_radial_check());
		Geom3D::Coordinate max = first[0].get_centroid() + ceil(first[0].get_radial_check());
		for (auto &kv : centroids) {
			const int bsite_id = kv.first;
			for (auto &centroid : kv.second) {
				Geom3D::Coordinate min2 = centroid.get_centroid() - ceil(centroid.get_radial_check());
				Geom3D::Coordinate max2 = centroid.get_centroid() + ceil(centroid.get_radial_check());
				if (min2.x() < min.x()) min.set_x(min2.x());
				if (min2.y() < min.y()) min.set_y(min2.y());
				if (min2.z() < min.z()) min.set_z(min2.z());
				if (max2.x() > max.x()) max.set_x(max2.x());
				if (max2.y() > max.y()) max.set_y(max2.y());
				if (max2.z() > max.z()) max.set_z(max2.z());
			}
		}
		dbgmsg("min point = " << min.pdb());
		dbgmsg("max point = " << max.pdb());
		const int total_gridpoints = 3*ceil((max.x()-min.x())/grid_spacing)
									*ceil((max.y()-min.y())/grid_spacing)
									*ceil((max.z()-min.z())/grid_spacing);
		cout <<  "approximately " << total_gridpoints << " gridpoints to evaluate\n\n";
		int model_number = 1;
		int points_kept = 0;
		int gridpoint_counter = 0;
		const double r = grid_spacing/2;
		const double max_d = min.distance(max); // distance between min and max
		const int last_column = ceil(max_d/r);
		const int last_row = ceil(max_d/(sqrt(3)*r));
		const int last_layer = ceil(max_d/(2*r*sqrt(6)/3));
	
		// initialize mapping between gridpoints and discretized 3D space
		__gmap.init(last_column + 1, last_row + 1, last_layer + 1);
	
		Geom3D::Coordinate eval;
		for(int column=0;column<=last_column;column++) {
			int even_column=(column%2==0) ? 1 : 0; // 1 if odd, 0 if even
			for(int row=0;row<=last_row;row++) {
				int even_row=(row%2==0) ? 1 : 0; // 1 if odd, 0 if even
				for(int layer=0;layer<=last_layer;layer++) {
					int even_layer=(layer%2==0) ? 1 : 0; // 1 if odd, 0 if even	
					if ((even_column==0 && even_row==0) || (even_column==1 && even_row==1)) {
						if (even_layer==1) {
							eval.set_x(min.x()+column*r);
							eval.set_y(min.y()+sqrt(3)*r*row);
							eval.set_z(min.z()+layer*2*r*sqrt(6)/3);
						}
						else {
							eval.set_x(min.x()+r+column*r);		 
							eval.set_y(min.y()+r/sqrt(3)+sqrt(3)*r*row); 
							eval.set_z(min.z()+layer*2*r*sqrt(6)/3);      
						}
						
						// mark that this point is not accessible
						//~ __gmap.data[column][row][layer] = false;
						//~ assert(column < __gmap.szi);
						//~ assert(row < __gmap.szj);
						//~ assert(layer < __gmap.szk);
						__gmap.data[column][row][layer] = nullptr;
						
						int okay_min=1;
						int okay_max=1;
						double closest = 10000.0;
						// if the point is within the radial_check of ANY of the centroids... 
						int bsite_id = -1;
						for (auto &kv : centroids) {
							for (auto &centroid : kv.second) {
								if (eval.distance(centroid.get_centroid()) <= centroid.get_radial_check()) {
									bsite_id = kv.first;
									goto END_LOOP;
								}
							}
						}
						END_LOOP:
						if (bsite_id != -1) {
							for (Molib::Atom *a : grid.get_neighbors(eval, dist_cutoff)) {
								Molib::Atom &atom = *a;
								//~ dbgmsg("before getting atom radius");
								const double vdW = atom.radius();
								//~ dbgmsg("vdW = " << vdW);
								const double eval_dist = excluded_radius+0.9*vdW;
								const double distance = atom.crd().distance(eval);
								if (distance <= eval_dist) {
									okay_min = 0;
									dbgmsg("distance = " << distance << " eval_dist = " << eval_dist);
									goto OUTER;
								}
								else {
									okay_min = 1;
									if (distance < closest) closest = distance;
								}
							}
							OUTER:
							if (closest>max_interatomic_distance) okay_max = 0;
							if (okay_min*okay_max > 0) {
								dbgmsg("before adding to __gridpoints");
								__gridpoints[bsite_id].push_back(
									Gpoint{
										eval, 
										IJK{column, row, layer}, 
										__score ? __score->compute_energy(grid, eval, *__ligand_idatm_types) : Array1d<double>()
									});
								dbgmsg("really out");
								//~ __gmap.data[column][row][layer] = true; 
								model_number++;
								points_kept++;
							}
						}
					}
					const int mod = gridpoint_counter % 10000;
					if (mod == 0) {
						double wall_secs = Benchmark::seconds_from_start();
						cout << "Processing gridpoint " << gridpoint_counter 
							<< " of approximately " << total_gridpoints 
							<< " (took " << wall_secs << " seconds)\n";
						Benchmark::reset();
					}
					gridpoint_counter++;
				}
				dbgmsg("column = " << column);
			}
		}
		
		// initialize gmap data here, because push_back can invalidate pointers...
		dbgmsg("initializing gmap");
		for (auto &kv : __gridpoints) {
			for (auto &gpoint : kv.second) {
				__gmap.data[gpoint.ijk().i][gpoint.ijk().j][gpoint.ijk().k] = &gpoint; 
			}
		}
		// the last ones that did not get to the next mod==0
		cout << points_kept << " points kept out of " << gridpoint_counter 
			<< " total gridpoints\n";
	}
	
	void Gpoints::__identify_gridpoints(const double &grid_spacing, const double &radial_check) {
	
		Benchmark::reset();
	
		// find the minimium and maximum coordinates
		Geom3D::Point center(0,0,0);
	
		Geom3D::Coordinate min = center - ceil(radial_check);
		Geom3D::Coordinate max = center + ceil(radial_check);
		dbgmsg("min point = " << min.pdb());
		dbgmsg("max point = " << max.pdb());
		const int total_gridpoints = 3*ceil((max.x()-min.x())/grid_spacing)
									*ceil((max.y()-min.y())/grid_spacing)
									*ceil((max.z()-min.z())/grid_spacing);
		cout <<  "approximately " << total_gridpoints << " gridpoints to evaluate\n\n";
		int points_kept = 0;
		int gridpoint_counter = 0;
		const double r = grid_spacing/2;
		const double max_d = min.distance(max); // distance between min and max
		const int last_column = ceil(max_d/r);
		const int last_row = ceil(max_d/(sqrt(3)*r));
		const int last_layer = ceil(max_d/(2*r*sqrt(6)/3));

		Geom3D::Coordinate eval;
		for(int column=0;column<=last_column;column++) {
			int even_column=(column%2==0) ? 1 : 0; // 1 if odd, 0 if even
			for(int row=0;row<=last_row;row++) {
				int even_row=(row%2==0) ? 1 : 0; // 1 if odd, 0 if even
				for(int layer=0;layer<=last_layer;layer++) {
					int even_layer=(layer%2==0) ? 1 : 0; // 1 if odd, 0 if even	
					if ((even_column==0 && even_row==0) || (even_column==1 && even_row==1)) {
						if (even_layer==1) {
							eval.set_x(min.x()+column*r);
							eval.set_y(min.y()+sqrt(3)*r*row);
							eval.set_z(min.z()+layer*2*r*sqrt(6)/3);
						}
						else {
							eval.set_x(min.x()+r+column*r);		 
							eval.set_y(min.y()+r/sqrt(3)+sqrt(3)*r*row); 
							eval.set_z(min.z()+layer*2*r*sqrt(6)/3);      
						}
						
						// if the point is within the radial_check of ANY of the centroids... 
						if (eval.distance(center) <= radial_check) {
							__gridpoints[0].push_back(Gpoint{eval, IJK{column, 
								row, layer}});
							dbgmsg("lastcolumn = " << last_column << " lastrow = " << last_row << " lastlayer = " << last_layer);
							dbgmsg("column = " << column << " row = " << row << " layer = " << layer);
							dbgmsg("gridpoint0 = " << __gridpoints[0].back().ijk());
							points_kept++;
						}
					}
					const int mod = gridpoint_counter % 10000;
					if (mod == 0) {
						double wall_secs = Benchmark::seconds_from_start();
						cout << "Processing gridpoint " << gridpoint_counter 
							<< " of approximately " << total_gridpoints 
							<< " (took " << wall_secs << " seconds)\n";
						Benchmark::reset();
					}
					gridpoint_counter++;
				}
				dbgmsg("column = " << column);
			}
		}

		// center the ijk coordinates of grid points around the center point (0,0,0)
		Gpoint cp = this->get_center_point(); // here we copy by value intentionally
		dbgmsg("center point = " << cp.ijk());
		for (auto &point : this->get_gridpoints0()) {
			dbgmsg("point = " << point.ijk() << " address = " << &point);
			point.ijk() = point.ijk() - cp.ijk();
			dbgmsg("centered point = " << point.ijk());
		}

		// the last ones that did not get to the next mod==0
		cout << points_kept << " points kept out of " << gridpoint_counter 
			<< " total gridpoints\n";
	}

};
