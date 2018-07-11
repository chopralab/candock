#include "grid.h"
#include "atom.h"
#include "debug.hpp"
#include "desc.h"
#include "eelement.h"
#include "ligand.h"
#include "probe.h"

//~ Grid::Grid(pair<Coor, Coor> pc) : minCrd(pc.first), maxCrd(pc.second),
//outsideGrid(false) {
//~
//~ dbgmsg("GRID> Allocating memory and determining size of grid" );
//~ clock_t start = clock();
//~
//~ /* najprej dolocimo stevilo celic, vzamemo nekoliko vecjo kocko okoli
//proteina zaradi ligandov */
//~ minCrd = minCrd - (2*MAXR + INTER_CHAIN_DIST) ;
//~ maxCrd = maxCrd + (2*MAXR + INTER_CHAIN_DIST) ;
//~
//~ cellX = (int) ((maxCrd.x - minCrd.x) / SPACING);
//~ cellY = (int) ((maxCrd.y - minCrd.y) / SPACING);
//~ cellZ = (int) ((maxCrd.z - minCrd.z) / SPACING);
//~
//~ //#ifndef NDEBUG
//~ dbgmsg("minCrd = " << minCrd.x << " " << minCrd.y << " " << minCrd.z );
//~ dbgmsg("maxCrd = " << maxCrd.x << " " << maxCrd.y << " " << maxCrd.z );
//~ dbgmsg("cells  = " << cellX << " " << cellY << " " << cellZ );
//~ //#endif
//~
//~
//~ /* nato alociramo pomnilnik */
//~ cell = new pElement** [cellX];
//~ for (int i = 0; i < cellX; i++) {
//~ cell[i] = new pElement* [cellY];
//~ for (int j = 0; j < cellY; j++) {
//~ cell[i][j] = new pElement [cellZ];
//~ for (int k = 0; k < cellZ; k++) {
//~
//~ cell[i][j][k] = NULL;
//~ }
//~ }
//~ }
//~ dbgmsg("GRID> Time = " << (double)(clock() - start) / CLOCKS_PER_SEC );
//~
//~ }
void Grid::init(pair<Coor, Coor> pc) {
    minCrd = pc.first;
    maxCrd = pc.second;

    dbgmsg("GRID> Allocating memory and determining size of grid");
#ifndef NDEBUG
    clock_t start = clock();
#endif
    /* najprej dolocimo stevilo celic, vzamemo nekoliko vecjo kocko okoli
     * proteina zaradi ligandov */
    minCrd = minCrd - (2 * MAXR + INTER_CHAIN_DIST);
    maxCrd = maxCrd + (2 * MAXR + INTER_CHAIN_DIST);

    cellX = (int)((maxCrd.x - minCrd.x) / SPACING);
    cellY = (int)((maxCrd.y - minCrd.y) / SPACING);
    cellZ = (int)((maxCrd.z - minCrd.z) / SPACING);

    //#ifndef NDEBUG
    dbgmsg("minCrd = " << minCrd.x << " " << minCrd.y << " " << minCrd.z);
    dbgmsg("maxCrd = " << maxCrd.x << " " << maxCrd.y << " " << maxCrd.z);
    dbgmsg("cells  = " << cellX << " " << cellY << " " << cellZ);
    //#endif

    /* nato alociramo pomnilnik */
    cell = new pElement**[cellX];
    for (int i = 0; i < cellX; i++) {
        cell[i] = new pElement*[cellY];
        for (int j = 0; j < cellY; j++) {
            cell[i][j] = new pElement[cellZ];
            for (int k = 0; k < cellZ; k++) {
                cell[i][j][k] = NULL;
            }
        }
    }
    dbgmsg("GRID> Time = " << (double)(clock() - start) / CLOCKS_PER_SEC);
}

Grid::~Grid() {
    if (outsideGrid)
        dbgmsg(
            "Warning (GRID) : Some spheres were outside grid (normal with "
            "--lig) ! ");

    dbgmsg("GRID> Deallocating memory");
#ifndef NDEBUG
    clock_t start = clock();
#endif
    Element* tmp;

    for (int i = 0; i < cellX; i++) {
        for (int j = 0; j < cellY; j++) {
            for (int k = 0; k < cellZ; k++) {
                while (cell[i][j][k] != NULL) {
                    tmp = cell[i][j][k];
                    cell[i][j][k] = cell[i][j][k]->next;
                    delete tmp;
                }
            }
            delete[] cell[i][j];
        }
        delete[] cell[i];
    }
    delete[] cell;

    dbgmsg("GRID> Time = " << (double)(clock() - start) / CLOCKS_PER_SEC);
}

void Grid::deallocate_content() {
    dbgmsg("GRID> Deallocate content only");
#ifndef NDEBUG
    clock_t start = clock();
#endif
    Element* tmp;

    for (int i = 0; i < cellX; i++) {
        for (int j = 0; j < cellY; j++) {
            for (int k = 0; k < cellZ; k++) {
                while (cell[i][j][k] != NULL) {
                    tmp = cell[i][j][k];
                    cell[i][j][k] = cell[i][j][k]->next;
                    delete tmp;
                    /* pomembno je, da je vsebina vsake celice pred uporabo NULL
                     */
                    cell[i][j][k] = NULL;
                }
            }
        }
    }
    dbgmsg("GRID> Time = " << (double)(clock() - start) / CLOCKS_PER_SEC);
}

void Grid::make_grid(Sphere* sphere) {
    int i, j, k;
    Element* tmp;
    i = (int)((sphere->crd.x - minCrd.x) / SPACING);
    j = (int)((sphere->crd.y - minCrd.y) / SPACING);
    k = (int)((sphere->crd.z - minCrd.z) / SPACING);

    /* dodamo samo, ce so i,j,k koordinate v mejah - zaradi prisotnih ligandov
     */
    if (i >= 0 && j >= 0 && k >= 0 && i < cellX && j < cellY && k < cellZ) {
        tmp = new Element();
        tmp->item = sphere;
        tmp->next = cell[i][j][k];
        cell[i][j][k] = tmp;

    } else {
        /* atom je izven grida - to se zgodi, ce imamo ligande */
        outsideGrid = true;
    }
}

void Grid::volume_slice_atom(Sphere* sphere, double max_dist, double distance) {
    /*
      Intented for use with Atom class. Uses EElement.
    */
    int x, y, z;
    int mx, my, mz, nx, ny, nz;
    Sphere* sphere2;
    Element* tmp;
    EElement* tmp2;

    mx = (int)((sphere->crd.x - minCrd.x - max_dist) / SPACING);
    if (mx < 0) mx = 0;
    my = (int)((sphere->crd.y - minCrd.y - max_dist) / SPACING);
    if (my < 0) my = 0;
    mz = (int)((sphere->crd.z - minCrd.z - max_dist) / SPACING);
    if (mz < 0) mz = 0;

    nx = (int)((sphere->crd.x - minCrd.x + max_dist) / SPACING);
    if (nx > cellX - 1) nx = cellX - 1;
    ny = (int)((sphere->crd.y - minCrd.y + max_dist) / SPACING);
    if (ny > cellY - 1) ny = cellY - 1;
    nz = (int)((sphere->crd.z - minCrd.z + max_dist) / SPACING);
    if (nz > cellZ - 1) nz = cellZ - 1;

    for (x = mx; x <= nx; x++)
        for (y = my; y <= ny; y++)
            for (z = mz; z <= nz; z++) {
                tmp = cell[x][y][z];
                while (tmp != NULL) {
                    sphere2 = tmp->item;
                    if (dist(sphere->crd, sphere2->crd) <
                        sphere->r + distance + sphere2->r)
                        if (sphere2->num != sphere->num) {
                            tmp2 = new EElement();
                            tmp2->item = sphere2;
                            tmp2->next = ((Atom*)sphere)->neighb;
                            ((Atom*)sphere)->neighb = tmp2;
                            //              dbgmsg("sphere->num = "<<
                            //              sphere->num << " x = " <<
                            //              sphere->crd.x << "sphere2->num = "<<
                            //              sphere2->num << " x = " <<
                            //              sphere2->crd.x);
                        }
                    tmp = tmp->next;
                }
            }
}

#ifdef CILE
void Grid::volume_slice_probe(Probe* probe, double max_dist, double distance) {
    int x, y, z;
    int mx, my, mz, nx, ny, nz;
    Atom* atom2;
    Element* tmp;
    EElement* tmp2;
    mx = (int)((probe->crd.x - minCrd.x - max_dist) / SPACING);
    if (mx < 0) mx = 0;
    my = (int)((probe->crd.y - minCrd.y - max_dist) / SPACING);
    if (my < 0) my = 0;
    mz = (int)((probe->crd.z - minCrd.z - max_dist) / SPACING);
    if (mz < 0) mz = 0;

    nx = (int)((probe->crd.x - minCrd.x + max_dist) / SPACING);
    if (nx > cellX - 1) nx = cellX - 1;
    ny = (int)((probe->crd.y - minCrd.y + max_dist) / SPACING);
    if (ny > cellY - 1) ny = cellY - 1;
    nz = (int)((probe->crd.z - minCrd.z + max_dist) / SPACING);
    if (nz > cellZ - 1) nz = cellZ - 1;
    for (x = mx; x <= nx; x++)
        for (y = my; y <= ny; y++)
            for (z = mz; z <= nz; z++) {
                tmp = cell[x][y][z];
                while (tmp != NULL) {
                    atom2 = (Atom*)tmp->item;
                    if (dist(probe->crd, atom2->crd) < distance + atom2->r) {
                        tmp2 = new EElement();
                        tmp2->item = atom2;
                        tmp2->next = probe->neighb;
                        probe->neighb = tmp2;
                    }
                    tmp = tmp->next;
                }
            }
}
#endif  // CILE

void Grid::volume_slice_desc(Sphere* sphere, double max_dist, double distance) {
    /*
      Intented for use with Descriptor class. Uses Element as neighb.
    */
    int x, y, z;
    int mx, my, mz, nx, ny, nz;
    Sphere* sphere2;
    Element *tmp, *tmp2;
    double dAB;

    mx = (int)((sphere->crd.x - minCrd.x - max_dist) / SPACING);
    if (mx < 0) mx = 0;
    my = (int)((sphere->crd.y - minCrd.y - max_dist) / SPACING);
    if (my < 0) my = 0;
    mz = (int)((sphere->crd.z - minCrd.z - max_dist) / SPACING);
    if (mz < 0) mz = 0;

    nx = (int)((sphere->crd.x - minCrd.x + max_dist) / SPACING);
    if (nx > cellX - 1) nx = cellX - 1;
    ny = (int)((sphere->crd.y - minCrd.y + max_dist) / SPACING);
    if (ny > cellY - 1) ny = cellY - 1;
    nz = (int)((sphere->crd.z - minCrd.z + max_dist) / SPACING);
    if (nz > cellZ - 1) nz = cellZ - 1;

    for (x = mx; x <= nx; x++)
        for (y = my; y <= ny; y++)
            for (z = mz; z <= nz; z++) {
                tmp = cell[x][y][z];
                while (tmp != NULL) {
                    sphere2 = tmp->item;
                    //          if
                    //          (((Descriptor*)sphere2)->colors.is_visited(color))
                    //          dbgmsg("sphere->r = " << sphere->r << "
                    //          sphere2->r = " << sphere2->r);

                    // pri deskriptorjih je sphere->r in sphere2->r je 0
                    if (sqrt(dAB = dist_fast(sphere->crd, sphere2->crd)) <
                        sphere->r + distance + sphere2->r)
                        if (sphere2->num != sphere->num) {
                            tmp2 = new Element();
                            tmp2->item = sphere2;
                            tmp2->dAB = dAB;
                            tmp2->next = ((Descriptor*)sphere)->neighb;
                            ((Descriptor*)sphere)->neighb = tmp2;
                        }
                    tmp = tmp->next;
                }
            }
}

void Grid::volume_slice_lig(Sphere* sphere, double max_dist, double distance) {
    /*
      Naredimo neighb listo za Ligand class. Sosedi so le ligandi istega tipa
      kot centralni ligand.
    */
    int x, y, z;
    int mx, my, mz, nx, ny, nz;
    Sphere* sphere2;
    Element *tmp, *tmp2;
    double dAB;

    mx = (int)((sphere->crd.x - minCrd.x - max_dist) / SPACING);
    if (mx < 0) mx = 0;
    my = (int)((sphere->crd.y - minCrd.y - max_dist) / SPACING);
    if (my < 0) my = 0;
    mz = (int)((sphere->crd.z - minCrd.z - max_dist) / SPACING);
    if (mz < 0) mz = 0;

    nx = (int)((sphere->crd.x - minCrd.x + max_dist) / SPACING);
    if (nx > cellX - 1) nx = cellX - 1;
    ny = (int)((sphere->crd.y - minCrd.y + max_dist) / SPACING);
    if (ny > cellY - 1) ny = cellY - 1;
    nz = (int)((sphere->crd.z - minCrd.z + max_dist) / SPACING);
    if (nz > cellZ - 1) nz = cellZ - 1;

    for (x = mx; x <= nx; x++)
        for (y = my; y <= ny; y++)
            for (z = mz; z <= nz; z++) {
                tmp = cell[x][y][z];
                while (tmp != NULL) {
                    sphere2 = tmp->item;
                    if (sqrt(dAB = dist_fast(sphere->crd, sphere2->crd)) <
                        distance)
                        // pointerji, ker num pri ligandih ni definiran,
                        // preverimo tip liganda, sosedi so istega tipa
                        if (sphere2 != sphere &&
                            ((Ligand*)sphere2)->type ==
                                ((Ligand*)sphere)->type) {
                            tmp2 = new Element();
                            tmp2->item = sphere2;
                            tmp2->dAB = dAB;
                            tmp2->next = ((Ligand*)sphere)->neighb;
                            ((Ligand*)sphere)->neighb = tmp2;
                            // stejemo sosede
                            ((Ligand*)sphere)->size_neighb++;
                        }
                    tmp = tmp->next;
                }
            }
}

void Grid::volume_slice_surf(Probe* probe, double max_dist, double distance) {
    int x, y, z;
    int mx, my, mz, nx, ny, nz;
    Atom* atom2;
    Element* tmp;
    mx = (int)((probe->crd.x - minCrd.x - max_dist) / SPACING);
    if (mx < 0) mx = 0;
    my = (int)((probe->crd.y - minCrd.y - max_dist) / SPACING);
    if (my < 0) my = 0;
    mz = (int)((probe->crd.z - minCrd.z - max_dist) / SPACING);
    if (mz < 0) mz = 0;

    nx = (int)((probe->crd.x - minCrd.x + max_dist) / SPACING);
    if (nx > cellX - 1) nx = cellX - 1;
    ny = (int)((probe->crd.y - minCrd.y + max_dist) / SPACING);
    if (ny > cellY - 1) ny = cellY - 1;
    nz = (int)((probe->crd.z - minCrd.z + max_dist) / SPACING);
    if (nz > cellZ - 1) nz = cellZ - 1;
    for (x = mx; x <= nx; x++)
        for (y = my; y <= ny; y++)
            for (z = mz; z <= nz; z++) {
                tmp = cell[x][y][z];
                while (tmp != NULL) {
                    atom2 = (Atom*)tmp->item;
                    if (dist(probe->crd, atom2->crd) < distance + atom2->r)
                        if (!atom2->colors.is_visited(probe->color)) {
                            atom2->colors.color[atom2->colors.size++] =
                                probe->color;
                        }
                    tmp = tmp->next;
                }
            }
}
