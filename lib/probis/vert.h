#ifndef _VERT_H
#define _VERT_H

class Row;


class Vert {
 public:
  Vert(int v, int c, int se, Row *r) { vertex = v;  color = c; size_edge = se; row = r; } // only for v[i][0]
  Vert(int l_id, int c_id, int v, int c, int se, Row *r) { vertex = v;  color = c; size_edge = se; row = r; lista_id = l_id; cluster_id = c_id; }
  Vert() {}
  int lista_id, cluster_id;   // potrebujemo pri iskanju max-weight klike, kjer so verteksi klastri aminokislin
  int vertex;
  int color;
  int size_edge;   // number of vertices that this vertex is connected to
  Row *row;
};




#endif // _VERT_H
