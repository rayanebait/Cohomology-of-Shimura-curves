#ifndef MAP_H
#define MAP_H

#include <pari.h>
#include <paripriv.h>

GEN rand_map(ulong n, int f);
GEN map_from_invol(GEN invol);
GEN map_dual(GEN M);
GEN map_numbers(GEN M);
GEN map_normalize(GEN M);
GEN map_graph(GEN M);


GEN map_face_index(GEN M);
GEN map_vertex_index(GEN M);
GEN map_dfsgen(GEN M);
GEN map_dfs(GEN M, GEN data);
GEN map_connected_components(GEN M, GEN CCS);

GEN map_is_reduced(GEN M);
GEN map_is_connected(GEN M);

GEN map_gluealongT(GEN M, int flag);
GEN findab(GEN Mone, GEN seed, GEN eindex);
GEN map_cutandpaste(GEN M, GEN seed, GEN eindex);
GEN map_genus0pres(int f, GEN M, GEN slpsgammai, GEN slpgis, GEN pointersgi);
GEN map_buildpres(GEN M, GEN seed, GEN slpsgammai, GEN slpgis, GEN pointersgi, int flag);
GEN map_get_presentation(GEN M, int flag);

/**/
GEN map_from_monodromy(GEN M, long d, GEN monodromy);

#endif
