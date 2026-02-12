#ifndef HMAP_H
#define HMAP_H

#include <pari.h>
#include <paripriv.h>

GEN rand_hmap(ulong n, int f);
GEN hmap_from_invol(GEN invol);
GEN hmap_dual(GEN M);
GEN hmap_numbers(GEN M);
GEN hmap_normalize(GEN M);
GEN hmap_graph(GEN M);


GEN hmap_face_index(GEN M);
GEN hmap_vertex_index(GEN M);
GEN hmap_dfsgen(GEN M);
GEN hmap_dfs(GEN M, GEN data);
GEN hmap_connected_components(GEN M, GEN CCS);

GEN hmap_is_reduced(GEN M);
GEN hmap_is_connected(GEN M);

GEN hmap_gluealongT(GEN M, int flag);
GEN findab(GEN Mone, GEN seed, GEN eindex);
GEN hmap_cutandpaste(GEN M, GEN seed, GEN eindex);
GEN hmap_genus0pres(int f, GEN M, GEN slpsgammai, GEN slpgis, GEN pointersgi);
GEN hmap_buildpres(GEN M, GEN seed, GEN slpsgammai, GEN slpgis, GEN pointersgi, int flag);
GEN hmap_get_presentation(GEN M, int flag);

/**/
GEN hmap_from_monodromy(GEN M, long d, GEN monodromy);

#endif
