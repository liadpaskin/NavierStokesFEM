function meshPart

    call ParMETIS V3 PartMeshKway

end


int ParMETIS V3 PartMeshKway (
idx t *elmdist, idx t *eptr, idx t *eind, idx t *elmwgt, idx t *wgtflag, idx t *numflag,
idx t *ncon, idx t *ncommonnodes, idx t *nparts, real t *tpwgts, real t *ubvec,
idx t *options, idx t *edgecut, idx t *part, MPI Comm *comm
)

idx t *elmdist:

This array describes how the elements of the mesh are distributed among the processors. It is anal-
ogous to the vtxdist array. Its contents are identical for every processor. (See discussion in
Section 4.2.3).

In addition to these four arrays, each processor also requires the array vtxdist[p + 1] that indicates the range of
vertices that are local to each processor. In particular, processor P i stores the vertices from vtxdist[i] up to (but
not including) vertex vtxdist[i + 1].
