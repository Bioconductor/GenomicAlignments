#include "GenomicAlignments.h"
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* cigar_utils.c */
	CALLMETHOD_DEF(valid_cigar, 2),
	CALLMETHOD_DEF(explode_cigar_ops, 2),
	CALLMETHOD_DEF(explode_cigar_op_lengths, 2),
	CALLMETHOD_DEF(cigar_op_table, 1),
	CALLMETHOD_DEF(cigar_ranges, 9),
	CALLMETHOD_DEF(cigar_width, 3),
	CALLMETHOD_DEF(cigar_narrow, 3),
	CALLMETHOD_DEF(cigar_qnarrow, 3),
	CALLMETHOD_DEF(ref_locs_to_query_locs, 4),
	CALLMETHOD_DEF(query_locs_to_ref_locs, 4),

/* encodeOverlaps_methods.c */
	CALLMETHOD_DEF(encode_overlaps1, 10),
	CALLMETHOD_DEF(RangesList_encode_overlaps, 7),
	CALLMETHOD_DEF(Hits_encode_overlaps, 10),

	{NULL, NULL, 0}
};

void R_init_GenomicAlignments(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

