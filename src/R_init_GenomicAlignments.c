#include "GenomicAlignments.h"
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

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

