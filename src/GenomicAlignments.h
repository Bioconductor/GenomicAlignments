#include <Rdefines.h>


/* encodeOverlaps_methods.c */

SEXP encode_overlaps1(
	SEXP query_start,
	SEXP query_width,
	SEXP query_space,
	SEXP query_break,
	SEXP flip_query,
	SEXP subject_start,
	SEXP subject_width,
	SEXP subject_space,
	SEXP as_matrix,
	SEXP as_raw
);

SEXP RangesList_encode_overlaps(
	SEXP query_starts,
	SEXP query_widths,
	SEXP query_spaces,
	SEXP query_breaks,
	SEXP subject_starts,
	SEXP subject_widths,
	SEXP subject_spaces
);

SEXP Hits_encode_overlaps(
	SEXP query_starts,
	SEXP query_widths,
	SEXP query_spaces,
	SEXP query_breaks,
	SEXP subject_starts,
	SEXP subject_widths,
	SEXP subject_spaces,
	SEXP query_hits,
	SEXP subject_hits,
	SEXP flip_query
);

