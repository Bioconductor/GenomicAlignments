#include <Rdefines.h>


/* cigar_utils.c */

const char *_get_cigar_parsing_error();

int _next_cigar_OP(
	const char *cigar_string,
	int offset,
	char *OP,
	int *OPL
);

SEXP valid_cigar(
	SEXP cigar,
	SEXP ans_type
);

SEXP explode_cigar_ops(
	SEXP cigar,
	SEXP ops
);

SEXP explode_cigar_op_lengths(
	SEXP cigar,
	SEXP ops
);

SEXP cigar_op_table(SEXP cigar);

SEXP cigar_ranges(
	SEXP cigar,
	SEXP flag,
	SEXP space,
	SEXP pos,
	SEXP f,
	SEXP ops,
	SEXP drop_empty_ranges,
	SEXP reduce_ranges,
	SEXP with_ops
);

SEXP cigar_width(
	SEXP cigar,
	SEXP flag,
	SEXP space
);

SEXP cigar_narrow(
	SEXP cigar,
	SEXP left_width,
	SEXP right_width
);

SEXP cigar_qnarrow(
	SEXP cigar,
	SEXP left_qwidth,
	SEXP right_qwidth
);


/* mapping_methods.c */

SEXP map_to_genome(
	SEXP start,
	SEXP end,
	SEXP cigar,
	SEXP pos
);

SEXP map_to_transcript(
	SEXP start,
	SEXP end,
	SEXP cigar,
	SEXP pos
);

SEXP ref_locs_to_query_locs(
	SEXP ref_locs,
	SEXP cigar,
	SEXP pos,
	SEXP narrow_left
);

SEXP query_locs_to_ref_locs(
	SEXP query_locs,
	SEXP cigar,
	SEXP pos,
	SEXP narrow_left
);


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

