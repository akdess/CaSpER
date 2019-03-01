#ifndef __GENOMIC_COORDINATES__
#define __GENOMIC_COORDINATES__

/*
Codebase indexing: 1 based and both start and end are inclusive.
*/
//#define CODEBASE_START_BASE (1)
//#define CODEBASE_END_BASE (1)
struct CODEBASE_COORDS
{
	static const int start_base = 1;
	static const int end_base = 1;
};

// VCF indexing, 1-based. The end is not set.
struct VCF_COORDS
{
	static const int start_base = 1;
	static const int end_base = 1;
};

/*
BED indexing: 0 based, start included, end excluded.
This applies to bedgraphs and also to narrowPeaks.
*/
//#define BED_START_BASE (0)
//#define BED_END_BASE (1)
struct BED_COORDS
{
	static const int start_base = 0;
	static const int end_base = 1;
};

/*
Interval indexing: 0 based, start included, end excluded.
*/
//#define INTERVAL_START_BASE (0)
//#define INTERVAL_END_BASE (1)
struct INTERVAL_COORDS
{
	static const int start_base = 0;
	static const int end_base = 1;
};

/*
Codebase indexing: 0 based, start included, end excluded.
This applies to bedgraphs and also to narrowPeaks.
//*/
//#define GFF_START_BASE (1)
//#define GFF_END_BASE (1)
struct GFF_COORDS
{
	static const int start_base = 1;
	static const int end_base = 1;
};

//#define ELAND_START_BASE (1)
//#define ELAND_START_BASE (1)
struct ELAND_COORDS
{
	static const int start_base = 1;
	static const int end_base = 1;
};

//#define tagAlign_START_BASE (1)
//#define tagAlign_START_BASE (1)
struct TAGALIGN_COORDS
{
	static const int start_base = 1;
	static const int end_base = 1;
};

//#define bowtie_START_BASE (0) // Default bowtie output.
//#define bowtie_END_BASE (0) // Default bowtie output.
struct BOWTIE_COORDS
{
	static const int start_base = 0;
	static const int end_base = 0;
};

struct dbSNP_COORDS
{
	static const int start_base = 1;
	static const int end_base = -1; // There is no end base for this.
};

struct SAM_COORDS
{
	static const int start_base = 1;
	static const int end_base = 1;
};

inline int translate_coord(int coord, int base_src, int base_dest)
{
	return(coord - base_src + base_dest);
}

#endif // __GENOMIC_COORDINATES__

