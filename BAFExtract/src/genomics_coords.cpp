#ifndef __INDEXING_PER_FORMAT__
#define __INDEXING_PER_FORMAT__

/*
This file lists the indexing bases for different sequence formats. This becomes necessary when output of one method is used as input for another 
method. Note that some
*/

//#define ELAND_BASE		(1)
//#define tagAlign_BASE	(1)
//#define bowtie_BASE		(0) // Default bowtie output.
//#define SAM_BASE		(0)

/*
Codebase indexing: 1 based and both start and end are inclusive.
*/
#define CODEBASE_START_BASE (1)
#define CODEBASE_END_BASE (1)

/*
BED indexing: 0 based, start included, end excluded.
This applies to bedgraphs and also to narrowPeaks.
*/
#define BED_START_BASE (0)
#define BED_END_BASE (1)

/*
Interval indexing: 0 based, start included, end excluded.
*/
#define INTERVAL_START_BASE (0)
#define INTERVAL_END_BASE (1)

/*
Codebase indexing: 0 based, start included, end excluded.
This applies to bedgraphs and also to narrowPeaks.
*/
#define GFF_START_BASE (1)
#define GFF_END_BASE (1)

#define ELAND_START_BASE (1)
#define ELAND_START_BASE (1)

#define tagAlign_START_BASE (1)
#define tagAlign_START_BASE (1)

#define bowtie_START_BASE (0) // Default bowtie output.
#define bowtie_END_BASE (0) // Default bowtie output.

#define SAM_START_BASE (0)
#define SAM_END_BASE (0)

#endif // __INDEXING_PER_FORMAT__

