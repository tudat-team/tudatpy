/* grpseg.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure GRPSEG ( Group DSK segments having shared attributes ) */
/* Subroutine */ int grpseg_(S_fp udcomp, logical *usetim, doublereal *timtol,
	 integer *n, integer *handls, integer *dlads, integer *mxpool, 
	integer *mxgrp, integer *ng, integer *sgptrs, integer *sgpool, 
	integer *seglst, logical *segtm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer head, node, tail, i__, j, k;
    logical match;
    extern /* Subroutine */ int chkin_(char *, ftnlen), lnkan_(integer *, 
	    integer *);
    extern integer lnktl_(integer *, integer *);
    extern logical failed_(void);
    extern /* Subroutine */ int lnkila_(integer *, integer *, integer *);
    logical tmatch;
    extern /* Subroutine */ int setmsg_(char *, ftnlen);
    extern logical return_(void);
    extern /* Subroutine */ int errint_(char *, integer *, ftnlen), sigerr_(
	    char *, ftnlen), chkout_(char *, ftnlen), lnkini_(integer *, 
	    integer *);

/* $ Abstract */

/*     Given arrays of handles and DLA descriptors of DSK segments, */
/*     and given a callback comparison function, group the segments */
/*     into subsets having matching attributes. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     DSKBRIEF */

/* $ Declarations */

/*     Include file dla.inc */

/*     This include file declares parameters for DLA format */
/*     version zero. */

/*        Version 3.0.1 17-OCT-2016 (NJB) */

/*           Corrected comment: VERIDX is now described as a DAS */
/*           integer address rather than a d.p. address. */

/*        Version 3.0.0 20-JUN-2006 (NJB) */

/*           Changed name of parameter DSCSIZ to DLADSZ. */

/*        Version 2.0.0 09-FEB-2005 (NJB) */

/*           Changed descriptor layout to make backward pointer */
/*           first element.  Updated DLA format version code to 1. */

/*           Added parameters for format version and number of bytes per */
/*           DAS comment record. */

/*        Version 1.0.0 28-JAN-2004 (NJB) */


/*     DAS integer address of DLA version code. */


/*     Linked list parameters */

/*     Logical arrays (aka "segments") in a DAS linked array (DLA) file */
/*     are organized as a doubly linked list.  Each logical array may */
/*     actually consist of character, double precision, and integer */
/*     components.  A component of a given data type occupies a */
/*     contiguous range of DAS addresses of that type.  Any or all */
/*     array components may be empty. */

/*     The segment descriptors in a SPICE DLA (DAS linked array) file */
/*     are connected by a doubly linked list.  Each node of the list is */
/*     represented by a pair of integers acting as forward and backward */
/*     pointers.  Each pointer pair occupies the first two integers of a */
/*     segment descriptor in DAS integer address space.  The DLA file */
/*     contains pointers to the first integers of both the first and */
/*     last segment descriptors. */

/*     At the DLA level of a file format implementation, there is */
/*     no knowledge of the data contents.  Hence segment descriptors */
/*     provide information only about file layout (in contrast with */
/*     the DAF system).  Metadata giving specifics of segment contents */
/*     are stored within the segments themselves in DLA-based file */
/*     formats. */


/*     Parameter declarations follow. */

/*     DAS integer addresses of first and last segment linked list */
/*     pointer pairs.  The contents of these pointers */
/*     are the DAS addresses of the first integers belonging */
/*     to the first and last link pairs, respectively. */

/*     The acronyms "LLB" and "LLE" denote "linked list begin" */
/*     and "linked list end" respectively. */


/*     Null pointer parameter. */


/*     Segment descriptor parameters */

/*     Each segment descriptor occupies a contiguous */
/*     range of DAS integer addresses. */

/*        The segment descriptor layout is: */

/*           +---------------+ */
/*           | BACKWARD PTR  | Linked list backward pointer */
/*           +---------------+ */
/*           | FORWARD PTR   | Linked list forward pointer */
/*           +---------------+ */
/*           | BASE INT ADDR | Base DAS integer address */
/*           +---------------+ */
/*           | INT COMP SIZE | Size of integer segment component */
/*           +---------------+ */
/*           | BASE DP ADDR  | Base DAS d.p. address */
/*           +---------------+ */
/*           | DP COMP SIZE  | Size of d.p. segment component */
/*           +---------------+ */
/*           | BASE CHR ADDR | Base DAS character address */
/*           +---------------+ */
/*           | CHR COMP SIZE | Size of character segment component */
/*           +---------------+ */

/*     Parameters defining offsets for segment descriptor elements */
/*     follow. */


/*     Descriptor size: */


/*     Other DLA parameters: */


/*     DLA format version.  (This number is expected to occur very */
/*     rarely at integer address VERIDX in uninitialized DLA files.) */


/*     Characters per DAS comment record. */


/*     End of include file dla.inc */


/*     Include file dskdsc.inc */

/*     This include file declares parameters for DSK segment descriptors. */

/* -       SPICELIB Version 1.0.0 08-FEB-2017 (NJB) */

/*           Updated version info. */

/*           22-JAN-2016 (NJB) */

/*              Added parameter for data class 2. Changed name of data */
/*              class 1 parameter. Corrected data class descriptions. */

/*           13-MAY-2010 (NJB) */

/*              Descriptor now contains two ID codes, one for the */
/*              surface, one for the associated ephemeris object. This */
/*              supports association of multiple surfaces with one */
/*              ephemeris object without creating file management */
/*              issues. */

/*              Room was added for coordinate system definition */
/*              parameters. */

/*               Flag arrays and model ID/component entries were deleted. */

/*            11-SEP-2008 (NJB) */


/*     DSK segment descriptors are implemented as an array of d.p. */
/*     numbers.  Note that each integer descriptor datum occupies one */
/*     d.p. value. */




/*     Segment descriptor parameters */

/*     Each segment descriptor occupies a contiguous */
/*     range of DAS d.p. addresses. */

/*        The DSK segment descriptor layout is: */

/*           +---------------------+ */
/*           | Surface ID code     | */
/*           +---------------------+ */
/*           | Center ID code      | */
/*           +---------------------+ */
/*           | Data class code     | */
/*           +---------------------+ */
/*           | Data type           | */
/*           +---------------------+ */
/*           | Ref frame code      | */
/*           +---------------------+ */
/*           | Coord sys code      | */
/*           +---------------------+ */
/*           | Coord sys parameters|  {10 elements} */
/*           +---------------------+ */
/*           | Min coord 1         | */
/*           +---------------------+ */
/*           | Max coord 1         | */
/*           +---------------------+ */
/*           | Min coord 2         | */
/*           +---------------------+ */
/*           | Max coord 2         | */
/*           +---------------------+ */
/*           | Min coord 3         | */
/*           +---------------------+ */
/*           | Max coord 3         | */
/*           +---------------------+ */
/*           | Start time          | */
/*           +---------------------+ */
/*           | Stop time           | */
/*           +---------------------+ */

/*     Parameters defining offsets for segment descriptor elements */
/*     follow. */


/*     Surface ID code: */


/*     Central ephemeris object NAIF ID: */


/*     Data class: */

/*     The "data class" is a code indicating the category of */
/*     data contained in the segment. */


/*     Data type: */


/*     Frame ID: */


/*     Coordinate system code: */


/*     Coordinate system parameter start index: */


/*     Number of coordinate system parameters: */


/*     Ranges for coordinate bounds: */


/*     Coverage time bounds: */


/*     Descriptor size (24): */


/*     Data class values: */

/*        Class 1 indicates a surface that can be represented as a */
/*                single-valued function of its domain coordinates. */

/*                An example is a surface defined by a function that */
/*                maps each planetodetic longitude and latitude pair to */
/*                a unique altitude. */


/*        Class 2 indicates a general surface. Surfaces that */
/*                have multiple points for a given pair of domain */
/*                coordinates---for example, multiple radii for a given */
/*                latitude and longitude---belong to class 2. */



/*     Coordinate system values: */

/*        The coordinate system code indicates the system to which the */
/*        tangential coordinate bounds belong. */

/*        Code 1 refers to the planetocentric latitudinal system. */

/*        In this system, the first tangential coordinate is longitude */
/*        and the second tangential coordinate is latitude. The third */
/*        coordinate is radius. */



/*        Code 2 refers to the cylindrical system. */

/*        In this system, the first tangential coordinate is radius and */
/*        the second tangential coordinate is longitude. The third, */
/*        orthogonal coordinate is Z. */



/*        Code 3 refers to the rectangular system. */

/*        In this system, the first tangential coordinate is X and */
/*        the second tangential coordinate is Y. The third, */
/*        orthogonal coordinate is Z. */



/*        Code 4 refers to the planetodetic/geodetic system. */

/*        In this system, the first tangential coordinate is longitude */
/*        and the second tangential coordinate is planetodetic */
/*        latitude. The third, orthogonal coordinate is altitude. */



/*     End of include file dskdsc.inc */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     UDCOMP     I   Comparison function for DSK segments. */
/*     USETIM     I   Time comparison flag. */
/*     TIMTOL     I   Time comparison tolerance. */
/*     N          I   Size of handle and DLA descriptor arrays. */
/*     HANDLS     I   Handles of DSK files. */
/*     DLADS      I   DLA descriptors of DSK segments. */
/*     MXPOOL     I   Maximum pool size. */
/*     MXGRP      I   Maximum group count. */
/*     NG         O   Number of segment groups. */
/*     SGPTRS     O   Pointers from groups to segment lists. */
/*     SGPOOL     O   Segment pool. */
/*     SEGLST     O   Map from pool nodes to input segment array indices. */
/*     SEGTM      O   Array of time match flags for segment groups. */

/* $ Detailed_Input */

/*     UDCOMP     is a callback function which is used by this */
/*                routine to compare attributes of DSK segments. */
/*                The calling sequence of UDCOMP is */

/*                   UDCOMP ( HAN1, DLADS1, HAN2,  DLADS2, */
/*                            TOL,  MATCH,  TMATCH        ) */

/*                Inputs: */

/*                   INTEGER           HAN1         {handle of segment 1} */
/*                   DOUBLE PRECISION  DLADS1(*)    {DLA descriptor of */
/*                                                   segment 1} */
/*                   INTEGER           HAN2         {handle of segment 2} */
/*                   DOUBLE PRECISION  DLADS2(*)    {DLA descriptor of */
/*                                                   segment 2} */
/*                   DOUBLE PRECISION  TOL          {Time tolerance} */

/*                Outputs: */

/*                   LOGICAL           MATCH        {flag indicating */
/*                                                   whether segments */
/*                                                   match} */

/*                   LOGICAL           TMATCH       {flag indicating */
/*                                                   whether segment */
/*                                                   time bounds match} */


/*     USETIM     is a logical flag that indicates whether to consider */
/*                segment time bounds when grouping segments. */

/*     TIMTOL     is a time tolerance to use for comparing time bounds. */
/*                Units are TDB seconds. */

/*     N          is the number of elements in the input array HANDLS */
/*                and the number of DLA descriptors in the input array */
/*                DLADS. */


/*     HANDLS */
/*     DLADS      are, respectively, arrays of handles of DSK files and */
/*                DLA descriptors. There are N entries in each array. */
/*                The Ith element of HANDLS and the Ith descriptor in */
/*                DLADS correspond to the Ith segment to be grouped. */

/*     MXPOOL     is the maximum number of entries that can be */
/*                accommodated in the linked list pool array SGPOOL. */

/*     MXGRP      is the maximum number of entries that can be */
/*                placed in the array SGPTRS. */


/* $ Detailed_Output */

/*     NG         is the number of groups comprised by the input */
/*                segments. */


/*     SGPTRS     is an array that maps group numbers to head nodes of */
/*                segment lists in SGPOOL. SGPTRS(I) is the index of the */
/*                head node of the segment list of the Ith group. */

/*                SGPTRS must be declared with size at least MXGRP. */


/*     SGPOOL     is a doubly linked list pool used to store segment */
/*                lists. Entries for a given segment group are linked */
/*                together. */

/*                SGPOOL must be declared with its second dimension */
/*                at least equal to MXPOOL. */


/*     SEGLST     is an array that maps nodes of SGPOOL to entries */
/*                in the input arrays HANDLS and DLADS. The Ith node */
/*                of SGPOOL */

/*                   SGPOOL(*,I) */

/*                corresponds to the segment designated by */

/*                   HANDLS(    SEGLST(I) ) */
/*                   DLADS ( *, SEGLST(I) ) */


/*     SEGTM      is an array of flags associated with segment groups. */
/*                The Ith element of SEGTIM is associated with the */
/*                group starting at node SGPTRS(I). Each element of */
/*                SEGTM that is associated with a group is .TRUE. if */
/*                and only if the segments in that group have matching */
/*                time bounds. The elements of SEGTM that are not */
/*                assocated with groups are undefined. */

/*                SEGTM must be declared with size at least MXGRP. */


/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the number of segments is not at least 1, the error */
/*         SPICE(INVALIDCOUNT) is signaled. */

/*     2)  If the pool size is not at least N, the error */
/*         SPICE(ARRAYTOOSMALL) is signaled. */

/*     3)  If the maximum group count is not at least 1, the error */
/*         SPICE(INVALIDSIZE) is signaled. Normally room for multiple */
/*         groups should be provided by the caller. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine groups segments for abbreviated summaries */
/*     created by DSKBRIEF. Segments that are considered by UDCOMP */
/*     to have matching attributes are grouped together. */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/*     DSKBRIEF Version 1.0.0, 15-FEB-2017 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("GRPSEG", (ftnlen)6);
    if (*n < 1) {
	setmsg_("Number of segments is #; must be at least 1.", (ftnlen)44);
	errint_("#", n, (ftnlen)1);
	sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
	chkout_("GRPSEG", (ftnlen)6);
	return 0;
    }
    if (*mxpool < *n) {
	setmsg_("Pool size is #; must be at least #.", (ftnlen)35);
	errint_("#", mxpool, (ftnlen)1);
	errint_("#", n, (ftnlen)1);
	sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
	chkout_("GRPSEG", (ftnlen)6);
	return 0;
    }
    if (*mxgrp < 1) {

/*        Normally MXGRP should allow room for multiple groups. However, */
/*        we can't know in advance how much room is needed, except that */
/*        the group count can't exceed N. We simply catch clearly */
/*        invalid values here. */

	setmsg_("Group array size is #; must be at least 1.", (ftnlen)42);
	errint_("#", mxgrp, (ftnlen)1);
	sigerr_("SPICE(INVALIDSIZE)", (ftnlen)18);
	chkout_("GRPSEG", (ftnlen)6);
	return 0;
    }

/*     Initialize the segment pool. */

    lnkini_(mxpool, sgpool);

/*     Initialize the time match flags. */

    i__1 = *mxgrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	segtm[i__ - 1] = TRUE_;
    }

/*     The first segment is automatically in the first group. */

    lnkan_(sgpool, &node);
    *ng = 1;
    sgptrs[*ng - 1] = node;
    seglst[node - 1] = 1;

/*     Now examine the other segments. If a segment matches one */
/*     we've already seen, the segment is added to the latter's */
/*     group. Otherwise, the segment becomes the first member of */
/*     a new group. */

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {

/*        Compare the Ith segment to the first segment */
/*        of each group, until we find a match or run out */
/*        of groups. */

	match = FALSE_;
	j = 0;
	while(! match && j < *ng) {
	    ++j;

/*           Compare the Ith segment against the first segment of the */
/*           Jth group. */

	    head = sgptrs[j - 1];
	    k = seglst[head - 1];
	    (*udcomp)(&handls[i__ - 1], &dlads[(i__ << 3) - 8], &handls[k - 1]
		    , &dlads[(k << 3) - 8], timtol, &match, &tmatch);
	    if (*usetim) {

/*              Time bounds are being considered; we have a match only */
/*              if the time bounds and the other attributes match. */

		match = match && tmatch;
	    }
	    if (match) {

/*              Indicate whether the group to which this segment belongs */
/*              has inconsistent time tags. */

		if (! tmatch) {
		    segtm[j - 1] = FALSE_;
		}
	    }
	    if (failed_()) {
		chkout_("GRPSEG", (ftnlen)6);
		return 0;
	    }
	}
	if (match) {

/*           The Ith segment belongs to group J. Link a node for the */
/*           segment to the tail of the list for this group. */

	    tail = lnktl_(&head, sgpool);
	    lnkan_(sgpool, &node);
	    if (failed_()) {
		chkout_("GRPSEG", (ftnlen)6);
		return 0;
	    }
	    lnkila_(&tail, &node, sgpool);
	    seglst[node - 1] = i__;
	} else {

/*           This segment is in a category of its own. Create a new */
/*           segment pointer for it. */

	    if (*ng == *mxgrp) {
		setmsg_("Size of group array is #; cannot add new element.", (
			ftnlen)49);
		errint_("#", mxgrp, (ftnlen)1);
		sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
		chkout_("GRPSEG", (ftnlen)6);
		return 0;
	    }
	    lnkan_(sgpool, &node);
	    if (failed_()) {
		chkout_("GRPSEG", (ftnlen)6);
		return 0;
	    }
	    ++(*ng);
	    sgptrs[*ng - 1] = node;

/*           Associate the Ith segment with the new node. */

	    seglst[node - 1] = i__;
	}
    }
    chkout_("GRPSEG", (ftnlen)6);
    return 0;
} /* grpseg_ */

