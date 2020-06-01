/* rc2grd.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure RC2GRD ( DSKBRIEF, rectangles to pixel grid ) */
/* Subroutine */ int rc2grd_(integer *nrec, doublereal *bnds1, doublereal *
	bnds2, integer *maxgrd, integer *maxord, logical *value, integer *
	ord1, integer *ord2, integer *civor1, integer *civor2, integer *
	pxmap1, integer *pxmap2, integer *nrows, integer *ncols, logical *
	grid)
{
    /* System generated locals */
    integer bnds1_dim2, bnds2_dim2, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer i__, j, k;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    static integer ngrid, bsize;
    extern /* Subroutine */ int errdp_(char *, doublereal *, ftnlen);
    extern logical failed_(void);
    extern /* Subroutine */ int sigerr_(char *, ftnlen);
    static integer rngmax;
    extern /* Subroutine */ int chkout_(char *, ftnlen), iovcmp_(doublereal *,
	     integer *, integer *, integer *, integer *), setmsg_(char *, 
	    ftnlen), errint_(char *, integer *, ftnlen);
    extern logical return_(void);
    static integer minpxx, minpxy, maxpxx, maxpxy, col, row;

/* $ Abstract */

/*     Map rectangle list to pixel grid. Each pixel is entirely covered */
/*     by at least one input rectangle or not covered by any rectangle, */
/*     except on its boundary. The pixels covered by at least one */
/*     rectangle are marked with a specified logical value. */

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
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     NREC       I   Number of rectangles. */
/*     BNDS1      I   Bounds of the first coordinates of rectangles. */
/*     BNDS2      I   Bounds of the second coordinates of rectangles. */
/*     MAXGRD     I   Maximum size of output grid. */
/*     MAXORD     I   Maximum size of output order vectors. */
/*     VALUE      I   Value identifying inclusion in rectangles. */
/*     ORD1       O   Order vector of first coordinates. */
/*     ORD2       O   Order vector of second coordinates. */
/*     CIVOR1     O   Compressed inverse order vector of first */
/*                    coordinates. */
/*     CIVOR2     O   Compressed inverse order vector of second */
/*                    coordinates. */
/*     PXMAP1     O   Mapping from pixel first coordinates to */
/*                    rectangle first coordinates. */
/*     PXMAP2     O   Mapping from pixel second coordinates to */
/*                    rectangle second coordinates. */
/*     NROWS      O   Number of rows in the output grid. */
/*     NCOLS      O   Number of columns in the output grid. */
/*     GRID       O   Output grid. */

/* $ Detailed_Input */

/*     NREC       is the number of rectangles in the input set. */
/*                These rectangles represent a spatial region to */
/*                be represented as a union of rectangular components */
/*                that overlap only at their boundaries. */

/*                In DSKBRIEF, the rectangles represent a region having */
/*                no spatial coverage, also called a "gap region." */

/*     BNDS1, */
/*     BNDS2      are, respectively, the bounds of the first and second */
/*                coordinates of the rectangles. The elements */

/*                   BNDS1(J,I), J = 1, 2 */
/*                   BNDS2(J,I), J = 1, 2 */

/*                are the lower and upper bounds of coordinates 1 and for */
/*                the Ithe rectangle. */

/*     MAXGRD     is the maximum size of the output grid. */

/*     MAXORD     is the maximum size of the output order vectors and */
/*                compressed inverse order vectors. */

/*     VALUE      is the logical value used in the pixel grid to mark */
/*                pixels that are in the region covered by the input */
/*                rectangles. */

/* $ Detailed_Output */

/*     ORD1, */
/*     ORD2       are, respectively, order vectors for 1-dimensional */
/*                arrays of all values of the first and second */
/*                coordinates. Upper and lower bounds are merged in */
/*                these arrays. */

/*     CIVOR1, */
/*     CIVOR2     are, respectively, compressed inverse order vectors */
/*                for 1-dimensional arrays of all values of the first */
/*                and second coordinates. */

/*     PXMAP1, */
/*     PXMAP2     are, respectively, arrays mapping pixel boundaries */
/*                to rectangle bounds. PXMAP1 maps pixel X-boundaries */
/*                to values in BNDS1; PXMAP2 maps pixel Y-boundaries */
/*                to values in BNDS2. */

/*                Pixel X-coordinates range from 1 to NCOLS+1; pixel */
/*                Y-coordinate range from 1 to NROWS+1. */

/*     NROWS      is the number of rows in the output pixel grid. */

/*     NCOLS      is number of columns in the output pixel grid. */

/*     GRID       is a pixel grid. GRID is an array of logical values; */
/*                the dimensions of GRID are NROWS x NCOLS. Each pixel */
/*                corresponds to a rectangle in the input coordinate */
/*                space that is either entirely in one of the input */
/*                rectangles, or is, with the exception of its edges, */
/*                outside of all of them. */

/*                The elements of GRID having the value VALUE are those */
/*                pixels corresponding to rectangles in the input */
/*                coordinate space that are entirely contained in an */
/*                input rectangle. */

/*                DSKBRIEF uses the marked pixels of GRID to denote */
/*                pixels belonging to regions covered by a DSK */
/*                segment. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If NREC is less than 1, the error SPICE(VALUEOUTOFRANGE) is */
/*         signaled. */

/*     2)  If MAXGRD is less than the required size of the output grid, */
/*         the error SPICE(VALUEOUTOFRANGE) is signaled. */

/*     3)  If MAXORD is less than 1, the error SPICE(VALUEOUTOFRANGE) is */
/*         signaled. */

/*     4)  No input rectangle may have an edge of zero length. If any */
/*         pair of lower and upper coordinate bounds is not strictly */
/*         increasing, the error SPICE(INVALIDBOUNDS) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine supports determination of DSK spatial coverage gaps. */

/*     The coverage "gap region" of a set of DSK segments is that */
/*     region, within the rectangle defined by the global extrema of the */
/*     extents of the segments' first and second coordinates, where */
/*     there is no spatial coverage. */

/*     The pixel grid created by this routine has the property that each */
/*     pixel corresponds to a rectangle in the input coordinate space */
/*     that is either entirely in one of the input rectangles, or is, */
/*     with the exception of its edges, outside of all of them. */

/*     The pixel grid simplifies the task of efficiently locating */
/*     "components" of the gap region, where components are rectangular */
/*     regions in the input coordinate space that overlap at most at */
/*     their edges. */

/*     Since the edges of pixels correspond to boundaries of the input */
/*     rectangles, it's easy to map a rectangular union of pixels */
/*     to a rectangle in the coordinate space of the input rectangles. */
/*     Each image rectangle is a component of the gap region. */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/*     DSKBRIEF Version 2.0.0, 06-OCT-2016 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local variables */


/*     Saved values */

    /* Parameter adjustments */
    bnds2_dim2 = *nrec;
    bnds1_dim2 = *nrec;

    /* Function Body */
    if (return_()) {
	return 0;
    }
    chkin_("RC2GRD", (ftnlen)6);

/*     Check input size arguments for obvious initialization errors. */

    if (*nrec < 1) {
	setmsg_("NREC is #; must be positive.", (ftnlen)28);
	errint_("#", nrec, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("RC2GRD", (ftnlen)6);
	return 0;
    }
    if (*maxgrd < 1) {
	setmsg_("MAXGRD is #; must be positive.", (ftnlen)30);
	errint_("#", maxgrd, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("RC2GRD", (ftnlen)6);
	return 0;
    }
    if (*maxord < 1) {
	setmsg_("MAXORD is #; must be positive.", (ftnlen)30);
	errint_("#", maxord, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("RC2GRD", (ftnlen)6);
	return 0;
    }

/*     All input rectangle heights and widths must be strictly */
/*     positive. */

    i__1 = *nrec;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (bnds1[(i__2 = (i__ << 1) - 1) < bnds1_dim2 << 1 && 0 <= i__2 ? 
		i__2 : s_rnge("bnds1", i__2, "rc2grd_", (ftnlen)304)] <= 
		bnds1[(i__3 = (i__ << 1) - 2) < bnds1_dim2 << 1 && 0 <= i__3 ?
		 i__3 : s_rnge("bnds1", i__3, "rc2grd_", (ftnlen)304)]) {
	    setmsg_("BNDS1(2,#) = #; BNDS1(1,#) = #. Rectangle widths (and h"
		    "eights) must be positive.", (ftnlen)80);
	    errint_("#", &i__, (ftnlen)1);
	    errdp_("#", &bnds1[(i__2 = (i__ << 1) - 1) < bnds1_dim2 << 1 && 0 
		    <= i__2 ? i__2 : s_rnge("bnds1", i__2, "rc2grd_", (ftnlen)
		    309)], (ftnlen)1);
	    errint_("#", &i__, (ftnlen)1);
	    errdp_("#", &bnds1[(i__2 = (i__ << 1) - 2) < bnds1_dim2 << 1 && 0 
		    <= i__2 ? i__2 : s_rnge("bnds1", i__2, "rc2grd_", (ftnlen)
		    311)], (ftnlen)1);
	    sigerr_("SPICE(INVALIDBOUNDS)", (ftnlen)20);
	    chkout_("RC2GRD", (ftnlen)6);
	    return 0;
	}
	if (bnds2[(i__2 = (i__ << 1) - 1) < bnds2_dim2 << 1 && 0 <= i__2 ? 
		i__2 : s_rnge("bnds2", i__2, "rc2grd_", (ftnlen)318)] <= 
		bnds2[(i__3 = (i__ << 1) - 2) < bnds2_dim2 << 1 && 0 <= i__3 ?
		 i__3 : s_rnge("bnds2", i__3, "rc2grd_", (ftnlen)318)]) {
	    setmsg_("BNDS2(2,#) = #; BNDS2(1,#) = #. Rectangle heights (and "
		    "widths) must be positive.", (ftnlen)80);
	    errint_("#", &i__, (ftnlen)1);
	    errdp_("#", &bnds2[(i__2 = (i__ << 1) - 1) < bnds2_dim2 << 1 && 0 
		    <= i__2 ? i__2 : s_rnge("bnds2", i__2, "rc2grd_", (ftnlen)
		    323)], (ftnlen)1);
	    errint_("#", &i__, (ftnlen)1);
	    errdp_("#", &bnds2[(i__2 = (i__ << 1) - 2) < bnds2_dim2 << 1 && 0 
		    <= i__2 ? i__2 : s_rnge("bnds2", i__2, "rc2grd_", (ftnlen)
		    325)], (ftnlen)1);
	    sigerr_("SPICE(INVALIDBOUNDS)", (ftnlen)20);
	    chkout_("RC2GRD", (ftnlen)6);
	    return 0;
	}
    }

/*     Find the order of the array of X bounds. We treat the array as a */
/*     one-dimensional array of length 2*NREC. */

/*     Produce the corresponding "compressed" inverse order vector. By */
/*     "compressed" we mean: suppose the set of input values were sorted */
/*     and compressed so that it contained no duplicates. For each */
/*     member of the original value array, map the member's index to the */
/*     index of the member in the compressed, sorted array. The */
/*     compressed inverse order vector contains this mapping. */

    bsize = *nrec << 1;
    iovcmp_(bnds1, &bsize, ord1, civor1, &rngmax);
    if (failed_()) {
	chkout_("RC2GRD", (ftnlen)6);
	return 0;
    }
    i__1 = bsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ord1[i__ - 1];
	k = (j - 1) / 2 + 1;
    }

/*     The width of the pixel grid is one less than the number of */
/*     distinct X bound values. */

    *ncols = rngmax - 1;

/*     Get the order vector and compressed inverse order vector of */
/*     the Y bounds. (Note we have the same number of X and Y */
/*     bounds.) */

    iovcmp_(bnds2, &bsize, ord2, civor2, &rngmax);
    if (failed_()) {
	chkout_("RC2GRD", (ftnlen)6);
	return 0;
    }

/*     The height of the pixel grid is one less than the number of */
/*     distinct Y bound values. */

    *nrows = rngmax - 1;

/*     Check the grid size again, now that we know how large it */
/*     needs to be. */

    ngrid = *nrows * *ncols;
    if (*maxgrd < ngrid) {
	setmsg_("MAXGRD is #; must be have size at least # in order to hold "
		"pixels for current set of rectangles.", (ftnlen)96);
	errint_("#", maxgrd, (ftnlen)1);
	errint_("#", &ngrid, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("RC2GRD", (ftnlen)6);
	return 0;
    }

/*     A program using this routine normally will need to map pixel */
/*     coordinates back to their corresponding d.p. values. Create */
/*     arrays to represent these mappings. */

/*     Note that, since the bounds arrays are generally larger than the */
/*     corresponding pixel grid dimensions, the mappings we're about to */
/*     perform (which map the integers 1:BIZE into the ranges of the */
/*     compressed inverse order vectors) are not 1-1. They're still */
/*     valid; the process is just a bit ungainly because it can involve */
/*     overwriting elements of the output array. Each time this happens, */
/*     the affected output array element gets overwritten with the same */
/*     value it already had. */

/*     We'll store the mappings in the arrays PXMAP1 and PXMAP2. */
/*     The pixel coordinates */

/*        ( I, J ) */

/*     correspond to the double precision coordinates */

/*        BNDS1( PXMAP1(I) ) */
/*        BNDS2( PXMAP2(J) ) */

/*     where we're treating BNDS1 and BNDS2 as one-dimensional */
/*     arrays of length BSIZE. */

    i__1 = bsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pxmap1[civor1[i__ - 1] - 1] = i__;
	pxmap2[civor2[i__ - 1] - 1] = i__;
    }

/*     Now map all rectangles to the integer indices of their */
/*     bounds in pixel space. Note that the pixel grid has */
/*     dimensions */

/*        ( NROWS, NCOLS ) */

/*     and the ranges of the integer coordinates of the */
/*     rectangle boundaries are */

/*        1 : NROWS + 1 */
/*        1 : NCOLS + 1 */

/*     We'll fill in the pixel grid to indicate which pixels are */
/*     covered by rectangles, and which ones lie in gaps. */

/*     Initialize the grid to indicate that it consists of one */
/*     large gap. */

    i__1 = ngrid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	grid[i__ - 1] = ! (*value);
    }

/*     For each input rectangle, mark the corresponding pixels */
/*     covered by the rectangle. Note that maximum pixel indices */
/*     are less by one than those of the corresponding rectangle */
/*     upper bound indices. */

    i__1 = *nrec;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Compute the bounds of the current rectangle in pixel */
/*        space. Recall that the all bounds for a given coordinate */
/*        (X or Y) are combined in a sequence of size 2*NREC. */

	j = (i__ - 1 << 1) + 1;
	minpxx = civor1[j - 1];
	minpxy = civor2[j - 1];
	k = i__ << 1;
	maxpxx = civor1[k - 1];
	maxpxy = civor2[k - 1];
	i__2 = maxpxx - 1;
	for (col = minpxx; col <= i__2; ++col) {
	    i__3 = maxpxy - 1;
	    for (row = minpxy; row <= i__3; ++row) {

/*              Mark the pixel at indices (ROW, COL) as */
/*              covered. */

		j = *nrows * (col - 1) + row;
		grid[j - 1] = *value;
	    }
	}
    }
    chkout_("RC2GRD", (ftnlen)6);
    return 0;
} /* rc2grd_ */

