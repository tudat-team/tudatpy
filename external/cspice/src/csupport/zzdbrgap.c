/* zzdbrgap.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c_b3 = 100000;
static integer c_b18 = 1000000;

/* $Procedure ZZDBRGAP ( DSKBRIEF, compute gaps in coverage ) */
/* Subroutine */ int zzdbrgap_(integer *corsys, integer *nrec, doublereal *
	bds1, doublereal *bds2, integer *maxn, integer *ncomp, doublereal *
	cbds1, doublereal *cbds2)
{
    /* System generated locals */
    integer bds2_dim2, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static logical grid[1000000];
    static integer srcs[200000], ordx[200000], ordy[200000], vset[200006], 
	    h__, i__, j, k;
    extern /* Subroutine */ int chkin_(char *, ftnlen), moved_(doublereal *, 
	    integer *, doublereal *);
    static integer ncols, nrows;
    extern /* Subroutine */ int rc2grd_(integer *, doublereal *, doublereal *,
	     integer *, integer *, logical *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, logical *)
	    ;
    extern logical failed_(void);
    static integer nr;
    extern /* Subroutine */ int fndcmp_(integer *, integer *, logical *, 
	    integer *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static logical gapval, coverd;
    extern /* Subroutine */ int reglon_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *);
    static doublereal minlon, maxlon;
    extern /* Subroutine */ int chkout_(char *, ftnlen);
    static doublereal outxbd[200000]	/* was [2][100000] */, outybd[200000]	
	    /* was [2][100000] */;
    static integer mrkset[200006], cmporx[200000], civorx[200000], civory[
	    200000], cmpory[200000], tmpset[200006];
    extern logical return_(void);
    static integer minpxx[200000], maxpxx[200000], maxpxy[200000], minpxy[
	    200000];
    extern /* Subroutine */ int ssizei_(integer *, integer *);

/* $ Abstract */

/*     Find spatial coverage gaps for a set of rectangles representing */
/*     the coverage of DSK segments. */

/*     The gap region is represented as a list of rectangles in a */
/*     specified coordinate system. The rectangles are disjoint except */
/*     for possible overlap at their edges. */

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
/*     CORSYS     I   Coordinate system code. */
/*     NREC       I   Number of input rectangles. */
/*     BDS1       I   Bounds of the rectangles' first coordinates. */
/*     BDS2       I   Bounds of the rectangles' second coordinates. */
/*     MAXN       I   Maximum number of output rectangles. */
/*     NCOMP      O   Number of output components. */
/*     CBDS1      O   Bounds of components' first coordinates. */
/*     CBDS2      O   Bounds of components' second coordinates. */
/*     MAXC       P   Maximum number of components. */
/*     MAXGRD     P   Maximum size of pixel grid. */

/* $ Detailed_Input */

/*     CORSYS     is the DSK code for the coordinate system in */
/*                which the input rectangles are represented. */

/*                The supported coordinate systems are */

/*                   Latitudinal */
/*                   Planetodetic */
/*                   Rectangular */

/*                The coordinate system plays little role in */
/*                determination of gaps, other than that longitude */
/*                requires special treatment so that longitude values */
/*                can be compared. */

/*     NREC       is the number of input rectangles. */

/*     BDS1, */
/*     BDS2       are, respectively, the bounds of the first and */
/*                second coordinates of the input rectangles. */
/*                The elements */

/*                   BDS1(J,I), J = 1, 2 */

/*                are the lower and upper bounds of the first */
/*                coordinate. The elements */

/*                   BDS2(J,I), J = 1, 2 */

/*                are the lower and upper bounds of the second */
/*                coordinate. */

/*                Units are radians for angular coordinates and */
/*                km for distance coordinates. */

/*     MAXN       is the maximum number of gap components that */
/*                can be returned. */

/* $ Detailed_Output */

/*     NCOMP      is the number of gap components found. */

/*     CBDS1, */
/*     CBDS2      are, respectively, the bounds of the first and second */
/*                coordinates of the rectangular components comprising */
/*                the gap region. */

/*                The bounds are expressed in the coordinate system */
/*                designated by CORSYS. */

/*                The components are disjoint except for possible */
/*                overlap at their edges. */

/*                Units of the bounds are those of the corresponding */
/*                input coordinates. */

/* $ Parameters */

/*     MAXC       is the maximum number of components that can be */
/*                accommodated by this routine's workspace arrays. */

/*     MAXGRD     is the maximum number of pixels in the workspace */
/*                grid used by this routine. */

/* $ Exceptions */

/*     1)  If any array used in the gap computation is too small, an */
/*         error will be signaled by a routine in the call tree of this */
/*         routine. */

/*     2)  Longitudes outside the range */

/*            -2*pi - ANGMRG : 2*pi + ANGMRG */

/*         are not accepted: if such a value is encountered, the */
/*         error will be diagnosed by a routine in the call tree */
/*         of this routine. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine supports determination of spatial coverage gaps of a */
/*     set of DSK segments. */

/*     The algorithms used in this routine and those it calls attempt */
/*     to operate efficiently: sorting operations are kept to a minimum. */

/*     The algorithms attempt to achieve efficiency in return for memory */
/*     consumption. */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    DSKBRIEF Version 1.0.0, 30-JAN-2017 (NJB) */

/*        Updated to accommodate FNDCMP's new argument list order. */

/*        Original version 1.0.0, 07-OCT-2016 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */

    /* Parameter adjustments */
    bds2_dim2 = *nrec;

    /* Function Body */
    if (return_()) {
	return 0;
    }
    chkin_("ZZDBRGAP", (ftnlen)8);

/*     Initialize the workspace sets. */

    ssizei_(&c_b3, mrkset);
    ssizei_(&c_b3, tmpset);
    ssizei_(&c_b3, vset);
    if (*corsys == 1 || *corsys == 4) {

/*        Adjust the input longitudes to make them usable by RC2GRD. */

	reglon_(nrec, bds1, maxn, &nr, &minlon, &maxlon, outxbd, srcs);

/*        Since REGLON may create new rectangles, the input set of "Y" */
/*        bounds (actually latitude) may not match up with the bounds in */
/*        OUTXBD. Create a new array of "Y" bounds parallel to the array */
/*        of X bounds. In this array, each rectangle has the Y bounds of */
/*        the source box from which it was created. */

	i__1 = nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    outybd[(i__2 = (i__ << 1) - 2) < 200000 && 0 <= i__2 ? i__2 : 
		    s_rnge("outybd", i__2, "zzdbrgap_", (ftnlen)287)] = bds2[(
		    i__4 = (srcs[(i__3 = i__ - 1) < 200000 && 0 <= i__3 ? 
		    i__3 : s_rnge("srcs", i__3, "zzdbrgap_", (ftnlen)287)] << 
		    1) - 2) < bds2_dim2 << 1 && 0 <= i__4 ? i__4 : s_rnge(
		    "bds2", i__4, "zzdbrgap_", (ftnlen)287)];
	    outybd[(i__2 = (i__ << 1) - 1) < 200000 && 0 <= i__2 ? i__2 : 
		    s_rnge("outybd", i__2, "zzdbrgap_", (ftnlen)288)] = bds2[(
		    i__4 = (srcs[(i__3 = i__ - 1) < 200000 && 0 <= i__3 ? 
		    i__3 : s_rnge("srcs", i__3, "zzdbrgap_", (ftnlen)288)] << 
		    1) - 1) < bds2_dim2 << 1 && 0 <= i__4 ? i__4 : s_rnge(
		    "bds2", i__4, "zzdbrgap_", (ftnlen)288)];
	}
    } else {

/*        Just transfer the input X and Y bounds. */

	i__1 = *nrec << 1;
	moved_(bds1, &i__1, outxbd);
	i__1 = *nrec << 1;
	moved_(bds2, &i__1, outybd);
	nr = *nrec;
    }

/*     Map the coordinate rectangles to a pixel grid. Mark */
/*     the coverage with the value .TRUE. */

    coverd = TRUE_;
    gapval = FALSE_;
    rc2grd_(&nr, outxbd, outybd, &c_b18, &c_b3, &coverd, ordx, ordy, civorx, 
	    civory, cmporx, cmpory, &nrows, &ncols, grid);

/*     Map the gaps in the pixel grid to a set of rectangles in pixel */
/*     space. */

    i__ = nrows * ncols;
    k = 0;
    i__1 = i__;
    for (j = 1; j <= i__1; ++j) {
	if (grid[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : s_rnge("grid",
		 i__2, "zzdbrgap_", (ftnlen)323)]) {
	    ++k;
	}
    }
    fndcmp_(&nrows, &ncols, &gapval, maxn, grid, vset, mrkset, tmpset, ncomp, 
	    minpxx, maxpxx, minpxy, maxpxy);
    if (failed_()) {
	chkout_("ZZDBRGAP", (ftnlen)8);
	return 0;
    }

/*     Map the gap rectangles, which are expressed in pixel coordinates, */
/*     to a set of rectangles in the input coordinate system. */

    i__1 = *ncomp;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        If the range of a pixel coordinate is */

/*           A:B */

/*        then the range of the indices of the corresponding bounds in */
/*        the compressed, ordered set of bounds is */

/*           A : B+1 */

/*        Map each bound index to the corresponding value using the */
/*        mappings output by RC2GRD. */

/*        We need to deal with the fact that the arrays CMPORX and */
/*        CMPORY treat the bounds arrays OUTXBD and OUTYBD as */
/*        one-dimensional. */

	j = cmporx[(i__3 = minpxx[(i__2 = i__ - 1) < 200000 && 0 <= i__2 ? 
		i__2 : s_rnge("minpxx", i__2, "zzdbrgap_", (ftnlen)361)] - 1) 
		< 200000 && 0 <= i__3 ? i__3 : s_rnge("cmporx", i__3, "zzdbr"
		"gap_", (ftnlen)361)];
	k = (j - 1) / 2 + 1;
	h__ = j - (k - 1 << 1);
	cbds1[(i__ << 1) - 2] = outxbd[(i__2 = h__ + (k << 1) - 3) < 200000 &&
		 0 <= i__2 ? i__2 : s_rnge("outxbd", i__2, "zzdbrgap_", (
		ftnlen)366)];
	j = cmporx[(i__3 = maxpxx[(i__2 = i__ - 1) < 200000 && 0 <= i__2 ? 
		i__2 : s_rnge("maxpxx", i__2, "zzdbrgap_", (ftnlen)369)]) < 
		200000 && 0 <= i__3 ? i__3 : s_rnge("cmporx", i__3, "zzdbrga"
		"p_", (ftnlen)369)];
	k = (j - 1) / 2 + 1;
	h__ = j - (k - 1 << 1);
	cbds1[(i__ << 1) - 1] = outxbd[(i__2 = h__ + (k << 1) - 3) < 200000 &&
		 0 <= i__2 ? i__2 : s_rnge("outxbd", i__2, "zzdbrgap_", (
		ftnlen)374)];
	j = cmpory[(i__3 = minpxy[(i__2 = i__ - 1) < 200000 && 0 <= i__2 ? 
		i__2 : s_rnge("minpxy", i__2, "zzdbrgap_", (ftnlen)377)] - 1) 
		< 200000 && 0 <= i__3 ? i__3 : s_rnge("cmpory", i__3, "zzdbr"
		"gap_", (ftnlen)377)];
	k = (j - 1) / 2 + 1;
	h__ = j - (k - 1 << 1);
	cbds2[(i__ << 1) - 2] = outybd[(i__2 = h__ + (k << 1) - 3) < 200000 &&
		 0 <= i__2 ? i__2 : s_rnge("outybd", i__2, "zzdbrgap_", (
		ftnlen)382)];
	j = cmpory[(i__3 = maxpxy[(i__2 = i__ - 1) < 200000 && 0 <= i__2 ? 
		i__2 : s_rnge("maxpxy", i__2, "zzdbrgap_", (ftnlen)385)]) < 
		200000 && 0 <= i__3 ? i__3 : s_rnge("cmpory", i__3, "zzdbrga"
		"p_", (ftnlen)385)];
	k = (j - 1) / 2 + 1;
	h__ = j - (k - 1 << 1);
	cbds2[(i__ << 1) - 1] = outybd[(i__2 = h__ + (k << 1) - 3) < 200000 &&
		 0 <= i__2 ? i__2 : s_rnge("outybd", i__2, "zzdbrgap_", (
		ftnlen)390)];
    }
    chkout_("ZZDBRGAP", (ftnlen)8);
    return 0;
} /* zzdbrgap_ */

