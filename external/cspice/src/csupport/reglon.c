/* reglon.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b17 = 1e-12;

/* $Procedure REGLON ( Regularize longitude intervals ) */
/* Subroutine */ int reglon_(integer *nivals, doublereal *bounds, integer *
	maxn, integer *nout, doublereal *minlon, doublereal *maxlon, 
	doublereal *outbds, integer *srcs)
{
    /* System generated locals */
    integer bounds_dim2, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    integer nreq;
    doublereal a, b;
    extern /* Subroutine */ int zznrmlon_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    integer i__;
    doublereal loclb;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    doublereal locub;
    extern doublereal dpmin_(void), dpmax_(void);
    extern /* Subroutine */ int errdp_(char *, doublereal *, ftnlen);
    extern doublereal twopi_(void);
    doublereal lb;
    extern logical failed_(void);
    doublereal ub;
    extern doublereal pi_(void);
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), errint_(char *, integer *, 
	    ftnlen);
    extern logical return_(void);
    extern doublereal dpr_(void);

/* $ Abstract */

/*     Regularize a set of longitude intervals. The output */
/*     intervals have their endpoints in ascending order, */
/*     and all output bounds are in the range MINLON:MAXLON. */

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

/*     File: dsktol.inc */


/*     This file contains declarations of tolerance and margin values */
/*     used by the DSK subsystem. */

/*     It is recommended that the default values defined in this file be */
/*     changed only by expert SPICE users. */

/*     The values declared in this file are accessible at run time */
/*     through the routines */

/*        DSKGTL  {DSK, get tolerance value} */
/*        DSKSTL  {DSK, set tolerance value} */

/*     These are entry points of the routine DSKTOL. */

/*        Version 1.0.0 27-FEB-2016 (NJB) */




/*     Parameter declarations */
/*     ====================== */

/*     DSK type 2 plate expansion factor */
/*     --------------------------------- */

/*     The factor XFRACT is used to slightly expand plates read from DSK */
/*     type 2 segments in order to perform ray-plate intercept */
/*     computations. */

/*     This expansion is performed to prevent rays from passing through */
/*     a target object without any intersection being detected. Such */
/*     "false miss" conditions can occur due to round-off errors. */

/*     Plate expansion is done by computing the difference vectors */
/*     between a plate's vertices and the plate's centroid, scaling */
/*     those differences by (1 + XFRACT), then producing new vertices by */
/*     adding the scaled differences to the centroid. This process */
/*     doesn't affect the stored DSK data. */

/*     Plate expansion is also performed when surface points are mapped */
/*     to plates on which they lie, as is done for illumination angle */
/*     computations. */

/*     This parameter is user-adjustable. */


/*     The keyword for setting or retrieving this factor is */


/*     Greedy segment selection factor */
/*     ------------------------------- */

/*     The factor SGREED is used to slightly expand DSK segment */
/*     boundaries in order to select segments to consider for */
/*     ray-surface intercept computations. The effect of this factor is */
/*     to make the multi-segment intercept algorithm consider all */
/*     segments that are sufficiently close to the ray of interest, even */
/*     if the ray misses those segments. */

/*     This expansion is performed to prevent rays from passing through */
/*     a target object without any intersection being detected. Such */
/*     "false miss" conditions can occur due to round-off errors. */

/*     The exact way this parameter is used is dependent on the */
/*     coordinate system of the segment to which it applies, and the DSK */
/*     software implementation. This parameter may be changed in a */
/*     future version of SPICE. */


/*     The keyword for setting or retrieving this factor is */


/*     Segment pad margin */
/*     ------------------ */

/*     The segment pad margin is a scale factor used to determine when a */
/*     point resulting from a ray-surface intercept computation, if */
/*     outside the segment's boundaries, is close enough to the segment */
/*     to be considered a valid result. */

/*     This margin is required in order to make DSK segment padding */
/*     (surface data extending slightly beyond the segment's coordinate */
/*     boundaries) usable: if a ray intersects the pad surface outside */
/*     the segment boundaries, the pad is useless if the intercept is */
/*     automatically rejected. */

/*     However, an excessively large value for this parameter is */
/*     detrimental, since a ray-surface intercept solution found "in" a */
/*     segment can supersede solutions in segments farther from the */
/*     ray's vertex. Solutions found outside of a segment thus can mask */
/*     solutions that are closer to the ray's vertex by as much as the */
/*     value of this margin, when applied to a segment's boundary */
/*     dimensions. */

/*     The keyword for setting or retrieving this factor is */


/*     Surface-point membership margin */
/*     ------------------------------- */

/*     The surface-point membership margin limits the distance */
/*     between a point and a surface to which the point is */
/*     considered to belong. The margin is a scale factor applied */
/*     to the size of the segment containing the surface. */

/*     This margin is used to map surface points to outward */
/*     normal vectors at those points. */

/*     If this margin is set to an excessively small value, */
/*     routines that make use of the surface-point mapping won't */
/*     work properly. */


/*     The keyword for setting or retrieving this factor is */


/*     Angular rounding margin */
/*     ----------------------- */

/*     This margin specifies an amount by which angular values */
/*     may deviate from their proper ranges without a SPICE error */
/*     condition being signaled. */

/*     For example, if an input latitude exceeds pi/2 radians by a */
/*     positive amount less than this margin, the value is treated as */
/*     though it were pi/2 radians. */

/*     Units are radians. */


/*     This parameter is not user-adjustable. */

/*     The keyword for retrieving this parameter is */


/*     Longitude alias margin */
/*     ---------------------- */

/*     This margin specifies an amount by which a longitude */
/*     value can be outside a given longitude range without */
/*     being considered eligible for transformation by */
/*     addition or subtraction of 2*pi radians. */

/*     A longitude value, when compared to the endpoints of */
/*     a longitude interval, will be considered to be equal */
/*     to an endpoint if the value is outside the interval */
/*     differs from that endpoint by a magnitude less than */
/*     the alias margin. */


/*     Units are radians. */


/*     This parameter is not user-adjustable. */

/*     The keyword for retrieving this parameter is */


/*     End of include file dsktol.inc */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     NIVALS     I   is the number of input longitude intervals. */
/*     BOUNDS     I   is an array of bounds of the input intervals. */
/*     MAXN       I   is the maximum number of output intervals. */
/*     NOUT       O   is the number of output intervals. */
/*     MINLON     O   is the minimum longitude of all output intervals. */
/*     MAXLON     O   is the maximum longitude of all output intervals. */
/*     OUTBDS     O   is an array of output longitude intervals. */
/*     SRCS       O   is an array mapping output to input intervals. */

/* $ Detailed_Input */

/*     NIVALS     is the number of input longitude intervals. */

/*     BOUNDS     is an array of upper and lower interval bounds. */
/*                Units are radians. */

/*                The elements */

/*                   BOUNDS(1,I) */
/*                   BOUNDS(2,I) */

/*                are, respectively, the lower and upper bounds of */
/*                the Ith interval. The upper bound may be less */
/*                than the lower bound. If the bounds are equal, they */
/*                are treated as though the upper bound exceeds the */
/*                lower bound by 2*pi. */

/*                 Bounds must be in the range */

/*                   [ -2*pi,  2*pi ] */

/*     MAXN       is the maximum number of output intervals that */
/*                can be returned. */

/* $ Detailed_Output */

/*     NOUT       is the number of output intervals. NOUT is greater */
/*                than or equal to NIVALS. */

/*     MINLON, */
/*     MAXLON     are, respectively, lower and upper bounds on */
/*                the range of the output longitudes. */

/*                   {MINLON, MAXLON} */

/*                is either */

/*                   {0, 2*pi}  or {-pi, pi} */

/*     OUTBDS     is an array of output longitude bounds. */
/*                Units are radians. */

/*                Each output interval represents a subset of some input */
/*                interval; the endpoints of an output interval may be */
/*                shifted by 2*pi relative to the endpoints of the */
/*                subset. */

/*                Input intervals for which the lower bound is greater */
/*                than the upper bound are broken up into two output */
/*                intervals. */

/*                The elements */

/*                   OUTBDS(1,I) */
/*                   OUTBDS(2,I) */

/*                are, respectively, the lower and upper bounds of */
/*                the Ith interval. The upper bound is always */
/*                greater than or equal to the lower bound. */
/*                Bounds are in the range */

/*                   [MINLON, MAXLON] */


/*     SRCS       is an array of indices that map the output intervals */
/*                to the source input intervals from which they were */
/*                derived. */

/* $ Parameters */

/*     ANGMRG     See the description in dsktol.inc. */

/* $ Exceptions */

/*     1)  If the output array doesn't have enough room to store */
/*         the output intervals, the error SPICE(ARRAYTOOSMALL) */
/*         is signaled. */

/*     2)  Longitudes outside the range */

/*            -2*pi - ANGMRG : 2*pi + ANGMRG */

/*         are not accepted: if such a value is encountered, the */
/*         error will be diagnosed by a routine in the call tree */
/*         of this routine. */

/*     3)  If the lower bound of a longitude interval matches the */
/*         upper bound, the error SPICE(ZEROBOUNDSEXTENT) is */
/*         signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine "regularizes" a set of longitude intervals: it maps */
/*     them to a set of output longitude intervals that has the same */
/*     coverage on the unit circle as the input set, such that each of */
/*     the output intervals has its endpoints in increasing order, and */
/*     all of the output intervals have their endpoints in a common */
/*     interval of length 2*pi. */

/*     The set of output intervals has the property that order */
/*     relationships between endpoints are valid. */

/*     This routine supports coverage gap determination. */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/*     DSKBRIEF Version 1.0.0, 21-FEB-2017 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local variables */

    /* Parameter adjustments */
    bounds_dim2 = *nivals;

    /* Function Body */
    if (return_()) {
	return 0;
    }
    chkin_("REGLON", (ftnlen)6);

/*     No output intervals have been found yet. */

    *nout = 0;
    if (*nivals == 0) {
	chkout_("REGLON", (ftnlen)6);
	return 0;
    }

/*     Get lower and upper bounds of input values. */

    *minlon = dpmax_();
    *maxlon = dpmin_();
    i__1 = *nivals;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lb = bounds[(i__2 = (i__ << 1) - 2) < bounds_dim2 << 1 && 0 <= i__2 ? 
		i__2 : s_rnge("bounds", i__2, "reglon_", (ftnlen)256)];
	ub = bounds[(i__2 = (i__ << 1) - 1) < bounds_dim2 << 1 && 0 <= i__2 ? 
		i__2 : s_rnge("bounds", i__2, "reglon_", (ftnlen)257)];

/*        Rectangles of zero longitude extent not allowed. */

	if (lb == ub) {
	    setmsg_("Longitude lower bound # (# degrees) equals upper bound.",
		     (ftnlen)55);
	    errdp_("#", &lb, (ftnlen)1);
	    d__1 = lb * dpr_();
	    errdp_("#", &d__1, (ftnlen)1);
	    sigerr_("SPICE(ZEROBOUNDSEXTENT)", (ftnlen)23);
	    chkout_("REGLON", (ftnlen)6);
	    return 0;
	}

/*        Adjust UB if necessary before deciding on the output */
/*        range. */

	if (ub < lb) {
	    ub += twopi_();
	}
/* Computing MIN */
	d__1 = min(lb,ub);
	*minlon = min(d__1,*minlon);
/* Computing MAX */
	d__1 = max(lb,ub);
	*maxlon = max(d__1,*maxlon);
    }

/*     If MAXLON and MINLON lie within the range */

/*        0 - ANGMRG : 2*pi + ANGMRG */

/*     we'll set the output longitudes to lie in the range */

/*        0 : 2*pi */


/*     If MAXLON and MINLON lie within the range */

/*        -pi - ANGMRG : pi + ANGMRG */

/*     we'll set the output longitudes to lie in the range */

/*        -pi : pi */


/*     We use the latter range if neither of the first two */
/*     conditions are met. */


    if (*minlon > -1e-12 && *maxlon < twopi_() + 1e-12) {
	a = 0.;
	b = twopi_();
    } else {

/*        We arbitrarily pick the output longitude range */

/*           -pi : pi */

	a = -pi_();
	b = pi_();
    }

/*     Set the output values of MINLON and MAXLON. */

    *minlon = a;
    *maxlon = b;

/*     Process each input interval. */

    i__1 = *nivals;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lb = bounds[(i__2 = (i__ << 1) - 2) < bounds_dim2 << 1 && 0 <= i__2 ? 
		i__2 : s_rnge("bounds", i__2, "reglon_", (ftnlen)337)];
	ub = bounds[(i__2 = (i__ << 1) - 1) < bounds_dim2 << 1 && 0 <= i__2 ? 
		i__2 : s_rnge("bounds", i__2, "reglon_", (ftnlen)338)];

/*        We'll adjust the inputs to ensure they're in range. */

/*        First, make sure we're starting with values in */
/*        the range [-2*pi, 2*pi]. */

	zznrmlon_(&lb, &ub, &c_b17, &loclb, &locub);
	if (failed_()) {
	    chkout_("REGLON", (ftnlen)6);
	    return 0;
	}

/*        Move each output into the range [A, B]. */

	if (loclb < a) {
	    loclb += twopi_();
	} else if (loclb > b) {
	    loclb -= twopi_();
	}
	if (locub < a) {
	    locub += twopi_();
	} else if (locub > b) {
	    locub -= twopi_();
	}

/*        Now the bounds are in range, but they may be */
/*        out of order. */

	if (loclb < locub) {

/*           The bounds are in order. Add the interval to */
/*           the list of output intervals. */

	    ++(*nout);
	    if (*nout > *maxn) {

/*              We're out of room. */

		setmsg_("Output arrays have room for # intervals we have fou"
			"nd # output intervals so far.", (ftnlen)80);
		errint_("#", maxn, (ftnlen)1);
		errint_("#", nout, (ftnlen)1);
		sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
		chkout_("REGLON", (ftnlen)6);
		return 0;
	    }
	    outbds[(*nout << 1) - 2] = loclb;
	    outbds[(*nout << 1) - 1] = locub;
	    srcs[*nout - 1] = i__;
	} else {

/*           The bounds are in range but out of order. */
/*           We'll split the input interval into two */
/*           output intervals. */

	    nreq = 0;
	    if (a < locub) {
		++nreq;
	    }
	    if (loclb < b) {
		++nreq;
	    }
	    if (*nout + nreq > *maxn) {

/*              We're out of room. */

		setmsg_("Output arrays have room for # intervals we have fou"
			"nd # output intervals so far.", (ftnlen)80);
		errint_("#", maxn, (ftnlen)1);
		i__2 = *nout + nreq;
		errint_("#", &i__2, (ftnlen)1);
		sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
		chkout_("REGLON", (ftnlen)6);
		return 0;
	    }

/*           The input interval "wraps around" the output boundaries. */

/*           The output intervals extend from A to the upper bound and */
/*           from the lower bound to B. */

	    if (a < locub) {
		++(*nout);
		outbds[(*nout << 1) - 2] = a;
		outbds[(*nout << 1) - 1] = locub;
		srcs[*nout - 1] = i__;
	    }
	    if (loclb < b) {
		++(*nout);
		outbds[(*nout << 1) - 2] = loclb;
		outbds[(*nout << 1) - 1] = b;
		srcs[*nout - 1] = i__;
	    }
	}
    }
    chkout_("REGLON", (ftnlen)6);
    return 0;
} /* reglon_ */

