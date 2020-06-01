/* zzellplt.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b27 = 0.;
static doublereal c_b31 = 1.;
static logical c_true = TRUE_;
static logical c_false = FALSE_;

/* $Procedure ZZELLPLT ( Tessellate an ellipsoid with triangular plates ) */
/* Subroutine */ int zzellplt_(doublereal *a, doublereal *b, doublereal *c__, 
	integer *nlon, integer *nlat, integer *maxv, integer *maxp, integer *
	nv, doublereal *verts, integer *np, integer *plates)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal dlat, dlon;
    extern /* Subroutine */ int vscl_(doublereal *, doublereal *, doublereal *
	    ), zzcapplt_(integer *, logical *, logical *, integer *, integer *
	    , integer *, integer *), zzgrdplt_(integer *, integer *, logical *
	    , integer *, integer *);
    integer i__, j, n;
    doublereal s;
    extern /* Subroutine */ int chkin_(char *, ftnlen), vpack_(doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    doublereal level;
    extern /* Subroutine */ int errdp_(char *, doublereal *, ftnlen);
    extern logical failed_(void);
    extern doublereal pi_(void);
    extern /* Subroutine */ int latrec_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), sigerr_(char *, ftnlen), chkout_(
	    char *, ftnlen), setmsg_(char *, ftnlen), errint_(char *, integer 
	    *, ftnlen);
    extern logical return_(void);
    doublereal dir[3], lat;
    integer bix;
    doublereal lon;
    integer nnp, nsp, pix, vix;

/* $ Abstract */

/*     Create a set of triangular plates covering a specified triaxial */
/*     ellipsoid. */

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

/*     DSK */

/* $ Keywords */

/*     ELLIPSOID */
/*     PLATE */
/*     TILE */
/*     TESSELLATE */

/* $ Declarations */
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     A          I   Length of ellipsoid semi-axis lying on the x-axis. */
/*     B          I   Length of ellipsoid semi-axis lying on the y-axis. */
/*     C          I   Length of ellipsoid semi-axis lying on the z-axis. */
/*     NLON       I   Number of longitude bands in plate set. */
/*     NLAT       I   Number of latitude bands in plate set. */
/*     MAXV       I   Maximum number of vertices to return. */
/*     MAXP       I   Maximum number of plates to return. */
/*     NV         O   Number of vertices in output array. */
/*     VERTS      O   Vertices. */
/*     NP         O   Number of plates in output array. */
/*     PLATES     O   Plates. */

/* $ Detailed_Input */

/*     A, */
/*     B, */
/*     C          are the lengths of the semi-axes of a triaxial */
/*                ellipsoid. The ellipsoid is centered at the origin and */
/*                oriented so that its axes lie on the x, y and z axes. */
/*                A, B, and C are the lengths of the semi-axes that */
/*                point in the x, y, and z directions respectively. */


/*     NLON       is the number of longitude bands in the output plate */
/*                set. Each longitude band is bounded by two meridians. */
/*                All longitude bands have equal angular extent in */
/*                longitude. The vertices of any plate lie on adjacent */
/*                longitude band boundaries. */


/*     NLAT       is the number of latitude bands in the output plate */
/*                set. The vertices of each band are bounded by two */
/*                cones of constant planetocentric latitude. All */
/*                latitude bands have equal angular extent in */
/*                planetocentric latitude. The vertices of any plate lie */
/*                on adjacent latitude band boundaries. */

/*                Each polar "cap" consists of one longitude band. */


/*     MAXV       is the maximum number of vertices to return. MAXV must */
/*                be at least */

/*                   ( NLON * ( NLAT - 1 ) )  +  2 */

/*                The array VERTS must have size at least 3*MAXV. */


/*     MAXP       is the maximum number of plates to return. MAXP must */
/*                be at least */

/*                   2 * NLON * ( NLAT - 1 ) */

/*                The array PLATES must have size at least 3*MAXP. */

/* $ Detailed_Output */

/*     NV         is the number of vertices in the output array VERTS. */

/*     VERTS      is an array containing the vertices of the output */
/*                plate set. There is a vertex at each intersection of a */
/*                latitude band boundary and a longitude band boundary. */
/*                The vertices at the north and south poles are at */
/*                indices NV and NV-1, respectively. The non-polar */
/*                vertex indices start at 1. Non-polar vertices are */
/*                indexed in top-down, left-to-right order, with */
/*                vertices of each latitude band stored contiguously. */

/*     NP         is the number of plates in the output array PLATES. */


/*     PLATES     is an array containing the tessellating plate set. */
/*                Each plate is an array of three vertex indices, where */
/*                the indices range from 1 to NV. */

/*                The vertices of any plate are ordered in the */
/*                right-handed sense about the outward normal direction */
/*                for that plate. */

/*                The non-polar plates---those not having a vertex at */
/*                the north or south pole---are indexed in top-down, */
/*                left-to-right order, with plates belonging to each */
/*                latitude band stored contiguously. Plates constituting */
/*                the north polar cap follow the non-polar plate set, */
/*                and plates constituting the south polar cap follow */
/*                those. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the length of any semi-axis of the ellipsoid is */
/*         non-positive, the error SPICE(INVALIDAXISLENGTH) is signaled. */

/*     2)  If NLAT is less than 2, the error SPICE(INVALIDCOUNT) is */
/*         signaled. */

/*     3)  If NLON is less than 3, the error SPICE(INVALIDCOUNT) is */
/*         signaled. */

/*     4)  If the number of vertices implied by the input values NLON */
/*         and NLAT exceeds MAXV, the error SPICE(ARRAYTOOSMALL) is */
/*         signaled. */

/*     5)  If the number of plates implied by the input values NLON */
/*         and NLAT exceeds MAXP, the error SPICE(ARRAYTOOSMALL) is */
/*         signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     The vertex and plate sets created by this routine are suitable */
/*     for use in a type 2 DSK segment. */

/*     While the primary purpose of this routine is to support testing, */
/*     there may be some user applications for which a tessellated plate */
/*     model is valuable. For example, computing an estimate of the area */
/*     of a specified surface region may be simplified by using a plate */
/*     model. */

/*     Note that, for ellipsoids having three distinct radii, the Z */
/*     components of the vertices on any latitude band boundary (except */
/*     the poles themselves) will vary with longitude. */

/*     Also note that the horizontal edge of a plate may extend beyond */
/*     the boundaries of the latitude band containing the plate. */

/* $ Examples */

/*     See use of C language version in tspice_c test families. */

/* $ Restrictions */

/*     1) For use only by TSPICE. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    TESTUTIL Version 1.0.0, 30-APR-2014 (NJB) */

/* -& */
/* $ Index_Entries */

/*     tessellate ellipsoid */

/* -& */

/*     SPICELIB functions */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("ZZELLPLT", (ftnlen)8);

/*     The semi-axes must have positive length. */

    if (*a <= 0. || *b <= 0. || *c__ <= 0.) {
	setmsg_("Semi-axis lengths:  A = #, B = #, C = #. ", (ftnlen)41);
	errdp_("#", a, (ftnlen)1);
	errdp_("#", b, (ftnlen)1);
	errdp_("#", c__, (ftnlen)1);
	sigerr_("SPICE(INVALIDAXISLENGTH)", (ftnlen)24);
	chkout_("ZZELLPLT", (ftnlen)8);
	return 0;
    }

/*     The longitude and latitude band counts must be realizable. */

    if (*nlat < 2) {
	setmsg_("The latitude band count must be at least 2 but was #.", (
		ftnlen)53);
	errint_("#", nlat, (ftnlen)1);
	sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
	chkout_("ZZELLPLT", (ftnlen)8);
	return 0;
    }
    if (*nlon < 3) {
	setmsg_("The longitude band count must be at least 3 but was #.", (
		ftnlen)54);
	errint_("#", nlon, (ftnlen)1);
	sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
	chkout_("ZZELLPLT", (ftnlen)8);
	return 0;
    }

/*     Compute the vertex and plate counts. Check against available */
/*     room. */

/*        Vertex count: there are NLAT-2 latitude bands, excluding */
/*                      the polar caps. These are bounded by NLAT-1 rows */
/*                      of vertices. Each row of vertices has NLON */
/*                      members. The caps add two vertices. */

/*        Plate count:  each latitude band, excluding the polar caps, */
/*                      contains 2*NLON plates. Each cap contains NLON */
/*                      plates. */


    *nv = *nlon * (*nlat - 1) + 2;
    *np = (*nlon << 1) * (*nlat - 1);
    if (*nv > *maxv) {
	setmsg_("The requested plate model requires # vertices but the maxim"
		"um vertex count is #.", (ftnlen)80);
	errint_("#", nv, (ftnlen)1);
	errint_("#", maxv, (ftnlen)1);
	sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
	chkout_("ZZELLPLT", (ftnlen)8);
	return 0;
    }
    if (*np > *maxp) {
	setmsg_("The requested plate model requires # plates but the maximum"
		" plate count is #.", (ftnlen)77);
	errint_("#", np, (ftnlen)1);
	errint_("#", maxp, (ftnlen)1);
	sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
	chkout_("ZZELLPLT", (ftnlen)8);
	return 0;
    }

/*     Create the vertex set. The north polar vertex is */
/*     at index 1; the south vertex is at index NV. It will */
/*     be convenient to make these the last two vertices. */

    vpack_(&c_b27, &c_b27, c__, &verts[(*nv - 1) * 3 - 3]);
    d__1 = -(*c__);
    vpack_(&c_b27, &c_b27, &d__1, &verts[*nv * 3 - 3]);

/*     The latitude bands are equally spaced in planetocentric */
/*     latitude. */

    dlat = pi_() / *nlat;
    dlon = pi_() * 2 / *nlon;
    vix = 1;
    i__1 = *nlat - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lat = pi_() / 2 - i__ * dlat;
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    lon = (j - 1) * dlon;

/*           Create a unit direction vector for the current */
/*           vertex. Scale this vector to make it lie on the */
/*           ellipsoid's surface; the scaled vector is the */
/*           current vertex. */

	    latrec_(&c_b31, &lon, &lat, dir);
/* Computing 2nd power */
	    d__1 = dir[0] / *a;
/* Computing 2nd power */
	    d__2 = dir[1] / *b;
/* Computing 2nd power */
	    d__3 = dir[2] / *c__;
	    level = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    s = 1. / sqrt(level);
	    vscl_(&s, dir, &verts[vix * 3 - 3]);

/*           Next vertex. */

	    ++vix;
	}
    }

/*     Create the plates for the latitude bounds other than */
/*     those belonging to the caps. */

/*     The first two inputs are the vertex row and column counts. */
/*     Next is a logical flag indicating whether longitude wrapping */
/*     should be used. */

    if (*nlat > 2) {
	i__1 = *nlat - 1;
	zzgrdplt_(&i__1, nlon, &c_true, &n, plates);
	if (failed_()) {
	    chkout_("ZZELLPLT", (ftnlen)8);
	    return 0;
	}
    }

/*     Add the north cap. This is a set of plates; the vertices */
/*     already have been computed. */

    pix = *np - (*nlon << 1) + 1;
    bix = 0;
    i__1 = *nv - 1;
    zzcapplt_(nlon, &c_true, &c_true, &bix, &i__1, &nnp, &plates[pix * 3 - 3])
	    ;
    if (failed_()) {
	chkout_("ZZELLPLT", (ftnlen)8);
	return 0;
    }

/*     Add the south cap. */

    pix += *nlon;
    bix = *nv - (*nlon + 2);
    zzcapplt_(nlon, &c_false, &c_true, &bix, nv, &nsp, &plates[pix * 3 - 3]);
    chkout_("ZZELLPLT", (ftnlen)8);
    return 0;
} /* zzellplt_ */

