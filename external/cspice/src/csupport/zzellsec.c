/* zzellsec.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b27 = 1.;
static doublereal c_b28 = 0.;
static logical c_true = TRUE_;
static logical c_false = FALSE_;

/* $Procedure ZZELLSEC ( Tessellate an ellipsoid section with plates ) */
/* Subroutine */ int zzellsec_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *minlon, doublereal *maxlon, doublereal *minlat, 
	doublereal *maxlat, integer *nlon, integer *nlat, integer *maxv, 
	integer *maxp, integer *nv, doublereal *verts, integer *np, integer *
	plates)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    logical ncap;
    doublereal dlat;
    logical scap;
    doublereal dlon;
    extern /* Subroutine */ int vscl_(doublereal *, doublereal *, doublereal *
	    );
    logical wrap;
    extern /* Subroutine */ int zzcapplt_(integer *, logical *, logical *, 
	    integer *, integer *, integer *, integer *), zzgrdplt_(integer *, 
	    integer *, logical *, integer *, integer *);
    integer i__, j, n;
    doublereal s;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer latlb;
    extern /* Subroutine */ int vpack_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);
    doublereal level;
    integer latub;
    extern /* Subroutine */ int errdp_(char *, doublereal *, ftnlen);
    integer ncols, nrows;
    extern logical failed_(void);
    extern doublereal pi_(void);
    extern /* Subroutine */ int latrec_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), sigerr_(char *, ftnlen), chkout_(
	    char *, ftnlen);
    integer polidx;
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen);
    doublereal lmxlon;
    extern logical return_(void);
    doublereal dir[3], lat;
    integer bix;
    doublereal lon;
    integer nnp, nsp, pix, vix;

/* $ Abstract */

/*     Create a set of triangular plates covering a specified section */
/*     of the surface of a triaxial ellipsoid. The boundaries of the */
/*     section are curves of constant planetocentric longitude and */
/*     latitude. */

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
/*     MINLON     I   Minimum longitude of section. */
/*     MAXLON     I   Minimum longitude of section. */
/*     MINLAT     I   Minimum latitude of section. */
/*     MAXLAT     I   Minimum latitude of section. */
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


/*     MINLON, */
/*     MAXLON, */
/*     MINLAT, */
/*     MAXLAT     are, respectively, the longitude and latitude bounds */
/*                of a section of the surface the triaxial ellipsoid. */
/*                The coordinate system is latitudinal. Units are */
/*                radians. */


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


/*     MAXV       is the maximum number of vertices to return. The */
/*                number of vertices created depends on which polar */
/*                caps are created. In all cases the number does not */
/*                exceed */

/*                   ( NLON + 1 ) * ( NLAT + 1 ) */

/*                The array VERTS must have size at least 3*MAXV. */


/*     MAXP       is the maximum number of plates to return. The */
/*                number of plates created depends on which polar */
/*                caps are created and whether longitude wrapping */
/*                is selected. In all cases the number does not */
/*                exceed */

/*                   2 * NLON * NLAT */

/*                The array PLATES must have size at least 3*MAXP. */

/* $ Detailed_Output */

/*     NV         is the number of vertices in the output array VERTS. */

/*     VERTS      is an array containing the vertices of the output */
/*                plate set. There is a vertex at each intersection of a */
/*                latitude band boundary and a longitude band boundary. */
/*                The vertices at the north and south poles are at */
/*                indices NV and NV-1, respectively, if both polar caps */
/*                are created. The non-polar vertex indices start at 1. */
/*                Non-polar vertices are indexed in top-down, */
/*                left-to-right order, with vertices of each latitude */
/*                band stored contiguously. */

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

/* -    TESTUTIL Version 1.0.0, 29-SEP-2014 (NJB) */

/*        Based on ZZELLPLT Version 1.0.0, 30-APR-2014 (NJB) */

/* -& */
/* $ Index_Entries */

/*     tessellate ellipsoid */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("ZZELLSEC", (ftnlen)8);

/*     The semi-axes must have positive length. */

    if (*a <= 0. || *b <= 0. || *c__ <= 0.) {
	setmsg_("Semi-axis lengths:  A = #, B = #, C = #. ", (ftnlen)41);
	errdp_("#", a, (ftnlen)1);
	errdp_("#", b, (ftnlen)1);
	errdp_("#", c__, (ftnlen)1);
	sigerr_("SPICE(INVALIDAXISLENGTH)", (ftnlen)24);
	chkout_("ZZELLSEC", (ftnlen)8);
	return 0;
    }

/*     The longitude and latitude band counts must be realizable. */

    if (*nlat < 2) {
	setmsg_("The latitude band count must be at least 2 but was #.", (
		ftnlen)53);
	errint_("#", nlat, (ftnlen)1);
	sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
	chkout_("ZZELLSEC", (ftnlen)8);
	return 0;
    }
    if (*nlon < 3) {
	setmsg_("The longitude band count must be at least 3 but was #.", (
		ftnlen)54);
	errint_("#", nlon, (ftnlen)1);
	sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
	chkout_("ZZELLSEC", (ftnlen)8);
	return 0;
    }

/*     Decide whether we have two distinct longitude boundaries. First */
/*     create a local maximum longitude that's greater than the minimum */
/*     longitude. The logical variable WRAP is .TRUE. if and only if we */
/*     have 2*pi - TOL radians of longitude coverage, where TOL is a */
/*     small value. */

    if (*maxlon > *minlon) {
	lmxlon = *maxlon;
    } else {
	lmxlon = *maxlon + pi_() * 2;
    }
    wrap = lmxlon - *minlon > pi_() * 2 - 1e-12;

/*     Decide whether we have north or south polar caps. */

    ncap = *maxlat > pi_() / 2 - 1e-12;
    scap = *minlat < -pi_() / 2 + 1e-12;

/*     Compute the vertex counts. */

    if (wrap) {

/*        Vertex count:  When both caps are present, there are NLAT-2 */
/*                       latitude bands, excluding the polar caps. These */
/*                       are bounded by NLAT-1 rows of vertices. Each */
/*                       row of vertices has NLON members. The caps add */
/*                       two vertices. */

	if (ncap && scap) {

/*           There are two polar caps. */

	    *nv = *nlon * (*nlat - 1) + 2;
	} else if (ncap || scap) {

/*           There's just one polar cap. Excluding the polar */
/*           vertex, there are NLAT rows of vertices. */

	    *nv = *nlon * *nlat + 1;
	} else {

/*           No polar caps. There are NLAT+1 rows of vertices. */

	    *nv = *nlon * (*nlat + 1);
	}
    } else {
	if (ncap && scap) {

/*           There are two polar caps. */

	    *nv = (*nlon + 1) * (*nlat - 1) + 2;
	} else if (ncap || scap) {

/*           There's just one polar cap. Excluding the polar */
/*           vertex, there are NLAT rows of vertices. */

	    *nv = (*nlon + 1) * *nlat + 1;
	} else {

/*           No polar caps. There are NLAT+1 rows of vertices. */

	    *nv = (*nlon + 1) * (*nlat + 1);
	}
    }

/*     Compute the plate counts. These depend on the set of */
/*     polar caps. */


/*        Plate count:   each latitude band, excluding the polar caps, */
/*                       contains 2*NLON plates. Each cap contains NLON */
/*                       plates. */

    if (ncap && scap) {

/*        There are two polar caps. */

	*np = (*nlon << 1) * (*nlat - 1);
    } else if (ncap || scap) {

/*        There's just one polar cap. Excluding the polar */
/*        vertex, there are NLAT rows of vertices. */

	*np = *nlon * ((*nlat << 1) - 1);
    } else {

/*        No polar caps. There are NLAT+1 rows of vertices. */

	*np = (*nlon << 1) * *nlat;
    }
    if (*nv > *maxv) {
	setmsg_("The requested plate model requires # vertices but the maxim"
		"um vertex count is #.", (ftnlen)80);
	errint_("#", nv, (ftnlen)1);
	errint_("#", maxv, (ftnlen)1);
	sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
	chkout_("ZZELLSEC", (ftnlen)8);
	return 0;
    }
    if (*np > *maxp) {
	setmsg_("The requested plate model requires # plates but the maximum"
		" plate count is #.", (ftnlen)77);
	errint_("#", np, (ftnlen)1);
	errint_("#", maxp, (ftnlen)1);
	sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
	chkout_("ZZELLSEC", (ftnlen)8);
	return 0;
    }

/*     Create the vertex set, excluding the polar caps. */

/*     LATLB will be the index of the first vertex row, excluding */
/*     the caps. LATUB will be the index of the last vertex row, */
/*     excluding the caps. */

    if (ncap && scap) {
	latlb = 2;
	latub = *nlat;
    } else if (ncap) {
	latlb = 2;
	latub = *nlat + 1;
    } else if (scap) {
	latlb = 1;
	latub = *nlat;
    } else {
	latlb = 1;
	latub = *nlat + 1;
    }

/*     NCOLS is the number of columns of vertices. */

    if (wrap) {
	ncols = *nlon;
    } else {
	ncols = *nlon + 1;
    }

/*     The latitude bands are equally spaced in planetocentric latitude. */

    dlat = (*maxlat - *minlat) / *nlat;
    dlon = (lmxlon - *minlon) / *nlon;
    vix = 1;
    i__1 = latub;
    for (i__ = latlb; i__ <= i__1; ++i__) {
	lat = *maxlat - (i__ - 1) * dlat;
	i__2 = ncols;
	for (j = 1; j <= i__2; ++j) {
	    lon = *minlon + (j - 1) * dlon;

/*           Create a unit direction vector for the current */
/*           vertex. Scale this vector to make it lie on the */
/*           ellipsoid's surface; the scaled vector is the */
/*           current vertex. */

	    latrec_(&c_b27, &lon, &lat, dir);
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

/*     Create the polar vertices if necessary. */

    if (ncap && scap) {
	vpack_(&c_b28, &c_b28, c__, &verts[(*nv - 1) * 3 - 3]);
	d__1 = -(*c__);
	vpack_(&c_b28, &c_b28, &d__1, &verts[*nv * 3 - 3]);
    } else if (ncap) {
	vpack_(&c_b28, &c_b28, c__, &verts[*nv * 3 - 3]);
    } else if (scap) {
	d__1 = -(*c__);
	vpack_(&c_b28, &c_b28, &d__1, &verts[*nv * 3 - 3]);
    }

/*     Create the plates for the latitude bounds other than */
/*     those belonging to the caps. */

/*     The first two inputs are the vertex row and column counts. */
/*     Next is a logical flag indicating whether longitude wrapping */
/*     should be used. */

    if (latub > latlb) {
	nrows = latub - latlb + 1;
	zzgrdplt_(&nrows, &ncols, &wrap, &n, plates);
	if (failed_()) {
	    chkout_("ZZELLSEC", (ftnlen)8);
	    return 0;
	}
    } else {
	n = 0;
    }
    if (ncap) {

/*        Add the north cap. This is a set of plates; the vertices */
/*        already have been computed. */

/*        PIX is the index of the first cap plate. BIX is the */
/*        base (predecessor) index of the first vertex in the */
/*        first vertex row. */

	pix = n + 1;
	bix = 0;

/*        POLIDX is the vertex index of the north polar vertex. */

	if (scap) {
	    polidx = *nv - 1;
	} else {
	    polidx = *nv;
	}
	zzcapplt_(&ncols, &c_true, &wrap, &bix, &polidx, &nnp, &plates[pix * 
		3 - 3]);
	if (failed_()) {
	    chkout_("ZZELLSEC", (ftnlen)8);
	    return 0;
	}
    }
    if (scap) {

/*        Add the south cap. */

	polidx = *nv;
	if (ncap) {
	    pix += nnp;
	    bix = *nv - (ncols + 2);
	} else {
	    pix = n + 1;
	    bix = *nv - (ncols + 1);
	}
	zzcapplt_(&ncols, &c_false, &wrap, &bix, &polidx, &nsp, &plates[pix * 
		3 - 3]);
    }
    chkout_("ZZELLSEC", (ftnlen)8);
    return 0;
} /* zzellsec_ */

