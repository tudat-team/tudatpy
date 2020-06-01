/* zzpsxtnt.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c_b6 = 16000002;
static integer c_b12 = 32000000;

/* $Procedure    ZZPSXTNT ( MKDSK, compute plate set extents ) */
/* Subroutine */ int zzpsxtnt_(integer *nv, doublereal *verts, integer *np, 
	integer *plates, doublereal *extent, doublereal *avplex)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    integer i__, j, k;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    extern doublereal dpmin_(void), dpmax_(void);
    doublereal pltbds[2];
    extern /* Subroutine */ int sigerr_(char *, ftnlen), setmsg_(char *, 
	    ftnlen), errint_(char *, integer *, ftnlen), chkout_(char *, 
	    ftnlen);
    extern logical return_(void);
    doublereal sum;

/* $ Abstract */

/*     Compute plate set extents and average plate extent. */

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

/*     DSK */
/*     FILES */

/* $ Declarations */

/*     Include file mkdsk02.inc */

/*     This include file declares parameters for DSK data type 2 */
/*     (plate model). On most platforms, these parameters are */
/*     indentical to those in the SPICELIB include file */

/*        dsk02.inc */

/*     However, on some platforms, the default parameters result */
/*     in excessive memory usage. For these platforms, the maximum */
/*     plate and vertex counts of type 2 segments created by MKDSK */
/*     have been reduced. */

/*     When support for the problematic platforms is discontinued, */
/*     references to this file may be replaced with references to */
/*     dsk02.inc. */


/* -       SPICELIB Version 1.0.0 17-FEB-2017 (NJB) */

/*          Based on SPICELIB include file dsk02.inc version */
/*          1.0.0 08-FEB-2017 (NJB) */






/*     DSK type 2 segment layout */
/*     ========================= */


/*     Each type 2 DSK segment has integer, d.p., and character */
/*     components.  The segment layout in DAS address space is as */
/*     follows: */


/*        Integer layout: */

/*           +-----------------+ */
/*           | NV              |  (# of vertices) */
/*           +-----------------+ */
/*           | NP              |  (# of plates ) */
/*           +-----------------+ */
/*           | NVXTOT          |  (total number of voxels) */
/*           +-----------------+ */
/*           | VGREXT          |  (voxel grid extents, 3 integers) */
/*           +-----------------+ */
/*           | CGRSCL          |  (coarse voxel grid scale, 1 integer) */
/*           +-----------------+ */
/*           | VOXNPT          |  (size of voxel-plate pointer list) */
/*           +-----------------+ */
/*           | VOXNPL          |  (size of voxel-plate list) */
/*           +-----------------+ */
/*           | VTXNPL          |  (size of vertex-plate list) */
/*           +-----------------+ */
/*           | PLATES          |  (NP 3-tuples of vertex IDs) */
/*           +-----------------+ */
/*           | VOXPTR          |  (voxel-plate pointer array) */
/*           +-----------------+ */
/*           | VOXPLT          |  (voxel-plate list) */
/*           +-----------------+ */
/*           | VTXPTR          |  (vertex-plate pointer array) */
/*           +-----------------+ */
/*           | VTXPLT          |  (vertex-plate list) */
/*           +-----------------+ */
/*           | CGRPTR          |  (coarse grid occupancy pointers) */
/*           +-----------------+ */



/*        D.p. layout: */

/*           +-----------------+ */
/*           | DSK descriptor  |  DSKDSZ elements */
/*           +-----------------+ */
/*           | Vertex bounds   |  6 values (min/max for each component) */
/*           +-----------------+ */
/*           | Voxel origin    |  3 elements */
/*           +-----------------+ */
/*           | Voxel size      |  1 element */
/*           +-----------------+ */
/*           | Vertices        |  3*NV elements */
/*           +-----------------+ */


/*     This local parameter MUST be kept consistent with */
/*     the parameter DSKDSZ which is declared in dskdsc.inc. */


/*     Integer item keyword parameters used by fetch routines: */


/*     Double precision item keyword parameters used by fetch routines: */


/*     The parameters below formerly were declared in pltmax.inc. */

/*     Limits on plate model capacity: */

/*     The maximum number of bodies, vertices and */
/*     plates in a plate model or collective thereof are */
/*     provided here. */

/*     These values can be used to dimension arrays, or to */
/*     use as limit checks. */

/*     The value of MAXPLT is determined from MAXVRT via */
/*     Euler's Formula for simple polyhedra having triangular */
/*     faces. */

/*     MAXVRT is the maximum number of vertices the triangular */
/*            plate model software will support. */


/*     MAXPLT is the maximum number of plates that the triangular */
/*            plate model software will support. */


/*     MAXNPV is the maximum allowed number of vertices, not taking into */
/*     account shared vertices. */

/*     Note that this value is not sufficient to create a vertex-plate */
/*     mapping for a model of maximum plate count. */


/*     MAXVOX is the maximum number of voxels. */


/*     MAXCGR is the maximum size of the coarse voxel grid. */


/*     MAXEDG is the maximum allowed number of vertex or plate */
/*     neighbors a vertex may have. */

/*     DSK type 2 spatial index parameters */
/*     =================================== */

/*        DSK type 2 spatial index integer component */
/*        ------------------------------------------ */

/*           +-----------------+ */
/*           | VGREXT          |  (voxel grid extents, 3 integers) */
/*           +-----------------+ */
/*           | CGRSCL          |  (coarse voxel grid scale, 1 integer) */
/*           +-----------------+ */
/*           | VOXNPT          |  (size of voxel-plate pointer list) */
/*           +-----------------+ */
/*           | VOXNPL          |  (size of voxel-plate list) */
/*           +-----------------+ */
/*           | VTXNPL          |  (size of vertex-plate list) */
/*           +-----------------+ */
/*           | CGRPTR          |  (coarse grid occupancy pointers) */
/*           +-----------------+ */
/*           | VOXPTR          |  (voxel-plate pointer array) */
/*           +-----------------+ */
/*           | VOXPLT          |  (voxel-plate list) */
/*           +-----------------+ */
/*           | VTXPTR          |  (vertex-plate pointer array) */
/*           +-----------------+ */
/*           | VTXPLT          |  (vertex-plate list) */
/*           +-----------------+ */


/*        Index parameters */


/*     Grid extent: */


/*     Coarse grid scale: */


/*     Voxel pointer count: */


/*     Voxel-plate list count: */


/*     Vertex-plate list count: */


/*     Coarse grid pointers: */


/*     Size of fixed-size portion of integer component: */


/*        DSK type 2 spatial index double precision component */
/*        --------------------------------------------------- */

/*           +-----------------+ */
/*           | Vertex bounds   |  6 values (min/max for each component) */
/*           +-----------------+ */
/*           | Voxel origin    |  3 elements */
/*           +-----------------+ */
/*           | Voxel size      |  1 element */
/*           +-----------------+ */



/*        Index parameters */

/*     Vertex bounds: */


/*     Voxel grid origin: */


/*     Voxel size: */


/*     Size of fixed-size portion of double precision component: */


/*     The limits below are used to define a suggested maximum */
/*     size for the integer component of the spatial index. */


/*     Maximum number of entries in voxel-plate pointer array: */


/*     Maximum cell size: */


/*     Maximum number of entries in voxel-plate list: */


/*     Spatial index integer component size: */


/*     End of include file mkdsk02.inc */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     NV         I   Vertex count. */
/*     VERTS      I   Vertices. */
/*     NP         I   Plate count. */
/*     PLATES     I   Plates. */
/*     EXTENT     O   Extents of the vertex set. */
/*     AVPLEX     O   Average plate extent. */

/* $ Detailed_Input */

/*     NV             is the number of vertices in the input vertex */
/*                    array. */

/*     VERTS          is an array of vertices associated with plates. */

/*     NP             is the number of plates in the input plate array. */

/*     PLATES         is an array of triangular plates. The element */

/*                       PLATES(J,I) */

/*                    is the index in the VERTS array of the Jth vertex */
/*                    of the Ith plate. The vertex indices are 1-based. */

/* $ Detailed_Output */

/*     EXTENT         is the extent of the input vertex set. See */
/*                    Particulars. */

/*     AVPLEX         is the average extent of the plates in the */
/*                    plate set. See Particulars. */

/* $ Parameters */

/*     See mkdsk02.inc */

/* $ Exceptions */

/*     This routine is meant to be operated in RETURN SPICE error */
/*     handling mode. The caller is expected to delete the DSK file if */
/*     an error occurs during file creation. */


/*     1)  If the vertex count is less than 3 or greater than MAXVRT, */
/*         the error SPICE(VALUEOUTOFRANGE) is signaled. */

/*     2)  If the plate count is less than 1 or greater than MAXPLT, */
/*         the error SPICE(VALUEOUTOFRANGE) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     The "extent" of a set of vertices is the set of extrema */
/*     of the coordinates of those vertices. */

/*     The extent of a plate along a coordinate axis is the length of */
/*     the orthogonal projection of the plate onto that axis. */

/*     The average "extent" of a set of plates is the average of the */
/*     extents of the plates along the coordinate axes. */

/* $ Examples */

/*     See usage in MKDSK. */

/* $ Restrictions */

/*     This routine should be called only from within MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 1.0.0, 18-FEB-2017 (NJB) */

/* -& */
/* $ Index_Entries */

/*     determine vertex set extents and average plate extent */

/* -& */
/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("ZZPSXTNT", (ftnlen)8);

/*     Reject invalid plate and vertex counts. */

    if (*nv < 3 || *nv > 16000002) {
	setmsg_("Vertex count NV = #; count must be in the range 3:#.", (
		ftnlen)52);
	errint_("#", nv, (ftnlen)1);
	errint_("#", &c_b6, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("ZZPSXTNT", (ftnlen)8);
	return 0;
    }
    if (*np < 1 || *np > 32000000) {
	setmsg_("Plate count NP = #; count must be in the range 1:#.", (
		ftnlen)51);
	errint_("#", np, (ftnlen)1);
	errint_("#", &c_b12, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("ZZPSXTNT", (ftnlen)8);
	return 0;
    }

/*     Compute extents of the plate set. These depend only on */
/*     the vertex set. */

    for (i__ = 1; i__ <= 3; ++i__) {
	extent[(i__1 = (i__ << 1) - 2) < 6 && 0 <= i__1 ? i__1 : s_rnge("ext"
		"ent", i__1, "zzpsxtnt_", (ftnlen)216)] = dpmax_();
	extent[(i__1 = (i__ << 1) - 1) < 6 && 0 <= i__1 ? i__1 : s_rnge("ext"
		"ent", i__1, "zzpsxtnt_", (ftnlen)217)] = dpmin_();
    }
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* Computing MIN */
	    d__1 = extent[(i__3 = (j << 1) - 2) < 6 && 0 <= i__3 ? i__3 : 
		    s_rnge("extent", i__3, "zzpsxtnt_", (ftnlen)226)], d__2 = 
		    verts[j + i__ * 3 - 4];
	    extent[(i__2 = (j << 1) - 2) < 6 && 0 <= i__2 ? i__2 : s_rnge(
		    "extent", i__2, "zzpsxtnt_", (ftnlen)226)] = min(d__1,
		    d__2);
/* Computing MAX */
	    d__1 = extent[(i__3 = (j << 1) - 1) < 6 && 0 <= i__3 ? i__3 : 
		    s_rnge("extent", i__3, "zzpsxtnt_", (ftnlen)227)], d__2 = 
		    verts[j + i__ * 3 - 4];
	    extent[(i__2 = (j << 1) - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge(
		    "extent", i__2, "zzpsxtnt_", (ftnlen)227)] = max(d__1,
		    d__2);
	}
    }

/*     Compute plate extents in all 3 dimensions; compute */
/*     the average of all extents. */

    sum = 0.;
    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        I is the current plate index. */

	for (j = 1; j <= 3; ++j) {

/*           J is the current coordinate index. */

	    pltbds[0] = dpmax_();
	    pltbds[1] = dpmin_();
	    for (k = 1; k <= 3; ++k) {

/*              K is the current vertex index for the current plate. */
/*              Account for this vertex's contributions to the plate's */
/*              extent in the direction of coordinate J. */

/* Computing MIN */
		d__1 = pltbds[0], d__2 = verts[j + plates[k + i__ * 3 - 4] * 
			3 - 4];
		pltbds[0] = min(d__1,d__2);
/* Computing MAX */
		d__1 = pltbds[1], d__2 = verts[j + plates[k + i__ * 3 - 4] * 
			3 - 4];
		pltbds[1] = max(d__1,d__2);
	    }

/*           Add the current plate's extent in the Jth coordinate. */

	    sum += pltbds[1] - pltbds[0];
	}
    }

/*     We've added three extents per plate, so we have 3*NP in all. */
/*     Note we've ensured NP is positive. */

    *avplex = sum / (*np * 3);
    chkout_("ZZPSXTNT", (ftnlen)8);
    return 0;
} /* zzpsxtnt_ */

