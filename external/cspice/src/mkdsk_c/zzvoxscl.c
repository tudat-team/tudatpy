/* zzvoxscl.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c_b6 = 100000;
static doublereal c_b18 = .33333333333333331;

/* $Procedure    ZZVOXSCL ( MKDSK, compute voxel scales ) */
/* Subroutine */ int zzvoxscl_(doublereal *extent, doublereal *avplex, 
	integer *trgcor, integer *trgfin, integer *corscl, doublereal *finscl)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    doublereal dpcs, q;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer maxnf;
    extern /* Subroutine */ int errdp_(char *, doublereal *, ftnlen);
    doublereal finsz, corsz, corsz0, lx, ly, lz;
    integer nx, ny;
    doublereal sx, sy, sz;
    integer nz, ncorse;
    extern /* Subroutine */ int sigerr_(char *, ftnlen);
    doublereal volume;
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), chkout_(char *, ftnlen);
    extern logical return_(void);

/* $ Abstract */

/*     Compute DSK type 2 coarse and fine voxel scales from plate set */
/*     extents and target voxel counts. */

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
/*     EXTENT     I   Extents of the vertex set. */
/*     AVPLEX     I   Average plate extent. */
/*     TRGCOR     I   Target coarse voxel count. */
/*     TRGFIN     I   Target fine voxel count. */
/*     CORSCL     O   Coarse voxel scale. */
/*     FINSCL     O   Fine voxel scale. */

/* $ Detailed_Input */

/*     EXTENT         is the extent of the input vertex set. Units */
/*                    are km. See Particulars. */

/*     AVPLEX         is the average extent of the plates in the */
/*                    plate set. Units are km. See Particulars. */

/*     TRGCOR         is the target coarse voxel count. This is an */
/*                    approximation to the final count; the actual count */
/*                    will be smaller than, but near, TRGCOR. */

/*     TRGFIN         is the target fine voxel count. This is an */
/*                    approximation to the final count; the actual count */
/*                    will be smaller than, but near, TRGFIN. */

/* $ Detailed_Output */

/*     CORSCL         is the coarse voxel scale compatible with the */
/*                    coarse grid determined by this routine. The grid */
/*                    is large enough to enclose the input vertex set. */
/*                    Its sides are aligned with the coordinate axes. */
/*                    CORSCL is the ratio of the edge length of coarse */
/*                    voxels to that of fine voxels. CORSCL is an */
/*                    integer. */

/*                    The edge length of coarse voxels is */

/*                       CORSCL * FINSCL * AVPLEX */


/*     FINSCL         is the fine voxel scale compatible with the */
/*                    coarse grid determined by this routine. */

/*                    The edge length of fine voxels is */

/*                       FINSCL * AVPLEX. */


/* $ Parameters */

/*     See mkdsk02.inc */

/* $ Exceptions */

/*     This routine is meant to be operated in RETURN SPICE error */
/*     handling mode. The caller is expected to delete the DSK file if */
/*     an error occurs during file creation. */


/*     1)  If the target coarse voxel count is less than 1 or greater */
/*         than MAXCGR, the error SPICE(VALUEOUTOFRANGE) is signaled. */

/*     2)  If the target fine voxel count is less than TRGCOR or greater */
/*         than MAXVOX, the error SPICE(VALUEOUTOFRANGE) is signaled. */

/*     3)  If the average plate extent is not strictly positive, the */
/*         error SPICE(VALUEOUTOFRANGE) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     The "extent" of a set of vertices is the set of extrema */
/*     of the coordinates of those vertices. */

/*     The extent of a plate along a coordinate axis is the length of */
/*     the orthogonal projection of the plate onto that axis. */

/*     The average "extent" of a set of plates is the average of the */
/*     extents of the plates along the coordinate axes. */

/*     The coarse voxel grid computed by this routine is large enough to */
/*     enclose the input vertex set. The grid has sides aligned with the */
/*     coordinate axes. */

/*     The voxel scales computed by this routine are such that, when the */
/*     coarse voxel grid is extended by 3 coarse voxels in each */
/*     dimension, the total number of coarse voxels does not exceed */
/*     MAXCGR, and the total number of fine voxels does not exceed */
/*     MAXVOX. */

/* $ Examples */

/*     See usage in MKDSK. */

/* $ Restrictions */

/*     This routine should be called only from within MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 1.0.0, 24-FEB-2017 (NJB) */

/* -& */
/* $ Index_Entries */

/*     determine type 2 dsk voxel scales */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("ZZVOXSCL", (ftnlen)8);

/*     Our target counts must be positive. */

    if (*trgcor < 1 || *trgcor > 100000) {
	setmsg_("Target coarse voxel count = #; count must be in the range 1"
		":#.", (ftnlen)62);
	errint_("#", trgcor, (ftnlen)1);
	errint_("#", &c_b6, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("ZZVOXSCL", (ftnlen)8);
	return 0;
    }
    if (*trgfin < *trgcor) {
	setmsg_("Target fine voxel count = #; target coarse count is #; targ"
		"et fine count must be greater than or equal to target coarse"
		" count.", (ftnlen)126);
	errint_("#", trgfin, (ftnlen)1);
	errint_("#", trgcor, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("ZZVOXSCL", (ftnlen)8);
	return 0;
    }

/*     Check the average plate extent. */

    if (*avplex <= 0.) {
	setmsg_("Average plate extent = #; extent must be strictly positive.",
		 (ftnlen)59);
	errdp_("#", avplex, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("ZZVOXSCL", (ftnlen)8);
	return 0;
    }

/*     Derive the coarse voxel size first. */


/*     Compute the lengths of the edges of the bounding box of the */
/*     vertex set. Compute the volume of the box. */

    lx = extent[1] - extent[0];
    ly = extent[3] - extent[2];
    lz = extent[5] - extent[4];

/*     In case we have zero extent in any dimension, assign an arbitrary */
/*     positive extent. */

    if (lx == 0.) {
	lx = 1.;
    }
    if (ly == 0.) {
	ly = 1.;
    }
    if (lz == 0.) {
	lz = 1.;
    }
    volume = lx * ly * lz;

/*     Compute the nominal coarse voxel size, given the target */
/*     coarse voxel count and the box volume. We know */

/*                                3 */
/*        volume = trgcor * corsz0 */


    d__1 = volume / *trgcor;
    corsz0 = pow_dd(&d__1, &c_b18);

/*     Compute the lengths of the bounding box edges as multiples of */
/*     CORSZ0. */

    sx = lx / corsz0;
    sy = ly / corsz0;
    sz = lz / corsz0;

/*     Since */

/*        SX * SY * SZ = VOLUME / CORSZ0**3 = TRGCOR */

/*     we know that */

/*        INT(SX) * INT(SY) * INT(SZ) <= TRGCOR */

/*     The integers above could serve as valid coarse voxel counts */
/*     in each dimension. However, we want to allow for as much */
/*     as 3 coarse voxels of padding in each dimension, so we'll */
/*     reduce the counts by 3 and compute the coarse voxel edge */
/*     length based on these reduced counts. */

/*     The portion of MKDSK code where the coarse voxel grid extension */
/*     takes place is inside the routine ZZMKSPIN. The code fragment */
/*     that carries out the extension is shown below: */


/*        C */
/*        C     Extend the coarse voxel grid by at least 1/2 */
/*        C     coarse voxel length along each degree of freedom. */
/*        C */
/*              XVMIN = DNINT ( XVMIN - 1.D0 ) */
/*              YVMIN = DNINT ( YVMIN - 1.D0 ) */
/*              ZVMIN = DNINT ( ZVMIN - 1.D0 ) */
/*              XVMAX = DNINT ( XVMAX + 1.D0 ) */
/*              YVMAX = DNINT ( YVMAX + 1.D0 ) */
/*              ZVMAX = DNINT ( ZVMAX + 1.D0 ) */


/* Computing MAX */
    i__1 = 1, i__2 = (integer) sx - 3;
    nx = max(i__1,i__2);
/* Computing MAX */
    i__1 = 1, i__2 = (integer) sy - 3;
    ny = max(i__1,i__2);
/* Computing MAX */
    i__1 = 1, i__2 = (integer) sz - 3;
    nz = max(i__1,i__2);

/*     Now we know the coarse voxel count, excluding padding voxels. */

    ncorse = nx * ny * nz;

/*     Compute the final coarse voxel size based on the bounding */
/*     box edge lengths and the coarse voxel counts along each */
/*     edge. */

/* Computing MAX */
    d__1 = lx / nx, d__2 = ly / ny, d__1 = max(d__1,d__2), d__2 = lz / nz;
    corsz = max(d__1,d__2);

/*     By construction, the coarse voxel grid created by ZZMKSPIN can */
/*     contain 3 extra coarse voxels in each dimension without exceeding */
/*     the target coarse voxel count. */

/*     Now compute the fine voxel count. We start with the target fine */
/*     voxel count and the computed coarse voxel count. From these we */
/*     can compute an estimated fractional coarse voxel scale. */

    d__1 = (doublereal) (*trgfin) / ncorse;
    dpcs = pow_dd(&d__1, &c_b18);

/*     The coarse scale is an integer. Round down the fractional scale; */
/*     this effectively enlarges the fine voxels. */

/*     Never set the coarse scale to a value below 1. */

/* Computing MAX */
    i__1 = (integer) dpcs;
    *corscl = max(i__1,1);

/*     Compute the maximum number of fine voxels that can occupy */
/*     the expanded fine grid. */

/* Computing 3rd power */
    i__1 = *corscl;
    maxnf = (nx + 3) * (ny + 3) * (nz + 3) * (i__1 * (i__1 * i__1));
    if (maxnf > 100000000) {

/*        The fine voxel count is too large. */

/*        We don't expect this to happen frequently, but if it does, */
/*        we'll reduce the coarse voxel scale, which effectively */
/*        increases the size of the fine voxels and reduces their */
/*        count. */

	q = 100000000 / (doublereal) ((nx + 3) * (ny + 3) * (nz + 3));
/* Computing MIN */
	i__1 = *corscl, i__2 = (integer) pow_dd(&q, &c_b18);
	*corscl = min(i__1,i__2);
    }

/*     The fine voxel edge length is the coarse voxel edge length, */
/*     scaled down by a factor of CORSCL. */

    finsz = corsz / *corscl;

/*     MKDSK requires the fine voxel scale, not the edge length. */

/*     The average plate extent times the fine voxel scale is the */
/*     fine voxel edge length, so we can derive the scale from */
/*     the edge length and extent. */

    *finscl = finsz / *avplex;
    chkout_("ZZVOXSCL", (ftnlen)8);
    return 0;
} /* zzvoxscl_ */

