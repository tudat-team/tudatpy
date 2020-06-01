/* zztrgnvx.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b3 = 1.2e7;
static integer c__1 = 1;
static integer c_b5 = 100000;
static integer c__8 = 8;
static integer c_b7 = 100000000;

/* $Procedure    ZZTRGNVX ( MKDSK, compute target voxel counts ) */
/* Subroutine */ int zztrgnvx_(integer *np, integer *trgcor, integer *trgfin)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);
    double sqrt(doublereal);

    /* Local variables */
    doublereal frac;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    extern integer brckti_(integer *, integer *, integer *);
    extern /* Subroutine */ int chkout_(char *, ftnlen);
    extern logical return_(void);

/* $ Abstract */


/*     Compute target coarse and fine voxel counts for a DSK type 2 */
/*     segment. */

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
/*     NP         O   Plate count. */
/*     TRGCOR     O   Target coarse voxel count. */
/*     TRGFIN     O   Target fine voxel count. */

/* $ Detailed_Input */

/*     NP             is the plate count of a type 2 DSK segment. */

/* $ Detailed_Output */

/*     TRGCOR         is the target coarse voxel count. This is an */
/*                    approximation to the final count; the actual count */
/*                    set by MKDSK will be smaller than, but near, */
/*                    TRGCOR. */

/*     TRGFIN         is the target fine voxel count. This is an */
/*                    approximation to the final count; the actual count */
/*                    set by MKDSK will be smaller than, but near, */
/*                    TRGFIN. */

/* $ Parameters */

/*     See mkdsk02.inc */

/* $ Exceptions */

/*     This routine is meant to be operated in RETURN SPICE error */
/*     handling mode. The caller is expected to delete the DSK file if */
/*     an error occurs during file creation. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     The target coarse and fine voxel counts produced by this routine */
/*     are set by heuristics. The counts are meant to yield reasonable */
/*     DSK type 2 ray-surface intercept computation speed. */

/*     Actual counts selected by MKDSK will be less than the respective */
/*     target counts. */

/*     See the source code for details. */

/* $ Examples */

/*     See usage in MKDSK. */

/* $ Restrictions */

/*     This routine should be called only from within MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 1.0.0, 19-FEB-2017 (NJB) */

/* -& */
/* $ Index_Entries */

/*     compute type 2 dsk target voxel counts */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("ZZTRGNVX", (ftnlen)8);

/*     All of the formulas below for target voxel counts are heuristics. */
/*     They may be updated if superior formulas are developed. */

    if ((doublereal) (*np) > 1e6) {

/*        For plate counts of NPLIM or more, we go for a full-size */
/*        coarse voxel grid. */

	*trgcor = 100000;

/*        The fine voxel count increases linearly as a function */
/*        of NP - NPLIM. The count at NP = MAXPLT is 2*MAXPLT. */

	d__1 = (*np - 1e6) * .64516129032258063 * 2;
	*trgfin = i_dnnt(&c_b3) + i_dnnt(&d__1);
    } else if (*np > 1000) {

/*        Scale down the coarse grid proportionally to */
/*        the square root of the ratio of NP to NPLIM. */

	frac = (doublereal) (*np) / 1e6;
	d__1 = sqrt(frac) * 100000;
	*trgcor = i_dnnt(&d__1);

/*        Scale down the fine voxel count as well. */

	d__1 = frac * 1.2e7;
	*trgfin = i_dnnt(&d__1);
    } else if (*np > 100) {
	*trgcor = 100;
	*trgfin = 1000;
    } else {
	*trgcor = 100;
	*trgfin = 100;
    }

/*     Bracket the estimates to ensure they're in range. */

    *trgcor = brckti_(trgcor, &c__1, &c_b5);
    *trgfin = brckti_(trgfin, &c__8, &c_b7);
    chkout_("ZZTRGNVX", (ftnlen)8);
    return 0;
} /* zztrgnvx_ */

