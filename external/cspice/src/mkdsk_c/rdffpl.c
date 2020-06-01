/* rdffpl.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c_b8 = 16000002;
static integer c_b22 = 32000000;
static integer c__1024 = 1024;
static integer c__2 = 2;
static integer c__3 = 3;

/* $Procedure   RDFFPL ( read triangular plate model flatfile ) */
/* Subroutine */ int rdffpl_(char *infile__, integer *plttyp, integer *nv, 
	doublereal *vrtces, integer *np, integer *plates, ftnlen infile_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer nrec;
    extern logical even_(integer *);
    static integer ntok, i__, j, n[6303750]	/* was [1025][1025][6] */, q;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    static char iline[255];
    extern /* Subroutine */ int rdnbl_(char *, char *, logical *, ftnlen, 
	    ftnlen), errch_(char *, char *, ftnlen, ftnlen), prsdp_(char *, 
	    doublereal *, ftnlen);
    static char error[160];
    extern logical eqstr_(char *, char *, ftnlen, ftnlen);
    static integer n0;
    extern doublereal vnorm_(doublereal *);
    extern /* Subroutine */ int vcrss_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal w1[3], w2[3];
    static integer nc, nd;
    extern logical failed_(void);
    static integer ni;
    extern /* Subroutine */ int rdffdi_(char *, integer *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, char *, logical *, 
	    ftnlen, ftnlen, ftnlen), nparsd_(char *, doublereal *, char *, 
	    integer *, ftnlen, ftnlen);
    static char format[10], tokens[80*3];
    extern /* Subroutine */ int chkout_(char *, ftnlen);
    extern logical return_(void);
    static integer ix1, ix2;
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), sigerr_(char *, ftnlen), extmsi_(char *, char 
	    *, integer *, ftnlen, ftnlen), nparsi_(char *, integer *, char *, 
	    integer *, ftnlen, ftnlen);
    static char arc[40*6];
    static doublereal ard[3];
    extern /* Subroutine */ int lparsm_(char *, char *, integer *, integer *, 
	    char *, ftnlen, ftnlen, ftnlen), prsint_(char *, integer *, 
	    ftnlen);
    static logical eof;
    extern /* Subroutine */ int tostdo_(char *, ftnlen);
    static integer ari[4];
    static doublereal vec[3];
    static integer ptr, face;

/* $ Abstract */

/*     Read a triangular plate model's flatfile's data file. */

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

/*     MKPLAT User's Guide */
/*     SURFACE */

/* $ Keywords */

/*     PLATE */
/*     FLATFILE */

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

/* $ Abstract */

/*     Include Section:  MKDSK Global Parameters */

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

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -    Version 4.0.0, 28-FEB-2017 (NJB) */

/*        Added declaration of version string VER. Previously */
/*        this declaration was located in prcinf.for. */

/*        Declarations of parameters */

/*           MAXCEL */
/*           MAXVXP */
/*           MAXNVLS */

/*        were moved to dsk02.inc. */

/*        Declarations of parameters */

/*           ZZMAXV */
/*           ZZMAXP */

/*        were deleted. */


/* -    Version 3.0.0, 20-OCT-2015 (NJB) */

/*        Parameter MAXQ was increased from 512 to 1024. */

/* -    Version 2.0.0, 26-MAR-2015 (NJB) */

/*        Declarations were added for the parameters */

/*           MAXCEL */
/*           MAXVXP */
/*           MXNVLS */


/* -    Version 1.0.0, 04-MAY-2010 (NJB) */

/* -& */

/*     MKDSK version: */


/*     Default time format: */


/*     Command line length: */


/*     SPICELIB cell lower bound: */


/*     Maximum file name length: */


/*     Output file line length: */


/*     Length of string for keyword value processing: */


/*     The maximum 'Q' value for Gaskell shape models. This value */
/*     is always a multiple of 4. */


/*     Input format parameters */


/*     The height grid "plate type" is the fifth format. */


/*     End Include Section:  MKDSK Global Parameters */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     INFILE     I   Name of plate data file. */
/*     PLTTYP     I   Integer code identifying type of data file. */
/*     NV         O   Number of vertices in model. */
/*     VRTCES     O   3 x NV array of vertices. */
/*     NP         O   Number of plates in model. */
/*     PLATES     O   3 x NP array of plates. */


/* $ Detailed_Input */

/*     INFILE     is the full pathname of the plate model flatfile. */

/*     PLTTYP     the integer code identifying the type of data */
/*                to read from INFILE. Two values are allowed: */

/*                    1    A standard plate-vertex data file */
/*                         containing vertex coordinates and */
/*                         the plate-vertex mappings. */

/*                    2    A Gaskel shape data file containing */
/*                         plate ordered vertex data. */

/*                    3    A vertex-facet table data file */
/*                         containing vertex coordinates and */
/*                         the facet listing. */

/*                    4    Rosetta/Osiris ".ver" format file */
/*                         containing vertex coordinates */
/*                         and a plate-vertex mapping. */

/* $ Detailed_Output */

/*     NV         is the number of vertices in the plate model. */


/*     VRTCES     is an array of the NV vertices given in a body-fixed */
/*                frame of reference. Elements */

/*                   VRTCES(J,I), J = 1 ... 3 */

/*                are, respectively, the X, Y, and Z coordinates of the */
/*                Ith vertex. */


/*     PLATES     is an array of the NP plates. Elements */

/*                   PLATES(J,I), J = 1 ... 3 */

/*                are, respectively, the indices of the vertices of the */
/*                Ith plate. The vertex indices range from 1 to NV. */

/*                The order of the vertices give the direction of the */
/*                plate's outward normal vector: the vector is parallel */
/*                to the cross product of the plate edge connecting */
/*                vertex 1 to vertex 2 and the plate edge connecting */
/*                vertex 2 to vertex 3, in that order. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the number of vertices (NV) exceeds the value of MAXVRT, */
/*        SPICE(TOOMANYVERTICES) signals. */

/*     2) If the number of plates (NP) exceeds the value of MAXPLT, */
/*        SPICE(TOOMANYPLATES) signals. */

/*     3) If the NPARSD routine fails while parsing a double string, */
/*        SPICE(BADDOUBLEPRECISION) signals. */

/*     4) If the NPARSI routine fails while parsing an integer string, */
/*        SPICE(BADINTEGER) signals. */

/*     5) If a data line lacks the expected format (data type 3), */
/*        SPICE(BADDATALINE) signals. */

/*     6) If the PLTTYP variable has an uncoded value, */
/*        SPICE(BADDATATYPE) signals. */

/*     7) If the Gaskell ICQ parameter "Q" is too large, the error */
/*        SPICE(QPARAMOUTOFRANGE) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     Input File Format */
/*     ----------------- */

/*     INFILE is assumed to be structured as follows (values */
/*            are space delimited within each input record). */

/*     -Plate/vertex data file (PLATE_TYPE 1): */

/*        Line  1: NV */
/*        NV = the number of vertices (INTEGER) */

/*        Lines 2 to (NV+1): VID X Y Z */
/*        VID         = the vertex ID (INTEGER) */
/*        X, Y, and Z = the vertex coordinates (DOUBLE PRECISION) */

/*        Line  (NV+2): NP */
/*        NP = number of triangular plates (INTEGER) */

/*        Lines (NV+3) to (NV+3+NP+1): PID V1 V2 V3 */
/*        PID            = the plate ID (INTEGER) */
/*        V1, V2, and V3 = the plate vertex ID's (INTEGER) */

/*     -Shape data file (PLATE_TYPE 2): */

/*        Line 1: Q */
/*        Q = number of vertex coordinates per side of the */
/*            cube face - not counting the origin. A */
/*            face has (Q+1)^2 vertex coordinates. */

/*        Lines 2 to (Q+1): X Y Z */
/*        X, Y, and Z = the vertex coordinates (DOUBLE PRECISION) */

/*     -Vertex/Facet table data file (PLATE_TYPE 3): */

/*        NV = the number of vertices (INTEGER) */

/*        Lines 1 to NV: 'V' X Y Z */
/*        'V'         = character flag indicating vertex data */
/*        X, Y, and Z = the vertex coordinates (DOUBLE PRECISION) */

/*        NP = number of triangular plates (INTEGER) */

/*        Lines (NV+1) to (NV+1+NP): 'F' V1 V2 V3 */
/*        'F'            = character flag indicating facet */
/*                         (plate-vertex) data */
/*        V1, V2, and V3 = the plate vertex ID's (INTEGER) */

/* $ Examples */

/*  C     The following include file defines MAXVRT and */
/*  C     MAXPLT, the maximum number of vertices and plates */
/*  C     used in the plate model software. */

/*        INCLUDE               'mkdsk02.inc' */

/*        INTEGER               NV */
/*        INTEGER               NP */
/*        INTEGER               PLATES ( 3, MAXPLT ) */
/*        INTEGER               PLTTYP */

/*        DOUBLE PRECISION      VRTCES ( 3, MAXVRT ) */

/*        CHARACTER*(80)        INFILE */

/*        CALL PROMPT ( 'Enter name of plate data file : ', INFILE ) */
/*        CALL PROMPT ( 'Enter type of plate data      : ', PLTTYP ) */

/*        CALL RDFFPL ( INFILE,  PLTTYP, NV, VERTCES, NP, PLATES ) */

/* $ Restrictions */

/*     It is the user's responsibility to properly dimension */
/*     arrays in the calling routine large enough to accept data */
/*     from the given file. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     J.A. Bytof      (JPL) */

/* $ Version */

/* -    MKDSK Version 5.1.0, 21-MAR-2017 (NJB) */

/*        Now references non-portable include file */

/*           mkdsk02.inc */

/*        instead of */

/*           dsk02.inc */

/*        This enables use of smaller maximum vertex and plate */
/*        counts on platforms where the default values consume */
/*        excessive memory. */

/*        Added checks for excessive vertex and plate counts. */

/*        Added FAILED checks following RDFFDI, RDNBL, and */
/*        PRSINT calls. */

/* -    MKDSK Version 5.0.0, 25-APR-2016 (NJB) */

/*        Bug fix: added check-out and return to input format 2 */
/*        branch, in the case where a Q value out of range is found. */



/*        07-APR-2015 (NJB) */

/*           Increased string length used to parse numeric tokens. */
/*           Updated message indicating completion of input file read. */
/*           Cleaned up debugging code. */

/*        05-AUG-2014 (NJB) */

/*           Argument list change: arguments PID and VID have been */
/*           removed; arguments X, Y, Z and V1, V2, V3 have been */
/*           replaced by arrays VRTCES and PLATES respectively. */

/*           Local array VEC is no longer used. */

/* -    MKDSK Version 4.0.0, 08-JUN-2010 (NJB) */

/*        Added capability of reading Rosetta Osiris ".ver" */
/*        format file. */

/* -    MKDSK Version 3.1.0, 04-MAY-2010 (NJB) */

/*        Changed INCLUDE file from platmax.inc to dsk02.inc. */
/*        Added INCLUDE statement referencing mkdsk.inc to */
/*        declare parameter MAXQ. */

/* -    MKDSK Version 3.0.1, 08-OCT-2009 (NJB) */

/*        Re-ordered header sections. */

/* -    MKDSK Version 3.0.0, 25-OCT-2004 (EDW) */

/*        Added capability to process Bob Gaskell shape files. */
/*        Added MAXQ parameter to 'pltmax.inc'. */

/* -    MKDSK Version 2.0.0, 22-OCT-1998 (JAB) */

/* -    MKDSK Version 1.0.0, 09-APR-1997 (JAB) */

/* -& */
/* $ Index_Entries */

/*     read triangular plate model flatfile */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    }
    chkin_("RDFFPL", (ftnlen)6);

/*     Initialize record counter. */

    nrec = 1;

/*     Read the data file corresponding to the PLTTYP ID. */

    if (*plttyp == 1) {

/*        Plate data type 1. */

/*        First line, a single integer, number of vertices. */

	s_copy(format, "I", (ftnlen)10, (ftnlen)1);
	rdffdi_(infile__, &nrec, format, &nd, &ni, &nc, ard, ari, arc, &eof, 
		infile_len, (ftnlen)10, (ftnlen)40);
	if (failed_()) {
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
	*nv = ari[0];

/*        Check the routine can handle the vertices. */

	if (*nv > 16000002) {
	    setmsg_("Number of vertices # exceeds limit #.", (ftnlen)37);
	    errint_("#", nv, (ftnlen)1);
	    errint_("#", &c_b8, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYVERTICES)", (ftnlen)22);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}

/*        Read data: integer double double double. */

	s_copy(format, "I D D D", (ftnlen)10, (ftnlen)7);
	i__1 = *nv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++nrec;
	    rdffdi_(infile__, &nrec, format, &nd, &ni, &nc, ard, ari, arc, &
		    eof, infile_len, (ftnlen)10, (ftnlen)40);
	    if (failed_()) {
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	    if (eof) {
		extmsi_("End of file after line #.", "#", &nrec, (ftnlen)25, (
			ftnlen)1);
	    }
	    vrtces[i__ * 3 - 3] = ard[0];
	    vrtces[i__ * 3 - 2] = ard[1];
	    vrtces[i__ * 3 - 1] = ard[2];
	}

/*        Plate data type 1 reads the plat-vertex information from */
/*        the plate file. */

/*        Read number of plates. */

	++nrec;
	s_copy(format, "I", (ftnlen)10, (ftnlen)1);
	rdffdi_(infile__, &nrec, format, &nd, &ni, &nc, ard, ari, arc, &eof, 
		infile_len, (ftnlen)10, (ftnlen)40);
	if (failed_()) {
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
	if (eof) {
	    extmsi_("End of file after line #.", "#", &nrec, (ftnlen)25, (
		    ftnlen)1);
	}

/*        Check we can process 'NP' plates. */

	*np = ari[0];
	if (*np > 32000000) {
	    setmsg_("Number of plates # exceeds limit #.", (ftnlen)35);
	    errint_("#", np, (ftnlen)1);
	    errint_("#", &c_b22, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYPLATES)", (ftnlen)20);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}

/*        Read each plate's ID and corresponding vertex set defining */
/*        the plate. */

	s_copy(format, "I I I I", (ftnlen)10, (ftnlen)7);
	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++nrec;
	    rdffdi_(infile__, &nrec, format, &nd, &ni, &nc, ard, ari, arc, &
		    eof, infile_len, (ftnlen)10, (ftnlen)40);
	    if (failed_()) {
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	    if (eof) {
		extmsi_("End of file after line #.", "#", &nrec, (ftnlen)25, (
			ftnlen)1);
	    }
	    plates[i__ * 3 - 3] = ari[1];
	    plates[i__ * 3 - 2] = ari[2];
	    plates[i__ * 3 - 1] = ari[3];
	}
    } else if (*plttyp == 2) {

/*        Plate data type 2. */

/*        Read a Gaskell shape file. Shape files contain only vertex */
/*        data, the plate-vertex mappings implicitly known from the */
/*        ordering of the data. The Gaskell model uses quadralaterals */
/*        plates, a plate defined by the vector set: */

/*        [ VEC( *, I, J  , FACE), VEC( *, I+1, J  , FACE), */
/*          VEC( *, I, J+1, FACE), VEC( *, I+1, J+1, FACE) ] */

/*        The file lists the vertex data in terms cube faces. Six */
/*        faces, each with Q+1 x Q+1 vertices. */

	s_copy(format, "I", (ftnlen)10, (ftnlen)1);
	rdffdi_(infile__, &nrec, format, &nd, &ni, &nc, ard, ari, arc, &eof, 
		infile_len, (ftnlen)10, (ftnlen)40);
	if (failed_()) {
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
	q = ari[0];
/* Computing 2nd power */
	i__1 = q + 1;
	*nv = i__1 * i__1 * 6;

/*        Does the read Q exceed the maximum value. */

	if (q > 1024) {
	    setmsg_("Shape parameter Q = #, exceeds maximum value #. This er"
		    "ror may indicate that the input file format is something"
		    " other than the Gaskell ICQ format.", (ftnlen)146);
	    errint_("#", &q, (ftnlen)1);
	    errint_("#", &c__1024, (ftnlen)1);
	    sigerr_("SPICE(QPARAMOUTOFRANGE)", (ftnlen)23);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}

/*        Check the routine can handle the vertices. */

	if (*nv > 16000002) {
	    setmsg_("Number of vertices # exceeds limit #.", (ftnlen)37);
	    errint_("#", nv, (ftnlen)1);
	    errint_("#", &c_b8, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYVERTICES)", (ftnlen)22);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}

/*        Read data: double double double */

	s_copy(format, "D D D", (ftnlen)10, (ftnlen)5);
	n0 = 0;
	for (face = 1; face <= 6; ++face) {
	    i__1 = q;
	    for (j = 0; j <= i__1; ++j) {
		i__2 = q;
		for (i__ = 0; i__ <= i__2; ++i__) {
		    ++nrec;

/*                 Read the vertex data. The 3-Vector */
/*                    _                 _ */
/*                   |  VEC(1,I,J,FACE)  | */
/*                   |  VEC(2,I,J,FACE)  | */
/*                   |_ VEC(3,I,J,FACE) _| */

/*                 contains the vertex coordinates. */

		    rdffdi_(infile__, &nrec, format, &nd, &ni, &nc, vec, ari, 
			    arc, &eof, infile_len, (ftnlen)10, (ftnlen)40);
		    if (failed_()) {
			chkout_("RDFFPL", (ftnlen)6);
			return 0;
		    }
		    if (eof) {
			extmsi_("End of file after line #.", "#", &nrec, (
				ftnlen)25, (ftnlen)1);
		    }

/*                Store the vertex ID to the 'N' array using the same */
/*                indexing as 'VEC'. */

		    ++n0;
		    n[(i__3 = i__ + (j + face * 1025) * 1025 - 1050625) < 
			    6303750 && 0 <= i__3 ? i__3 : s_rnge("n", i__3, 
			    "rdffpl_", (ftnlen)608)] = n0;

/*                Copy the vertex coordinates to the corresponding */
/*                plate model arrays. This operation is somewhat */
/*                repetitive, but clarity is most important. */

		    vrtces[n0 * 3 - 3] = vec[0];
		    vrtces[n0 * 3 - 2] = vec[1];
		    vrtces[n0 * 3 - 1] = vec[2];
		}
	    }
	}

/*        Plate data type 2 calculates the plate-vertex mappings */
/*        from the ordering of the vertex data. */

/*        Number of plates based on the 'Q' value. */

/* Computing 2nd power */
	i__1 = q;
	*np = i__1 * i__1 * 12;

/*        Check we can process 'NP' plates. */

	if (*np > 32000000) {
	    setmsg_("Number of plates # exceeds limit #.", (ftnlen)35);
	    errint_("#", np, (ftnlen)1);
	    errint_("#", &c_b22, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYPLATES)", (ftnlen)20);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}

/*        Expand a Gaskell shape model vertex ordering to */
/*        form a plate-vertex map. */


/*        Cube edges. Associate later IDs with */
/*        the first occurrence. */

	i__1 = q - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    n[(i__2 = i__ + (q + 6150) * 1025 - 1050625) < 6303750 && 0 <= 
		    i__2 ? i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)657)] =
		     n[(i__3 = q - i__ + (q + 4100) * 1025 - 1050625) < 
		    6303750 && 0 <= i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_",
		     (ftnlen)657)];
	    n[(i__2 = i__ + 5253125) < 6303750 && 0 <= i__2 ? i__2 : s_rnge(
		    "n", i__2, "rdffpl_", (ftnlen)658)] = n[(i__3 = i__ + (q 
		    + 2050) * 1025 - 1050625) < 6303750 && 0 <= i__3 ? i__3 : 
		    s_rnge("n", i__3, "rdffpl_", (ftnlen)658)];
	    n[(i__2 = i__ + 4202500) < 6303750 && 0 <= i__2 ? i__2 : s_rnge(
		    "n", i__2, "rdffpl_", (ftnlen)659)] = n[(i__3 = q + (q - 
		    i__ + 1025) * 1025 - 1050625) < 6303750 && 0 <= i__3 ? 
		    i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)659)];
	    n[(i__2 = i__ + 3151875) < 6303750 && 0 <= i__2 ? i__2 : s_rnge(
		    "n", i__2, "rdffpl_", (ftnlen)660)] = n[(i__3 = q - i__) <
		     6303750 && 0 <= i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_"
		    , (ftnlen)660)];
	    n[(i__2 = i__ + 2101250) < 6303750 && 0 <= i__2 ? i__2 : s_rnge(
		    "n", i__2, "rdffpl_", (ftnlen)661)] = n[(i__3 = (i__ + 
		    1025) * 1025 - 1050625) < 6303750 && 0 <= i__3 ? i__3 : 
		    s_rnge("n", i__3, "rdffpl_", (ftnlen)661)];
	    n[(i__2 = i__ + 1050625) < 6303750 && 0 <= i__2 ? i__2 : s_rnge(
		    "n", i__2, "rdffpl_", (ftnlen)662)] = n[(i__3 = i__ + (q 
		    + 1025) * 1025 - 1050625) < 6303750 && 0 <= i__3 ? i__3 : 
		    s_rnge("n", i__3, "rdffpl_", (ftnlen)662)];
	}
	i__1 = q - 1;
	for (j = 1; j <= i__1; ++j) {
	    n[(i__2 = q + (j + 6150) * 1025 - 1050625) < 6303750 && 0 <= i__2 
		    ? i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)668)] = n[(
		    i__3 = j + (q + 5125) * 1025 - 1050625) < 6303750 && 0 <= 
		    i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)668)];
	    n[(i__2 = q + (j + 5125) * 1025 - 1050625) < 6303750 && 0 <= i__2 
		    ? i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)669)] = n[(
		    i__3 = (j + 4100) * 1025 - 1050625) < 6303750 && 0 <= 
		    i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)669)];
	    n[(i__2 = q + (j + 4100) * 1025 - 1050625) < 6303750 && 0 <= i__2 
		    ? i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)670)] = n[(
		    i__3 = (j + 3075) * 1025 - 1050625) < 6303750 && 0 <= 
		    i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)670)];
	    n[(i__2 = q + (j + 3075) * 1025 - 1050625) < 6303750 && 0 <= i__2 
		    ? i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)671)] = n[(
		    i__3 = (j + 2050) * 1025 - 1050625) < 6303750 && 0 <= 
		    i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)671)];
	    n[(i__2 = (j + 6150) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? 
		    i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)672)] = n[(
		    i__3 = q - j + (q + 3075) * 1025 - 1050625) < 6303750 && 
		    0 <= i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)
		    672)];
	    n[(i__2 = (j + 5125) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? 
		    i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)673)] = n[(
		    i__3 = q + (j + 2050) * 1025 - 1050625) < 6303750 && 0 <= 
		    i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)673)];
	}

/*        Cube corners. Associate later IDs with */
/*        the first occurrence. */

	n[2101250] = n[0];
	n[(i__1 = q + 3151875) < 6303750 && 0 <= i__1 ? i__1 : s_rnge("n", 
		i__1, "rdffpl_", (ftnlen)682)] = n[0];
	n[1050625] = n[(i__1 = (q + 1025) * 1025 - 1050625) < 6303750 && 0 <= 
		i__1 ? i__1 : s_rnge("n", i__1, "rdffpl_", (ftnlen)683)];
	n[(i__1 = q + 2101250) < 6303750 && 0 <= i__1 ? i__1 : s_rnge("n", 
		i__1, "rdffpl_", (ftnlen)684)] = n[(i__2 = (q + 1025) * 1025 
		- 1050625) < 6303750 && 0 <= i__2 ? i__2 : s_rnge("n", i__2, 
		"rdffpl_", (ftnlen)684)];
	n[3151875] = n[(i__1 = q) < 6303750 && 0 <= i__1 ? i__1 : s_rnge(
		"n", i__1, "rdffpl_", (ftnlen)685)];
	n[(i__1 = q + 4202500) < 6303750 && 0 <= i__1 ? i__1 : s_rnge("n", 
		i__1, "rdffpl_", (ftnlen)686)] = n[(i__2 = q) < 6303750 && 0 
		<= i__2 ? i__2 : s_rnge("n", i__2, "rdffpl_", (ftnlen)686)];
	n[4202500] = n[(i__1 = q + (q + 1025) * 1025 - 1050625) < 6303750 && 
		0 <= i__1 ? i__1 : s_rnge("n", i__1, "rdffpl_", (ftnlen)687)];
	n[(i__1 = q + 1050625) < 6303750 && 0 <= i__1 ? i__1 : s_rnge("n", 
		i__1, "rdffpl_", (ftnlen)688)] = n[(i__2 = q + (q + 1025) * 
		1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : s_rnge("n", 
		i__2, "rdffpl_", (ftnlen)688)];
	n[5253125] = n[(i__1 = (q + 2050) * 1025 - 1050625) < 6303750 && 0 <= 
		i__1 ? i__1 : s_rnge("n", i__1, "rdffpl_", (ftnlen)689)];
	n[(i__1 = q + (q + 3075) * 1025 - 1050625) < 6303750 && 0 <= i__1 ? 
		i__1 : s_rnge("n", i__1, "rdffpl_", (ftnlen)690)] = n[(i__2 = 
		(q + 2050) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : 
		s_rnge("n", i__2, "rdffpl_", (ftnlen)690)];
	n[(i__1 = (q + 5125) * 1025 - 1050625) < 6303750 && 0 <= i__1 ? i__1 :
		 s_rnge("n", i__1, "rdffpl_", (ftnlen)691)] = n[(i__2 = q + (
		q + 2050) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : 
		s_rnge("n", i__2, "rdffpl_", (ftnlen)691)];
	n[(i__1 = q + 5253125) < 6303750 && 0 <= i__1 ? i__1 : s_rnge("n", 
		i__1, "rdffpl_", (ftnlen)692)] = n[(i__2 = q + (q + 2050) * 
		1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : s_rnge("n", 
		i__2, "rdffpl_", (ftnlen)692)];
	n[(i__1 = q + (q + 4100) * 1025 - 1050625) < 6303750 && 0 <= i__1 ? 
		i__1 : s_rnge("n", i__1, "rdffpl_", (ftnlen)693)] = n[(i__2 = 
		(q + 3075) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : 
		s_rnge("n", i__2, "rdffpl_", (ftnlen)693)];
	n[(i__1 = (q + 6150) * 1025 - 1050625) < 6303750 && 0 <= i__1 ? i__1 :
		 s_rnge("n", i__1, "rdffpl_", (ftnlen)694)] = n[(i__2 = (q + 
		3075) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : 
		s_rnge("n", i__2, "rdffpl_", (ftnlen)694)];
	n[(i__1 = q + (q + 5125) * 1025 - 1050625) < 6303750 && 0 <= i__1 ? 
		i__1 : s_rnge("n", i__1, "rdffpl_", (ftnlen)695)] = n[(i__2 = 
		(q + 4100) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : 
		s_rnge("n", i__2, "rdffpl_", (ftnlen)695)];
	n[(i__1 = q + (q + 6150) * 1025 - 1050625) < 6303750 && 0 <= i__1 ? 
		i__1 : s_rnge("n", i__1, "rdffpl_", (ftnlen)696)] = n[(i__2 = 
		(q + 4100) * 1025 - 1050625) < 6303750 && 0 <= i__2 ? i__2 : 
		s_rnge("n", i__2, "rdffpl_", (ftnlen)696)];
	n0 = 0;
	for (face = 1; face <= 6; ++face) {
	    i__1 = q - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		i__2 = q - 1;
		for (j = 0; j <= i__2; ++j) {

/*                 As mentioned, the Gaskell shape model uses */
/*                 quadrilateral plates as opposed to triangular used by */
/*                 the NAIF plate system. We can reduce a quadrilateral */
/*                 into two triangles by adding a vector connecting */
/*                 opposite facing vertices creating two triangles */
/*                 (1 and 2): */

/*                               V1 */
/*                               /\ */
/*                              / 1\ */
/*                          V4 /____\ V2 */
/*                             \  2 / */
/*                              \  / */
/*                               \/ */
/*                               V3 */

/*                 or */
/*                                V1 */
/*                               /|\ */
/*                              / | \ */
/*                         V4  /  |2 \ V2 */
/*                             \ 1|  / */
/*                              \ | / */
/*                               \|/ */
/*                                V3 */

/*                 We chose the connecting vector to minimize the cross */
/*                 product magnitude of the connected vertex vectors. */

/*                    connection = MIN( V1 x V3, V2 X V4) */


/*                 Connection 1. */

/* Computing 2nd power */
		    i__3 = q + 1;
		    ix1 = (face - 1) * (i__3 * i__3) + (q + 1) * j + i__ + 1;
/* Computing 2nd power */
		    i__3 = q + 1;
		    ix2 = (face - 1) * (i__3 * i__3) + (q + 1) * (j + 1) + (
			    i__ + 1) + 1;
		    vcrss_(&vrtces[ix1 * 3 - 3], &vrtces[ix2 * 3 - 3], w1);

/*                 Connection 2. */

/* Computing 2nd power */
		    i__3 = q + 1;
		    ix1 = (face - 1) * (i__3 * i__3) + (q + 1) * j + (i__ + 1)
			     + 1;
/* Computing 2nd power */
		    i__3 = q + 1;
		    ix2 = (face - 1) * (i__3 * i__3) + (q + 1) * (j + 1) + 
			    i__ + 1;
		    vcrss_(&vrtces[ix1 * 3 - 3], &vrtces[ix2 * 3 - 3], w2);

/*                 Calculate the magnitudes of the cross products; */
/*                 branch based on the minimum value. */

/*                 Fill two entries of the plate-vertex mapping */
/*                 from the divided quadralateral. */

		    if (vnorm_(w1) <= vnorm_(w2)) {

/*                    Apply connection 1. */

			++n0;
			plates[n0 * 3 - 3] = n[(i__3 = i__ + (j + face * 1025)
				 * 1025 - 1050625) < 6303750 && 0 <= i__3 ? 
				i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)
				770)];
			plates[n0 * 3 - 2] = n[(i__3 = i__ + 1 + (j + 1 + 
				face * 1025) * 1025 - 1050625) < 6303750 && 0 
				<= i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", 
				(ftnlen)771)];
			plates[n0 * 3 - 1] = n[(i__3 = i__ + 1 + (j + face * 
				1025) * 1025 - 1050625) < 6303750 && 0 <= 
				i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (
				ftnlen)772)];
			++n0;
			plates[n0 * 3 - 3] = n[(i__3 = i__ + (j + face * 1025)
				 * 1025 - 1050625) < 6303750 && 0 <= i__3 ? 
				i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)
				775)];
			plates[n0 * 3 - 2] = n[(i__3 = i__ + (j + 1 + face * 
				1025) * 1025 - 1050625) < 6303750 && 0 <= 
				i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (
				ftnlen)776)];
			plates[n0 * 3 - 1] = n[(i__3 = i__ + 1 + (j + 1 + 
				face * 1025) * 1025 - 1050625) < 6303750 && 0 
				<= i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", 
				(ftnlen)777)];
		    } else {

/*                    Apply connection 2. */

			++n0;
			plates[n0 * 3 - 3] = n[(i__3 = i__ + (j + face * 1025)
				 * 1025 - 1050625) < 6303750 && 0 <= i__3 ? 
				i__3 : s_rnge("n", i__3, "rdffpl_", (ftnlen)
				784)];
			plates[n0 * 3 - 2] = n[(i__3 = i__ + (j + 1 + face * 
				1025) * 1025 - 1050625) < 6303750 && 0 <= 
				i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (
				ftnlen)785)];
			plates[n0 * 3 - 1] = n[(i__3 = i__ + 1 + (j + face * 
				1025) * 1025 - 1050625) < 6303750 && 0 <= 
				i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (
				ftnlen)786)];
			++n0;
			plates[n0 * 3 - 3] = n[(i__3 = i__ + 1 + (j + face * 
				1025) * 1025 - 1050625) < 6303750 && 0 <= 
				i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (
				ftnlen)789)];
			plates[n0 * 3 - 2] = n[(i__3 = i__ + (j + 1 + face * 
				1025) * 1025 - 1050625) < 6303750 && 0 <= 
				i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", (
				ftnlen)790)];
			plates[n0 * 3 - 1] = n[(i__3 = i__ + 1 + (j + 1 + 
				face * 1025) * 1025 - 1050625) < 6303750 && 0 
				<= i__3 ? i__3 : s_rnge("n", i__3, "rdffpl_", 
				(ftnlen)791)];
		    }
		}
	    }
	}
    } else if (*plttyp == 3) {

/*        Plate data type 3. */

/*        Read data. The data format depends on the */
/*        prefix marker either 'v' (vertex) or 'f' (facet). */
/*        Since we don't know the prefix a priori, */
/*        extract the four data elements, then parse. */

	*np = 0;
	*nv = 0;
	s_copy(format, "C C C C", (ftnlen)10, (ftnlen)7);
	eof = FALSE_;
	while(! eof) {
	    rdffdi_(infile__, &nrec, format, &nd, &ni, &nc, ard, ari, arc, &
		    eof, infile_len, (ftnlen)10, (ftnlen)40);
	    if (failed_()) {
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	    if (! eof) {

/*              Branch based on the type of data, 'v' or 'f'. */

		if (eqstr_(arc, "v", (ftnlen)40, (ftnlen)1)) {

/*                 Vertex data. Process a line of three doubles. */
/*                 The algorithm includes an implicit assumption that */
/*                 the data order defines the vertex IDs.  If not, */
/*                 the output plate files will be useless. */

		    ++(*nv);
		    if (*nv > 16000002) {
			setmsg_("Number of vertices # exceeds limit #.", (
				ftnlen)37);
			errint_("#", nv, (ftnlen)1);
			errint_("#", &c_b8, (ftnlen)1);
			sigerr_("SPICE(TOOMANYVERTICES)", (ftnlen)22);
			chkout_("RDFFPL", (ftnlen)6);
			return 0;
		    }
		    for (j = 2; j <= 4; ++j) {
			nparsd_(arc + ((i__1 = j - 1) < 6 && 0 <= i__1 ? i__1 
				: s_rnge("arc", i__1, "rdffpl_", (ftnlen)857))
				 * 40, &ard[(i__2 = j - 2) < 3 && 0 <= i__2 ? 
				i__2 : s_rnge("ard", i__2, "rdffpl_", (ftnlen)
				857)], error, &ptr, (ftnlen)40, (ftnlen)160);
			if (ptr != 0) {
			    setmsg_("D.P. error (#) in line #.", (ftnlen)25);
			    errch_("#", error, (ftnlen)1, (ftnlen)160);
			    errint_("#", &nrec, (ftnlen)1);
			    sigerr_("SPICE(BADDOUBLEPRECISION)", (ftnlen)25);
			    chkout_("RDFFPL", (ftnlen)6);
			    return 0;
			}
		    }
		    vrtces[*nv * 3 - 3] = ard[0];
		    vrtces[*nv * 3 - 2] = ard[1];
		    vrtces[*nv * 3 - 1] = ard[2];
		} else if (eqstr_(arc, "f", (ftnlen)40, (ftnlen)1)) {

/*                 Plate-vertex data. Process a line of three integers. */
/*                 The algorithm includes an implicit assumption that */
/*                 the data order defines the plate IDs. If not, */
/*                 the output plate files will be useless. */

		    ++(*np);
		    if (*np > 32000000) {
			setmsg_("Number of plates # exceeds limit #.", (
				ftnlen)35);
			errint_("#", np, (ftnlen)1);
			errint_("#", &c_b22, (ftnlen)1);
			sigerr_("SPICE(TOOMANYPLATES)", (ftnlen)20);
			chkout_("RDFFPL", (ftnlen)6);
			return 0;
		    }
		    for (j = 2; j <= 4; ++j) {
			nparsi_(arc + ((i__1 = j - 1) < 6 && 0 <= i__1 ? i__1 
				: s_rnge("arc", i__1, "rdffpl_", (ftnlen)901))
				 * 40, &ari[(i__2 = j - 2) < 4 && 0 <= i__2 ? 
				i__2 : s_rnge("ari", i__2, "rdffpl_", (ftnlen)
				901)], error, &ptr, (ftnlen)40, (ftnlen)160);
			if (ptr != 0) {
			    setmsg_("Integer error (#) in line #.", (ftnlen)
				    28);
			    errch_("#", error, (ftnlen)1, (ftnlen)160);
			    errint_("#", &nrec, (ftnlen)1);
			    sigerr_("SPICE(BADINTEGER)", (ftnlen)17);
			    chkout_("RDFFPL", (ftnlen)6);
			    return 0;
			}
		    }
		    plates[*np * 3 - 3] = ari[0];
		    plates[*np * 3 - 2] = ari[1];
		    plates[*np * 3 - 1] = ari[2];
		} else {

/*                 This block executes if the first characters of */
/*                 a type 3 data lacks a 'v' or 'f'. */

		    setmsg_("Bad data line at record #.", (ftnlen)26);
		    errint_("#", &nrec, (ftnlen)1);
		    sigerr_("SPICE(BADDATALINE)", (ftnlen)18);
		    chkout_("RDFFPL", (ftnlen)6);
		    return 0;
		}
		++nrec;
	    }
	}

/*        All type 3 data read. We should have a non zero NP and NV, */
/*        if not, something failed. */

	if (*np == 0 || *nv == 0) {
	    setmsg_("Read error, type 3 data file. Num plates = #1, num vert"
		    "s = #2.Both should be non-zero", (ftnlen)85);
	    errint_("#1", np, (ftnlen)2);
	    errint_("#2", nv, (ftnlen)2);
	    sigerr_("SPICE(DATAREADFAILED)", (ftnlen)21);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
    } else if (*plttyp == 4) {

/*        We have a Rosetta Osiris style ".ver" file. */

/*        Get the vertex and plate count from the first non-blank */
/*        input line. */

	rdnbl_(infile__, iline, &eof, infile_len, (ftnlen)255);
	if (failed_()) {
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
	lparsm_(iline, " ,", &c__2, &ntok, tokens, (ftnlen)255, (ftnlen)2, (
		ftnlen)80);
	if (ntok != 2) {
	    setmsg_("Vertex and plate count were expected on first line of i"
		    "nput file. Line was #.", (ftnlen)77);
	    errch_("#", iline, (ftnlen)1, (ftnlen)255);
	    sigerr_("SPICE(UNRECOGNIZEDFORMAT)", (ftnlen)25);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
	prsint_(tokens, nv, (ftnlen)80);
	prsint_(tokens + 80, np, (ftnlen)80);
	if (failed_()) {
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
	if (*nv > 16000002) {
	    setmsg_("Number of vertices # exceeds limit #.", (ftnlen)37);
	    errint_("#", nv, (ftnlen)1);
	    errint_("#", &c_b8, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYVERTICES)", (ftnlen)22);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}
	if (*np > 32000000) {
	    setmsg_("Number of plates # exceeds limit #.", (ftnlen)35);
	    errint_("#", np, (ftnlen)1);
	    errint_("#", &c_b22, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYPLATES)", (ftnlen)20);
	    chkout_("RDFFPL", (ftnlen)6);
	    return 0;
	}

/*        Read the vertex data and store it in the output vertex array. */

	i__1 = *nv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rdnbl_(infile__, iline, &eof, infile_len, (ftnlen)255);
	    if (failed_()) {
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	    if (eof) {
		setmsg_("Expected to find # vertices in input file # but ran"
			" out of data after reading # lines of vertex data.", (
			ftnlen)101);
		errint_("#", nv, (ftnlen)1);
		errch_("#", infile__, (ftnlen)1, infile_len);
		errint_("#", &i__, (ftnlen)1);
		sigerr_("SPICE(FILETRUNCATED)", (ftnlen)20);
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }

/*           Parse the input line; we expect it to contain the */
/*           components of the Ith vertex. */

	    lparsm_(iline, " ,", &c__3, &ntok, tokens, (ftnlen)255, (ftnlen)2,
		     (ftnlen)80);
	    if (ntok != 3) {
		setmsg_("Three vertex components were expected on current li"
			"ne of input file. Line was #.", (ftnlen)80);
		errch_("#", iline, (ftnlen)1, (ftnlen)255);
		sigerr_("SPICE(UNRECOGNIZEDFORMAT)", (ftnlen)25);
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	    prsdp_(tokens, &vrtces[i__ * 3 - 3], (ftnlen)80);
	    prsdp_(tokens + 80, &vrtces[i__ * 3 - 2], (ftnlen)80);
	    prsdp_(tokens + 160, &vrtces[i__ * 3 - 1], (ftnlen)80);
	    if (failed_()) {
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	}

/*        Read the plate indices. Discard the first */
/*        and every second non-blank line of plate data; these */
/*        lines all contain the same vertex count (3). */

	j = 0;
	i__1 = *np << 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rdnbl_(infile__, iline, &eof, infile_len, (ftnlen)255);
	    if (failed_()) {
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	    if (eof) {
		setmsg_("Expected to find # plates in input file # but ran o"
			"ut of data after reading # lines of plate data.", (
			ftnlen)98);
		errint_("#", np, (ftnlen)1);
		errch_("#", infile__, (ftnlen)1, infile_len);
		errint_("#", &i__, (ftnlen)1);
		sigerr_("SPICE(FILETRUNCATED)", (ftnlen)20);
		chkout_("RDFFPL", (ftnlen)6);
		return 0;
	    }
	    if (even_(&i__)) {

/*              This line should contain the components of the Jth plate. */

		++j;

/*              Parse the current line and store the plate components. */

		lparsm_(iline, " ,", &c__3, &ntok, tokens, (ftnlen)255, (
			ftnlen)2, (ftnlen)80);
		if (ntok != 3) {
		    setmsg_("Three plate components were expected on current"
			    " line of input file. Line was #.", (ftnlen)79);
		    errch_("#", iline, (ftnlen)1, (ftnlen)255);
		    sigerr_("SPICE(UNRECOGNIZEDFORMAT)", (ftnlen)25);
		    chkout_("RDFFPL", (ftnlen)6);
		    return 0;
		}
		prsint_(tokens, &plates[j * 3 - 3], (ftnlen)80);
		prsint_(tokens + 80, &plates[j * 3 - 2], (ftnlen)80);
		prsint_(tokens + 160, &plates[j * 3 - 1], (ftnlen)80);
		if (failed_()) {
		    chkout_("RDFFPL", (ftnlen)6);
		    return 0;
		}
	    }
	}
    } else {
	setmsg_("Unkown plate data type: #.", (ftnlen)26);
	errint_("#", plttyp, (ftnlen)1);
	sigerr_("SPICE(BADDATATYPE)", (ftnlen)18);
	chkout_("RDFFPL", (ftnlen)6);
	return 0;
    }
    tostdo_("...Done reading plate model input file.", (ftnlen)39);
    tostdo_(" ", (ftnlen)1);
    chkout_("RDFFPL", (ftnlen)6);
    return 0;
} /* rdffpl_ */

