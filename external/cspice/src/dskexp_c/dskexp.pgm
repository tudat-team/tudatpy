/* dskexp.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__8 = 8;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__3 = 3;

/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static char keys[40*8] = "-dsk                                    " "-te"
	    "xt                                   " "-format                 "
	    "                " "-prec                                   " 
	    "-usage                                  " "-u                  "
	    "                    " "-help                                   " 
	    "-h                                      ";

    /* System generated locals */
    address a__1[2], a__2[3];
    integer i__1, i__2[2], i__3[3], i__4;
    doublereal d__1, d__2;
    char ch__1[273], ch__2[68], ch__3[279], ch__4[60], ch__5[59], ch__6[57], 
	    ch__7[288];
    cilist ci__1;
    icilist ici__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer i_dnnt(doublereal *);
    double d_lg10(doublereal *);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_wsfe(cilist *), e_wsfe(void), f_clos(cllist *);

    /* Local variables */
    static integer ndig;
    static char dpfm[50];
    static integer nfil, prec, nseg, unit, i__, j, n;
    extern /* Subroutine */ int chkin_(char *, ftnlen), dskgd_(integer *, 
	    integer *, doublereal *), errch_(char *, char *, ftnlen, ftnlen), 
	    usage_(void), dskp02_(integer *, integer *, integer *, integer *, 
	    integer *, integer *), dskv02_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *);
    static integer segno;
    static logical found;
    extern /* Subroutine */ int repmi_(char *, char *, integer *, char *, 
	    ftnlen, ftnlen, ftnlen), movei_(integer *, integer *, integer *), 
	    dskz02_(integer *, integer *, integer *, integer *);
    static integer dtype;
    extern logical eqstr_(char *, char *, ftnlen, ftnlen);
    static doublereal verts[48000006]	/* was [3][16000002] */;
    static integer dladsc[8], handle;
    extern /* Subroutine */ int dlabfs_(integer *, integer *, logical *);
    static integer np, fmtcde, nv;
    extern /* Subroutine */ int dlafns_(integer *, integer *, integer *, 
	    logical *);
    static char cmd[2000], numfmt[50], dsk[255], numstr[50], outfil[255], 
	    outstr[80], stprec[50], txtfil[255], txtfmt[50];
    static doublereal dskdsc[24];
    static integer nxtdsc[8], plates[96000000]	/* was [3][32000000] */;
    extern /* Subroutine */ int getcml_(char *, ftnlen), kxtrct_(char *, char 
	    *, integer *, char *, logical *, char *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), byebye_(char *, ftnlen), prsint_(char *, integer *, 
	    ftnlen), setmsg_(char *, ftnlen), sigerr_(char *, ftnlen), 
	    tostdo_(char *, ftnlen), dasopr_(char *, integer *, ftnlen), 
	    errint_(char *, integer *, ftnlen), suffix_(char *, integer *, 
	    char *, ftnlen, ftnlen), txtopn_(char *, integer *, ftnlen), 
	    chkout_(char *, ftnlen);

/* $ Abstract */

/*     Export a SPICE DSK file to a set of one or more text */
/*     data files. The output files are suitable for input to */
/*     MKDSK. */

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
/*     PCK */
/*     TIME */

/* $ Keywords */

/*     FILES */
/*     TOPOGRAPHY */

/* $ Files */

/*     Inputs:   DSK file. The file should contain only type 2 */
/*               segments. */

/*     Outputs:  A set of text files, one for each type 2 segment */
/*               in the input file. */


/* $ Particulars */

/*     This is a command-line program. The command syntax is */

/*         dskexp    -dsk <dsk> */
/*                   -text <output name> */
/*                   -format <MKDSK format code/name> */
/*                 [ -prec <# of vertex mantissa digits (1:17)> ] */

/*     The supported format codes are */

/*        MKDSK format code 1: "plate-vertex" */
/*        Gaskell format (code 2) is not supported */
/*        MKDSK format code 3: "vertex-facet" or "obj" */
/*        MKDSK format code 4: "ver" */

/*     The MKDSK user's guide documents these formats. */


/*     A syntax summary can be dumped by typing any of the */
/*     commands */

/*        dskexp -usage */
/*        dskexp -u */
/*        dskexp -help */
/*        dskexp -h */

/* $ Examples */

/*     The commands shown below create a text output file from */
/*     an input DSK called */

/*         phobos.bds */


/*     1)  Create a vertex-facet format output file. This example uses */
/*         default precision for the output vertices. */

/*           dskexp -dsk phobos.bds -text phobos.obj -format vertex-facet */

/*     2)  Create a vertex-facet format output file. Use 9-digit */
/*         mantissas for the vertices. The format name "obj" can be used */
/*         to indicate vertex-facet format. */

/*           dskexp -dsk phobos.bds -text phobos.obj -format obj -prec 9 */

/*     3)  Create a vertex-facet format output file. Use 9-digit */
/*         mantissas for the vertices. The format code 3 can be used to */
/*         indicate vertex-facet format. */

/*           dskexp -dsk phobos.bds -text phobos.obj -format 3 -prec 9 */

/*     4)  Create a plate-vertex format output file. This example uses */
/*         default precision for the output vertices. */

/*           dskexp -dsk phobos.bds -text phobos.obj -format plate-vertex */

/*     5)  Create a plate-vertex format output file. Use the integer */
/*         code for this format. */

/*           dskexp -dsk phobos.bds -text phobos.obj -format 1 */

/*     6)  Create a Rosetta "ver" format output file. */

/*           dskexp -dsk phobos.bds -text phobos.obj -format ver */

/*     7)  Create a Rosetta "ver" format output file. Use the integer */
/*         code for this format. */

/*           dskexp -dsk phobos.bds -text phobos.obj -format 4 */


/* $ Restrictions */

/*     Currently only type 2 DSK segments can be exported. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -    DSKEXP Version 1.0.0, 07-MAR-2017 (NJB) */

/*        Previous verfsoin 07-APR-2015 (NJB) */

/* -& */

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


/*     SPICELIB functions */


/*     Include file dsk02.inc */

/*     This include file declares parameters for DSK data type 2 */
/*     (plate model). */

/* -       SPICELIB Version 1.0.0 08-FEB-2017 (NJB) */

/*          Updated version info. */

/*           22-JAN-2016 (NJB) */

/*              Now includes spatial index parameters. */

/*           26-MAR-2015 (NJB) */

/*              Updated to increase MAXVRT to 16000002. MAXNPV */
/*              has been changed to (3/2)*MAXPLT. Set MAXVOX */
/*              to 100000000. */

/*           13-MAY-2010 (NJB) */

/*              Updated to reflect new no-record design. */

/*           04-MAY-2010 (NJB) */

/*              Updated for new type 2 segment design. Now uses */
/*              a local parameter to represent DSK descriptor */
/*              size (NB). */

/*           13-SEP-2008 (NJB) */

/*              Updated to remove albedo information. */
/*              Updated to use parameter for DSK descriptor size. */

/*           27-DEC-2006 (NJB) */

/*              Updated to remove minimum and maximum radius information */
/*              from segment layout.  These bounds are now included */
/*              in the segment descriptor. */

/*           26-OCT-2006 (NJB) */

/*              Updated to remove normal, center, longest side, albedo, */
/*              and area keyword parameters. */

/*           04-AUG-2006 (NJB) */

/*              Updated to support coarse voxel grid.  Area data */
/*              have now been removed. */

/*           10-JUL-2006 (NJB) */


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


/*     End of include file dsk02.inc */


/*     Local parameters */


/*     Output format codes */


/*     Local variables */


/*     Saved variables */


/*     Initial values */

    chkin_("DSKEXP", (ftnlen)6);
    getcml_(cmd, (ftnlen)2000);

/*     Look for a help request first. */

    for (i__ = 5; i__ <= 8; ++i__) {
	kxtrct_(keys + ((i__1 = i__ - 1) < 8 && 0 <= i__1 ? i__1 : s_rnge(
		"keys", i__1, "dskexp_", (ftnlen)262)) * 40, keys, &c__8, cmd,
		 &found, outstr, (ftnlen)40, (ftnlen)40, (ftnlen)2000, (
		ftnlen)80);
	if (found) {
	    usage_();
	    byebye_("FAILURE", (ftnlen)7);
	}
    }

/*     Look for normal command parameters. */

    kxtrct_(keys, keys, &c__8, cmd, &found, dsk, (ftnlen)40, (ftnlen)40, (
	    ftnlen)2000, (ftnlen)255);
    if (! found) {
	usage_();
	byebye_("FAILURE", (ftnlen)7);
    }
    kxtrct_(keys + 40, keys, &c__8, cmd, &found, txtfil, (ftnlen)40, (ftnlen)
	    40, (ftnlen)2000, (ftnlen)255);
    if (! found) {
	usage_();
	byebye_("FAILURE", (ftnlen)7);
    }
    kxtrct_(keys + 80, keys, &c__8, cmd, &found, txtfmt, (ftnlen)40, (ftnlen)
	    40, (ftnlen)2000, (ftnlen)50);
    if (! found) {
	usage_();
	byebye_("FAILURE", (ftnlen)7);
    }
    kxtrct_(keys + 120, keys, &c__8, cmd, &found, stprec, (ftnlen)40, (ftnlen)
	    40, (ftnlen)2000, (ftnlen)50);
    if (found) {
	prsint_(stprec, &prec, (ftnlen)50);
	if (prec < 1 || prec > 17) {
	    setmsg_("Precision must be in the range 1:17.", (ftnlen)36);
	    sigerr_("SPICE(BADPRECVALUE)", (ftnlen)19);
	}
	s_copy(dpfm, "(1PE#.#)", (ftnlen)50, (ftnlen)8);
	i__1 = prec + 9;
	repmi_(dpfm, "#", &i__1, dpfm, (ftnlen)50, (ftnlen)1, (ftnlen)50);
	i__1 = prec - 1;
	repmi_(dpfm, "#", &i__1, dpfm, (ftnlen)50, (ftnlen)1, (ftnlen)50);
    } else {
	s_copy(dpfm, "(1PE25.16)", (ftnlen)50, (ftnlen)10);
    }
    tostdo_(" ", (ftnlen)1);
/* Writing concatenation */
    i__2[0] = 18, a__1[0] = "Input DSK:        ";
    i__2[1] = 255, a__1[1] = dsk;
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)273);
    tostdo_(ch__1, (ftnlen)273);
/* Writing concatenation */
    i__2[0] = 18, a__1[0] = "Output text file: ";
    i__2[1] = 255, a__1[1] = txtfil;
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)273);
    tostdo_(ch__1, (ftnlen)273);
/* Writing concatenation */
    i__2[0] = 18, a__1[0] = "Output format:    ";
    i__2[1] = 50, a__1[1] = txtfmt;
    s_cat(ch__2, a__1, i__2, &c__2, (ftnlen)68);
    tostdo_(ch__2, (ftnlen)68);
    tostdo_(" ", (ftnlen)1);

/*     Identify output format code. */

    if (eqstr_(txtfmt, "1", (ftnlen)50, (ftnlen)1) || eqstr_(txtfmt, "plate-"
	    "vertex", (ftnlen)50, (ftnlen)12)) {
	fmtcde = 1;
    } else if (eqstr_(txtfmt, "3", (ftnlen)50, (ftnlen)1) || eqstr_(txtfmt, 
	    "vertex-facet", (ftnlen)50, (ftnlen)12) || eqstr_(txtfmt, "obj", (
	    ftnlen)50, (ftnlen)3)) {
	fmtcde = 3;
    } else if (eqstr_(txtfmt, "4", (ftnlen)50, (ftnlen)1) || eqstr_(txtfmt, 
	    "ver", (ftnlen)50, (ftnlen)3)) {
	fmtcde = 4;
    } else {
	setmsg_("Output format <#> is not recognized.", (ftnlen)36);
	errch_("#", txtfmt, (ftnlen)1, (ftnlen)50);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
    }

/*     Scan DSK file; make sure all segments have */
/*     supported data types. */

    dasopr_(dsk, &handle, (ftnlen)255);
    segno = 0;
    nfil = 0;
    dlabfs_(&handle, nxtdsc, &found);
    if (! found) {
	setmsg_("No segments were found in DSK file #.", (ftnlen)37);
	errch_("#", dsk, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NODSKSEGMENT)", (ftnlen)19);
    }

/*     Process the segments of the file. */

    while(found) {
	movei_(nxtdsc, &c__8, dladsc);
	++segno;
	dskgd_(&handle, dladsc, dskdsc);
	dtype = i_dnnt(&dskdsc[3]);

/*        Right now type 2 is all we do. */

	if (dtype != 2) {
	    setmsg_("Segment # in DSK file # has data type #.", (ftnlen)40);
	    errint_("#", &segno, (ftnlen)1);
	    errch_("#", dsk, (ftnlen)1, (ftnlen)255);
	    errint_("#", &dtype, (ftnlen)1);
	    sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	}

/*        Find the next segment. */

	dlafns_(&handle, dladsc, nxtdsc, &found);
    }

/*     Create the output format for the file number suffix, if */
/*     it's needed. */

    nseg = segno;
    if (nseg > 1) {
	d__2 = (doublereal) nseg;
	d__1 = d_lg10(&d__2);
	ndig = i_dnnt(&d__1) + 1;
	s_copy(numfmt, "(I@.@)", (ftnlen)50, (ftnlen)6);
	for (i__ = 1; i__ <= 2; ++i__) {
	    repmi_(numfmt, "@", &ndig, numfmt, (ftnlen)50, (ftnlen)1, (ftnlen)
		    50);
	}
    }

/*     Scan the input file; write data to the output file. */

    segno = 0;
    dlabfs_(&handle, nxtdsc, &found);
    while(found) {
	++segno;
	movei_(nxtdsc, &c__8, dladsc);
	tostdo_(" ", (ftnlen)1);
	s_copy(outstr, "Processing segment #", (ftnlen)80, (ftnlen)20);
	repmi_(outstr, "#", &segno, outstr, (ftnlen)80, (ftnlen)1, (ftnlen)80)
		;
	tostdo_(outstr, (ftnlen)80);
	dskgd_(&handle, dladsc, dskdsc);
	dtype = i_dnnt(&dskdsc[3]);
	if (dtype == 2) {

/*           Making it here means we can try to open the output file. */

	    ++nfil;
	    if (nseg > 1) {
		s_copy(outfil, txtfil, (ftnlen)255, (ftnlen)255);
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 50;
		ici__1.iciunit = numstr;
		ici__1.icifmt = numfmt;
		s_wsfi(&ici__1);
		do_fio(&c__1, (char *)&nfil, (ftnlen)sizeof(integer));
		e_wsfi();
		suffix_("_", &c__0, outfil, (ftnlen)1, (ftnlen)255);
		suffix_(numstr, &c__0, outfil, (ftnlen)50, (ftnlen)255);
	    } else {
		s_copy(outfil, txtfil, (ftnlen)255, (ftnlen)255);
	    }
	    txtopn_(outfil, &unit, (ftnlen)255);
	    tostdo_("  Segment is type 2", (ftnlen)19);
	    dskz02_(&handle, dladsc, &nv, &np);
	    s_copy(outstr, "   Number of vertices: #", (ftnlen)80, (ftnlen)24)
		    ;
	    repmi_(outstr, "#", &nv, outstr, (ftnlen)80, (ftnlen)1, (ftnlen)
		    80);
	    tostdo_(outstr, (ftnlen)80);
	    s_copy(outstr, "   Number of plates:   #", (ftnlen)80, (ftnlen)24)
		    ;
	    repmi_(outstr, "#", &np, outstr, (ftnlen)80, (ftnlen)1, (ftnlen)
		    80);
	    tostdo_(outstr, (ftnlen)80);
	    tostdo_("   Reading vertices...", (ftnlen)22);
	    dskv02_(&handle, dladsc, &c__1, &nv, &n, verts);
	    tostdo_("    Reading plates...", (ftnlen)21);
	    dskp02_(&handle, dladsc, &c__1, &np, &n, plates);
/* Writing concatenation */
	    i__2[0] = 24, a__1[0] = "    Writing output file ";
	    i__2[1] = 255, a__1[1] = outfil;
	    s_cat(ch__3, a__1, i__2, &c__2, (ftnlen)279);
	    tostdo_(ch__3, (ftnlen)279);

/*           Format code 1: plate-vertex table. */

	    if (fmtcde == 1) {
		tostdo_("    Writing vertices...", (ftnlen)23);
		ci__1.cierr = 0;
		ci__1.ciunit = unit;
		ci__1.cifmt = "(1X,I10)";
		s_wsfe(&ci__1);
		do_fio(&c__1, (char *)&nv, (ftnlen)sizeof(integer));
		e_wsfe();
		i__1 = nv;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ci__1.cierr = 0;
		    ci__1.ciunit = unit;
/* Writing concatenation */
		    i__3[0] = 9, a__2[0] = "(1X,I11,3";
		    i__3[1] = 50, a__2[1] = dpfm;
		    i__3[2] = 1, a__2[2] = ")";
		    ci__1.cifmt = (s_cat(ch__4, a__2, i__3, &c__3, (ftnlen)60)
			    , ch__4);
		    s_wsfe(&ci__1);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&verts[(i__4 = j + i__ * 3 - 4) 
				< 48000006 && 0 <= i__4 ? i__4 : s_rnge("ver"
				"ts", i__4, "dskexp_", (ftnlen)505)], (ftnlen)
				sizeof(doublereal));
		    }
		    e_wsfe();
		}
		tostdo_("    Writing plates...", (ftnlen)21);
		ci__1.cierr = 0;
		ci__1.ciunit = unit;
		ci__1.cifmt = "(1X,I10)";
		s_wsfe(&ci__1);
		do_fio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
		e_wsfe();
		i__1 = np;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ci__1.cierr = 0;
		    ci__1.ciunit = unit;
		    ci__1.cifmt = "(1X,4I10)";
		    s_wsfe(&ci__1);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&plates[(i__4 = j + i__ * 3 - 4)
				 < 96000000 && 0 <= i__4 ? i__4 : s_rnge(
				"plates", i__4, "dskexp_", (ftnlen)515)], (
				ftnlen)sizeof(integer));
		    }
		    e_wsfe();
		}

/*           Format code 3: vertex-facet table. Also called */
/*           "obj" format. */

	    } else if (fmtcde == 3) {
		tostdo_("    Writing vertices...", (ftnlen)23);
		i__1 = nv;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ci__1.cierr = 0;
		    ci__1.ciunit = unit;
/* Writing concatenation */
		    i__3[0] = 8, a__2[0] = "(1X,A3,3";
		    i__3[1] = 50, a__2[1] = dpfm;
		    i__3[2] = 1, a__2[2] = ")";
		    ci__1.cifmt = (s_cat(ch__5, a__2, i__3, &c__3, (ftnlen)59)
			    , ch__5);
		    s_wsfe(&ci__1);
		    do_fio(&c__1, " v ", (ftnlen)3);
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&verts[(i__4 = j + i__ * 3 - 4) 
				< 48000006 && 0 <= i__4 ? i__4 : s_rnge("ver"
				"ts", i__4, "dskexp_", (ftnlen)530)], (ftnlen)
				sizeof(doublereal));
		    }
		    e_wsfe();
		}
		tostdo_("    Writing plates...", (ftnlen)21);
		i__1 = np;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ci__1.cierr = 0;
		    ci__1.ciunit = unit;
		    ci__1.cifmt = "(1X,A3,3I10)";
		    s_wsfe(&ci__1);
		    do_fio(&c__1, " f ", (ftnlen)3);
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&plates[(i__4 = j + i__ * 3 - 4)
				 < 96000000 && 0 <= i__4 ? i__4 : s_rnge(
				"plates", i__4, "dskexp_", (ftnlen)538)], (
				ftnlen)sizeof(integer));
		    }
		    e_wsfe();
		}

/*           Format code 4: Rosetta "ver" format. */

	    } else if (fmtcde == 4) {
		ci__1.cierr = 0;
		ci__1.ciunit = unit;
		ci__1.cifmt = "(1X,2I10)";
		s_wsfe(&ci__1);
		do_fio(&c__1, (char *)&nv, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
		e_wsfe();
		tostdo_("    Writing vertices...", (ftnlen)23);
		i__1 = nv;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ci__1.cierr = 0;
		    ci__1.ciunit = unit;
/* Writing concatenation */
		    i__3[0] = 6, a__2[0] = "(1X, 3";
		    i__3[1] = 50, a__2[1] = dpfm;
		    i__3[2] = 1, a__2[2] = ")";
		    ci__1.cifmt = (s_cat(ch__6, a__2, i__3, &c__3, (ftnlen)57)
			    , ch__6);
		    s_wsfe(&ci__1);
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&verts[(i__4 = j + i__ * 3 - 4) 
				< 48000006 && 0 <= i__4 ? i__4 : s_rnge("ver"
				"ts", i__4, "dskexp_", (ftnlen)553)], (ftnlen)
				sizeof(doublereal));
		    }
		    e_wsfe();
		}
		tostdo_("    Writing plates...", (ftnlen)21);
		i__1 = np;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ci__1.cierr = 0;
		    ci__1.ciunit = unit;
		    ci__1.cifmt = "(1X,A1)";
		    s_wsfe(&ci__1);
		    do_fio(&c__1, "3", (ftnlen)1);
		    e_wsfe();
		    ci__1.cierr = 0;
		    ci__1.ciunit = unit;
		    ci__1.cifmt = "(1X,A2,3I10)";
		    s_wsfe(&ci__1);
		    do_fio(&c__1, "  ", (ftnlen)2);
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&plates[(i__4 = j + i__ * 3 - 4)
				 < 96000000 && 0 <= i__4 ? i__4 : s_rnge(
				"plates", i__4, "dskexp_", (ftnlen)563)], (
				ftnlen)sizeof(integer));
		    }
		    e_wsfe();
		}
	    } else {
		setmsg_("Output format <#> is not recognized.", (ftnlen)36);
		errch_("#", txtfmt, (ftnlen)1, (ftnlen)50);
		sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = unit;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
/* Writing concatenation */
	    i__2[0] = 33, a__1[0] = "    Finished writing output file ";
	    i__2[1] = 255, a__1[1] = outfil;
	    s_cat(ch__7, a__1, i__2, &c__2, (ftnlen)288);
	    tostdo_(ch__7, (ftnlen)288);
	}

/*        Find the next segment. */

	dlafns_(&handle, dladsc, nxtdsc, &found);
    }
    tostdo_(" ", (ftnlen)1);
    chkout_("DSKEXP", (ftnlen)6);
    return 0;
} /* MAIN__ */

/* Subroutine */ int usage_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int repmc_(char *, char *, char *, char *, ftnlen,
	     ftnlen, ftnlen, ftnlen);
    char versn[80];
    extern /* Subroutine */ int tostdo_(char *, ftnlen);
    char usestr[80*12];
    extern /* Subroutine */ int tkvrsn_(char *, char *, ftnlen, ftnlen);
    char outstr[80];

    tostdo_(" ", (ftnlen)1);
    s_copy(outstr, "DSKEXP --- Version 1.0.0, March 7, 2017 --- Toolkit Vers"
	    "ion #", (ftnlen)80, (ftnlen)61);
    tkvrsn_("TOOLKIT", versn, (ftnlen)7, (ftnlen)80);
    repmc_(outstr, "#", versn, outstr, (ftnlen)80, (ftnlen)1, (ftnlen)80, (
	    ftnlen)80);
    tostdo_(outstr, (ftnlen)80);
    s_copy(usestr, " ", (ftnlen)80, (ftnlen)1);
    s_copy(usestr + 80, "  Usage: dskexp -dsk <dsk> -text <output name> -for"
	    "mat <MKDSK format code/name>", (ftnlen)80, (ftnlen)79);
    s_copy(usestr + 160, "                 [ -prec <# of vertex mantissa dig"
	    "its (1:17)> ]", (ftnlen)80, (ftnlen)63);
    s_copy(usestr + 240, " ", (ftnlen)80, (ftnlen)1);
    s_copy(usestr + 320, "     MKDSK format code 1: \"plate-vertex\"", (
	    ftnlen)80, (ftnlen)40);
    s_copy(usestr + 400, "     Gaskell format (code 2) is not supported", (
	    ftnlen)80, (ftnlen)45);
    s_copy(usestr + 480, "     MKDSK format code 3: \"vertex-facet\" or \""
	    "obj\"", (ftnlen)80, (ftnlen)49);
    s_copy(usestr + 560, "     MKDSK format code 4: \"ver\"", (ftnlen)80, (
	    ftnlen)31);
    s_copy(usestr + 640, "     MKDSK height grid format (code 5) is not supp"
	    "orted", (ftnlen)80, (ftnlen)55);
    s_copy(usestr + 720, " ", (ftnlen)80, (ftnlen)1);
    s_copy(usestr + 800, "  Example: dskexp -dsk phobos.bds -text phobos.txt"
	    " -format plate-vertex -prec 10", (ftnlen)80, (ftnlen)80);
    s_copy(usestr + 880, " ", (ftnlen)80, (ftnlen)1);
    for (i__ = 1; i__ <= 12; ++i__) {
	tostdo_(usestr + ((i__1 = i__ - 1) < 12 && 0 <= i__1 ? i__1 : s_rnge(
		"usestr", i__1, "usage_", (ftnlen)640)) * 80, (ftnlen)80);
    }
    return 0;
} /* usage_ */

/* Main program alias */ int dskexp_ () { MAIN__ (); return 0; }
