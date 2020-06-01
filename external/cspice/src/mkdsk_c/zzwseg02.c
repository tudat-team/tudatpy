/* zzwseg02.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b17 = 1.;
static integer c_b22 = 16000002;
static integer c_b23 = 32000000;
static integer c_b29 = 60000000;
static integer c_b30 = 16000000;
static integer c_b31 = 68000000;
static integer c_b32 = 164100012;
static integer c__3 = 3;
static integer c_b38 = 100000000;

/* $Procedure    ZZWSEG02 ( MKDSK, write type 2 DSK segment ) */
/* Subroutine */ int zzwseg02_(char *infile__, integer *handle, ftnlen 
	infile_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal dval, last;
    static integer i__;
    extern /* Subroutine */ int zzvoxscl_(doublereal *, doublereal *, integer 
	    *, integer *, integer *, doublereal *);
    static doublereal scale;
    static char frame[32];
    extern /* Subroutine */ int zztrgnvx_(integer *, integer *, integer *), 
	    chkin_(char *, ftnlen), getp02_(integer *, doublereal *, integer *
	    );
    static integer cells[120000000]	/* was [2][60000000] */;
    extern /* Subroutine */ int zzpsxtnt_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *), dskw02_(integer *, 
	    integer *, integer *, integer *, char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    ftnlen), movei_(integer *, integer *, integer *);
    static integer dtype;
    static doublereal first;
    extern logical eqstr_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dskrb2_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *), 
	    dskmi2_(integer *, doublereal *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *, integer *, logical *, integer 
	    *, integer *, doublereal *, integer *);
    static doublereal mncor1, mncor2, mncor3, mxcor1, mxcor2, mxcor3;
    extern logical failed_(void);
    static integer np, nv, centid, dclass, cgrscl;
    extern /* Subroutine */ int getgen_(integer *, integer *, char *, 
	    doublereal *, doublereal *, integer *, integer *, char *, char *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, logical *, ftnlen, ftnlen, ftnlen), rdffpl_(char *,
	     integer *, integer *, doublereal *, integer *, integer *, ftnlen)
	    , mkgrid_(char *, integer *, char *, char *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern logical return_(void);
    static char aunits[255], dunits[255];
    static doublereal avplex, corpar[10], extent[6]	/* was [2][3] */, 
	    spaixd[10], voxscl, vrtces[48000006]	/* was [3][16000002] 
	    */;
    static integer corsys, nvxptr, nvxtot, plates[96000000]	/* was [3][
	    32000000] */, plttyp, spaixi[164100012], surfid, trgcor, trgfin, 
	    vgrext[3];
    static logical makvpm;
    extern /* Subroutine */ int chkout_(char *, ftnlen), setmsg_(char *, 
	    ftnlen), errint_(char *, integer *, ftnlen), sigerr_(char *, 
	    ftnlen), tostdo_(char *, ftnlen), convrt_(doublereal *, char *, 
	    char *, doublereal *, ftnlen, ftnlen), vsclip_(doublereal *, 
	    doublereal *);

/* $ Abstract */

/*     Write a type 2 segment to a DSK file. */

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
/*     INFILE     I   Name of plate model input file. */
/*     HANDLE     I   DAS file handle of DSK. */

/* $ Detailed_Input */

/*     INFILE         is the name of the MKDSK input plate data file */
/*                    from which a type 2 DSK segment is to be created. */
/*                    See the MKDSK User's Guide for details. */

/*     HANDLE         is the handle of a DSK that is open for writing. */
/*                    The file is left open by this routine and must be */
/*                    closed by the calling application. */

/* $ Detailed_Output */

/*     None. */

/*     This routine operates by side effects. See Particulars */
/*     for details. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     This routine is meant to be operated in RETURN SPICE error */
/*     handling mode. The caller is expected to delete the DSK file if */
/*     an error occurs during file creation. */


/*     1) If the setup file requests creation of a segment having */
/*        an unrecognized data type, the error SPICE(NOTSUPPORTED) */
/*        is signaled. */

/*     2) If a new DSK having the specified name cannot be created, */
/*        the error will be diagnosed by routines in the call tree */
/*        of this routine. */

/*     3) If an error occurs while writing comments to the DSK file, */
/*        the error will be diagnosed by routines in the call tree */
/*        of this routine. */

/*     4) If an error occurs while writing data to the DSK file, */
/*        the error will be diagnosed by routines in the call tree */
/*        of this routine. */

/*     5) If an error is present in the setup file, the error will */
/*        be diagnosed, if possible, by routines in the call tree */
/*        of this routine. */

/* $ Files */

/*     Input */
/*     ----- */

/*     1) This routine expects the input plate data file designated */
/*        by INFILE to conform to an expected format. The supported */
/*        formats are described in the MKDSK user's guide. */

/*     2) This routine assumes that a MKDSK setup file has already */
/*        been loaded into the kernel pool. */

/*     Output */
/*     ------ */

/*     This routine writes a type 2 DSK segment to a DSK file that */
/*     has been opened for write access. */

/* $ Particulars */

/*     This routine parses an input plate data file, creates */
/*     a spatial index for the plate model, and writes a */
/*     type 2 DSK segment to the file designated by HANDLE. */

/* $ Examples */

/*     See usage in MKDSK. */

/* $ Restrictions */

/*     This routine should be called only from within MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 5.0.0, 24-FEB-2017 (NJB) */

/*        Added automatic voxel scale determination, for the */
/*        case where scales are not specified in the setup file. */

/*        Now includes mkdsk02.inc rather than dsk02.inc. */

/*        Last update 17-AUG-2016 (NJB) */

/*           Now supports planetodetic and rectangular coordinates. */
/*           Now supports unit conversion for vertices. */

/* -    MKDSK Version 4.0.0, 30-DEC-2015 (NJB) */

/*        Re-written to make use of the SPICELIB DSK type 2 */
/*        spatial index creation routine DSKMI2. */

/*        The minimum radius of the plate set is now computed */
/*        rather than estimated. */

/*        Updated to use new ZZMKSPIN interface. This routine */
/*        no longer fills in the coarse voxel grid pointer array, */
/*        since ZZMKSPIN performs that task. */

/*        Now imports parameter declarations for pointer and cell */
/*        array bounds from mkdsk.inc. */


/* -    MKDSK Version 3.0.0, 06-AUG-2014 (NJB) */

/*        No longer passes PID and VID arrays to RDFFPL. */
/*        No longer uses arrays */

/*          V1, V2, V3, X, Y, Z */

/*        No longer declares workspace array PNTRS, since this */
/*        array is no longer needed by ZZMKSPIN. */

/* -    MKDSK Version 2.0.0, 05-MAY-2014 (NJB) */

/*        Now does not compute the vertex plate map. The call */
/*        to DSKW02 sets the vertex-plate map array size to 0. */

/*        Last update was Version 2.0.0, 03-MAY-2014 (NJB) */

/*        Now has improved error checking. Calls ZZ* routines rather */
/*        than original versions. */

/* -    SPICELIB Version 1.0.0, 29-JUN-2010 (NJB) */

/*        Removed variable VERTS, which occupied 9*MAXPLT d.p. */
/*        numbers; in other words, over 0.5Gb of memory. */
/*        Re-wrote calls to PLCALC and MKSPIN accordingly. */

/*        Changed dimension of VTXPTR from MAXNPV to MAXVRT; */
/*        given the current values of these parameters, namely */
/*        ~24e6 and ~4e6, this saves about 80Mb of memory. */

/*        Declared the arrays CELLS and PNTRS here rather than */
/*        separately in MKSPIN and VRTCOM; this saves about */
/*        3*MAXNPV integers, or about 290Mb of memory. Re-wrote */
/*        calls to MKSPIN and VRTCOM accordingly. */


/* -& */
/* $ Index_Entries */

/*     write type 2 segment to dsk file */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Maximum size of vertex-plate list: */


/*     Size of integer part of spatial index: */


/*     Local variables */


/*     Saved variables */

/*     Save all variables to minimize stack overflow problems. */

    if (return_()) {
	return 0;
    }
    chkin_("ZZWSEG02", (ftnlen)8);

/*     Extract general DSK parameters from the setup file. */

    getgen_(&surfid, &centid, frame, &first, &last, &dclass, &dtype, aunits, 
	    dunits, &corsys, corpar, &mncor1, &mxcor1, &mncor2, &mxcor2, &
	    makvpm, (ftnlen)32, (ftnlen)255, (ftnlen)255);
    if (failed_()) {
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }

/*     Check the data class. */

    if (dclass < 1 || dclass > 2) {
	setmsg_("Data class was #. The only supported values are 1 and 2.", (
		ftnlen)56);
	errint_("#", &dclass, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }

/*     Verify the segment type. */

    if (dtype != 2) {
	setmsg_("Segment type must be 2 but actually was #.", (ftnlen)42);
	errint_("#", &dtype, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }

/*     Get type 2 specific parameters. If the coarse voxel scale */
/*     is set to zero, compute the scales. */

    getp02_(&plttyp, &voxscl, &cgrscl);
    if (failed_()) {
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }

/*     Read the plate model input file. */

    tostdo_("Reading plate model input file...", (ftnlen)33);
    if (plttyp < 5) {

/*        The input file contains plates and vertices. */

	rdffpl_(infile__, &plttyp, &nv, vrtces, &np, plates, infile_len);
	if (failed_()) {
	    chkout_("ZZWSEG02", (ftnlen)8);
	    return 0;
	}

/*        Convert vertices to have units of km, if necessary. */

	if (! eqstr_(dunits, "KM", (ftnlen)255, (ftnlen)2) && ! eqstr_(dunits,
		 "KILOMETERS", (ftnlen)255, (ftnlen)10)) {

/*           Let SCALE be the number of km equivalent to one input unit. */

	    convrt_(&c_b17, dunits, "KM", &scale, (ftnlen)255, (ftnlen)2);
	    if (failed_()) {
		chkout_("ZZWSEG02", (ftnlen)8);
		return 0;
	    }
	    i__1 = nv;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vsclip_(&scale, &vrtces[(i__2 = i__ * 3 - 3) < 48000006 && 0 
			<= i__2 ? i__2 : s_rnge("vrtces", i__2, "zzwseg02_", (
			ftnlen)394)]);
	    }
	}
    } else {

/*        The input file contains a height grid. */

	mkgrid_(infile__, &plttyp, aunits, dunits, &corsys, corpar, &c_b22, &
		c_b23, &nv, vrtces, &np, plates, infile_len, (ftnlen)255, (
		ftnlen)255);

/*        The vertices output from MKGRID always have units of km. */

    }
    if (failed_()) {
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }
    if (cgrscl == 0) {

/*        Compute vertex and plate extents. */

	zzpsxtnt_(&nv, vrtces, &np, plates, extent, &avplex);

/*        Compute target coarse and fine voxel counts. */

	zztrgnvx_(&np, &trgcor, &trgfin);
	if (failed_()) {
	    chkout_("ZZWSEG02", (ftnlen)8);
	    return 0;
	}

/*        Compute coarse and fine voxel scales. */

	zzvoxscl_(extent, &avplex, &trgcor, &trgfin, &cgrscl, &voxscl);
	if (failed_()) {
	    tostdo_("Could not generate voxel scales. Set voxel scales in se"
		    "tup file.", (ftnlen)64);
	    chkout_("ZZWSEG02", (ftnlen)8);
	    return 0;
	}
    }

/*     Generate the spatial index. */

    tostdo_("Generating Spatial Index...", (ftnlen)27);
    dskmi2_(&nv, vrtces, &np, plates, &voxscl, &cgrscl, &c_b29, &c_b30, &
	    c_b31, &makvpm, &c_b32, cells, spaixd, spaixi);
    if (failed_()) {
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }
    movei_(spaixi, &c__3, vgrext);
    nvxtot = vgrext[0] * vgrext[1] * vgrext[2];
    if (nvxtot > 100000000) {
	setmsg_("NVXTOT (#) should be smaller than MAXVOX (#).", (ftnlen)45);
	errint_("#", &nvxtot, (ftnlen)1);
	errint_("#", &c_b38, (ftnlen)1);
	sigerr_("SPICE(VOXELGRIDTOOBIG)", (ftnlen)22);
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }
    nvxptr = spaixi[4];
    if (nvxptr > 100000000) {
	setmsg_("NVXPTR (#) should be smaller than MAXVXP (#).", (ftnlen)45);
	errint_("#", &nvxptr, (ftnlen)1);
	errint_("#", &c_b30, (ftnlen)1);
	sigerr_("SPICE(POINTERSETTOOBIG)", (ftnlen)23);
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }

/*     Convert units of coordinate parameters, if necessary. */

    if (corsys == 1 || corsys == 4) {

/*        Convert input coordinate bounds to radians and km. */

	convrt_(&mncor1, aunits, "RADIANS", &dval, (ftnlen)255, (ftnlen)7);
	mncor1 = dval;
	convrt_(&mxcor1, aunits, "RADIANS", &dval, (ftnlen)255, (ftnlen)7);
	mxcor1 = dval;
	convrt_(&mncor2, aunits, "RADIANS", &dval, (ftnlen)255, (ftnlen)7);
	mncor2 = dval;
	convrt_(&mxcor2, aunits, "RADIANS", &dval, (ftnlen)255, (ftnlen)7);
	mxcor2 = dval;
	if (corsys == 4) {

/*           Convert equatorial radius to km. */

	    convrt_(&corpar[6], dunits, "KM", &dval, (ftnlen)255, (ftnlen)2);
	    corpar[6] = dval;
	}
    } else if (corsys == 3) {
	convrt_(&mncor1, dunits, "KM", &dval, (ftnlen)255, (ftnlen)2);
	mncor1 = dval;
	convrt_(&mxcor1, dunits, "KM", &dval, (ftnlen)255, (ftnlen)2);
	mxcor1 = dval;
	convrt_(&mncor2, dunits, "KM", &dval, (ftnlen)255, (ftnlen)2);
	mncor2 = dval;
	convrt_(&mxcor2, dunits, "KM", &dval, (ftnlen)255, (ftnlen)2);
	mxcor2 = dval;
    } else {
	setmsg_("Coordinate system # is not supported.", (ftnlen)37);
	errint_("#", &corsys, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("ZZWSEG02", (ftnlen)8);
	return 0;
    }

/*     Compute bounds on the third coordinate. */

    dskrb2_(&nv, vrtces, &np, plates, &corsys, corpar, &mncor3, &mxcor3);

/*     Write a type 2 (plate model) segment to the DSK file. */

    dskw02_(handle, &centid, &surfid, &dclass, frame, &corsys, corpar, &
	    mncor1, &mxcor1, &mncor2, &mxcor2, &mncor3, &mxcor3, &first, &
	    last, &nv, vrtces, &np, plates, spaixd, spaixi, (ftnlen)32);
    chkout_("ZZWSEG02", (ftnlen)8);
    return 0;
} /* zzwseg02_ */

