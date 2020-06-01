/* sum04.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__7 = 7;
static integer c__12 = 12;
static integer c__26 = 26;
static integer c__25 = 25;
static integer c__9 = 9;
static integer c__11 = 11;
static integer c__27 = 27;

/* $Procedure   SUM04 ( DSKBRIEF, summarize type 4 segment ) */
/* Subroutine */ int sum04_(integer *handle, integer *dladsc, integer *nsig)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_dnnt(doublereal *), s_rnge(char *, integer, char *, integer);

    /* Local variables */
    doublereal xdsc[7];
    char xstr[132];
    doublereal cent1, cent2;
    integer b, e, i__, n, ibase;
    extern /* Subroutine */ int dskb04_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, doublereal *,
	     doublereal *, doublereal *, doublereal *), dskd04_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *), chkin_(char *, ftnlen), dskgd_(integer *, integer *, 
	    doublereal *);
    integer ngdim;
    extern /* Subroutine */ int repmc_(char *, char *, char *, char *, ftnlen,
	     ftnlen, ftnlen, ftnlen), repmd_(char *, char *, doublereal *, 
	    integer *, char *, ftnlen, ftnlen, ftnlen), repmf_(char *, char *,
	     doublereal *, integer *, char *, char *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    doublereal nvdsc[7];
    extern /* Subroutine */ int repmi_(char *, char *, integer *, char *, 
	    ftnlen, ftnlen, ftnlen);
    doublereal accdsc[7], co1psz, co2psz;
    integer nc, nf;
    doublereal hscale, dpbuff[1], refdsc[7];
    char accstr[132], itpstr[132], outlin[132];
    doublereal dpnfmt, dskdsc[24];
    extern doublereal dpr_(void);
    doublereal itpdsc[7], nulval, prjdsc[7];
    integer corsys, dskfmt, grdims[20]	/* was [2][10] */, nr, numfmt, nw;
    logical nullok;
    extern /* Subroutine */ int tostdo_(char *, ftnlen), setmsg_(char *, 
	    ftnlen), errint_(char *, integer *, ftnlen), sigerr_(char *, 
	    ftnlen), chkout_(char *, ftnlen), dasrdi_(integer *, integer *, 
	    integer *, integer *);
    char fmt1[132];

/* $ Abstract */

/*     Display type 4-specific summary of contents of a SPICE */
/*     Digital Shape Kernel (DSK). */

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
/*     DSKBRIEF */

/* $ Declarations */

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


/*     Include file dsk04.inc */

/*     This include file declares parameters for DSK data type 4 */
/*     (double precision raster values). */

/*        Version 4.1.0 07-OCT-2016 (NJB) */

/*           Corrected some long lines for delivery to DSKBRIEF. */

/*        Version 4.0.0 12-NOV-2012 (NJB) */

/*           Added parameter MARGIN, which is used to handle minor */
/*           round-off discrepancies. */

/*           Changed name of parameter BILCPT to NVBCPT. */

/*        Version 3.0.0 18-SEP-2012 (NJB) */

/*     Each type 4 DSK segment has only a double precision DLA component. */
/*     The segment layout in DAS address space is as follows: */

/*        D.p. layout: */


/*   +---------------------+ */
/*   | DSK descriptor      |  DSKDSZ elements */
/*   +---------------------+ */
/*   | Format version ID   |  Format version ID code */
/*   +---------------------+ */
/*   | NR                  |  number of rows */
/*   +---------------------+ */
/*   | NC                  |  number of pixel columns */
/*   +---------------------+ */
/*   | Coord 1 pixel size  |  size in units of coordinate 1 */
/*   +---------------------+ */
/*   | Coord 2 pixel size  |  size in units of coordinate 2 */
/*   +---------------------+ */
/*   | Center coord 1      |  Coordinates of center of */
/*   +---------------------+  pixel grid */
/*   | Center coord 2      | */
/*   +---------------------+ */
/*   | Nulls ok flag       |  0/1 for nulls not allowed/allowed */
/*   +---------------------+ */
/*   | Null value          |  Integer value indicating "no data" */
/*   +---------------------+ */
/*   | Height scale        |  Height units in km. */
/*   +---------------------+ */
/*   | Numeric format      |  Numeric format, e.g. packed 16-bit integer */
/*   +---------------------+ */
/*   |Interpolation dsc ptr|  Pointer to interpolation descriptor */
/*   +---------------------+ */
/*   | Intercept descr ptr |  Pointer to intercept descriptor */
/*   +---------------------+ */
/*   | Normal descr ptr    |  Pointer to normal vector descriptor */
/*   +---------------------+ */
/*   | Accel descr ptr     |  Pointer to acceleration descriptor */
/*   +---------------------+ */
/*   | Ref surface dsc ptr |  Reference surface descriptor pointer */
/*   +---------------------+ */
/*   | Projection dsc ptr  |  Map projection descriptor pointer */
/*   +---------------------+ */
/*   | Coarse grid ptr     |  Pointer to coarse grid structure */
/*   +---------------------+ */
/*   | Data start ptr      |  Pointer to start of pixel data */
/*   +---------------------+ */
/*   | Interp descr        |  Interpolation method descriptor */
/*   +---------------------+ */
/*   | Intercept descr     |  Intercept method descriptor */
/*   +---------------------+ */
/*   | Normal vec descr    |  Normal vector calc descriptor */
/*   +---------------------+ */
/*   | Accel descr         |  Intercept acceleration descriptor */
/*   +---------------------+ */
/*   | Ref surface descr   |  Reference surface descriptor */
/*   +---------------------+ */
/*   | Map projection descr|  Map projection descriptor */
/*   +---------------------+ */
/*   | Coarse grid struct  |  Coarse grid structure */
/*   +---------------------+ */


/*        Integer layout: */

/*   +---------------------+ */
/*   | Grid dimension count|  Depth of nested grids */
/*   +---------------------+ */
/*   | Grid dimensions     |  Row and column sizes of nested grids */
/*   +---------------------+ */
/*   | Pixel data          |  Pixels containing height values */
/*   +---------------------+ */



/*     Indices of descriptor items at fixed locations */
/*     ============================================== */

/*     The following indices are relative to the base address */
/*     of the double precision DAS address range of the DLA */
/*     segment containing the described type 4 DSK segment. */
/*     Indices are 1-based: index 1 is the first DAS double */
/*     precision number of the segment. */


/*     Index of DSK descriptor: */


/*     DSK descriptor size: */

/*     This local parameter MUST be kept consistent with */
/*     the parameter DSKDSZ which is declared in dskdsc.inc. */


/*     Index of format version: */


/*     Index of row count: */


/*     Index of column count: */


/*     Index of pixel size in direction of coordinate 1: */


/*     Index of pixel size in direction of coordinate 2: */


/*     Index of coordinate 1 of pixel grid center: */


/*     Index of coordinate 2 of pixel grid center: */


/*     Index of flag indicating whether null values are allowed: */


/*     Index of null value indicator: */


/*     Index of height scale: */


/*     Index of numeric format: */


/*     Index of pointer to interpolation algorithm descriptor: */


/*     Index of pointer to intercept algorithm descriptor: */


/*     Index of pointer to normal vector computation descriptor: */


/*     Index of pointer to intercept acceleration descriptor: */


/*     Index of reference surface decriptor pointer: */


/*     Index of projection decriptor pointer: */


/*     Index of pointer to coarse grid structure: */


/*     Index of pointer to pixel data: */


/*     Indices of integer items */
/*     ======================== */

/*     Index of grid dimension count: */


/*     Index of grid dimensions: */


/*     Maximum number of grid dimensions: */


/*     The grid dimension array has fixed size of 2*DMAX integers. */


/*     Double precision item keyword parameters used by fetch routines */
/*     =============================================================== */

/*     DSK descriptor: */


/*     Format version: */


/*     Keyword for row count: */


/*     Keyword for column count: */


/*     Keyword for pixel size in direction of coordinate 1: */


/*     Keyword for pixel size in direction of coordinate 2: */


/*     Keyword for coordinate 1 of pixel grid center: */


/*     Keyword for coordinate 2 of pixel grid center: */


/*     Keyword for flag indicating whether null values are allowed: */


/*     Keyword for null value indicator: */


/*     Keyword for height scale: */


/*     Keyword for numeric format: */


/*     Keyword for pointer to reference surface descriptor: */


/*     Keyword for pointer to map projection descriptor: */


/*     Keyword for pointer to interpolation algorithm descriptor: */


/*     Keyword for pointer to intercept algorithm descriptor: */


/*     Keyword for pointer to normal vector computation descriptor: */


/*     Keyword for pointer to intercept acceleration descriptor: */


/*     Keyword for pointer to coarse grid structure: */


/*     Keyword for pointer to pixel data: */


/*     Keyword for interpolation algorithm descriptor: */


/*     Keyword for intercept algorithm descriptor: */


/*     Keyword for normal vector algorithm descriptor: */


/*     Keyword for coarse grid descriptor: */


/*     Keyword for reference surface descriptor: */


/*     Keyword for map projection descriptor: */


/*     Keyword for intercept acceleration descriptor: */


/*     Integer item keyword parameters used by fetch routines */
/*     =============================================================== */

/*     Note: keyword codes for integer items occupy a range */
/*     distinct from that for keywords for d.p. items. The */
/*     first keyword code for an integer item follows the */
/*     last d.p. code. */


/*     Keyword for nested grid dimension count: */


/*     Keyword for grid dimensions: */

/*     Keyword for pixel data: */

/*     This keyword refers to unpacked data: each 16-bit value is */
/*     returned in an individual integer. The pixel grid has NR rows and */
/*     NC columns. There are two 16-bit pixels per stored integer. The */
/*     data are stored in row-major order. */


/*     Keyword for raw data: */

/*     This keyword refers to stored integers, each containing two */
/*     16-bit data values.  The raw grid has NR rows and */
/*     NC/2 columns. */



/*     Constants */
/*     ========= */

/*     The constant NULL is used to indicate a null pointer or */
/*     undefined value. Since this value is stored in double precision */
/*     DAS words, its type is double precision. */


/*     The constants TRUE and FALSE are numeric codes having */
/*     the meanings associated with their names. */


/*     The constant MAXDSZ is the maximum size of any algorithm */
/*     descriptor used in the type 4 DSK segment structure. */


/*     Numeric format */
/*     ============== */

/*        Codes for 32-bit, packed 16-bit, and packed 8-bit integer */
/*        data: */


/*     Code for equivalenced 32-bit real data: */


/*     Reference surface descriptor */
/*     ============================ */


/*         Reference surface codes */
/*         ----------------------- */

/*            Code REFSPH indicates that the reference surface is a */
/*            sphere. */

/*            For spheres, a single radius is stored in the parameter */
/*            section of the descriptor. */


/*            Code REFELL indicates that the reference surface is a */
/*            triaxial ellipsoid. */

/*            For ellipsoids, three radii are stored in the parameter */
/*            section of the descriptor. */


/*            Code REFINH indicates that the reference surface is a */
/*            triaxial ellipsoid having radii inherited from the segment */
/*            descriptor. Any radii supplied in the reference surface */
/*            descriptor are ignored. */

/*            The reference surface subtype is also ignored when the */
/*            reference surface code is REFINH. */


/*         Reference surface subtype codes */
/*         ------------------------------- */

/*            Code REFRAD indicates that surface heights are */
/*            measured relative to the reference surface */
/*            in the radial direction. */

/*            REFRAD is the usual choice for segments */
/*            using planetocentric coordinates. */


/*            Code REFNRM indicates that surface heights are */
/*            measured relative to the reference surface */
/*            in the outward normal direction. */

/*            REFNRM is the usual choice for segments */
/*            using planetodetic coordinates. */


/*         Reference surface index parameters */
/*         ---------------------------------- */

/*         Reference surface code within descriptor: */


/*         Index of descriptor size within descriptor: */


/*         Index of subtype within descriptor: */


/*         Index of parameters within descriptor: */


/*     Map projection descriptor */
/*     ========================= */


/*         Map projection codes */
/*         -------------------- */

/*            Code PRJEQR indicates that the map projection is */
/*            equirectangular. */


/*            Code PRJSTE indicates that the map projection is */
/*            stereographic. */


/*            Code PRJMOL indicates that the map projection is */
/*            MGS MOLA-style stereographic. In this projection, */
/*            radial distance on the tangent plane is proportional */
/*            to twice the sine of half of the colatitude of the */
/*            projected point. In a conventional stereographic */
/*            projection, the sine function is replaced by the */
/*            tangent function. */


/*         Map projection descriptor index parameters */
/*         ------------------------------------------ */

/*         Index of projection code within descriptor: */


/*         Index of descriptor size within descriptor: */


/*         For stereographic projections (including the MOLA-style */
/*         projection), the Euler angle matrix product */

/*            [TWIST]  [Pi/2 - COLAT]  [LON] */
/*                   3               2      3 */

/*         describes the orientation of the projection axes relative to */
/*         the base frame. Angular units are radians. */

/*         The tangent point of the projection plane corresponds to the */
/*         intersection of the transformed Z-axis with the unit sphere. */
/*         X and Y pixel scales are used to transform direction vectors */
/*         in the base frame to points on the tangent plane. */

/*         Indices of longitude, latitude, and twist within the */
/*         descriptor: */


/*        Indices of projection scales in the X and Y directions: */


/*     Interpolation algorithm descriptor */
/*     ================================== */

/*     The interpolation algorithm descriptor is a variable-size set of */
/*     parameters that specify the interpolation algorithm used to */
/*     provide values of the third coordinate for arbitrary pairs of */
/*     first and second coordinate values. */

/*     The descriptor starts with a code identifying the class of */
/*     algorithm. This is followed by the size of the descriptor. If */
/*     applicable, additional parameters follow the size. */

/*         Interpolation algorithm codes */
/*         ----------------------------- */

/*            Code ITPNO indicates that no interpolation */
/*            is performed. */


/*            Code ITPBIL denotes bilinear interpolation. */


/*            Code ITPNN3 denotes planar interpolation using */
/*            the three nearest neighboring data values. */


/*         Interpolation descriptor index parameters */
/*         ----------------------------------------- */

/*         Index of interpolation code within descriptor: */


/*         Index of descriptor size within descriptor: */


/*         Index of first parameter within descriptor: */


/*     Intercept algorithm descriptor */
/*     ============================== */

/*     The ray-surface intercept algorithm descriptor is a variable-size */
/*     set of parameters that specify the intercept algorithm used for */
/*     the associated segment. */

/*     The descriptor starts with a code identifying the class of */
/*     algorithm. This is followed by the size of the descriptor. If */
/*     applicable, additional parameters follow the size. */

/*         Intercept algorithm codes */
/*         ------------------------- */

/*            Code XCASTP indicates that constant angular */
/*            size stepping is performed. This algorithm is */
/*            applicable only to class 1 surfaces. */


/*            Code XCPSTP indicates that constant orthogonal projection */
/*            length stepping is performed. This algorithm is applicable */
/*            only to surfaces associated with rectangular coordinates. */


/*            Code XCELBD denotes stepping from one cell boundary */
/*            to the next, where a "cell" is region over which */
/*            the same set of data are used for interpolation. */

/*            Examples of cells: */

/*               1) When no interpolation is used, cell boundaries */
/*                  coincide with pixel boundaries. */

/*               2) When bilinear interpolation is used, cell boundaries */
/*                  are curves, with one of the first two coordinates */
/*                  held constant, connecting centers of adjacent pixels. */
/*                  For example, if the coordinate system is latitudinal, */
/*                  such curves have either constant longitude or */
/*                  latitude. */

/*               3) When 3-nearest-neighbor interpolation is used, */
/*                  cells are triangular regions connecting centers */
/*                  of the 3 pixels used for interpolation. */



/*         Intercept descriptor index parameters */
/*         ------------------------------------- */

/*         Index of intercept algorithm code within descriptor: */


/*         Index of descriptor size within descriptor: */


/*         Index of step parameter within descriptor: */


/*         Index of tolerance parameter within descriptor: */


/*     Intercept acceleration algorithm descriptor */
/*     =========================================== */


/*         Intercept acceleration algorithm codes */
/*         -------------------------------------- */


/*         Code ACCNO indicates no acceleration algorithm */
/*         is used. */


/*         Code ACCBIG indicatest that bounds on height */
/*         differences in the COORD1 and COORD2 directions */
/*         are used to calculate a "big" initial step */
/*         during the root bracketing process. These bounds */
/*         are determined during the segment writing process. */



/*         Intercept acceleration descriptor index parameters */
/*         -------------------------------------------------- */

/*         Index of acceleration algorithm code within descriptor: */


/*         Index of descriptor size within descriptor: */


/*         Index of maximum magnitude of height difference */
/*         in the COORD1 direction: */


/*         Index of maximum magnitude of height difference */
/*         in the COORD2 direction: */


/*     Normal vector calculation algorithm descriptor */
/*     ============================================== */


/*     The normal vector algorithm descriptor is a variable-size set of */
/*     parameters that specify the normal vector computation algorithm */
/*     used for the associated segment. */

/*     The descriptor starts with a code identifying the class of */
/*     algorithm. This is followed by the size of the descriptor. If */
/*     applicable, additional parameters follow the size. */

/*         Normal vector algorithm codes */
/*         ----------------------------- */

/*            Code NVBCPT indicates that normal vectors are compatible */
/*            with binlinear interpolation. This algorithm may be used */
/*            only when the interpolation method is bilinear. */

/*            Normal vector directions at boundaries between */
/*            interpolation regions---that is, where at least one pixel */
/*            coordinate is integral---are guaranteed to be compatible */
/*            with one of the adjacent regions, but are otherwise */
/*            unspecified. Note that discontinuity of normal vectors */
/*            should be expected at these boundaries. */



/*            Code NV1PAR indicates that estimates of the 1st order */
/*            partial derivatives of the third coordinate in the coord 1 */
/*            and coord 2 directions are used to estimate the normal */
/*            vector direction. */


/*         Normal vector descriptor index parameters */
/*         ------------------------------------- */

/*         Index of normal vector algorithm code within descriptor: */


/*         Index of descriptor size within descriptor: */


/*         Index of first parameter within descriptor: */


/*     Additional DSK type 4 parameters */
/*     ================================ */


/*     MARGIN is used to determine when computed coordinates are */
/*     close enough to the boundaries of a coverage region to be */
/*     considered to be located within that region. */


/*     End of include file dsk04.inc */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE     I   Handle of DSK file. */
/*     DLADSC     I   DLA descriptor of segment. */
/*     NSIG       I   Number of significant digits in floating point */
/*                    output. */

/* $ Detailed_Input */

/*     HANDLE     is the handle of a DSK file containing a segment */
/*                to be summarized. */

/*     DLADSC     is the DLA descriptor of a segment to be summarized. */

/*     NSIG       is the number of significant digits in floating point */
/*                numeric output. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If an unrecognized coordinate system is encountered, this */
/*         routine signals the error SPICE(NOTSUPPORTED). */

/*     2)  If a map projection not compatible with the input coordinate */
/*         system is encountered, this routine signals the error */
/*         SPICE(NOTSUPPORTED). */

/* $ Files */

/*     See the input HANDLE. */

/* $ Particulars */

/*     This routine displays detailed summary information for a */
/*     specified type 4 DSK segment. The display is written to */
/*     standard output. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     1) The expected range of NSIG is 6:17. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -    DSKBRIEF Version 1.0.0, 06-OCT-2016 (NJB) */

/*        Version 4.0.0 04-OCT-2016 (NJB) */
/*        Version 4.0.0 04-OCT-2013 (NJB) */
/*        Version 3.2.0 14-NOV-2012 (NJB) */
/*        Version 3.1.0 11-NOV-2012 (NJB) */
/*        Version 3.0.0 04-OCT-2012 (NJB) */
/*        Version 2.0.0 22-SEP-2012 (NJB) */
/*        Version 1.0.0 17-AUG-2012 (NJB) */

/* -& */
/* $ Index_Entries */

/*     summarize type 4 dsk segment */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Initial values */

    chkin_("SUM04", (ftnlen)5);

/*     Create d.p. format string. */

    nw = *nsig + 7;
    nf = *nsig - 1;
    s_copy(fmt1, "(1PE@W.@F)", (ftnlen)132, (ftnlen)10);
    repmi_(fmt1, "@W", &nw, fmt1, (ftnlen)132, (ftnlen)2, (ftnlen)132);
    repmi_(fmt1, "@F", &nf, fmt1, (ftnlen)132, (ftnlen)2, (ftnlen)132);

/*     Display type 4 parameters. */

    tostdo_(" ", (ftnlen)1);
    tostdo_("Type 4 parameters", (ftnlen)17);
    tostdo_("-----------------", (ftnlen)17);

/*     Check the DSK format version. */

    dskd04_(handle, dladsc, &c__2, &c__1, &c__1, &n, dpbuff);
    dskfmt = i_dnnt(dpbuff);
    if (dskfmt != 3) {
	setmsg_("DSK format version was expected to be 3 but was #. Only ver"
		"sion 3 is supported by this software version.", (ftnlen)104);
	errint_("#", &dskfmt, (ftnlen)1);
	sigerr_("SPICE(VERSIONMISMATCH)", (ftnlen)22);
	chkout_("SUM04", (ftnlen)5);
	return 0;
    }
    dskb04_(handle, dladsc, &c__7, &c__7, &c__7, &nr, &nc, &co1psz, &co2psz, &
	    cent1, &cent2, &nullok, &nulval, itpdsc, xdsc, nvdsc);

/*     Get and display the numeric data format. */

    s_copy(outlin, "   Height data numeric format:         #", (ftnlen)132, (
	    ftnlen)40);
    dskd04_(handle, dladsc, &c__12, &c__1, &c__1, &n, &dpnfmt);
    numfmt = i_dnnt(&dpnfmt);
    if (numfmt == 2) {
	repmc_(outlin, "#", "16-bit integer", outlin, (ftnlen)132, (ftnlen)1, 
		(ftnlen)14, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else if (numfmt == 1) {
	repmc_(outlin, "#", "32-bit integer", outlin, (ftnlen)132, (ftnlen)1, 
		(ftnlen)14, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else {
	repmc_(outlin, "#", "Not supported", outlin, (ftnlen)132, (ftnlen)1, (
		ftnlen)13, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    }

/*     Get the map projection descriptor; display the projection. */

    dskd04_(handle, dladsc, &c__26, &c__1, &c__7, &n, prjdsc);
    s_copy(outlin, "   Map projection:                     #", (ftnlen)132, (
	    ftnlen)40);
    if (prjdsc[0] == 1.) {
	repmc_(outlin, "#", "Equirectangular", outlin, (ftnlen)132, (ftnlen)1,
		 (ftnlen)15, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else if (prjdsc[0] == 2.) {
	repmc_(outlin, "#", "Stereographic", outlin, (ftnlen)132, (ftnlen)1, (
		ftnlen)13, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else {
	repmc_(outlin, "#", "Not supported", outlin, (ftnlen)132, (ftnlen)1, (
		ftnlen)13, (ftnlen)132);
    }

/*     Get the reference surface; display the surface parameters. */

    dskd04_(handle, dladsc, &c__25, &c__1, &c__7, &n, refdsc);
    s_copy(outlin, "   Reference surface:                  #", (ftnlen)132, (
	    ftnlen)40);
    if (refdsc[0] == 3.) {
	repmc_(outlin, "#", "From coordinate system (above)", outlin, (ftnlen)
		132, (ftnlen)1, (ftnlen)30, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else if (refdsc[0] == 1.) {
	repmc_(outlin, "#", "Sphere", outlin, (ftnlen)132, (ftnlen)1, (ftnlen)
		6, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
	s_copy(outlin, "   Radius (km):                        #", (ftnlen)
		132, (ftnlen)40);
	repmd_(outlin, "#", &refdsc[3], &c__9, outlin, (ftnlen)132, (ftnlen)1,
		 (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else {
	repmc_(outlin, "#", "Not supported", outlin, (ftnlen)132, (ftnlen)1, (
		ftnlen)13, (ftnlen)132);
    }

/*     Get the DSK descriptor. */

    dskgd_(handle, dladsc, dskdsc);

/*     Get the coordinate system. */

    corsys = i_dnnt(&dskdsc[5]);

/*     Show row and column counts. */

    s_copy(outlin, "   Number of pixel grid rows:          #", (ftnlen)132, (
	    ftnlen)40);
    repmi_(outlin, "#", &nr, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132);
    tostdo_(outlin, (ftnlen)132);
    s_copy(outlin, "   Number of pixel grid columns:       #", (ftnlen)132, (
	    ftnlen)40);
    repmi_(outlin, "#", &nc, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132);
    tostdo_(outlin, (ftnlen)132);

/*     Show nested grid dimensions if there are at least two levels. */


/*     Fetch the number of nested grid dimensions and the */
/*     dimensions themselves. */

    ibase = dladsc[2];
    b = ibase + 1;
    dasrdi_(handle, &b, &b, &ngdim);
    s_copy(outlin, "   Integer grid nesting levels:        #", (ftnlen)132, (
	    ftnlen)40);
    repmi_(outlin, "#", &ngdim, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132);
    tostdo_(outlin, (ftnlen)132);
    if (ngdim >= 2) {
	b = ibase + 2;
	e = b + (ngdim << 1) - 1;
	dasrdi_(handle, &b, &e, grdims);
	i__1 = ngdim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(outlin, "      Grid dimensions at level #:      # rows x "
		    "# columns", (ftnlen)132, (ftnlen)57);
	    repmi_(outlin, "#", &i__, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)
		    132);
	    repmi_(outlin, "#", &grdims[(i__2 = (i__ << 1) - 2) < 20 && 0 <= 
		    i__2 ? i__2 : s_rnge("grdims", i__2, "sum04_", (ftnlen)
		    372)], outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132);
	    repmi_(outlin, "#", &grdims[(i__2 = (i__ << 1) - 1) < 20 && 0 <= 
		    i__2 ? i__2 : s_rnge("grdims", i__2, "sum04_", (ftnlen)
		    373)], outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	}
    }

/*     Show pixel dimensions and grid center coordinates. */

/*     Note: the code below is applicable only to the */
/*     equirectangular projection. */

    if (prjdsc[0] == 1.) {
	if (corsys == 1) {

/*           Show coordinate 1 pixel dimension. */

	    s_copy(outlin, "   Longitude pixel dimension (deg):    #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = co1psz * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "E", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);

/*           Show coordinate 2 pixel dimension. */

	    s_copy(outlin, "   Latitude pixel dimension (deg):     #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = co2psz * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "E", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);

/*           Show grid center coordinates. */

	    s_copy(outlin, "   Pixel grid center longitude (deg):  #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = cent1 * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "F", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	    s_copy(outlin, "   Pixel grid center latitude (deg):   #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = cent2 * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "F", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	} else if (corsys == 4) {

/*           Show coordinate 1 pixel dimension. */

	    s_copy(outlin, "   Longitude pixel dimension   (deg):  #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = co1psz * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "E", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);

/*           Show coordinate 2 pixel dimension. */

	    s_copy(outlin, "   Latitude pixel dimension    (deg):  #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = co2psz * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "E", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);

/*           Show grid center coordinates. */

	    s_copy(outlin, "   Pixel grid center longitude (deg):  #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = cent1 * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "F", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	    s_copy(outlin, "   Pixel grid center latitude  (deg):  #", (
		    ftnlen)132, (ftnlen)40);
	    d__1 = cent2 * dpr_();
	    repmf_(outlin, "#", &d__1, nsig, "F", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	} else {
	    setmsg_("This coordinate system is not supported for the equirec"
		    "tangular projection. Coordinate system code: #", (ftnlen)
		    101);
	    errint_("#", &corsys, (ftnlen)1);
	    sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	}
    } else if (prjdsc[0] == 2.) {
	if (corsys == 3) {

/*           Show coordinate 1 pixel dimension. */

	    s_copy(outlin, "   X pixel dimension (km):             #", (
		    ftnlen)132, (ftnlen)40);
	    repmf_(outlin, "#", &co1psz, nsig, "E", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);

/*           Show coordinate 2 pixel dimension. */

	    s_copy(outlin, "   Y pixel dimension (km):             #", (
		    ftnlen)132, (ftnlen)40);
	    repmf_(outlin, "#", &co2psz, nsig, "E", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);

/*           Show grid center coordinates. */

	    s_copy(outlin, "   Pixel grid center X (km):           #", (
		    ftnlen)132, (ftnlen)40);
	    repmf_(outlin, "#", &cent1, nsig, "F", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	    s_copy(outlin, "   Pixel grid center Y (km):           #", (
		    ftnlen)132, (ftnlen)40);
	    repmf_(outlin, "#", &cent2, nsig, "F", outlin, (ftnlen)132, (
		    ftnlen)1, (ftnlen)1, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	} else {
	    setmsg_("This coordinate system is not supported for the project"
		    "ion having code #. Coordinate system code: #", (ftnlen)99)
		    ;
	    i__1 = i_dnnt(prjdsc);
	    errint_("#", &i__1, (ftnlen)1);
	    errint_("#", &corsys, (ftnlen)1);
	    sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	    chkout_("SUM04", (ftnlen)5);
	    return 0;
	}
    } else {
	setmsg_("Unrecognized coordinate system code: #", (ftnlen)38);
	errint_("#", &corsys, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("SUM04", (ftnlen)5);
	return 0;
    }
    s_copy(outlin, "   Null values allowed:                #", (ftnlen)132, (
	    ftnlen)40);
    if (nullok) {
	repmc_(outlin, "#", "Yes", outlin, (ftnlen)132, (ftnlen)1, (ftnlen)3, 
		(ftnlen)132);
    } else {
	repmc_(outlin, "#", "No", outlin, (ftnlen)132, (ftnlen)1, (ftnlen)2, (
		ftnlen)132);
    }
    tostdo_(outlin, (ftnlen)132);
    if (nullok) {
	s_copy(outlin, "   Null value parameter:               #", (ftnlen)
		132, (ftnlen)40);
	i__1 = i_dnnt(&nulval);
	repmi_(outlin, "#", &i__1, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)
		132);
	tostdo_(outlin, (ftnlen)132);
    }

/*     Show the height scale. */

    dskd04_(handle, dladsc, &c__11, &c__1, &c__1, &n, &hscale);
    s_copy(outlin, "   Height units in km:                 #", (ftnlen)132, (
	    ftnlen)40);
    repmd_(outlin, "#", &hscale, nsig, outlin, (ftnlen)132, (ftnlen)1, (
	    ftnlen)132);
    tostdo_(outlin, (ftnlen)132);
    s_copy(outlin, "   Interpolation method:               #", (ftnlen)132, (
	    ftnlen)40);
    if (itpdsc[0] == 1.) {
	s_copy(itpstr, "None: return raw height data", (ftnlen)132, (ftnlen)
		28);
    } else if (itpdsc[0] == 2.) {
	s_copy(itpstr, "Bilinear", (ftnlen)132, (ftnlen)8);
    } else if (itpdsc[0] == 3.) {
	s_copy(itpstr, "3 nearest neighbor linear", (ftnlen)132, (ftnlen)25);
    } else {
	setmsg_("Bad interpolation code: #", (ftnlen)25);
	i__1 = i_dnnt(itpdsc);
	errint_("#", &i__1, (ftnlen)1);
	sigerr_("BUG", (ftnlen)3);
	chkout_("SUM04", (ftnlen)5);
	return 0;
    }
    repmc_(outlin, "#", itpstr, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132, (
	    ftnlen)132);
    tostdo_(outlin, (ftnlen)132);
    s_copy(outlin, "   Surface intercept method:           #", (ftnlen)132, (
	    ftnlen)40);
    if (xdsc[0] == 1.) {
	s_copy(xstr, "Constant angular step", (ftnlen)132, (ftnlen)21);
    } else if (xdsc[0] == 2.) {
	s_copy(xstr, "Constant projected step", (ftnlen)132, (ftnlen)23);
    } else if (xdsc[0] == 3.) {
	s_copy(xstr, "Cell boundary step", (ftnlen)132, (ftnlen)18);
    } else {
	setmsg_("Bad intercept code: #", (ftnlen)21);
	i__1 = i_dnnt(xdsc);
	errint_("#", &i__1, (ftnlen)1);
	sigerr_("BUG", (ftnlen)3);
	chkout_("SUM04", (ftnlen)5);
	return 0;
    }
    repmc_(outlin, "#", xstr, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132, (
	    ftnlen)132);
    tostdo_(outlin, (ftnlen)132);
    if (xdsc[0] == 1.) {
	s_copy(outlin, "      Step size                (deg):  #", (ftnlen)
		132, (ftnlen)40);
	d__1 = xdsc[2] * dpr_();
	repmf_(outlin, "#", &d__1, nsig, "E", outlin, (ftnlen)132, (ftnlen)1, 
		(ftnlen)1, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
	s_copy(outlin, "      Convergence tolerance      (m):  #", (ftnlen)
		132, (ftnlen)40);
	d__1 = xdsc[3] * 1e3;
	repmf_(outlin, "#", &d__1, nsig, "E", outlin, (ftnlen)132, (ftnlen)1, 
		(ftnlen)1, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else if (xdsc[0] == 2.) {
	s_copy(outlin, "      Step size(km):                #", (ftnlen)132, (
		ftnlen)37);
	d__1 = xdsc[2] * dpr_();
	repmf_(outlin, "#", &d__1, nsig, "E", outlin, (ftnlen)132, (ftnlen)1, 
		(ftnlen)1, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
	s_copy(outlin, "      Convergence tolerance (km):   #", (ftnlen)132, (
		ftnlen)37);
	repmf_(outlin, "#", &xdsc[3], nsig, "E", outlin, (ftnlen)132, (ftnlen)
		1, (ftnlen)1, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else {
	setmsg_("Bad intercept code: #", (ftnlen)21);
	i__1 = i_dnnt(xdsc);
	errint_("#", &i__1, (ftnlen)1);
	sigerr_("BUG", (ftnlen)3);
	chkout_("SUM04", (ftnlen)5);
	return 0;
    }

/*     Display surface intercept acceleration information: */

    dskd04_(handle, dladsc, &c__27, &c__1, &c__7, &n, accdsc);
    s_copy(outlin, "      Acceleration algorithm:          #", (ftnlen)132, (
	    ftnlen)40);
    if (accdsc[0] == 1.) {
	s_copy(accstr, "Not enabled", (ftnlen)132, (ftnlen)11);
	repmc_(outlin, "#", accstr, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)
		132, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else if (accdsc[0] == 2.) {
	s_copy(accstr, "Large initial step", (ftnlen)132, (ftnlen)18);
	repmc_(outlin, "#", accstr, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)
		132, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
	s_copy(outlin, "         Maximum magnitude of height", (ftnlen)132, (
		ftnlen)36);
	tostdo_(outlin, (ftnlen)132);
	s_copy(outlin, "         deltas between adjacent pixels:", (ftnlen)
		132, (ftnlen)40);
	tostdo_(outlin, (ftnlen)132);
	s_copy(outlin, "            In COORD1 direction (km):  #", (ftnlen)
		132, (ftnlen)40);
	if (corsys == 4) {
	    repmc_(outlin, "COORD1", "longitude", outlin, (ftnlen)132, (
		    ftnlen)6, (ftnlen)9, (ftnlen)132);
	}
	d__1 = accdsc[2] * hscale;
	repmd_(outlin, "#", &d__1, nsig, outlin, (ftnlen)132, (ftnlen)1, (
		ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
	s_copy(outlin, "            In COORD2 direction  (km):  #", (ftnlen)
		132, (ftnlen)41);
	if (corsys == 4) {
	    repmc_(outlin, "COORD2", "latitude", outlin, (ftnlen)132, (ftnlen)
		    6, (ftnlen)8, (ftnlen)132);
	}
	d__1 = accdsc[3] * hscale;
	repmd_(outlin, "#", &d__1, nsig, outlin, (ftnlen)132, (ftnlen)1, (
		ftnlen)132);
	tostdo_(outlin, (ftnlen)132);
    } else {
	setmsg_("Bad accleration code: #", (ftnlen)23);
	i__1 = i_dnnt(accdsc);
	errint_("#", &i__1, (ftnlen)1);
	sigerr_("BUG", (ftnlen)3);
	chkout_("SUM04", (ftnlen)5);
	return 0;
    }
    s_copy(outlin, "   Normal vector calculation method:   #", (ftnlen)132, (
	    ftnlen)40);
    if (nvdsc[0] == 2.) {
	s_copy(xstr, "Defined by partial derivatives", (ftnlen)132, (ftnlen)
		30);
    } else if (nvdsc[0] == 1.) {
	s_copy(xstr, "Bilinear compatible", (ftnlen)132, (ftnlen)19);
    } else {
	setmsg_("Bad normal vector code: #", (ftnlen)25);
	i__1 = i_dnnt(nvdsc);
	errint_("#", &i__1, (ftnlen)1);
	sigerr_("BUG", (ftnlen)3);
	chkout_("SUM04", (ftnlen)5);
	return 0;
    }
    repmc_(outlin, "#", xstr, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)132, (
	    ftnlen)132);
    tostdo_(outlin, (ftnlen)132);
    chkout_("SUM04", (ftnlen)5);
    return 0;
} /* sum04_ */

