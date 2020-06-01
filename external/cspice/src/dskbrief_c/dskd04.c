/* dskd04.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure DSKD04 ( DSK, fetch d.p. type 4 data ) */
/* Subroutine */ int dskd04_(integer *handle, integer *dladsc, integer *item, 
	integer *start, integer *room, integer *n, doublereal *values)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer i_dnnt(doublereal *), s_rnge(char *, integer, char *, integer);

    /* Local variables */
    integer addr__, size, b, e;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer dsize;
    extern /* Subroutine */ int moved_(doublereal *, integer *, doublereal *);
    integer dpbase;
    doublereal dpbuff[100];
    integer maxadd, dpsize;
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), dasrdd_(integer *, integer *, integer *, doublereal *);

/* $ Abstract */

/*     Fetch double precision data from a type 4 DSK segment. */

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

/*     DAS */
/*     DSK */

/* $ Keywords */

/*     DAS */
/*     DSK */
/*     FILES */

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
/*     HANDLE     I   DSK file handle. */
/*     DLADSC     I   DLA descriptor. */
/*     ITEM       I   Keyword identifying item to fetch. */
/*     START      I   Start index. */
/*     ROOM       I   Amount of room in output array. */
/*     N          O   Number of values returned. */
/*     VALUES     O   Array containing requested item. */

/* $ Detailed_Input */

/*     HANDLE         is the handle of a DSK file containing a type 4 */
/*                    segment from which data are to be fetched. */

/*     DLADSC         is the DLA descriptor associated with the segment */
/*                    from which data are to be fetched. */

/*     ITEM           is an integer "keyword" parameter designating the */
/*                    item to fetch. In the descriptions below, note */
/*                    that "model" refers to the model represented by */
/*                    the designated segment.  This model may be a */
/*                    subset of a larger model. */

/*                    Names and meanings of parameters supported by this */
/*                    routine are: */

/*                       KWDSC      DSK descriptor of segment. See the */
/*                                  INCLUDE file dskdsc.inc for a */
/*                                  discussion of the contents of DSK */
/*                                  descriptors. Note that DSK */
/*                                  descriptors are not to be confused */
/*                                  with DLA descriptors, which contain */
/*                                  segment component base address and */
/*                                  size information. */


/*     START          is the start index within specified data item from */
/*                    which data are to be fetched. The index of the */
/*                    first element of each data item is 1. */

/*     ROOM           is the amount of room in the output array.  It */
/*                    is permissible to provided an output array */
/*                    that has too little room to fetch an item in */
/*                    one call. */

/* $ Detailed_Output */

/*     N              is the number of elements fetched to the output */
/*                    array VALUES.  N is normally in the range */
/*                    1:ROOM; if an error occurs on the call, N is */
/*                    undefined. */

/*     VALUES         is a contiguous set of elements of the item */
/*                    designated by ITEM.  The correspondence of */
/*                    VALUES at the elements of the data item is: */

/*                       VALUES(1)      ITEM(START) */
/*                         ...             ... */
/*                       VALUES(N)      ITEM(START+N-1) */

/*                    If an error occurs on the call, VALUES is */
/*                    undefined. */

/* $ Parameters */

/*     See the INCLUDE files */

/*         dla.inc */
/*         dsk04.inc */
/*         dskdsc.inc */

/* $ Exceptions */

/*     1) If the input handle is invalid, the error will be diagnosed by */
/*        routines in the call tree of this routine. */

/*     2) If a file read error occurs, the error will be diagnosed by */
/*        routines in the call tree of this routine. */

/*     3) If the input DLA descriptor is invalid, the effect of this */
/*        routine is undefined. The error *may* be diagnosed by routines */
/*        in the call  tree of this routine, but there are no */
/*        guarantees. */

/*     4) If ROOM is non-positive, the error SPICE(VALUEOUTOFRANGE) */
/*        is signaled. */

/*     5) If the input keyword parameter is not recognized, the error */
/*        SPICE(NOTSUPPORTED) is signaled. */

/*     6) If START is less than 1 or greater than the size of the */
/*        item to be fetched, the error SPICE(INDEXOUTOFRANGE) is */
/*        signaled. */

/* $ Files */

/*     See input argument HANDLE. */

/* $ Particulars */

/*     DSK files are built using the DLA low-level format and */
/*     the DAS architecture; DLA files are a specialized type of DAS */
/*     file in which data are organized as a doubly linked list of */
/*     segments.  Each segment's data belong to contiguous components of */
/*     character, double precision, and integer type. Depending on */
/*     the segment type, some of these components may be empty. */

/*     Note that the DSK descriptor for the segment is not needed by */
/*     this routine; the DLA descriptor contains the base address and */
/*     size information for the integer, double precision, and character */
/*     components of the segment, and these suffice for the purpose of */
/*     fetching data. */

/* $ Examples */

/*     1) */



/* $ Restrictions */

/*     1) This is a prototype routine.  The interface is not */
/*        expected to change, but there are no guarantees. */

/*     2) This routine uses discovery check-in to boost */
/*        execution speed.  However, this routine is in */
/*        violation of NAIF standards for use of discovery */
/*        check-in:  routines called from this routine may */
/*        signal errors.  If errors are signaled in called */
/*        routines, this routine's name will be missing */
/*        from the traceback message. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.0, 19-SEP-2012 (NJB) */

/*        Updated to handle new type 4 format. New items */
/*        are: number format, reference surface descriptor, */
/*        projection descriptor. */

/* -    SPICELIB Version 1.0.0, 17-AUG-2012 (NJB) */

/* -& */
/* $ Index_Entries */

/*     fetch d.p. data from a type_4_dsk_segment */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Use discovery check-in.  This is done for efficiency; note */
/*     however that this routine does not meet SPICE standards for */
/*     discovery check-in eligibility. */

/*     Type 1 DSK segments contain only d.p. information, */
/*     so there's only one DLA component to consider. */

    dpbase = dladsc[4];
    dpsize = dladsc[5];
    maxadd = dpbase + dpsize;
/*      WRITE (*,*) 'DP BASE = ', DPBASE */
/*      WRITE (*,*) 'DP SIZE = ', DPSIZE */
/*      WRITE (*,*) 'ITEM    = ', ITEM */
/*      WRITE (*,*) ' ' */
    if (*room <= 0) {
	chkin_("DSKD04", (ftnlen)6);
	setmsg_("ROOM was #; must be positive.", (ftnlen)29);
	errint_("#", room, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("DSKD04", (ftnlen)6);
	return 0;
    }

/*     Perform the best check on START that we can without */
/*     having determined the size of the item to fetch. */

    if (*start < 1 || *start > dpsize) {
	chkin_("DSKD04", (ftnlen)6);
	setmsg_("START must be in the range defined by the size of the data "
		"associated with the keyword parameter #, namely 1:#.  Actual"
		" value of START was #.", (ftnlen)141);
	errint_("#", item, (ftnlen)1);
	errint_("#", &dpsize, (ftnlen)1);
	errint_("#", start, (ftnlen)1);
	sigerr_("SPICE(INDEXOUTOFRANGE)", (ftnlen)22);
	chkout_("DSKD04", (ftnlen)6);
	return 0;
    }

/*     Branch based on the item to be returned. In each case, set */
/*     the start address of the item and the size of the item. */

    if (*item == 1) {

/*        Return the DSK descriptor. */
/*        First locate the data. */

	b = dpbase + 1;
	size = 24;
    } else if (*item == 2) {

/*        Return the format version. */

	b = dpbase + 25;
	size = 1;
    } else if (*item == 3) {

/*        Return the row count. */

	b = dpbase + 26;
	size = 1;
    } else if (*item == 4) {

/*        Return the column count. */

	b = dpbase + 27;
	size = 1;
    } else if (*item == 7) {

/*        Return the first coordinate of the pixel grid center. */

	b = dpbase + 30;
	size = 1;
    } else if (*item == 8) {

/*        Return the second coordinate of the pixel grid center. */

	b = dpbase + 31;
	size = 1;
    } else if (*item == 5) {

/*        Return the pixel size in the direction of coordinate 1. */

	b = dpbase + 28;
	size = 1;
    } else if (*item == 6) {

/*        Return the pixel size in the direction of coordinate 2. */

	b = dpbase + 29;
	size = 1;
    } else if (*item == 9) {

/*        Return the flag indicating whether nulls are allowed. */
/*        This is a d.p. parameter that maps to a Boolean value. */

	b = dpbase + 32;
	size = 1;
    } else if (*item == 10) {

/*        Return the null value parameter. */

	b = dpbase + 33;
	size = 1;
    } else if (*item == 11) {

/*        Return the height scale parameter. */

	b = dpbase + 34;
	size = 1;
    } else if (*item == 12) {

/*        Return the number format parameter. */

	b = dpbase + 35;
	size = 1;
    } else if (*item == 13) {

/*        Return the pointer to the reference surface */
/*        descriptor. */

	b = dpbase + 40;
	size = 1;
    } else if (*item == 14) {

/*        Return the pointer to the map projection */
/*        descriptor. */

	b = dpbase + 41;
	size = 1;
    } else if (*item == 15) {

/*        Return the pointer to the interpolation algorithm */
/*        descriptor. */

	b = dpbase + 36;
	size = 1;
    } else if (*item == 16) {

/*        Return the pointer to the intercept algorithm */
/*        descriptor. */

	b = dpbase + 37;
	size = 1;
    } else if (*item == 17) {

/*        Return the pointer to the normal vector algorithm */
/*        descriptor. */

	b = dpbase + 38;
	size = 1;
    } else if (*item == 18) {

/*        Return the pointer to the intercept acceleration */
/*        parameters. */

	b = dpbase + 39;
	size = 1;
    } else if (*item == 19) {

/*        Return the pointer to the coarse grid structure. */

	b = dpbase + 42;
	size = 1;
    } else if (*item == 20) {

/*        Return the pointer to the pixel data. */

	b = dpbase + 43;
	size = 1;
    } else if (*item == 25) {

/*        Return the reference surface descriptor. */

/*        We'll need the pointer to this descriptor first. */

	addr__ = dpbase + 40;
	dasrdd_(handle, &addr__, &addr__, dpbuff);

/*        We won't know the descriptor's size until we */
/*        actually read the descriptor, so prepare to */
/*        read the descriptor of maximum size. We'll */
/*        adjust the size before returning. */

	b = dpbase + i_dnnt(dpbuff);
	size = 7;
/* Computing MIN */
	i__1 = b + size - 1;
	e = min(i__1,maxadd);
	dasrdd_(handle, &b, &e, dpbuff);

/*        Extract the actual descriptor size from the descriptor. */

/* Computing MIN */
	i__1 = *room, i__2 = i_dnnt(&dpbuff[1]);
	dsize = min(i__1,i__2);

/*        Move as much of the descriptor as is requested and */
/*        will fit to the output array. */

/* Computing MIN */
	i__1 = *room, i__2 = dsize - *start + 1;
	*n = min(i__1,i__2);
	moved_(&dpbuff[(i__1 = *start - 1) < 100 && 0 <= i__1 ? i__1 : s_rnge(
		"dpbuff", i__1, "dskd04_", (ftnlen)500)], n, values);

/*        The outputs are set. */

	return 0;
    } else if (*item == 26) {

/*        Return the map projection descriptor. */

/*        We'll need the pointer to this descriptor first. */

	addr__ = dpbase + 41;
	dasrdd_(handle, &addr__, &addr__, dpbuff);

/*        We won't know the descriptor's size until we */
/*        actually read the descriptor, so prepare to */
/*        read the descriptor of maximum size. We'll */
/*        adjust the size before returning. */

	b = dpbase + i_dnnt(dpbuff);
	size = 7;
/* Computing MIN */
	i__1 = b + size - 1;
	e = min(i__1,maxadd);
	dasrdd_(handle, &b, &e, dpbuff);

/*        Extract the actual descriptor size from the descriptor. */

/* Computing MIN */
	i__1 = *room, i__2 = i_dnnt(&dpbuff[1]);
	dsize = min(i__1,i__2);

/*        Move as much of the descriptor as is requested and */
/*        will fit to the output array. */

/* Computing MIN */
	i__1 = *room, i__2 = dsize - *start + 1;
	*n = min(i__1,i__2);
	moved_(&dpbuff[(i__1 = *start - 1) < 100 && 0 <= i__1 ? i__1 : s_rnge(
		"dpbuff", i__1, "dskd04_", (ftnlen)543)], n, values);

/*        The outputs are set. */

	return 0;
    } else if (*item == 21) {

/*        Return the interpolation algorithm descriptor. */

/*        We'll need the pointer to this descriptor first. */

	addr__ = dpbase + 36;
	dasrdd_(handle, &addr__, &addr__, dpbuff);

/*        We won't know the descriptor's size until we */
/*        actually read the descriptor, so prepare to */
/*        read the descriptor of maximum size. We'll */
/*        adjust the size before returning. */

	b = dpbase + i_dnnt(dpbuff);
	size = 7;
/* Computing MIN */
	i__1 = b + size - 1;
	e = min(i__1,maxadd);
	dasrdd_(handle, &b, &e, dpbuff);

/*        Extract the actual descriptor size from the descriptor. */

/* Computing MIN */
	i__1 = *room, i__2 = i_dnnt(&dpbuff[1]);
	dsize = min(i__1,i__2);

/*        Move as much of the descriptor as is requested and */
/*        will fit to the output array. */

/* Computing MIN */
	i__1 = *room, i__2 = dsize - *start + 1;
	*n = min(i__1,i__2);
	moved_(&dpbuff[(i__1 = *start - 1) < 100 && 0 <= i__1 ? i__1 : s_rnge(
		"dpbuff", i__1, "dskd04_", (ftnlen)586)], n, values);

/*        The outputs are set. */

	return 0;
    } else if (*item == 22) {

/*        Return the intercept algorithm descriptor. */

/*        We'll need the pointer to this descriptor */
/*        first. */

	addr__ = dpbase + 37;
	dasrdd_(handle, &addr__, &addr__, dpbuff);

/*        We won't know the descriptor's size until we */
/*        actually read the descriptor, so prepare to */
/*        read the descriptor of maximum size. We'll */
/*        adjust the size before returning. */

	b = dpbase + i_dnnt(dpbuff);
	size = 7;
/* Computing MIN */
	i__1 = b + size - 1;
	e = min(i__1,maxadd);
	dasrdd_(handle, &b, &e, dpbuff);

/*        Extract the actual descriptor size from the descriptor. */

/* Computing MIN */
	i__1 = *room, i__2 = i_dnnt(&dpbuff[1]);
	dsize = min(i__1,i__2);

/*        Move as much of the descriptor as is requested and */
/*        will fit to the output array. */

/* Computing MIN */
	i__1 = *room, i__2 = dsize - *start + 1;
	*n = min(i__1,i__2);
	moved_(&dpbuff[(i__1 = *start - 1) < 100 && 0 <= i__1 ? i__1 : s_rnge(
		"dpbuff", i__1, "dskd04_", (ftnlen)629)], n, values);

/*        The outputs are set. */

	return 0;
    } else if (*item == 23) {

/*        Return the normal vector algorithm descriptor. */

/*        We'll need the pointer to this descriptor first. */

	addr__ = dpbase + 38;
	dasrdd_(handle, &addr__, &addr__, dpbuff);

/*        We won't know the descriptor's size until we */
/*        actually read the descriptor, so prepare to */
/*        read the descriptor of maximum size. We'll */
/*        adjust the size before returning. */

	b = dpbase + i_dnnt(dpbuff);
	size = 7;
/* Computing MIN */
	i__1 = b + size - 1;
	e = min(i__1,maxadd);
	dasrdd_(handle, &b, &e, dpbuff);

/*        Extract the actual descriptor size from the descriptor. */

/* Computing MIN */
	i__1 = *room, i__2 = i_dnnt(&dpbuff[1]);
	dsize = min(i__1,i__2);

/*        Move as much of the descriptor as is requested and */
/*        will fit to the output array. */

/* Computing MIN */
	i__1 = *room, i__2 = dsize - *start + 1;
	*n = min(i__1,i__2);
	moved_(&dpbuff[(i__1 = *start - 1) < 100 && 0 <= i__1 ? i__1 : s_rnge(
		"dpbuff", i__1, "dskd04_", (ftnlen)671)], n, values);

/*        The outputs are set. */

	return 0;
    } else if (*item == 27) {

/*        Return the intercept acceleration descriptor. */

/*        We'll need the pointer to this descriptor first. */

	addr__ = dpbase + 39;
	dasrdd_(handle, &addr__, &addr__, dpbuff);

/*        We won't know the descriptor's size until we */
/*        actually read the descriptor, so prepare to */
/*        read the descriptor of maximum size. We'll */
/*        adjust the size before returning. */

	b = dpbase + i_dnnt(dpbuff);
	size = 7;
/* Computing MIN */
	i__1 = b + size - 1;
	e = min(i__1,maxadd);
	dasrdd_(handle, &b, &e, dpbuff);

/*        Extract the actual descriptor size from the descriptor. */

/* Computing MIN */
	i__1 = *room, i__2 = i_dnnt(&dpbuff[1]);
	dsize = min(i__1,i__2);

/*        Move as much of the descriptor as is requested and */
/*        will fit to the output array. */

/* Computing MIN */
	i__1 = *room, i__2 = dsize - *start + 1;
	*n = min(i__1,i__2);
	moved_(&dpbuff[(i__1 = *start - 1) < 100 && 0 <= i__1 ? i__1 : s_rnge(
		"dpbuff", i__1, "dskd04_", (ftnlen)713)], n, values);

/*        The outputs are set. */

	return 0;
    } else if (*item == 24) {
	chkin_("DSKD04", (ftnlen)6);
	setmsg_("Keyword parameter #, which denotes the coarse grid data str"
		"ucture, is not yet supported.", (ftnlen)88);
	errint_("#", item, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("DSKD04", (ftnlen)6);
	return 0;
    } else {
	chkin_("DSKD04", (ftnlen)6);
	setmsg_("Keyword parameter # was not recognized.", (ftnlen)39);
	errint_("#", item, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("DSKD04", (ftnlen)6);
	return 0;
    }

/*     The valid range for START is 1:SIZE. */

    if (*start < 1 || *start > size) {
	chkin_("DSKD04", (ftnlen)6);
	setmsg_("START must be in the range defined by the size of the data "
		"associated with the keyword parameter #, namely 1:#.  Actual"
		" value of START was #.", (ftnlen)141);
	errint_("#", item, (ftnlen)1);
	errint_("#", &size, (ftnlen)1);
	errint_("#", start, (ftnlen)1);
	sigerr_("SPICE(INDEXOUTOFRANGE)", (ftnlen)22);
	chkout_("DSKD04", (ftnlen)6);
	return 0;
    }

/*     Read the requested data.  We already have the start address B. */

/* Computing MIN */
    i__1 = *room, i__2 = size - *start + 1;
    *n = min(i__1,i__2);
    e = b + *n - 1;
    dasrdd_(handle, &b, &e, values);
    return 0;
} /* dskd04_ */

