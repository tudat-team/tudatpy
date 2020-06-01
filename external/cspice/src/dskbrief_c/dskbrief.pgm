/* dskbrief.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__25000 = 25000;
static integer c__1 = 1;
static integer c__10 = 10;
static integer c__3 = 3;
static integer c__17 = 17;
static integer c__2 = 2;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__6 = 6;
static integer c__19 = 19;
static integer c__21 = 21;
static integer c__23 = 23;
static integer c__8 = 8;
static integer c_b91 = 100000;
static integer c__10000 = 10000;
static integer c__24 = 24;

/* $Program    DSKBRIEF ( BRIEF DSK summary ) */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static char optlst[32*10] = "-h                              " "-v      "
	    "                        " "-ext                            " 
	    "-full                           " "-a                          "
	    "    " "-seg                            " "-gaps                 "
	    "          " "-d                              " "-u              "
	    "                " "-tg                             ";
    static doublereal timtol = 0.;
    static logical all = FALSE_;
    static logical ext = FALSE_;
    static logical full = FALSE_;
    static logical gap = FALSE_;
    static logical seg = FALSE_;
    static logical ust = FALSE_;

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__3, i__4, i__5;
    doublereal d__1, d__2;
    char ch__1[268];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer i_dnnt(doublereal *), s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer node, nseg, ndsk;
    static doublereal xbds[200000]	/* was [2][100000] */, ybds[200000]	
	    /* was [2][100000] */;
    static integer nsig;
    extern /* Subroutine */ int sum02_(integer *, integer *, integer *);
    static integer srcs[4];
    extern /* Subroutine */ int sum04_(integer *, integer *, integer *);
    static integer nout;
    extern /* Subroutine */ int zzdbrgap_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *);
    static integer b, e, i__, j, k;
    extern /* Subroutine */ int kdata_(integer *, char *, char *, char *, 
	    char *, integer *, logical *, ftnlen, ftnlen, ftnlen, ftnlen);
    static char fname[255];
    static integer dlads[800000]	/* was [8][100000] */;
    extern /* Subroutine */ int chkin_(char *, ftnlen), dskgd_(integer *, 
	    integer *, doublereal *);
    extern logical beint_(char *, ftnlen);
    extern /* Subroutine */ int errch_(char *, char *, ftnlen, ftnlen), 
	    repmc_(char *, char *, char *, char *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static doublereal gxbds[200000]	/* was [2][100000] */, gybds[200000]	
	    /* was [2][100000] */;
    extern /* Subroutine */ int moved_(doublereal *, integer *, doublereal *);
    static integer segno, ncomp;
    static logical found;
    static integer dskno;
    static logical segtm[100000];
    extern /* Subroutine */ int repmi_(char *, char *, integer *, char *, 
	    ftnlen, ftnlen, ftnlen), movei_(integer *, integer *, integer *);
    static integer b1, e1;
    static doublereal minco1, minco2, maxco1, maxco2, maxco3, minco3;
    static integer ng, dladsc[8], handle;
    static logical didfil;
    extern /* Subroutine */ int dlabfs_(integer *, integer *, logical *);
    static logical diddsk;
    static char banner[80];
    extern /* Subroutine */ int dlafns_(integer *, integer *, integer *, 
	    logical *);
    static doublereal comdsc[24];
    extern integer brckti_(integer *, integer *, integer *), isrchc_(char *, 
	    integer *, char *, ftnlen, ftnlen), lnknxt_(integer *, integer *);
    extern logical exists_(char *, ftnlen);
    extern /* Subroutine */ int attcmp_();
    static char errmsg[255], filtyp[10], numstr[32], option[32], outlin[255], 
	    source[255];
    static doublereal dskdsc[24], grpdsc[24], lonbds[4]	/* was [2][2] */, 
	    maxtdb, maxlon, minlon, mintdb, outbds[8]	/* was [2][4] */;
    static char cmd[25500];
    static integer corsys, currnt[8], grpsiz, hanlst[100000], nivals, nxtdsc[
	    8], seglst[100000], sgpool[200012]	/* was [2][100006] */, sgptrs[
	    10000], typcde;
    static logical opterr;
    extern /* Subroutine */ int prcinf_(char *, ftnlen), getcml_(char *, 
	    ftnlen), byebye_(char *, ftnlen), setmsg_(char *, ftnlen), 
	    errint_(char *, integer *, ftnlen), sigerr_(char *, ftnlen), 
	    fndnwd_(char *, integer *, integer *, integer *, ftnlen), prsint_(
	    char *, integer *, ftnlen), tostdo_(char *, ftnlen), furnsh_(char 
	    *, ftnlen), ktotal_(char *, integer *, ftnlen), dasopr_(char *, 
	    integer *, ftnlen), dspdsc_(doublereal *, integer *, integer *, 
	    integer *), grpseg_(U_fp, logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, logical *), reglon_(integer *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *), dspgap_(logical *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), chkout_(char *, ftnlen);

/* $ Abstract */

/*     DSKBRIEF is a utility program that provides brief summaries of */
/*     the contents of one or more DSK files. */

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

/*     DSKBRIEF User's Guide. */

/* $ Keywords */

/*     FILES */
/*     UTILITY */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     None. */

/* $ Files */

/*     1)  All summary tasks performed require a name for at least one */
/*         DSK file to be provided. */

/*     2)  To display names of surfaces, and appropriate text kernel */
/*         files containing surface name-ID mappings must be provided. */

/* $ Particulars */

/*     For usage details see the DSKBRIEF User's Guide. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -       DSKBRIEF Version 6.0.0 21-FEB-2017 (NJB) */

/*           Gap summary capability has been added. */

/*           Now, by default, summarizes groups of segments */
/*           having matching attributes, when those segments */
/*           belong to the same DSK file. Optionally summarizes */
/*           groups of segments from multiple DSK files. */

/*           Summary of type 4 segments is now supported. */

/* -       DSKBRIEF Version 5.0.0 05-MAY-2010 (NJB) */
/* -       DSKBRIEF Version 4.0.0 14-SEP-2008 (NJB) */
/* -       DSKBRIEF Version 3.0.0 28-DEC-2006 (NJB) */
/* -       DSKBRIEF Version 2.0.0 04-AUG-2006 (NJB) */
/* -       DSKBRIEF Version 1.0.0 10-JUL-2006 (NJB) */

/* -& */

/*     Global parameters */


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


/*     File: dsktol.inc */


/*     This file contains declarations of tolerance and margin values */
/*     used by the DSK subsystem. */

/*     It is recommended that the default values defined in this file be */
/*     changed only by expert SPICE users. */

/*     The values declared in this file are accessible at run time */
/*     through the routines */

/*        DSKGTL  {DSK, get tolerance value} */
/*        DSKSTL  {DSK, set tolerance value} */

/*     These are entry points of the routine DSKTOL. */

/*        Version 1.0.0 27-FEB-2016 (NJB) */




/*     Parameter declarations */
/*     ====================== */

/*     DSK type 2 plate expansion factor */
/*     --------------------------------- */

/*     The factor XFRACT is used to slightly expand plates read from DSK */
/*     type 2 segments in order to perform ray-plate intercept */
/*     computations. */

/*     This expansion is performed to prevent rays from passing through */
/*     a target object without any intersection being detected. Such */
/*     "false miss" conditions can occur due to round-off errors. */

/*     Plate expansion is done by computing the difference vectors */
/*     between a plate's vertices and the plate's centroid, scaling */
/*     those differences by (1 + XFRACT), then producing new vertices by */
/*     adding the scaled differences to the centroid. This process */
/*     doesn't affect the stored DSK data. */

/*     Plate expansion is also performed when surface points are mapped */
/*     to plates on which they lie, as is done for illumination angle */
/*     computations. */

/*     This parameter is user-adjustable. */


/*     The keyword for setting or retrieving this factor is */


/*     Greedy segment selection factor */
/*     ------------------------------- */

/*     The factor SGREED is used to slightly expand DSK segment */
/*     boundaries in order to select segments to consider for */
/*     ray-surface intercept computations. The effect of this factor is */
/*     to make the multi-segment intercept algorithm consider all */
/*     segments that are sufficiently close to the ray of interest, even */
/*     if the ray misses those segments. */

/*     This expansion is performed to prevent rays from passing through */
/*     a target object without any intersection being detected. Such */
/*     "false miss" conditions can occur due to round-off errors. */

/*     The exact way this parameter is used is dependent on the */
/*     coordinate system of the segment to which it applies, and the DSK */
/*     software implementation. This parameter may be changed in a */
/*     future version of SPICE. */


/*     The keyword for setting or retrieving this factor is */


/*     Segment pad margin */
/*     ------------------ */

/*     The segment pad margin is a scale factor used to determine when a */
/*     point resulting from a ray-surface intercept computation, if */
/*     outside the segment's boundaries, is close enough to the segment */
/*     to be considered a valid result. */

/*     This margin is required in order to make DSK segment padding */
/*     (surface data extending slightly beyond the segment's coordinate */
/*     boundaries) usable: if a ray intersects the pad surface outside */
/*     the segment boundaries, the pad is useless if the intercept is */
/*     automatically rejected. */

/*     However, an excessively large value for this parameter is */
/*     detrimental, since a ray-surface intercept solution found "in" a */
/*     segment can supersede solutions in segments farther from the */
/*     ray's vertex. Solutions found outside of a segment thus can mask */
/*     solutions that are closer to the ray's vertex by as much as the */
/*     value of this margin, when applied to a segment's boundary */
/*     dimensions. */

/*     The keyword for setting or retrieving this factor is */


/*     Surface-point membership margin */
/*     ------------------------------- */

/*     The surface-point membership margin limits the distance */
/*     between a point and a surface to which the point is */
/*     considered to belong. The margin is a scale factor applied */
/*     to the size of the segment containing the surface. */

/*     This margin is used to map surface points to outward */
/*     normal vectors at those points. */

/*     If this margin is set to an excessively small value, */
/*     routines that make use of the surface-point mapping won't */
/*     work properly. */


/*     The keyword for setting or retrieving this factor is */


/*     Angular rounding margin */
/*     ----------------------- */

/*     This margin specifies an amount by which angular values */
/*     may deviate from their proper ranges without a SPICE error */
/*     condition being signaled. */

/*     For example, if an input latitude exceeds pi/2 radians by a */
/*     positive amount less than this margin, the value is treated as */
/*     though it were pi/2 radians. */

/*     Units are radians. */


/*     This parameter is not user-adjustable. */

/*     The keyword for retrieving this parameter is */


/*     Longitude alias margin */
/*     ---------------------- */

/*     This margin specifies an amount by which a longitude */
/*     value can be outside a given longitude range without */
/*     being considered eligible for transformation by */
/*     addition or subtraction of 2*pi radians. */

/*     A longitude value, when compared to the endpoints of */
/*     a longitude interval, will be considered to be equal */
/*     to an endpoint if the value is outside the interval */
/*     differs from that endpoint by a magnitude less than */
/*     the alias margin. */


/*     Units are radians. */


/*     This parameter is not user-adjustable. */

/*     The keyword for retrieving this parameter is */


/*     End of include file dsktol.inc */


/*     External routines */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Initial values */


/*     Option flags */

    chkin_("DSKBRIEF", (ftnlen)8);

/*     Set the default value of the number of significant */
/*     digits in floating point outputs. */

    nsig = 6;

/*     Initialize banner. */

    s_copy(banner, " ", (ftnlen)80, (ftnlen)1);
    for (i__ = 1; i__ <= 78; ++i__) {
	*(unsigned char *)&banner[i__ - 1] = '=';
    }

/*     Display the DSKBRIEF version. */

    prcinf_("VERSION", (ftnlen)7);

/*     Get the DSK file name, and the names of any */
/*     additional files, from the command line. */

    s_copy(cmd, " ", (ftnlen)25500, (ftnlen)1);
    getcml_(cmd, (ftnlen)25500);
    if (s_cmp(cmd, " ", (ftnlen)25500, (ftnlen)1) == 0) {
	prcinf_("USAGE", (ftnlen)5);
	byebye_("SUCCESS", (ftnlen)7);
    }
    if (s_cmp(cmd + 25000, " ", (ftnlen)500, (ftnlen)1) != 0) {
	setmsg_("Input command is too long: non-blank characters were found "
		"past index #.", (ftnlen)72);
	errint_("#", &c__25000, (ftnlen)1);
	sigerr_("SPICE(COMMANDTOOLONG)", (ftnlen)21);
    }

/*     Identify all options specified on the command line. */

    fndnwd_(cmd, &c__1, &b, &e, (ftnlen)25500);
    while(b > 0) {
	s_copy(option, cmd + (b - 1), (ftnlen)32, e - (b - 1));
	i__ = isrchc_(option, &c__10, optlst, (ftnlen)32, (ftnlen)32);
	if (i__ > 0) {

/*           We have a recognized option. */

	    if (s_cmp(option, "-h", (ftnlen)32, (ftnlen)2) == 0) {
		prcinf_("HELP", (ftnlen)4);
		byebye_("SUCCESS", (ftnlen)7);
	    } else if (s_cmp(option, "-u", (ftnlen)32, (ftnlen)2) == 0) {
		prcinf_("USAGE", (ftnlen)5);
		byebye_("SUCCESS", (ftnlen)7);
	    } else if (s_cmp(option, "-v", (ftnlen)32, (ftnlen)2) == 0) {
		byebye_("SUCCESS", (ftnlen)7);
	    } else if (s_cmp(option, "-a", (ftnlen)32, (ftnlen)2) == 0) {
		all = TRUE_;
	    } else if (s_cmp(option, "-full", (ftnlen)32, (ftnlen)5) == 0) {

/*              The -full option implies the -seg option. */

		full = TRUE_;
		seg = TRUE_;
	    } else if (s_cmp(option, "-seg", (ftnlen)32, (ftnlen)4) == 0) {
		seg = TRUE_;
	    } else if (s_cmp(option, "-ext", (ftnlen)32, (ftnlen)4) == 0) {
		ext = TRUE_;
	    } else if (s_cmp(option, "-tg", (ftnlen)32, (ftnlen)3) == 0) {
		ust = TRUE_;
	    } else if (s_cmp(option, "-gaps", (ftnlen)32, (ftnlen)5) == 0) {
		gap = TRUE_;
	    } else if (s_cmp(option, "-d", (ftnlen)32, (ftnlen)2) == 0) {

/*              The next token should be an integer. */

		i__1 = e + 1;
		fndnwd_(cmd, &i__1, &b1, &e1, (ftnlen)25500);
		if (b1 > 0) {
		    s_copy(numstr, cmd + (b1 - 1), (ftnlen)32, e1 - (b1 - 1));
		    if (! beint_(numstr, (ftnlen)32)) {
			setmsg_("String <#> following -d option was not an i"
				"nteger.", (ftnlen)50);
			errch_("#", numstr, (ftnlen)1, (ftnlen)32);
			sigerr_("SPICE(INVALIDINTEGER)", (ftnlen)21);
		    } else {
			prsint_(numstr, &nsig, (ftnlen)32);
			nsig = brckti_(&nsig, &c__3, &c__17);
		    }

/*                 Erase the number from the command. */

		    s_copy(cmd + (b1 - 1), " ", e1 - (b1 - 1), (ftnlen)1);
		} else {
		    setmsg_("An integer in the range 6:17 must follow the -d"
			    " option.", (ftnlen)55);
		    sigerr_("SPICE(SYNTAXERROR)", (ftnlen)18);
		}
	    } else {
		setmsg_("BUG: unrecognized option <#>.", (ftnlen)29);
		errch_("#", option, (ftnlen)1, (ftnlen)32);
		sigerr_("SPICE(BUG)", (ftnlen)10);
	    }

/*           Erase the option from the command. */

	    s_copy(cmd + (b - 1), " ", e - (b - 1), (ftnlen)1);
	}
	i__1 = e + 1;
	fndnwd_(cmd, &i__1, &b, &e, (ftnlen)25500);
    }

/*     Check for option combinations that involve overrides. */

    opterr = FALSE_;
    if (full || seg) {
	all = FALSE_;
	ext = FALSE_;
	ust = FALSE_;
    }
    if (ust) {
	ext = TRUE_;
    }
    if (opterr) {
	tostdo_(errmsg, (ftnlen)255);
	tostdo_(" ", (ftnlen)1);
	tostdo_("Run DSKBRIEF without command line options to see program us"
		"age.", (ftnlen)63);
	tostdo_(" ", (ftnlen)1);
	byebye_("FAILURE", (ftnlen)7);
    }

/*     Load all files specified on the command line. */

    didfil = FALSE_;
    fndnwd_(cmd, &c__1, &b, &e, (ftnlen)25500);
    while(b > 0) {

/*        The current word is not an option; assume it's a file */
/*        name. */

	s_copy(fname, cmd + (b - 1), (ftnlen)255, e - (b - 1));
	if (! exists_(fname, (ftnlen)255)) {
	    setmsg_("Token <#> is neither a recognized option nor the name o"
		    "f an existing file.", (ftnlen)74);
	    errch_("#", fname, (ftnlen)1, (ftnlen)255);
	    sigerr_("SPICE(SYNTAXERROR)", (ftnlen)18);
	}
	furnsh_(fname, (ftnlen)255);
	didfil = TRUE_;
	i__1 = e + 1;
	fndnwd_(cmd, &i__1, &b, &e, (ftnlen)25500);
    }

/*     Initialize variables used to group segments. */

    nseg = 0;
    ng = 0;

/*     Now summarize all DSK files. */

    ktotal_("DSK", &ndsk, (ftnlen)3);
    diddsk = ndsk > 0;
    i__1 = ndsk;
    for (dskno = 1; dskno <= i__1; ++dskno) {

/*        Get the name and handle of the current DSK. */

	kdata_(&dskno, "DSK", fname, filtyp, source, &handle, &found, (ftnlen)
		3, (ftnlen)255, (ftnlen)10, (ftnlen)255);
	if (! found) {
	    setmsg_("KDATA did not found the DSK at index #. This is a bug.", 
		    (ftnlen)54);
	    errint_("#", &dskno, (ftnlen)1);
	    sigerr_("SPICE(BUG)", (ftnlen)10);
	}
	if (seg) {

/*           Summarize the current file segment by segment. */

	    dasopr_(fname, &handle, (ftnlen)255);
	    tostdo_(" ", (ftnlen)1);
/* Writing concatenation */
	    i__2[0] = 13, a__1[0] = "Summary for: ";
	    i__2[1] = 255, a__1[1] = fname;
	    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)268);
	    tostdo_(ch__1, (ftnlen)268);
	    tostdo_(" ", (ftnlen)1);

/*           Search forward through the file, extracting DLA segment */
/*           descriptors as we go. */

	    segno = 0;
	    dlabfs_(&handle, dladsc, &found);
	    while(found) {
		++segno;
		tostdo_(banner, (ftnlen)80);
		s_copy(outlin, "Segment number # of file #", (ftnlen)255, (
			ftnlen)26);
		repmi_(outlin, "#", &segno, outlin, (ftnlen)255, (ftnlen)1, (
			ftnlen)255);
		repmc_(outlin, "#", fname, outlin, (ftnlen)255, (ftnlen)1, (
			ftnlen)255, (ftnlen)255);
		tostdo_(outlin, (ftnlen)255);
		tostdo_(" ", (ftnlen)1);

/*              Get the DSK descriptor from the current segment. */

		dskgd_(&handle, dladsc, dskdsc);
		typcde = i_dnnt(&dskdsc[3]);

/*              Show: */

/*                 body ID */
/*                 surface ID */
/*                 reference frame */
/*                 data type */
/*                 class */
/*                 coordinate system */
/*                 surface coverage */
/*                 time bounds */

		dspdsc_(dskdsc, &c__1, &c__2, &nsig);
		dspdsc_(dskdsc, &c__1, &c__1, &nsig);
		dspdsc_(dskdsc, &c__1, &c__5, &nsig);
		dspdsc_(dskdsc, &c__1, &c__4, &nsig);
		dspdsc_(dskdsc, &c__1, &c__3, &nsig);
		dspdsc_(dskdsc, &c__1, &c__6, &nsig);
		dspdsc_(dskdsc, &c__1, &c__17, &nsig);
		dspdsc_(dskdsc, &c__1, &c__19, &nsig);
		dspdsc_(dskdsc, &c__1, &c__21, &nsig);
		dspdsc_(dskdsc, &c__1, &c__23, &nsig);
		if (full) {

/*                 Display type-specific parameters: */

		    if (typcde == 2) {

/*                    Display type 2 parameters. */

			sum02_(&handle, dladsc, &nsig);
		    } else if (typcde == 4) {

/*                    Display type 4 parameters. */

			sum04_(&handle, dladsc, &nsig);
		    } else {
			setmsg_("Segment # has data type #.", (ftnlen)26);
			errint_("#", &segno, (ftnlen)1);
			errint_("#", &typcde, (ftnlen)1);
			sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
		    }
		}

/*              Copy the last read DLA descriptor to CURRNT.  Fetch */
/*              the next descriptor. */

		movei_(dladsc, &c__8, currnt);
		dlafns_(&handle, currnt, dladsc, &found);
	    }
	    tostdo_(banner, (ftnlen)80);
	} else {

/*           This branch supports the default and -a cases. */

/*           Just store information about this file's segments. We'll */
/*           create a summary for the file, or a global summary, later. */

	    if (! all) {

/*              The list we're going to build applies only to the */
/*              current file. */

		nseg = 0;
	    }
	    dasopr_(fname, &handle, (ftnlen)255);
	    dlabfs_(&handle, dladsc, &found);
	    while(found) {

/*              We found another segment. */

		++nseg;
		if (nseg == 100000) {
		    setmsg_("Size of segment array is #; cannot add a new el"
			    "ement.", (ftnlen)53);
		    errint_("#", &c_b91, (ftnlen)1);
		    sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
		}

/*              Buffer the handle and descriptor. */

		hanlst[(i__3 = nseg - 1) < 100000 && 0 <= i__3 ? i__3 : 
			s_rnge("hanlst", i__3, "dskbrief_", (ftnlen)666)] = 
			handle;
		movei_(dladsc, &c__8, &dlads[(i__3 = (nseg << 3) - 8) < 
			800000 && 0 <= i__3 ? i__3 : s_rnge("dlads", i__3, 
			"dskbrief_", (ftnlen)668)]);

/*              Fetch the next descriptor. */

		dlafns_(&handle, dladsc, nxtdsc, &found);
		if (found) {
		    movei_(nxtdsc, &c__8, dladsc);
		}
	    }
	    if (! all) {

/*              Group the segments according to their attributes. */

		grpseg_((U_fp)attcmp_, &ust, &timtol, &nseg, hanlst, dlads, &
			c_b91, &c__10000, &ng, sgptrs, sgpool, seglst, segtm);
/* Writing concatenation */
		i__2[0] = 13, a__1[0] = "Summary for: ";
		i__2[1] = 255, a__1[1] = fname;
		s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)268);
		tostdo_(ch__1, (ftnlen)268);
		tostdo_(" ", (ftnlen)1);
		i__3 = ng;
		for (i__ = 1; i__ <= i__3; ++i__) {

/*                 Determine coordinate bounds for the Ith group. We'll */
/*                 need the coordinate bounds for each segment in the */
/*                 group. */

/*                 The first segment in the group has attributes, other */
/*                 than coordinate bounds, that are representative of */
/*                 the group. Extract these from the DSK descriptor of */
/*                 this segment. */

/*                 J is the index in the segment list of the first */
/*                 segment of the Ith segment group. */

		    j = seglst[(i__5 = sgptrs[(i__4 = i__ - 1) < 10000 && 0 <=
			     i__4 ? i__4 : s_rnge("sgptrs", i__4, "dskbrief_",
			     (ftnlen)707)] - 1) < 100000 && 0 <= i__5 ? i__5 :
			     s_rnge("seglst", i__5, "dskbrief_", (ftnlen)707)]
			    ;
		    grpsiz = 1;

/*                 Get the DSK descriptor for the first segment of the */
/*                 Ith group. This is the source of the "common" */
/*                 parameters. */

		    dskgd_(&hanlst[(i__4 = j - 1) < 100000 && 0 <= i__4 ? 
			    i__4 : s_rnge("hanlst", i__4, "dskbrief_", (
			    ftnlen)715)], &dlads[(i__5 = (j << 3) - 8) < 
			    800000 && 0 <= i__5 ? i__5 : s_rnge("dlads", i__5,
			     "dskbrief_", (ftnlen)715)], comdsc);

/*                 Initialize the coordinate bounds. */

		    minco1 = comdsc[16];
		    maxco1 = comdsc[17];
		    minco2 = comdsc[18];
		    maxco2 = comdsc[19];
		    minco3 = comdsc[20];
		    maxco3 = comdsc[21];

/*                 Store the domain coordinate bounds for use in gap */
/*                 detection. */

		    xbds[(i__4 = (grpsiz << 1) - 2) < 200000 && 0 <= i__4 ? 
			    i__4 : s_rnge("xbds", i__4, "dskbrief_", (ftnlen)
			    730)] = minco1;
		    xbds[(i__4 = (grpsiz << 1) - 1) < 200000 && 0 <= i__4 ? 
			    i__4 : s_rnge("xbds", i__4, "dskbrief_", (ftnlen)
			    731)] = maxco1;
		    ybds[(i__4 = (grpsiz << 1) - 2) < 200000 && 0 <= i__4 ? 
			    i__4 : s_rnge("ybds", i__4, "dskbrief_", (ftnlen)
			    733)] = minco2;
		    ybds[(i__4 = (grpsiz << 1) - 1) < 200000 && 0 <= i__4 ? 
			    i__4 : s_rnge("ybds", i__4, "dskbrief_", (ftnlen)
			    734)] = maxco2;

/*                 Initialize the time bounds. */

		    mintdb = comdsc[22];
		    maxtdb = comdsc[23];

/*                 Now update the bounds using those of each segment in */
/*                 the group. We can skip the first segment. */

		    node = lnknxt_(&sgptrs[(i__4 = i__ - 1) < 10000 && 0 <= 
			    i__4 ? i__4 : s_rnge("sgptrs", i__4, "dskbrief_", 
			    (ftnlen)746)], sgpool);
		    while(node > 0) {
			j = seglst[(i__4 = node - 1) < 100000 && 0 <= i__4 ? 
				i__4 : s_rnge("seglst", i__4, "dskbrief_", (
				ftnlen)750)];
			dskgd_(&hanlst[(i__4 = j - 1) < 100000 && 0 <= i__4 ? 
				i__4 : s_rnge("hanlst", i__4, "dskbrief_", (
				ftnlen)752)], &dlads[(i__5 = (j << 3) - 8) < 
				800000 && 0 <= i__5 ? i__5 : s_rnge("dlads", 
				i__5, "dskbrief_", (ftnlen)752)], dskdsc);

/*                    Longitude is a special case. As usual. */

			if (comdsc[5] == 1. || comdsc[5] == 4.) {

/*                       The segment longitude ranges must be compatible */
/*                       in order to take meaningful minima and maxima */
/*                       of the longitude bounds. We'll regularize the */
/*                       combination of the longitude bounds we have so */
/*                       far and those of the current segment. */

			    nivals = 2;
			    lonbds[0] = minco1;
			    lonbds[1] = maxco1;
			    lonbds[2] = dskdsc[16];
			    lonbds[3] = dskdsc[17];
			    reglon_(&nivals, lonbds, &c__4, &nout, &minlon, &
				    maxlon, outbds, srcs);

/*                       The output longitude bounds belong to a */
/*                       consistent range: either -180:180 or 0:360. */


/*                       Note that we can't use the extrema from the */
/*                       first segment of the group that we obtained */
/*                       before the loop start, since these bounds */
/*                       haven't been regularized. */

			    minco1 = outbds[0];
			    maxco1 = outbds[1];
			    i__4 = nout;
			    for (k = 2; k <= i__4; ++k) {
/* Computing MIN */
				d__1 = minco1, d__2 = outbds[(i__5 = (k << 1) 
					- 2) < 8 && 0 <= i__5 ? i__5 : s_rnge(
					"outbds", i__5, "dskbrief_", (ftnlen)
					789)];
				minco1 = min(d__1,d__2);
/* Computing MAX */
				d__1 = maxco1, d__2 = outbds[(i__5 = (k << 1) 
					- 1) < 8 && 0 <= i__5 ? i__5 : s_rnge(
					"outbds", i__5, "dskbrief_", (ftnlen)
					790)];
				maxco1 = max(d__1,d__2);
			    }
			} else {
			    minco1 = min(minco1,dskdsc[16]);
			    maxco1 = max(maxco1,dskdsc[17]);
			}
			minco2 = min(minco2,dskdsc[18]);
			maxco2 = max(maxco2,dskdsc[19]);
			minco3 = min(minco3,dskdsc[20]);
			maxco3 = max(maxco3,dskdsc[21]);

/*                    Store the domain coordinate bounds for use in gap */
/*                    detection. Note that we must use the segment */
/*                    bounds, not the group bounds. */

			++grpsiz;
			xbds[(i__4 = (grpsiz << 1) - 2) < 200000 && 0 <= i__4 
				? i__4 : s_rnge("xbds", i__4, "dskbrief_", (
				ftnlen)812)] = dskdsc[16];
			xbds[(i__4 = (grpsiz << 1) - 1) < 200000 && 0 <= i__4 
				? i__4 : s_rnge("xbds", i__4, "dskbrief_", (
				ftnlen)813)] = dskdsc[17];
			ybds[(i__4 = (grpsiz << 1) - 2) < 200000 && 0 <= i__4 
				? i__4 : s_rnge("ybds", i__4, "dskbrief_", (
				ftnlen)815)] = dskdsc[18];
			ybds[(i__4 = (grpsiz << 1) - 1) < 200000 && 0 <= i__4 
				? i__4 : s_rnge("ybds", i__4, "dskbrief_", (
				ftnlen)816)] = dskdsc[19];

/*                    Update the time bounds. */

			mintdb = min(mintdb,dskdsc[22]);
			maxtdb = max(maxtdb,dskdsc[23]);

/*                    Look at the next segment of the group. */

			node = lnknxt_(&node, sgpool);
		    }

/*                 Create a DSK descriptor that contains the group's */
/*                 coordinate bounds. */

		    moved_(comdsc, &c__24, grpdsc);
		    grpdsc[16] = minco1;
		    grpdsc[17] = maxco1;
		    grpdsc[18] = minco2;
		    grpdsc[19] = maxco2;
		    grpdsc[20] = minco3;
		    grpdsc[21] = maxco3;

/*                 Show: */

/*                    body ID */
/*                    surface ID */
/*                    reference frame */
/*                    coordinate system */
/*                    surface coverage */

		    dspdsc_(grpdsc, &c__1, &c__2, &nsig);
		    dspdsc_(grpdsc, &c__1, &c__1, &nsig);
		    dspdsc_(grpdsc, &c__1, &c__5, &nsig);

/*                 Display data type and class if -ext was */
/*                 specified. */

		    if (ext) {
			dspdsc_(grpdsc, &c__1, &c__4, &nsig);
			dspdsc_(grpdsc, &c__1, &c__3, &nsig);
		    }
		    dspdsc_(grpdsc, &c__1, &c__6, &nsig);
		    dspdsc_(grpdsc, &c__1, &c__17, &nsig);
		    dspdsc_(grpdsc, &c__1, &c__19, &nsig);
		    dspdsc_(grpdsc, &c__1, &c__21, &nsig);
		    if (ext) {

/*                    Display time range. */

			dspdsc_(grpdsc, &c__1, &c__23, &nsig);
		    }

/*                 If time is not used for segment grouping and */
/*                 the members of the group don't have compatible */
/*                 times, issue a warning. */

		    if (! ust && ! segtm[(i__4 = i__ - 1) < 100000 && 0 <= 
			    i__4 ? i__4 : s_rnge("segtm", i__4, "dskbrief_", (
			    ftnlen)882)]) {
			tostdo_("    ***Segments have inconsistent time cove"
				"rage.***", (ftnlen)51);
		    }

/*                 Determine whether this set of segments has */
/*                 coverage gaps. */

		    corsys = i_dnnt(&grpdsc[5]);
		    zzdbrgap_(&corsys, &grpsiz, xbds, ybds, &c_b91, &ncomp, 
			    gxbds, gybds);

/*                 If the user has commanded gap display, and if */
/*                 there are gaps, here's where we display them. */

		    dspgap_(&gap, &corsys, &nsig, &ncomp, gxbds, gybds);
		    tostdo_(" ", (ftnlen)1);
		}

/*              End of group loop. */

	    }

/*           End of default case. */

	}

/*        End of the case conditional block. */

    }

/*     We've processed all DSK files listed on the command line. */

    if (all && diddsk && ! seg) {

/*        Group the segments according to their attributes. */

	grpseg_((U_fp)attcmp_, &ust, &timtol, &nseg, hanlst, dlads, &c_b91, &
		c__10000, &ng, sgptrs, sgpool, seglst, segtm);
	tostdo_("Summary for: all DSK files", (ftnlen)26);
	tostdo_(" ", (ftnlen)1);
	i__1 = ng;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Determine coordinate bounds for the Ith group. */
/*           We'll need the coordinate bounds for each segment */
/*           in the group. */

/*           The first segment in the group has attributes, other */
/*           than coordinate bounds, that are representative of */
/*           the group. Extract these from the DSK descriptor of */
/*           this segment. */

/*           J is the index in the segment list of the first */
/*           segment of the Ith segment group. */

	    j = seglst[(i__4 = sgptrs[(i__3 = i__ - 1) < 10000 && 0 <= i__3 ? 
		    i__3 : s_rnge("sgptrs", i__3, "dskbrief_", (ftnlen)947)] 
		    - 1) < 100000 && 0 <= i__4 ? i__4 : s_rnge("seglst", i__4,
		     "dskbrief_", (ftnlen)947)];

/*           Get the DSK descriptor for the first segment of */
/*           the Ith group. This is the source of the "common" */
/*           parameters. */

	    dskgd_(&hanlst[(i__3 = j - 1) < 100000 && 0 <= i__3 ? i__3 : 
		    s_rnge("hanlst", i__3, "dskbrief_", (ftnlen)954)], &dlads[
		    (i__4 = (j << 3) - 8) < 800000 && 0 <= i__4 ? i__4 : 
		    s_rnge("dlads", i__4, "dskbrief_", (ftnlen)954)], comdsc);

/*           Initialize the coordinate bounds. */

	    minco1 = comdsc[16];
	    maxco1 = comdsc[17];
	    minco2 = comdsc[18];
	    maxco2 = comdsc[19];
	    minco3 = comdsc[20];
	    maxco3 = comdsc[21];
	    grpsiz = 1;

/*           Store the domain coordinate bounds for use in gap */
/*           detection. */

	    xbds[(i__3 = (grpsiz << 1) - 2) < 200000 && 0 <= i__3 ? i__3 : 
		    s_rnge("xbds", i__3, "dskbrief_", (ftnlen)971)] = minco1;
	    xbds[(i__3 = (grpsiz << 1) - 1) < 200000 && 0 <= i__3 ? i__3 : 
		    s_rnge("xbds", i__3, "dskbrief_", (ftnlen)972)] = maxco1;
	    ybds[(i__3 = (grpsiz << 1) - 2) < 200000 && 0 <= i__3 ? i__3 : 
		    s_rnge("ybds", i__3, "dskbrief_", (ftnlen)974)] = minco2;
	    ybds[(i__3 = (grpsiz << 1) - 1) < 200000 && 0 <= i__3 ? i__3 : 
		    s_rnge("ybds", i__3, "dskbrief_", (ftnlen)975)] = maxco2;

/*           Initialize the time bounds. */

	    mintdb = comdsc[22];
	    maxtdb = comdsc[23];

/*           Now update the bounds using those of each segment in */
/*           the group. We can skip the first segment. */

	    node = lnknxt_(&sgptrs[(i__3 = i__ - 1) < 10000 && 0 <= i__3 ? 
		    i__3 : s_rnge("sgptrs", i__3, "dskbrief_", (ftnlen)987)], 
		    sgpool);
	    while(node > 0) {
		j = seglst[(i__3 = node - 1) < 100000 && 0 <= i__3 ? i__3 : 
			s_rnge("seglst", i__3, "dskbrief_", (ftnlen)991)];
		dskgd_(&hanlst[(i__3 = j - 1) < 100000 && 0 <= i__3 ? i__3 : 
			s_rnge("hanlst", i__3, "dskbrief_", (ftnlen)993)], &
			dlads[(i__4 = (j << 3) - 8) < 800000 && 0 <= i__4 ? 
			i__4 : s_rnge("dlads", i__4, "dskbrief_", (ftnlen)993)
			], dskdsc);

/*              Longitude is a special case. As usual. */

		if (comdsc[5] == 1. || comdsc[5] == 4.) {

/*                 The segment longitude ranges must be compatible in */
/*                 order to take meaningful minima and maxima of the */
/*                 longitude bounds. We'll regularize the combination of */
/*                 the longitude bounds we have so far and those of the */
/*                 current segment. */

		    nivals = 2;
		    lonbds[0] = minco1;
		    lonbds[1] = maxco1;
		    lonbds[2] = dskdsc[16];
		    lonbds[3] = dskdsc[17];
		    reglon_(&nivals, lonbds, &c__4, &nout, &minlon, &maxlon, 
			    outbds, srcs);

/*                 The output longitude bounds belong to a consistent */
/*                 range: either -180:180 or 0:360. */

/*                 Note that we can't use the extrema from the first */
/*                 segment of the group that we obtained before the */
/*                 loop start, since these bounds haven't been */
/*                 regularized. */

		    minco1 = outbds[0];
		    maxco1 = outbds[1];
		    i__3 = nout;
		    for (k = 2; k <= i__3; ++k) {
/* Computing MIN */
			d__1 = minco1, d__2 = outbds[(i__4 = (k << 1) - 2) < 
				8 && 0 <= i__4 ? i__4 : s_rnge("outbds", i__4,
				 "dskbrief_", (ftnlen)1029)];
			minco1 = min(d__1,d__2);
/* Computing MAX */
			d__1 = maxco1, d__2 = outbds[(i__4 = (k << 1) - 1) < 
				8 && 0 <= i__4 ? i__4 : s_rnge("outbds", i__4,
				 "dskbrief_", (ftnlen)1030)];
			maxco1 = max(d__1,d__2);
		    }
		} else {
		    minco1 = min(minco1,dskdsc[16]);
		    maxco1 = max(maxco1,dskdsc[17]);
		}
		minco2 = min(minco2,dskdsc[18]);
		maxco2 = max(maxco2,dskdsc[19]);
		minco3 = min(minco3,dskdsc[20]);
		maxco3 = max(maxco3,dskdsc[21]);

/*              Store the domain coordinate bounds for use in gap */
/*              detection. Note that we must use the segment bounds, */
/*              not the group bounds. */

		++grpsiz;
		xbds[(i__3 = (grpsiz << 1) - 2) < 200000 && 0 <= i__3 ? i__3 :
			 s_rnge("xbds", i__3, "dskbrief_", (ftnlen)1051)] = 
			dskdsc[16];
		xbds[(i__3 = (grpsiz << 1) - 1) < 200000 && 0 <= i__3 ? i__3 :
			 s_rnge("xbds", i__3, "dskbrief_", (ftnlen)1052)] = 
			dskdsc[17];
		ybds[(i__3 = (grpsiz << 1) - 2) < 200000 && 0 <= i__3 ? i__3 :
			 s_rnge("ybds", i__3, "dskbrief_", (ftnlen)1054)] = 
			dskdsc[18];
		ybds[(i__3 = (grpsiz << 1) - 1) < 200000 && 0 <= i__3 ? i__3 :
			 s_rnge("ybds", i__3, "dskbrief_", (ftnlen)1055)] = 
			dskdsc[19];

/*              Update the group time bounds. */

		mintdb = min(mintdb,dskdsc[22]);
		maxtdb = max(maxtdb,dskdsc[23]);

/*              Look at the next segment of the group. */

		node = lnknxt_(&node, sgpool);
	    }

/*           Now display information about the group. */

/*           Create a DSK descriptor that contains the */
/*           group's coordinate bounds. */

	    moved_(comdsc, &c__24, grpdsc);
	    grpdsc[16] = minco1;
	    grpdsc[17] = maxco1;
	    grpdsc[18] = minco2;
	    grpdsc[19] = maxco2;
	    grpdsc[20] = minco3;
	    grpdsc[21] = maxco3;

/*           Show: */

/*              body ID */
/*              surface ID */
/*              reference frame */
/*              coordinate system */
/*              surface coverage */

	    dspdsc_(grpdsc, &c__1, &c__2, &nsig);
	    dspdsc_(grpdsc, &c__1, &c__1, &nsig);
	    dspdsc_(grpdsc, &c__1, &c__5, &nsig);
	    if (ext) {
		dspdsc_(grpdsc, &c__1, &c__4, &nsig);
		dspdsc_(grpdsc, &c__1, &c__3, &nsig);
	    }
	    dspdsc_(grpdsc, &c__1, &c__6, &nsig);
	    dspdsc_(grpdsc, &c__1, &c__17, &nsig);
	    dspdsc_(grpdsc, &c__1, &c__19, &nsig);
	    dspdsc_(grpdsc, &c__1, &c__21, &nsig);
	    if (ext) {
		dspdsc_(grpdsc, &c__1, &c__23, &nsig);
	    }
	    if (! ust && ! segtm[(i__3 = i__ - 1) < 100000 && 0 <= i__3 ? 
		    i__3 : s_rnge("segtm", i__3, "dskbrief_", (ftnlen)1111)]) 
		    {
		tostdo_("    ***Segments have inconsistent time coverage.***",
			 (ftnlen)51);
	    }

/*           Determine whether this set of segments has coverage gaps. */

	    corsys = i_dnnt(&grpdsc[5]);
	    zzdbrgap_(&corsys, &grpsiz, xbds, ybds, &c_b91, &ncomp, gxbds, 
		    gybds);

/*           If the user has commanded gap display, and if */
/*           there are gaps, here's where we display them. */

	    dspgap_(&gap, &corsys, &nsig, &ncomp, gxbds, gybds);
	    tostdo_(" ", (ftnlen)1);
	}

/*        End of group loop. */

    }
    if (didfil && ! diddsk) {
	tostdo_("No DSK files were provided -- no summary will be displayed.",
		 (ftnlen)59);
	tostdo_(" ", (ftnlen)1);
	tostdo_("Run DSKBRIEF without command line options to see program us"
		"age.", (ftnlen)63);
	tostdo_(" ", (ftnlen)1);
    }
    chkout_("DSKBRIEF", (ftnlen)8);
    return 0;
} /* MAIN__ */

/* Main program alias */ int dskbrief_ () { MAIN__ (); return 0; }
