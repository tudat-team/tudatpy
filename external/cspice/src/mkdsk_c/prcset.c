/* prcset.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__10 = 10;
static integer c__100 = 100;

/* $Procedure PRCSET ( Process setup file for MKDSK---umbrella routine ) */
/* Subroutine */ int prcset_0_(int n__, char *setup, char *input, char *
	output, char *cmtfil, integer *surfid, integer *centid, char *frame, 
	doublereal *first, doublereal *last, integer *dclass, integer *dtype, 
	char *aunits, char *dunits, integer *corsys, doublereal *corpar, 
	doublereal *mncor1, doublereal *mxcor1, doublereal *mncor2, 
	doublereal *mxcor2, integer *pltype, doublereal *voxscl, integer *
	cgrscl, logical *wrap, logical *mkncap, logical *mkscap, logical *
	rowmaj, logical *topdwn, logical *leftrt, doublereal *refval, 
	doublereal *hscale, integer *ncols, integer *nrows, doublereal *
	lftcor, doublereal *topcor, doublereal *colstp, doublereal *rowstp, 
	logical *makvpm, ftnlen setup_len, ftnlen input_len, ftnlen 
	output_len, ftnlen cmtfil_len, ftnlen frame_len, ftnlen aunits_len, 
	ftnlen dunits_len)
{
    /* Initialized data */

    static char csynms[255*4] = "LATITUDINAL                                "
	    "                                                                "
	    "                                                                "
	    "                                                                "
	    "                    " "CYLINDRICAL                              "
	    "                                                                "
	    "                                                                "
	    "                                                                "
	    "                      " "RECTANGULAR                            "
	    "                                                                "
	    "                                                                "
	    "                                                                "
	    "                        " "PLANETODETIC                         "
	    "                                                                "
	    "                                                                "
	    "                                                                "
	    "                          ";
    static logical init = FALSE_;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    char cval[255];
    doublereal f;
    integer i__, n;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    logical csfnd, fsfnd;
    extern /* Subroutine */ int ucase_(char *, char *, ftnlen, ftnlen), 
	    errch_(char *, char *, ftnlen, ftnlen);
    logical found;
    extern /* Subroutine */ int errdp_(char *, doublereal *, ftnlen);
    extern logical eqstr_(char *, char *, ftnlen, ftnlen);
    char words[255*2];
    extern /* Subroutine */ int ljust_(char *, char *, ftnlen, ftnlen);
    char vtype[1];
    extern /* Subroutine */ int bods2c_(char *, integer *, logical *, ftnlen),
	     str2et_(char *, doublereal *, ftnlen);
    extern logical failed_(void);
    doublereal re;
    extern /* Subroutine */ int cleard_(integer *, doublereal *);
    doublereal rp;
    integer framid;
    extern integer esrchc_(char *, integer *, char *, ftnlen, ftnlen);
    char centnm[36];
    extern logical exists_(char *, ftnlen), return_(void);
    char surfnm[36], sysnam[255], timstr[50], unistr[255*2];
    integer ntoken;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), furnsh_(char *, ftnlen), 
	    gcpool_(char *, integer *, integer *, integer *, char *, logical *
	    , ftnlen, ftnlen), dtpool_(char *, logical *, integer *, char *, 
	    ftnlen, ftnlen);
    char lsk[255];
    extern /* Subroutine */ int gipool_(char *, integer *, integer *, integer 
	    *, integer *, logical *, ftnlen), srfscc_(char *, integer *, 
	    integer *, logical *, ftnlen), errint_(char *, integer *, ftnlen),
	     lparse_(char *, char *, integer *, integer *, char *, ftnlen, 
	    ftnlen, ftnlen), namfrm_(char *, integer *, ftnlen), gdpool_(char 
	    *, integer *, integer *, integer *, doublereal *, logical *, 
	    ftnlen);

/* $ Abstract */

/*     Umbrella routine for MKDSK setup file processing. */

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

/*     None. */

/* $ Keywords */

/*     MKDSK */

/* $ Declarations */
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

/* $ Brief_I/O */

/*     VARIABLE  I/O  Entry points */
/*     --------  ---  -------------------------------------------------- */
/*     SETUP      I   GETSET */
/*     INPUT     I/O  GETSET */
/*     OUTPUT    I/O  GETSET */
/*     CMTFIL     O   GETSET */


/* $ Detailed_Input */

/*     See the entry points for a description of their inputs. */

/* $ Detailed_Output */

/*     See the entry points for a description of their outputs. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     If this routine is called directly, the error SPICE(BOGUSENTRY) */
/*     is signaled.  See the entry points for descriptions of */
/*     exceptions specific to those routines. */

/* $ Files */

/*     This suite of routines reads and returns information from an */
/*     MKDSK setup file. See the MKDSK setup template for a list of */
/*     supported setup parameters. */

/* $ Particulars */

/*     The entry points in this package are */

/*        GETSET  {Get setup file information for MKDSK} */
/*        GETTYP  {Get segment data type} */
/*        GETGEN  {Get general DSK parameters} */
/*        GETP02  {Get type 2 parameters} */

/*     GETSET must be called before the other entry points may be called. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     This routine is intended for use only within the program */
/*     MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     B.V. Semenov    (JPL) */

/* $ Version */

/* -    MKDSK Version 3.0.0, 08-MAR-2017 (NJB) (BVS) */

/*        Updated to support automatic voxel scale setting. */
/*        Updated to support plate type 5 (height grid) input */
/*        format. */

/*        Some error handling bugs were corrected. */

/*        Last update 19-JAN-2016 (NJB) */

/*           Updated to support surface name-ID translation */
/*           and new coordinate systems. Updated header. */

/* -    MKDSK Version 2.0.0, 29-JUN-2010 (NJB) */

/*        Updated shape and DSK file keywords. */

/* -    MKDSK Version 1.0.0, 15-APR-2010 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Upper bound on coarse voxel scale. This is used to weed out */
/*     nonsense input values. */


/*     Local variables */


/*     Saved variables */


/*     Initial values */

    /* Parameter adjustments */
    if (corpar) {
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_getset;
	case 2: goto L_gettyp;
	case 3: goto L_getgen;
	case 4: goto L_getp02;
	case 5: goto L_getg05;
	}

    if (return_()) {
	return 0;
    }
    chkin_("PRCSET", (ftnlen)6);
    sigerr_("SPICE(BOGUSENTRY)", (ftnlen)17);
    chkout_("PRCSET", (ftnlen)6);
    return 0;
/* $Procedure GETSET ( Get setup file information for MKDSK ) */

L_getset:
/* $ Abstract */

/*     Get the names of the input shape file, the output DSK file, and */
/*     the comment file from an MKDSK setup file. */

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

/*     None. */

/* $ Keywords */

/*     MKDSK */

/* $ Declarations */

/*     CHARACTER*(*)         SETUP */
/*     CHARACTER*(*)         INPUT */
/*     CHARACTER*(*)         OUTPUT */
/*     CHARACTER*(*)         CMTFIL */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     SETUP      I   Name of setup file. */
/*     INPUT     I/O  Name of input shape file. */
/*     OUTPUT    I/O  Name of output DSK file. */
/*     CMTFIL     O   Name of comment file. */

/* $ Detailed_Input */

/*     SETUP          is the name of an MKDSK setup file. */

/*     INPUT          is the name of a shape file to be converted */
/*                    to DSK format. This file conforms to the */
/*                    format specification given by the MKDSK */
/*                    User's Guide. */

/*     OUTPUT         is the name of a DSK file to be written. */
/*                    OUTPUT must be a new file. */

/* $ Detailed_Output */

/*     CMTFIL         is the name of a comment file whose contents */
/*                    are to be added to the comment area of */
/*                    the DSK file created by this program.  The */
/*                    comment file contents precede the default */
/*                    comments added by MKDSK. */

/*                    If no comment file is specified, CMTFIL is */
/*                    returned blank. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the setup file name is blank, the error */
/*         SPICE(BLANKFILENAME) is signaled. */

/*     2)  If the setup file doesn't exist, the error */
/*         SPICE(FILENOTFOUND) is signaled. */

/*     3)  If the name of the surface to be represented by the DSK is */
/*         not specified in the setup file, the error */
/*         SPICE(NOSURFACENAME) is signaled. */

/*     4)  If the surface name is present but cannot be mapped to an */
/*         integer ID code, the error SPICE(NOTRANSLATION) is signaled. */

/*     5)  If the input shape file doesn't exist, the error */
/*         SPICE(FILENOTFOUND) is signaled. */

/*     6)  If the DSK file name is not specified on the command line */
/*         and doesn't appear in the setup file, the error */
/*         SPICE(NOFILESPEC) is signaled. */

/*     7)  If the output file name matches that of an existing file, */
/*         the error SPICE(FILEEXISTS) is signaled. */

/*     8)  If a comment file keyword is present in the setup file */
/*         but the associated value does not parse as a quoted string, */
/*         the error SPICE(TYPEMISMATCH) is signaled. */

/*     9)  If an DSK start time is present in the setup file */
/*         but the associated value does not parse as a quoted */
/*         string, the error SPICE(TYPEMISMATCH) is signaled. */

/*     10) If an DSK stop time is present in the setup file */
/*         but the associated value does not parse as a quoted */
/*         string, the error SPICE(TYPEMISMATCH) is signaled. */


/* $ Files */

/*     See the descriptions of INPUT and OUTPUT above. */

/* $ Particulars */

/*     None. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     This routine is intended for use only within the program */
/*     MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     B.V. Semenov    (JPL) */

/* $ Version */

/* -    MKDSK Version 3.0.0, 24-FEB-2017 (NJB) (BVS) */

/*        Leapseconds kernel assignment is now optional. */


/*        19-JAN-2016 (NJB) */

/*        Updated to support surface name-ID translation */
/*        and new coordinate systems. */

/*        Alpha DSK Version 1.0.0, 15-APR-2010 (NJB) */

/* -& */
    if (return_()) {
	return 0;
    }
    chkin_("GETSET", (ftnlen)6);

/*     Check the setup file name. */

    if (s_cmp(setup, " ", setup_len, (ftnlen)1) == 0) {
	setmsg_("Setup file name may not be blank.", (ftnlen)33);
	sigerr_("SPICE(BLANKFILENAME)", (ftnlen)20);
	chkout_("GETSET", (ftnlen)6);
	return 0;
    }
    if (! exists_(setup, setup_len)) {
	setmsg_("Setup file <#> was not found.", (ftnlen)29);
	errch_("#", setup, (ftnlen)1, setup_len);
	sigerr_("SPICE(FILENOTFOUND)", (ftnlen)19);
	chkout_("GETSET", (ftnlen)6);
	return 0;
    }

/*     Load the setup file. */

    furnsh_(setup, setup_len);
    if (failed_()) {
	chkout_("GETSET", (ftnlen)6);
	return 0;
    }

/*     Check the input file name.  If the name is available, it */
/*     supersedes an input file name supplied in the setup file. */
/*     If the name is blank, the setup file must supply the name. */

    if (s_cmp(input, " ", input_len, (ftnlen)1) == 0) {

/*        Extract the input file name from the kernel pool. */
/*        Set the INPUT argument to the specified file name. */

	gcpool_("INPUT_SHAPE_FILE", &c__1, &c__1, &n, input, &found, (ftnlen)
		16, input_len);
	if (! found) {
	    setmsg_("Input file was not specified on the command line or in "
		    "the setup file.", (ftnlen)70);
	    sigerr_("SPICE(NOFILESPEC)", (ftnlen)17);
	    chkout_("GETSET", (ftnlen)6);
	    return 0;
	}
    }
    if (! exists_(input, input_len)) {
	setmsg_("Input file <#> was not found.", (ftnlen)29);
	errch_("#", input, (ftnlen)1, input_len);
	sigerr_("SPICE(FILENOTFOUND)", (ftnlen)19);
	chkout_("GETSET", (ftnlen)6);
	return 0;
    }

/*     Check the output file name.  If the name is available, it */
/*     supersedes an output file name supplied in the setup file. */
/*     If the name is blank, the setup file must supply the name. */

    if (s_cmp(output, " ", output_len, (ftnlen)1) == 0) {

/*        Extract the output file name from the kernel pool. */
/*        Set the INPUT argument to the specified file name. */

	gcpool_("OUTPUT_DSK_FILE", &c__1, &c__1, &n, output, &found, (ftnlen)
		15, output_len);
	if (! found) {
	    setmsg_("Output file was not specified on the command line or in"
		    " the setup file.", (ftnlen)71);
	    sigerr_("SPICE(NOFILESPEC)", (ftnlen)17);
	    chkout_("GETSET", (ftnlen)6);
	    return 0;
	}
    }
    if (exists_(output, output_len)) {
	setmsg_("Output file <#> already exists.", (ftnlen)31);
	errch_("#", output, (ftnlen)1, output_len);
	sigerr_("SPICE(FILEEXISTS)", (ftnlen)17);
	chkout_("GETSET", (ftnlen)6);
	return 0;
    }

/*     Obtain the name of the leapseconds kernel and load the kernel. */

    gcpool_("LEAPSECONDS_FILE", &c__1, &c__1, &n, lsk, &found, (ftnlen)16, (
	    ftnlen)255);
    if (found) {
	furnsh_(lsk, (ftnlen)255);
	if (failed_()) {
	    chkout_("GETSET", (ftnlen)6);
	    return 0;
	}
    }

/*     See whether a comment file specification was given. */
/*     Capture the file name if so. */

    dtpool_("COMMENT_FILE", &found, &n, vtype, (ftnlen)12, (ftnlen)1);
    if (found) {
	if (*(unsigned char *)vtype != 'C') {
	    setmsg_("Comment file name was not given a character string valu"
		    "e in the setup file.", (ftnlen)75);
	    sigerr_("SPICE(TYPEMISMATCH)", (ftnlen)19);
	    chkout_("GETSET", (ftnlen)6);
	    return 0;
	}
	gcpool_("COMMENT_FILE", &c__1, &c__1, &n, cmtfil, &found, (ftnlen)12, 
		cmtfil_len);
    }
    if (! found) {
	s_copy(cmtfil, " ", cmtfil_len, (ftnlen)1);
    }
    init = TRUE_;
    chkout_("GETSET", (ftnlen)6);
    return 0;

/*     Get segment data type. */


L_gettyp:
    chkin_("GETTYP", (ftnlen)6);
    gipool_("DATA_TYPE", &c__1, &c__1, &n, dtype, &found, (ftnlen)9);
    if (! found) {
	setmsg_("No segment data type was provided in the setup file.", (
		ftnlen)52);
	sigerr_("SPICE(MISSINGDATATYPE)", (ftnlen)22);
	chkout_("GETTYP", (ftnlen)6);
	return 0;
    }
    chkout_("GETTYP", (ftnlen)6);
    return 0;

/*     Get general DSK parameters. */


L_getgen:

/*     Fetch the surface name; convert the name to an ID code. */

    chkin_("GETGEN", (ftnlen)6);
    gcpool_("SURFACE_NAME", &c__1, &c__1, &n, surfnm, &found, (ftnlen)12, (
	    ftnlen)36);
    if (! found) {
	setmsg_("No surface name was provided in the setup file. Note that t"
		"he surface must be specified by a string.", (ftnlen)100);
	sigerr_("SPICE(MISSINGSURFACE)", (ftnlen)21);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Fetch the central body name; convert the name to an ID code. */

    gcpool_("CENTER_NAME", &c__1, &c__1, &n, centnm, &found, (ftnlen)11, (
	    ftnlen)36);
    if (! found) {
	setmsg_("No central body name was provided in the setup file. Note t"
		"hat the body must be specified by a string.", (ftnlen)102);
	sigerr_("SPICE(MISSINGCENTER)", (ftnlen)20);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Convert the central body name to an ID code. */

    bods2c_(centnm, centid, &found, (ftnlen)36);
    if (! found) {
	setmsg_("The central body name # could not be mapped to an integer I"
		"D code.", (ftnlen)66);
	errch_("#", centnm, (ftnlen)1, (ftnlen)36);
	sigerr_("SPICE(NOTRANSLATION)", (ftnlen)20);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Convert the surface name to an ID code. */

    srfscc_(surfnm, centid, surfid, &found, (ftnlen)36);
    if (! found) {
	setmsg_("The surface name # could not be mapped to an integer ID cod"
		"e.", (ftnlen)61);
	errch_("#", surfnm, (ftnlen)1, (ftnlen)36);
	sigerr_("SPICE(NOTRANSLATION)", (ftnlen)20);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     See whether an output DSK start time was given. */
/*     Capture the value if so. */

    dtpool_("START_TIME", &found, &n, vtype, (ftnlen)10, (ftnlen)1);
    if (found) {
	if (*(unsigned char *)vtype != 'C') {
	    setmsg_("DSK start time was not given a character string value i"
		    "n the setup file.", (ftnlen)72);
	    sigerr_("SPICE(TYPEMISMATCH)", (ftnlen)19);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gcpool_("START_TIME", &c__1, &c__1, &n, timstr, &found, (ftnlen)10, (
		ftnlen)50);
	if (found) {
	    str2et_(timstr, first, (ftnlen)50);
	    if (failed_()) {
		chkout_("GETGEN", (ftnlen)6);
		return 0;
	    }
	} else {
	    setmsg_("DTPOOL says a start time was provided in the setup file"
		    ", but GCPOOL can't find it (?)", (ftnlen)85);
	    sigerr_("SPICE(BUG)", (ftnlen)10);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
    } else {
	setmsg_("No start time was provided in the setup file.", (ftnlen)45);
	sigerr_("SPICE(NOSTARTTIME)", (ftnlen)18);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     See whether an output DSK stop time was given. */
/*     Capture the value if so. */

    dtpool_("STOP_TIME", &found, &n, vtype, (ftnlen)9, (ftnlen)1);
    if (found) {
	if (*(unsigned char *)vtype != 'C') {
	    setmsg_("DSK stop time was not given a character string value in"
		    " the setup file.", (ftnlen)71);
	    sigerr_("SPICE(TYPEMISMATCH)", (ftnlen)19);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gcpool_("STOP_TIME", &c__1, &c__1, &n, timstr, &found, (ftnlen)9, (
		ftnlen)50);
	if (found) {
	    str2et_(timstr, last, (ftnlen)50);
	    if (failed_()) {
		chkout_("GETGEN", (ftnlen)6);
		return 0;
	    }
	} else {
	    setmsg_("DTPOOL says a stop time was provided in the setup file,"
		    " but GCPOOL can't find it (?)", (ftnlen)84);
	    sigerr_("SPICE(BUG)", (ftnlen)10);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
    } else {
	setmsg_("No stop time was provided in the setup file.", (ftnlen)44);
	sigerr_("SPICE(NOSTOPTIME)", (ftnlen)17);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Get segment data class. */

    gipool_("DATA_CLASS", &c__1, &c__1, &n, dclass, &found, (ftnlen)10);
    if (! found) {
	setmsg_("No segment data class was provided in the setup file.", (
		ftnlen)53);
	sigerr_("SPICE(MISSINGDATACLASS)", (ftnlen)23);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Get segment data type. */

    gipool_("DATA_TYPE", &c__1, &c__1, &n, dtype, &found, (ftnlen)9);
    if (! found) {
	setmsg_("No segment data type was provided in the setup file.", (
		ftnlen)52);
	sigerr_("SPICE(MISSINGDATATYPE)", (ftnlen)22);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     See whether a units string was provided. */

    dtpool_("INPUT_DATA_UNITS", &found, &n, vtype, (ftnlen)16, (ftnlen)1);
    if (found) {
	if (*(unsigned char *)vtype != 'C') {
	    setmsg_("Units specification was not given a character string va"
		    "lue in the setup file.", (ftnlen)77);
	    sigerr_("SPICE(TYPEMISMATCH)", (ftnlen)19);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gcpool_("INPUT_DATA_UNITS", &c__1, &c__2, &n, unistr, &found, (ftnlen)
		16, (ftnlen)255);
	if (found) {
	    if (n != 2) {

/*              We need both distance and angular units. */

		setmsg_("Improperly formatted unit specification in setup fi"
			"le: number of strings found was #. Both angular unit"
			"s and distance units must be specified.", (ftnlen)142)
			;
		errint_("#", &n, (ftnlen)1);
		sigerr_("SPICE(SYNTAXERROR)", (ftnlen)18);
		chkout_("GETGEN", (ftnlen)6);
		return 0;
	    }

/*           Parse the unit specifications. */

	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		lparse_(unistr + ((i__2 = i__ - 1) < 2 && 0 <= i__2 ? i__2 : 
			s_rnge("unistr", i__2, "prcset_", (ftnlen)1015)) * 
			255, "=", &c__2, &ntoken, words, (ftnlen)255, (ftnlen)
			1, (ftnlen)255);
		if (ntoken != 2) {
		    setmsg_("Improperly formatted unit specification in setu"
			    "p file: #", (ftnlen)56);
		    errch_("#", unistr + ((i__2 = i__ - 1) < 2 && 0 <= i__2 ? 
			    i__2 : s_rnge("unistr", i__2, "prcset_", (ftnlen)
			    1021)) * 255, (ftnlen)1, (ftnlen)255);
		    sigerr_("SPICE(SYNTAXERROR)", (ftnlen)18);
		    chkout_("GETGEN", (ftnlen)6);
		    return 0;
		}
		if (eqstr_(words, "ANGLES", (ftnlen)255, (ftnlen)6)) {
		    ljust_(words + 255, aunits, (ftnlen)255, aunits_len);
		    ucase_(aunits, aunits, aunits_len, aunits_len);
		} else if (eqstr_(words, "DISTANCES", (ftnlen)255, (ftnlen)9))
			 {
		    ljust_(words + 255, dunits, (ftnlen)255, dunits_len);
		    ucase_(dunits, dunits, dunits_len, dunits_len);

/*                 Map "KILOMETERS" to "KM"; the latter is */
/*                 recognized by CONVRT. */

		    if (s_cmp(dunits, "KILOMETERS", dunits_len, (ftnlen)10) ==
			     0) {
			s_copy(dunits, "KM", dunits_len, (ftnlen)2);
		    }
		} else {
		    setmsg_("Unrecognized dimension # in unit specification "
			    "in setup file: #", (ftnlen)63);
		    errch_("#", words, (ftnlen)1, (ftnlen)255);
		    errch_("#", unistr + ((i__2 = i__ - 1) < 2 && 0 <= i__2 ? 
			    i__2 : s_rnge("unistr", i__2, "prcset_", (ftnlen)
			    1050)) * 255, (ftnlen)1, (ftnlen)255);
		    sigerr_("SPICE(SYNTAXERROR)", (ftnlen)18);
		    chkout_("GETGEN", (ftnlen)6);
		    return 0;
		}
	    }
	} else {
	    setmsg_("DTPOOL says a units assignment was provided in the setu"
		    "p file, but GCPOOL can't find it (?)", (ftnlen)91);
	    sigerr_("SPICE(BUG)", (ftnlen)10);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
    } else {
	setmsg_("No unit specification was provided in the setup file.", (
		ftnlen)53);
	sigerr_("SPICE(NOUNITSPEC)", (ftnlen)17);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Get the reference frame name. */

    gcpool_("REF_FRAME_NAME", &c__1, &c__1, &n, frame, &found, (ftnlen)14, 
	    frame_len);
    if (! found) {
	setmsg_("No reference frame name was provided in the setup file.", (
		ftnlen)55);
	sigerr_("SPICE(MISSINGFRAME)", (ftnlen)19);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Verify that the frame can be mapped to an ID code. */

    namfrm_(frame, &framid, frame_len);
    if (framid == 0) {
	setmsg_("Reference frame name # could not be mapped to a frame ID co"
		"de.", (ftnlen)62);
	errch_("#", frame, (ftnlen)1, frame_len);
	sigerr_("SPICE(NOTRANSLATION)", (ftnlen)20);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Get the coordinate system name. */

    gcpool_("COORDINATE_SYSTEM", &c__1, &c__1, &n, sysnam, &found, (ftnlen)17,
	     (ftnlen)255);
    if (! found) {
	setmsg_("No coordinate system name was provided in the setup file.", (
		ftnlen)57);
	sigerr_("SPICE(MISSINGCOORDSYS)", (ftnlen)22);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Map the coordinate system name to an ID. */

    *corsys = esrchc_(sysnam, &c__4, csynms, (ftnlen)255, (ftnlen)255);
    if (*corsys == 0) {
	setmsg_("Coordinate system name # was not recognized. ", (ftnlen)45);
	errch_("#", sysnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    } else if (*corsys == 2) {
	setmsg_("Cylindrical coordinates are not supported by this program.", 
		(ftnlen)58);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("GETGEN", (ftnlen)6);
	return 0;
    }

/*     Get coordinate bounds. */

    if (*corsys == 1 || *corsys == 4) {
	gdpool_("MINIMUM_LONGITUDE", &c__1, &c__1, &n, mncor1, &found, (
		ftnlen)17);
	if (! found) {
	    setmsg_("No minimum longitude was provided in the setup file.", (
		    ftnlen)52);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gdpool_("MAXIMUM_LONGITUDE", &c__1, &c__1, &n, mxcor1, &found, (
		ftnlen)17);
	if (! found) {
	    setmsg_("No maximum longitude was provided in the setup file.", (
		    ftnlen)52);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gdpool_("MINIMUM_LATITUDE", &c__1, &c__1, &n, mncor2, &found, (ftnlen)
		16);
	if (! found) {
	    setmsg_("No minimum latitude was provided in the setup file.", (
		    ftnlen)51);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gdpool_("MAXIMUM_LATITUDE", &c__1, &c__1, &n, mxcor2, &found, (ftnlen)
		16);
	if (! found) {
	    setmsg_("No maximum latitude was provided in the setup file.", (
		    ftnlen)51);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
    }
    if (*corsys == 3) {
	gdpool_("MINIMUM_X", &c__1, &c__1, &n, mncor1, &found, (ftnlen)9);
	if (! found) {
	    setmsg_("No minimum X-coordinate was provided in the setup file.",
		     (ftnlen)55);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gdpool_("MAXIMUM_X", &c__1, &c__1, &n, mxcor1, &found, (ftnlen)9);
	if (! found) {
	    setmsg_("No maximum X-coordinate was provided in the setup file.",
		     (ftnlen)55);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gdpool_("MINIMUM_Y", &c__1, &c__1, &n, mncor2, &found, (ftnlen)9);
	if (! found) {
	    setmsg_("No minimum Y-coordinate was provided in the setup file.",
		     (ftnlen)55);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gdpool_("MAXIMUM_Y", &c__1, &c__1, &n, mxcor2, &found, (ftnlen)9);
	if (! found) {
	    setmsg_("No maximum Y-coordinate was provided in the setup file.",
		     (ftnlen)55);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
    }

/*     Get coordinate parameters, if necessary. */

    cleard_(&c__10, corpar);
    if (*corsys == 4) {
	gdpool_("EQUATORIAL_RADIUS", &c__1, &c__1, &n, &re, &found, (ftnlen)
		17);
	if (! found) {
	    setmsg_("No equatorial radius for the planetodetic coordinate sy"
		    "stem's reference spheroid was provided in the setup file."
		    , (ftnlen)112);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	gdpool_("POLAR_RADIUS", &c__1, &c__1, &n, &rp, &found, (ftnlen)12);
	if (! found) {
	    setmsg_("No polar radius for the planetodetic coordinate system'"
		    "s reference spheroid was provided in the setup file.", (
		    ftnlen)107);
	    sigerr_("SPICE(MISSINGCOORDBOUND)", (ftnlen)24);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
	if (re <= 0. || rp <= 0.) {
	    setmsg_("In the setup file, the equatorial radius = #; the polar"
		    " radius = #. Both radii must be positive.", (ftnlen)96);
	    errdp_("#", &re, (ftnlen)1);
	    errdp_("#", &rp, (ftnlen)1);
	    sigerr_("SPICE(INVALIDRADII)", (ftnlen)19);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}

/*        Map the radii to a flattening coefficient. */

	corpar[0] = re;
	f = (re - rp) / re;
	corpar[1] = f;
    }

/*     Fetch the "make vertex-plate mapping" flag, if it's present. */

    gcpool_("MAKE_VERTEX_PLATE_MAP", &c__1, &c__1, &n, cval, &found, (ftnlen)
	    21, (ftnlen)255);
    if (! found) {

/*        The flag is not required to be present. By default, */
/*        no map is created. */

	*makvpm = FALSE_;
    } else {
	if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
	    *makvpm = TRUE_;
	} else if (eqstr_(cval, "NO", (ftnlen)255, (ftnlen)2)) {
	    *makvpm = FALSE_;
	} else {
	    setmsg_("\"Make vertex-plate map\" flag must be YES or NO but wa"
		    "s #.", (ftnlen)57);
	    errch_("#", cval, (ftnlen)1, (ftnlen)255);
	    sigerr_("SPICE(INVALIDFLAG)", (ftnlen)18);
	    chkout_("GETGEN", (ftnlen)6);
	    return 0;
	}
    }
    chkout_("GETGEN", (ftnlen)6);
    return 0;

/*     Get type 2 parameters. */


L_getp02:
    chkin_("GETP02", (ftnlen)6);

/*     Get the plate model type. */

    gipool_("PLATE_TYPE", &c__1, &c__1, &n, pltype, &found, (ftnlen)10);
    if (! found) {
	setmsg_("No plate model type was provided in the setup file.", (
		ftnlen)51);
	sigerr_("SPICE(MISSINGPLATETYPE)", (ftnlen)23);
	chkout_("GETP02", (ftnlen)6);
	return 0;
    }

/*     Get the voxel scale. */

    gdpool_("FINE_VOXEL_SCALE", &c__1, &c__1, &n, voxscl, &fsfnd, (ftnlen)16);

/*     Get the coarse voxel scale. */

    gipool_("COARSE_VOXEL_SCALE", &c__1, &c__1, &n, cgrscl, &csfnd, (ftnlen)
	    18);

/*     It's ok if no scales were provided; otherwise both */
/*     must be provided. */

    if (! fsfnd && ! csfnd) {

/*        Return scales set to zero. The scales will be */
/*        determined automatically. */

	*voxscl = 0.;
	*cgrscl = 0;
    } else if (fsfnd && csfnd) {

/*        Both scales were provided; check them. */

	if (*voxscl <= 0.) {
	    setmsg_("Fine voxel scale must be strictly positive but was #. ("
		    "This scale normally should be greater than or equal to 1"
		    ".0.)", (ftnlen)115);
	    errdp_("#", voxscl, (ftnlen)1);
	    sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	    chkout_("GETP02", (ftnlen)6);
	    return 0;
	}
	if (*cgrscl < 1) {
	    setmsg_("Coarse voxel scale must be greater than or equal to 1, "
		    "but was #.", (ftnlen)65);
	    errint_("#", cgrscl, (ftnlen)1);
	    sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	    chkout_("GETP02", (ftnlen)6);
	    return 0;
	} else if (*cgrscl > 100) {
	    setmsg_("Coarse voxel scale must be less than or equal to #, but"
		    " was #. (Normally this scale should not exceed 20.)", (
		    ftnlen)106);
	    errint_("#", &c__100, (ftnlen)1);
	    errint_("#", cgrscl, (ftnlen)1);
	    sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	    chkout_("GETP02", (ftnlen)6);
	    return 0;
	}
    } else if (! fsfnd) {
	setmsg_("No fine voxel scale was provided in the setup file, but a c"
		"oarse voxel scale was provided. Either add an assignment for"
		" the fine voxel scale, or provide neither scale, in which ca"
		"se the scales will be set automatically.", (ftnlen)219);
	sigerr_("SPICE(MISSINGVOXELSCALE)", (ftnlen)24);
	chkout_("GETP02", (ftnlen)6);
	return 0;
    } else if (! csfnd) {
	setmsg_("No coarse voxel scale was provided in the setup file, but a"
		" fine voxel scale was provided. Either add an assignment for"
		" the coarse voxel scale, or provide neither scale, in which "
		"case the scales will be set automatically.", (ftnlen)221);
	sigerr_("SPICE(MISSINGVOXELSCALE)", (ftnlen)24);
	chkout_("GETP02", (ftnlen)6);
	return 0;
    }
    chkout_("GETP02", (ftnlen)6);
    return 0;

/*     Get plate type 5 (height grid) parameters. */


L_getg05:
    if (return_()) {
	return 0;
    }
    chkin_("GETG05", (ftnlen)6);

/*     Initialize flags that are set conditionally. */

    *wrap = FALSE_;
    *mkncap = FALSE_;
    *mkscap = FALSE_;

/*     Initialize REFVAL. */

    *refval = 0.;
    if (*corsys == 3) {
	gcpool_("WRAP_LONGITUDE", &c__1, &c__1, &n, cval, &found, (ftnlen)14, 
		(ftnlen)255);
	if (found) {
	    if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
		setmsg_("The longitude wrap flag does not apply to rectangul"
			"ar coordinates. Set the flag value to 'NO' or delete"
			" the longitude wrap assignment from the setup file.", 
			(ftnlen)154);
		sigerr_("SPICE(SPURIOUSKEYWORD)", (ftnlen)22);
		chkout_("GETG05", (ftnlen)6);
		return 0;
	    }
	}
	gcpool_("MAKE_NORTH_POLAR_CAP", &c__1, &c__1, &n, cval, &found, (
		ftnlen)20, (ftnlen)255);
	if (found) {
	    if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
		setmsg_("The north polar cap flag does not apply to rectangu"
			"lar coordinates. Set the flag value to 'NO' or delet"
			"e the north polar cap flag assignment from the setup"
			" file.", (ftnlen)161);
		sigerr_("SPICE(SPURIOUSKEYWORD)", (ftnlen)22);
		chkout_("GETG05", (ftnlen)6);
		return 0;
	    }
	}
	gcpool_("MAKE_SOUTH_POLAR_CAP", &c__1, &c__1, &n, cval, &found, (
		ftnlen)20, (ftnlen)255);
	if (found) {
	    if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
		setmsg_("The south polar cap flag does not apply to rectangu"
			"lar coordinates.  Set the flag value to 'NO' or dele"
			"te the south polar cap flag assignment from the setu"
			"p file.", (ftnlen)162);
		sigerr_("SPICE(SPURIOUSKEYWORD)", (ftnlen)22);
		chkout_("GETG05", (ftnlen)6);
		return 0;
	    }
	}

/*        Get the coordinate of the top row. */

	gdpool_("TOP_ROW_Y_COORDINATE", &c__1, &c__1, &n, topcor, &found, (
		ftnlen)20);
	if (! found) {
	    setmsg_("No Y-coordinate of the top row was  provided in the set"
		    "up file.", (ftnlen)63);
	    sigerr_("SPICE(MISSINGTOPCOR)", (ftnlen)20);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}

/*        Get the coordinate of the left column. */

	gdpool_("LEFT_COLUMN_X_COORDINATE", &c__1, &c__1, &n, lftcor, &found, 
		(ftnlen)24);
	if (! found) {
	    setmsg_("No X-coordinate of the left column was  provided in the"
		    " setup file.", (ftnlen)67);
	    sigerr_("SPICE(MISSINGLEFTCOR)", (ftnlen)21);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}
    } else {

/*        This is a lon/lat coordinate system. */

/*        Get the longitude wrap flag. */

	gcpool_("WRAP_LONGITUDE", &c__1, &c__1, &n, cval, &found, (ftnlen)14, 
		(ftnlen)255);
	if (! found) {
	    setmsg_("No longitude wrap flag was provided in the setup file.", 
		    (ftnlen)54);
	    sigerr_("SPICE(MISSINGWRAPFLAG)", (ftnlen)22);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}
	if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
	    *wrap = TRUE_;
	} else if (eqstr_(cval, "NO", (ftnlen)255, (ftnlen)2)) {
	    *wrap = FALSE_;
	} else {
	    setmsg_("Longitude wrap flag must be YES or NO but was #.", (
		    ftnlen)48);
	    errch_("#", cval, (ftnlen)1, (ftnlen)255);
	    sigerr_("SPICE(INVALIDFLAG)", (ftnlen)18);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}

/*        Get the north polar cap flag. */

	gcpool_("MAKE_NORTH_POLAR_CAP", &c__1, &c__1, &n, cval, &found, (
		ftnlen)20, (ftnlen)255);
	if (! found) {
	    setmsg_("No north polar cap flag was provided in the setup file.",
		     (ftnlen)55);
	    sigerr_("SPICE(MISSINGNCAPFLAG)", (ftnlen)22);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}
	if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
	    *mkncap = TRUE_;
	} else if (eqstr_(cval, "NO", (ftnlen)255, (ftnlen)2)) {
	    *mkncap = FALSE_;
	} else {
	    setmsg_("North polar cap flag must be YES or NO but was #.", (
		    ftnlen)49);
	    errch_("#", cval, (ftnlen)1, (ftnlen)255);
	    sigerr_("SPICE(INVALIDFLAG)", (ftnlen)18);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}

/*        Get the south polar cap flag. */

	gcpool_("MAKE_SOUTH_POLAR_CAP", &c__1, &c__1, &n, cval, &found, (
		ftnlen)20, (ftnlen)255);
	if (! found) {
	    setmsg_("No south polar cap flag was provided in the setup file.",
		     (ftnlen)55);
	    sigerr_("SPICE(MISSINGSCAPFLAG)", (ftnlen)22);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}
	if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
	    *mkscap = TRUE_;
	} else if (eqstr_(cval, "NO", (ftnlen)255, (ftnlen)2)) {
	    *mkscap = FALSE_;
	} else {
	    setmsg_("South polar cap flag must be YES or NO but was #.", (
		    ftnlen)49);
	    errch_("#", cval, (ftnlen)1, (ftnlen)255);
	    sigerr_("SPICE(INVALIDFLAG)", (ftnlen)18);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}

/*        Get the coordinate of the top row. */

	gdpool_("TOP_ROW_LATITUDE", &c__1, &c__1, &n, topcor, &found, (ftnlen)
		16);
	if (! found) {
	    setmsg_("No latitude of the top row was  provided in the setup f"
		    "ile.", (ftnlen)59);
	    sigerr_("SPICE(MISSINGTOPCOR)", (ftnlen)20);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}

/*        Get the coordinate of the left column. */

	gdpool_("LEFT_COLUMN_LONGITUDE", &c__1, &c__1, &n, lftcor, &found, (
		ftnlen)21);
	if (! found) {
	    setmsg_("No longitude of the left column was  provided in the se"
		    "tup file.", (ftnlen)64);
	    sigerr_("SPICE(MISSINGLEFTCOR)", (ftnlen)21);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}
    }

/*     Get the row major grid order flag. */

    gcpool_("INPUT_GRID_ORDER_ROW_MAJOR", &c__1, &c__1, &n, cval, &found, (
	    ftnlen)26, (ftnlen)255);
    if (! found) {
	setmsg_("No row major flag was provided in the setup file.", (ftnlen)
		49);
	sigerr_("SPICE(MISSINGROWMAJFLAG)", (ftnlen)24);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }
    if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
	*rowmaj = TRUE_;
    } else if (eqstr_(cval, "NO", (ftnlen)255, (ftnlen)2)) {
	*rowmaj = FALSE_;
    } else {
	setmsg_("Row major flag must be YES or NO but was #.", (ftnlen)43);
	errch_("#", cval, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(INVALIDFLAG)", (ftnlen)18);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }

/*     Get the top-down grid order flag. */

    gcpool_("COLUMN_VALUE_ORDER_TOP_DOWN", &c__1, &c__1, &n, cval, &found, (
	    ftnlen)27, (ftnlen)255);
    if (! found) {
	setmsg_("No top-down flag was provided in the setup file.", (ftnlen)
		48);
	sigerr_("SPICE(MISSINGTOPDOWNFLAG)", (ftnlen)25);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }
    if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
	*topdwn = TRUE_;
    } else if (eqstr_(cval, "NO", (ftnlen)255, (ftnlen)2)) {
	*topdwn = FALSE_;
    } else {
	setmsg_("Top-down flag must be YES or NO but was #.", (ftnlen)42);
	errch_("#", cval, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(INVALIDFLAG)", (ftnlen)18);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }

/*     Get the left-right grid order flag. */

    gcpool_("ROW_VALUE_ORDER_LEFT_RIGHT", &c__1, &c__1, &n, cval, &found, (
	    ftnlen)26, (ftnlen)255);
    if (! found) {
	setmsg_("No left-right flag was provided in the setup file.", (ftnlen)
		50);
	sigerr_("SPICE(MISSINGLEFTRTFLAG)", (ftnlen)24);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }
    if (eqstr_(cval, "YES", (ftnlen)255, (ftnlen)3)) {
	*leftrt = TRUE_;
    } else if (eqstr_(cval, "NO", (ftnlen)255, (ftnlen)2)) {
	*leftrt = FALSE_;
    } else {
	setmsg_("Left-right flag must be YES or NO but was #.", (ftnlen)44);
	errch_("#", cval, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(INVALIDFLAG)", (ftnlen)18);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }
    if (*corsys != 4) {

/*        Get the height reference value. */

	gdpool_("HEIGHT_REFERENCE", &c__1, &c__1, &n, refval, &found, (ftnlen)
		16);
	if (! found) {
	    setmsg_("No height reference was provided in the setup file.", (
		    ftnlen)51);
	    sigerr_("SPICE(MISSINGHEIGHTREF)", (ftnlen)23);
	    chkout_("GETG05", (ftnlen)6);
	    return 0;
	}
    }

/*     Get the height scale value. */

    gdpool_("HEIGHT_SCALE", &c__1, &c__1, &n, hscale, &found, (ftnlen)12);
    if (! found) {
	setmsg_("No height scale factor was provided in the setup file.", (
		ftnlen)54);
	sigerr_("SPICE(MISSINGHSCALE)", (ftnlen)20);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }

/*     Get the column count. */

    gipool_("NUMBER_OF_COLUMNS", &c__1, &c__1, &n, ncols, &found, (ftnlen)17);
    if (! found) {
	setmsg_("No column count was provided in the setup file.", (ftnlen)47)
		;
	sigerr_("SPICE(MISSINGNCOLS)", (ftnlen)19);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }

/*     Get the row count. */

    gipool_("NUMBER_OF_ROWS", &c__1, &c__1, &n, nrows, &found, (ftnlen)14);
    if (! found) {
	setmsg_("No row count was provided in the setup file.", (ftnlen)44);
	sigerr_("SPICE(MISSINGNROWS)", (ftnlen)19);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }

/*     Get the column step size. */

    gdpool_("COLUMN_STEP_SIZE", &c__1, &c__1, &n, colstp, &found, (ftnlen)16);
    if (! found) {
	setmsg_("No column step size was  provided in the setup file.", (
		ftnlen)52);
	sigerr_("SPICE(MISSINGCOLSTEP)", (ftnlen)21);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }

/*     Get the row step size. */

    gdpool_("ROW_STEP_SIZE", &c__1, &c__1, &n, rowstp, &found, (ftnlen)13);
    if (! found) {
	setmsg_("No row step size was  provided in the setup file.", (ftnlen)
		49);
	sigerr_("SPICE(MISSINGROWSTEP)", (ftnlen)21);
	chkout_("GETG05", (ftnlen)6);
	return 0;
    }
    chkout_("GETG05", (ftnlen)6);
    return 0;
} /* prcset_ */

/* Subroutine */ int prcset_(char *setup, char *input, char *output, char *
	cmtfil, integer *surfid, integer *centid, char *frame, doublereal *
	first, doublereal *last, integer *dclass, integer *dtype, char *
	aunits, char *dunits, integer *corsys, doublereal *corpar, doublereal 
	*mncor1, doublereal *mxcor1, doublereal *mncor2, doublereal *mxcor2, 
	integer *pltype, doublereal *voxscl, integer *cgrscl, logical *wrap, 
	logical *mkncap, logical *mkscap, logical *rowmaj, logical *topdwn, 
	logical *leftrt, doublereal *refval, doublereal *hscale, integer *
	ncols, integer *nrows, doublereal *lftcor, doublereal *topcor, 
	doublereal *colstp, doublereal *rowstp, logical *makvpm, ftnlen 
	setup_len, ftnlen input_len, ftnlen output_len, ftnlen cmtfil_len, 
	ftnlen frame_len, ftnlen aunits_len, ftnlen dunits_len)
{
    return prcset_0_(0, setup, input, output, cmtfil, surfid, centid, frame, 
	    first, last, dclass, dtype, aunits, dunits, corsys, corpar, 
	    mncor1, mxcor1, mncor2, mxcor2, pltype, voxscl, cgrscl, wrap, 
	    mkncap, mkscap, rowmaj, topdwn, leftrt, refval, hscale, ncols, 
	    nrows, lftcor, topcor, colstp, rowstp, makvpm, setup_len, 
	    input_len, output_len, cmtfil_len, frame_len, aunits_len, 
	    dunits_len);
    }

/* Subroutine */ int getset_(char *setup, char *input, char *output, char *
	cmtfil, ftnlen setup_len, ftnlen input_len, ftnlen output_len, ftnlen 
	cmtfil_len)
{
    return prcset_0_(1, setup, input, output, cmtfil, (integer *)0, (integer *
	    )0, (char *)0, (doublereal *)0, (doublereal *)0, (integer *)0, (
	    integer *)0, (char *)0, (char *)0, (integer *)0, (doublereal *)0, 
	    (doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)
	    0, (integer *)0, (doublereal *)0, (integer *)0, (logical *)0, (
	    logical *)0, (logical *)0, (logical *)0, (logical *)0, (logical *)
	    0, (doublereal *)0, (doublereal *)0, (integer *)0, (integer *)0, (
	    doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0,
	     (logical *)0, setup_len, input_len, output_len, cmtfil_len, (
	    ftnint)0, (ftnint)0, (ftnint)0);
    }

/* Subroutine */ int gettyp_(integer *dtype)
{
    return prcset_0_(2, (char *)0, (char *)0, (char *)0, (char *)0, (integer *
	    )0, (integer *)0, (char *)0, (doublereal *)0, (doublereal *)0, (
	    integer *)0, dtype, (char *)0, (char *)0, (integer *)0, (
	    doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0,
	     (doublereal *)0, (integer *)0, (doublereal *)0, (integer *)0, (
	    logical *)0, (logical *)0, (logical *)0, (logical *)0, (logical *)
	    0, (logical *)0, (doublereal *)0, (doublereal *)0, (integer *)0, (
	    integer *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0, (
	    doublereal *)0, (logical *)0, (ftnint)0, (ftnint)0, (ftnint)0, (
	    ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0);
    }

/* Subroutine */ int getgen_(integer *surfid, integer *centid, char *frame, 
	doublereal *first, doublereal *last, integer *dclass, integer *dtype, 
	char *aunits, char *dunits, integer *corsys, doublereal *corpar, 
	doublereal *mncor1, doublereal *mxcor1, doublereal *mncor2, 
	doublereal *mxcor2, logical *makvpm, ftnlen frame_len, ftnlen 
	aunits_len, ftnlen dunits_len)
{
    return prcset_0_(3, (char *)0, (char *)0, (char *)0, (char *)0, surfid, 
	    centid, frame, first, last, dclass, dtype, aunits, dunits, corsys,
	     corpar, mncor1, mxcor1, mncor2, mxcor2, (integer *)0, (
	    doublereal *)0, (integer *)0, (logical *)0, (logical *)0, (
	    logical *)0, (logical *)0, (logical *)0, (logical *)0, (
	    doublereal *)0, (doublereal *)0, (integer *)0, (integer *)0, (
	    doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0,
	     makvpm, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, frame_len, 
	    aunits_len, dunits_len);
    }

/* Subroutine */ int getp02_(integer *pltype, doublereal *voxscl, integer *
	cgrscl)
{
    return prcset_0_(4, (char *)0, (char *)0, (char *)0, (char *)0, (integer *
	    )0, (integer *)0, (char *)0, (doublereal *)0, (doublereal *)0, (
	    integer *)0, (integer *)0, (char *)0, (char *)0, (integer *)0, (
	    doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0,
	     (doublereal *)0, pltype, voxscl, cgrscl, (logical *)0, (logical *
	    )0, (logical *)0, (logical *)0, (logical *)0, (logical *)0, (
	    doublereal *)0, (doublereal *)0, (integer *)0, (integer *)0, (
	    doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0,
	     (logical *)0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (
	    ftnint)0, (ftnint)0, (ftnint)0);
    }

/* Subroutine */ int getg05_(integer *corsys, logical *wrap, logical *mkncap, 
	logical *mkscap, logical *rowmaj, logical *topdwn, logical *leftrt, 
	doublereal *refval, doublereal *hscale, integer *ncols, integer *
	nrows, doublereal *lftcor, doublereal *topcor, doublereal *colstp, 
	doublereal *rowstp)
{
    return prcset_0_(5, (char *)0, (char *)0, (char *)0, (char *)0, (integer *
	    )0, (integer *)0, (char *)0, (doublereal *)0, (doublereal *)0, (
	    integer *)0, (integer *)0, (char *)0, (char *)0, corsys, (
	    doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0,
	     (doublereal *)0, (integer *)0, (doublereal *)0, (integer *)0, 
	    wrap, mkncap, mkscap, rowmaj, topdwn, leftrt, refval, hscale, 
	    ncols, nrows, lftcor, topcor, colstp, rowstp, (logical *)0, (
	    ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (
	    ftnint)0);
    }

