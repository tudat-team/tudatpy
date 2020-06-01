/* prcinf.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/* $Procedure PRCINF ( Process an information request ) */
/* Subroutine */ int prcinf_(char *inftyp, ftnlen inftyp_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    char ch__1[133];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rnge(char *, integer, char *, integer);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen), errch_(char *, char *,
	     ftnlen, ftnlen);
    extern integer rtrim_(char *, ftnlen);
    extern logical eqstr_(char *, char *, ftnlen, ftnlen);
    char begmrk[80], endmrk[80];
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), tostdo_(char *, ftnlen);
    extern logical return_(void);
    static char hlptxt[80*18], verstr[80];
    extern /* Subroutine */ int tkvrsn_(char *, char *, ftnlen, ftnlen);
    static char usgtxt[80*17], tmptxt[80*132];

/* $ Abstract */

/*     Process an information request:  display "help," "usage," */
/*     "template, or program version information. */

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

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     INFTYP     I   Type of information to display. */

/* $ Detailed_Input */

/*     INFTYP         is a character string indicating the type */
/*                    of information to display.  The options are: */

/*                       'HELP'        Dump the introductory */
/*                                     paragraph of the user's guide. */

/*                       'TEMPLATE'    Display a setup file template. */

/*                       'USAGE'       Display a terse description */
/*                                     of the program's invocation */
/*                                     syntax. */

/*                       'VERSION'     Display the program version */
/*                                     and creation date. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the value of INFTYP is not recognized, the error */
/*        SPICE(NOTSUPPORTED) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine centralizes message display functions for */
/*     MKDSK. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     1) For use only within program MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 4.0.0, 04-APR-2017 (NJB) */

/*        Updated setup file template. Moved declaration of */
/*        version string to the include file mkdsk.inc. */

/* -    MKDSK Version 3.0.0, 30-JUN-2014 (NJB) */

/*        Updated version string. */

/* -    MKDSK Version 2.0.0, 29-JUN-2010 (NJB) */

/*        Updated template to reflect keyword changes. */

/* -    MKDSK Version 1.0.0, 08-JUN-2010 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Initial values */

    if (return_()) {
	return 0;
    }
    chkin_("PRCINF", (ftnlen)6);
    if (first) {

/*        This lovely mess was created using Bill Taber's ImportText */
/*        pipe. */

	*(unsigned char *)begmrk = '\\';
	s_copy(begmrk + 1, "begindata", (ftnlen)79, (ftnlen)9);
	*(unsigned char *)endmrk = '\\';
	s_copy(endmrk + 1, "begintext", (ftnlen)79, (ftnlen)9);
	s_copy(hlptxt, "   MKDSK is a SPICE Toolkit utility program that con"
		"verts a shape", (ftnlen)80, (ftnlen)65);
	s_copy(hlptxt + 80, "   data file having a recognized format to a SP"
		"ICE DSK (\"Digital", (ftnlen)80, (ftnlen)64);
	s_copy(hlptxt + 160, "   Shape Kernel\") file.", (ftnlen)80, (ftnlen)
		23);
	s_copy(hlptxt + 240, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 320, "   MKDSK requires as inputs a shape data file,"
		" a leapseconds file, and", (ftnlen)80, (ftnlen)70);
	s_copy(hlptxt + 400, "   a setup file containing commands that contr"
		"ol MKDSK's operation.", (ftnlen)80, (ftnlen)67);
	s_copy(hlptxt + 480, "   Execute MKDSK with the -t option to see a s"
		"etup file template:", (ftnlen)80, (ftnlen)65);
	s_copy(hlptxt + 560, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 640, "       % mkdsk -t", (ftnlen)80, (ftnlen)17);
	s_copy(hlptxt + 720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 800, "   The user may optionally specify a text file"
		" containing descriptive", (ftnlen)80, (ftnlen)69);
	s_copy(hlptxt + 880, "   text to be placed in the comment area of th"
		"e DSK (doing this is", (ftnlen)80, (ftnlen)66);
	s_copy(hlptxt + 960, "   highly recommended).", (ftnlen)80, (ftnlen)
		23);
	s_copy(hlptxt + 1040, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 1120, "   For documentation purposes the contents of"
		" the MKDSK setup file are", (ftnlen)80, (ftnlen)70);
	s_copy(hlptxt + 1200, "   automatically placed at the end of the com"
		"ment area of the DSK file.", (ftnlen)80, (ftnlen)71);
	s_copy(hlptxt + 1280, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 1360, "   See the MKDSK User's Guide for further inf"
		"ormation.", (ftnlen)80, (ftnlen)54);
	s_copy(tmptxt, "  Complete MKDSK Setup File Template:", (ftnlen)80, (
		ftnlen)37);
	s_copy(tmptxt + 80, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 160, begmrk, (ftnlen)80, (ftnlen)80);
	s_copy(tmptxt + 240, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 320, "   INPUT_SHAPE_FILE       = 'Name of input sha"
		"pe data file'", (ftnlen)80, (ftnlen)59);
	s_copy(tmptxt + 400, "   OUTPUT_DSK_FILE        = 'Name of output DS"
		"K file'", (ftnlen)80, (ftnlen)53);
	s_copy(tmptxt + 480, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 560, "   Optional assignment:", (ftnlen)80, (ftnlen)
		23);
	s_copy(tmptxt + 640, "   COMMENT_FILE           = 'Name of optional "
		"comment file'", (ftnlen)80, (ftnlen)59);
	s_copy(tmptxt + 720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 800, "   Optional assignment:", (ftnlen)80, (ftnlen)
		23);
	s_copy(tmptxt + 880, "   LEAPSECONDS_FILE       = 'Name of leapsecon"
		"ds file'", (ftnlen)80, (ftnlen)54);
	s_copy(tmptxt + 960, "                            A leapseconds kern"
		"el is required;", (ftnlen)80, (ftnlen)61);
	s_copy(tmptxt + 1040, "                            it can be named u"
		"sing the", (ftnlen)80, (ftnlen)53);
	s_copy(tmptxt + 1120, "                            KERNELS_TO_LOAD a"
		"ssignment.", (ftnlen)80, (ftnlen)55);
	s_copy(tmptxt + 1200, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 1280, "   Optional assignment:", (ftnlen)80, (ftnlen)
		23);
	s_copy(tmptxt + 1360, "   KERNELS_TO_LOAD        = ( 'Kernel_1' 'Ker"
		"nel_2' ... )", (ftnlen)80, (ftnlen)57);
	s_copy(tmptxt + 1440, "                            List any addition"
		"al kernels needed.", (ftnlen)80, (ftnlen)63);
	s_copy(tmptxt + 1520, "                            Note that a leaps"
		"econds kernel can be", (ftnlen)80, (ftnlen)65);
	s_copy(tmptxt + 1600, "                            supplied using th"
		"is assignment.", (ftnlen)80, (ftnlen)59);
	s_copy(tmptxt + 1680, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 1760, "   CENTER_NAME            = 'Central body nam"
		"e'", (ftnlen)80, (ftnlen)47);
	s_copy(tmptxt + 1840, "   SURFACE_NAME           = 'Surface name'", (
		ftnlen)80, (ftnlen)42);
	s_copy(tmptxt + 1920, "   REF_FRAME_NAME         = 'Reference frame "
		"name'", (ftnlen)80, (ftnlen)50);
	s_copy(tmptxt + 2000, "   START_TIME             = 'Start time'", (
		ftnlen)80, (ftnlen)40);
	s_copy(tmptxt + 2080, "   STOP_TIME              = 'Stop time'", (
		ftnlen)80, (ftnlen)39);
	s_copy(tmptxt + 2160, "   DATA_CLASS             = 1 for single-valu"
		"ed surface", (ftnlen)80, (ftnlen)55);
	s_copy(tmptxt + 2240, "                              topography (for"
		" latitudinal", (ftnlen)80, (ftnlen)57);
	s_copy(tmptxt + 2320, "                              coordinates, th"
		"is implies each", (ftnlen)80, (ftnlen)60);
	s_copy(tmptxt + 2400, "                              ray emanating f"
		"rom the reference", (ftnlen)80, (ftnlen)62);
	s_copy(tmptxt + 2480, "                              frame's origin "
		"intersects the", (ftnlen)80, (ftnlen)59);
	s_copy(tmptxt + 2560, "                              surface once)", (
		ftnlen)80, (ftnlen)43);
	s_copy(tmptxt + 2640, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 2720, "                              or", (ftnlen)80, 
		(ftnlen)32);
	s_copy(tmptxt + 2800, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 2880, "                            2 for arbitrary t"
		"opography (e.g.", (ftnlen)80, (ftnlen)60);
	s_copy(tmptxt + 2960, "                              that of a dumbb"
		"ell-shaped asteroid)", (ftnlen)80, (ftnlen)65);
	s_copy(tmptxt + 3040, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 3120, "   INPUT_DATA_UNITS       = ( 'ANGLES    = an"
		"gular unit'", (ftnlen)80, (ftnlen)56);
	s_copy(tmptxt + 3200, "                              'DISTANCES = di"
		"stance unit' )", (ftnlen)80, (ftnlen)59);
	s_copy(tmptxt + 3280, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 3360, "   COORDINATE_SYSTEM      = 'LATITUDINAL'  or",
		 (ftnlen)80, (ftnlen)45);
	s_copy(tmptxt + 3440, "                            'RECTANGULAR'  or",
		 (ftnlen)80, (ftnlen)45);
	s_copy(tmptxt + 3520, "                            'PLANETODETIC'", (
		ftnlen)80, (ftnlen)42);
	s_copy(tmptxt + 3600, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 3680, "      For latitudinal coordinates:", (ftnlen)
		80, (ftnlen)34);
	s_copy(tmptxt + 3760, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 3840, "      MINIMUM_LATITUDE    = lower latitude bo"
		"und in selected units", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 3920, "      MAXIMUM_LATITUDE    = upper latitude bo"
		"und in selected units", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 4000, "      MINIMUM_LONGITUDE   = lower longitude b"
		"ound in selected units", (ftnlen)80, (ftnlen)67);
	s_copy(tmptxt + 4080, "      MAXIMUM_LONGITUDE   = upper longitude b"
		"ound in selected units", (ftnlen)80, (ftnlen)67);
	s_copy(tmptxt + 4160, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 4240, "      For rectangular coordinates:", (ftnlen)
		80, (ftnlen)34);
	s_copy(tmptxt + 4320, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 4400, "      MINIMUM_X           = lower X coordinat"
		"e bound in selected units", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 4480, "      MAXIMUM_X           = upper X coordinat"
		"e bound in selected units", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 4560, "      MINIMUM_Y           = lower Y coordinat"
		"e bound in selected units", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 4640, "      MAXIMUM_Y           = upper Y coordinat"
		"e bound in selected units", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 4720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 4800, "      For planetodetic coordinates:", (ftnlen)
		80, (ftnlen)35);
	s_copy(tmptxt + 4880, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 4960, "      MINIMUM_LATITUDE    = lower latitude bo"
		"und in selected units", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 5040, "      MAXIMUM_LATITUDE    = upper latitude bo"
		"und in selected units", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 5120, "      MINIMUM_LONGITUDE   = lower longitude b"
		"ound in selected units", (ftnlen)80, (ftnlen)67);
	s_copy(tmptxt + 5200, "      MAXIMUM_LONGITUDE   = upper longitude b"
		"ound in selected units", (ftnlen)80, (ftnlen)67);
	s_copy(tmptxt + 5280, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 5360, "      EQUATORIAL_RADIUS   = equatorial sphero"
		"id radius in selected units", (ftnlen)80, (ftnlen)72);
	s_copy(tmptxt + 5440, "      POLAR_RADIUS        = polar spheroid ra"
		"dius in selected units", (ftnlen)80, (ftnlen)67);
	s_copy(tmptxt + 5520, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 5600, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 5680, "   DATA_TYPE              = 2 (triangular pla"
		"tes)", (ftnlen)80, (ftnlen)49);
	s_copy(tmptxt + 5760, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 5840, "      For data type 2:", (ftnlen)80, (ftnlen)
		22);
	s_copy(tmptxt + 5920, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 6000, "      PLATE_TYPE          = 1  for plate-vert"
		"ex table", (ftnlen)80, (ftnlen)53);
	s_copy(tmptxt + 6080, "                            2  for Gaskell sh"
		"ape model", (ftnlen)80, (ftnlen)54);
	s_copy(tmptxt + 6160, "                            3  for vertex-fac"
		"et table", (ftnlen)80, (ftnlen)53);
	s_copy(tmptxt + 6240, "                            4  for Rosetta/OS"
		"IRIS \"ver\" table", (ftnlen)80, (ftnlen)61);
	s_copy(tmptxt + 6320, "                            5  for rectangula"
		"r height grid", (ftnlen)80, (ftnlen)58);
	s_copy(tmptxt + 6400, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 6480, "      The following two assignments are optio"
		"nal.", (ftnlen)80, (ftnlen)49);
	s_copy(tmptxt + 6560, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 6640, "      Optional assignments:", (ftnlen)80, (
		ftnlen)27);
	s_copy(tmptxt + 6720, "      FINE_VOXEL_SCALE    = Double precision "
		"value > 0.0", (ftnlen)80, (ftnlen)56);
	s_copy(tmptxt + 6800, "      COARSE_VOXEL_SCALE  = Integer >= 1", (
		ftnlen)80, (ftnlen)40);
	s_copy(tmptxt + 6880, "                            If these assignme"
		"nts are not provided,", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 6960, "                            MKDSK will set th"
		"e voxel scales automatically.", (ftnlen)80, (ftnlen)74);
	s_copy(tmptxt + 7040, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 7120, "      Optional assignment:", (ftnlen)80, (
		ftnlen)26);
	s_copy(tmptxt + 7200, "      MAKE_VERTEX_PLATE_MAP =  'YES' or 'NO'", 
		(ftnlen)80, (ftnlen)44);
	s_copy(tmptxt + 7280, "                               If this assign"
		"ment is not provided,", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 7360, "                               MKDSK will not"
		" create a vertex-plate mapping.", (ftnlen)80, (ftnlen)76);
	s_copy(tmptxt + 7440, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 7520, "         For plate type 5, the following addi"
		"tional assignments", (ftnlen)80, (ftnlen)63);
	s_copy(tmptxt + 7600, "         are required:", (ftnlen)80, (ftnlen)
		22);
	s_copy(tmptxt + 7680, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 7760, "         WRAP_LONGITUDE              = connec"
		"t leftmost column to", (ftnlen)80, (ftnlen)65);
	s_copy(tmptxt + 7840, "                                       rightm"
		"ost column: 'YES' or 'NO'", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 7920, "         MAKE_NORTH_POLAR_CAP        = extend"
		" plate set to north pole:", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 8000, "                                       'YES' "
		"or 'NO'", (ftnlen)80, (ftnlen)52);
	s_copy(tmptxt + 8080, "         MAKE_SOUTH_POLAR_CAP        = extend"
		" plate set to south pole:", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 8160, "                                       'YES' "
		"or 'NO'", (ftnlen)80, (ftnlen)52);
	s_copy(tmptxt + 8240, "         INPUT_GRID_ORDER_ROW_MAJOR  = input "
		"data set is row-major:", (ftnlen)80, (ftnlen)67);
	s_copy(tmptxt + 8320, "                                       'YES' "
		"or 'NO'", (ftnlen)80, (ftnlen)52);
	s_copy(tmptxt + 8400, "         COLUMN_VALUE_ORDER_TOP_DOWN = input "
		"data set is top-down:", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 8480, "                                       'YES' "
		"or 'NO'", (ftnlen)80, (ftnlen)52);
	s_copy(tmptxt + 8560, "         ROW_VALUE_ORDER_LEFT_RIGHT  = input "
		"data set is left-right:", (ftnlen)80, (ftnlen)68);
	s_copy(tmptxt + 8640, "                                       'YES' "
		"or 'NO'", (ftnlen)80, (ftnlen)52);
	s_copy(tmptxt + 8720, "         HEIGHT_SCALE                = value "
		"by which to multiply the", (ftnlen)80, (ftnlen)69);
	s_copy(tmptxt + 8800, "                                       height"
		" data to convert to km", (ftnlen)80, (ftnlen)67);
	s_copy(tmptxt + 8880, "         NUMBER_OF_ROWS              = number"
		" of rows in grid", (ftnlen)80, (ftnlen)61);
	s_copy(tmptxt + 8960, "         NUMBER_OF_COLUMNS           = number"
		" of columns in grid", (ftnlen)80, (ftnlen)64);
	s_copy(tmptxt + 9040, "         COLUMN_STEP_SIZE            = column"
		" separation: longitude or X step", (ftnlen)80, (ftnlen)77);
	s_copy(tmptxt + 9120, "         ROW_STEP_SIZE               = row se"
		"paration: latitude or Y step", (ftnlen)80, (ftnlen)73);
	s_copy(tmptxt + 9200, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 9280, "            For plate type 5, latitudinal or "
		"planetodetic coordinates:", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 9360, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 9440, "            LEFT_COLUMN_LONGITUDE    = longit"
		"ude of leftmost column of grid", (ftnlen)80, (ftnlen)75);
	s_copy(tmptxt + 9520, "            TOP_ROW_LATITUDE         = latitu"
		"de of top row of grid", (ftnlen)80, (ftnlen)66);
	s_copy(tmptxt + 9600, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 9680, "            For plate type 5, rectangular coo"
		"rdinates:", (ftnlen)80, (ftnlen)54);
	s_copy(tmptxt + 9760, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 9840, "            LEFT_COLUMN_X_COORDINATE = X-coor"
		"dinate of leftmost column of grid", (ftnlen)80, (ftnlen)78);
	s_copy(tmptxt + 9920, "            TOP_ROW_Y_COORDINATE     = Y-coor"
		"dinate of top row of grid", (ftnlen)80, (ftnlen)70);
	s_copy(tmptxt + 10000, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 10080, "            For plate type 5, latitudinal or"
		" rectangular coordinates:", (ftnlen)80, (ftnlen)69);
	s_copy(tmptxt + 10160, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 10240, "            HEIGHT_REFERENCE         = value"
		" to add to input heights; units", (ftnlen)80, (ftnlen)75);
	s_copy(tmptxt + 10320, "                                       are g"
		"iven by INPUT_DATA_UNITS", (ftnlen)80, (ftnlen)68);
	s_copy(tmptxt + 10400, " ", (ftnlen)80, (ftnlen)1);
	s_copy(tmptxt + 10480, endmrk, (ftnlen)80, (ftnlen)80);
	s_copy(usgtxt, "     Program usage:", (ftnlen)80, (ftnlen)19);
	s_copy(usgtxt + 80, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 160, "              > mkdsk   [-setup <setup file na"
		"me>]", (ftnlen)80, (ftnlen)50);
	s_copy(usgtxt + 240, "                        [-input <input shape d"
		"ata file name>]", (ftnlen)80, (ftnlen)61);
	s_copy(usgtxt + 320, "                        [-output <output DSK f"
		"ile name>]", (ftnlen)80, (ftnlen)56);
	s_copy(usgtxt + 400, "                        [-h|-help]", (ftnlen)80,
		 (ftnlen)34);
	s_copy(usgtxt + 480, "                        [-t|-template]", (
		ftnlen)80, (ftnlen)38);
	s_copy(usgtxt + 560, "                        [-u|-usage]", (ftnlen)
		80, (ftnlen)35);
	s_copy(usgtxt + 640, "                        [-v|-version]", (ftnlen)
		80, (ftnlen)37);
	s_copy(usgtxt + 720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 800, "     If a setup file name isn't provided on th"
		"e command line, the", (ftnlen)80, (ftnlen)65);
	s_copy(usgtxt + 880, "     program will prompt for it. It will not p"
		"rompt for the input", (ftnlen)80, (ftnlen)65);
	s_copy(usgtxt + 960, "     or output file names; these file names mu"
		"st be provided on", (ftnlen)80, (ftnlen)63);
	s_copy(usgtxt + 1040, "     the command line or in the setup file. I"
		"f input and output", (ftnlen)80, (ftnlen)63);
	s_copy(usgtxt + 1120, "     file names are provided on the command l"
		"ine, any file names", (ftnlen)80, (ftnlen)64);
	s_copy(usgtxt + 1200, "     assigned using setup keywords are ignore"
		"d. The input file", (ftnlen)80, (ftnlen)62);
	s_copy(usgtxt + 1280, "     must already exist and the output file m"
		"ust be a new file.", (ftnlen)80, (ftnlen)63);
	first = FALSE_;
    }
    if (eqstr_(inftyp, "TEMPLATE", inftyp_len, (ftnlen)8)) {

/*        Display the template text. */

	for (i__ = 1; i__ <= 132; ++i__) {
	    tostdo_(tmptxt + ((i__1 = i__ - 1) < 132 && 0 <= i__1 ? i__1 : 
		    s_rnge("tmptxt", i__1, "prcinf_", (ftnlen)510)) * 80, (
		    ftnlen)80);
	}
	tostdo_(" ", (ftnlen)1);
    } else if (eqstr_(inftyp, "USAGE", inftyp_len, (ftnlen)5)) {

/*        Display the usage text. */

	for (i__ = 1; i__ <= 17; ++i__) {
	    tostdo_(usgtxt + ((i__1 = i__ - 1) < 17 && 0 <= i__1 ? i__1 : 
		    s_rnge("usgtxt", i__1, "prcinf_", (ftnlen)521)) * 80, (
		    ftnlen)80);
	}
	tostdo_(" ", (ftnlen)1);
    } else if (eqstr_(inftyp, "VERSION", inftyp_len, (ftnlen)7)) {

/*        Create and display "version" message. */

	tkvrsn_("TOOLKIT", verstr, (ftnlen)7, (ftnlen)80);
	tostdo_(" ", (ftnlen)1);
/* Writing concatenation */
	i__2[0] = 53, a__1[0] = "MKDSK Program; Ver. 2.0.0, 28-FEB-2017; Too"
		"lkit Ver. ";
	i__2[1] = rtrim_(verstr, (ftnlen)80), a__1[1] = verstr;
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)133);
	tostdo_(ch__1, rtrim_(verstr, (ftnlen)80) + 53);
	tostdo_(" ", (ftnlen)1);
    } else if (eqstr_(inftyp, "HELP", inftyp_len, (ftnlen)4)) {
	for (i__ = 1; i__ <= 18; ++i__) {
	    tostdo_(hlptxt + ((i__1 = i__ - 1) < 18 && 0 <= i__1 ? i__1 : 
		    s_rnge("hlptxt", i__1, "prcinf_", (ftnlen)542)) * 80, (
		    ftnlen)80);
	}
	tostdo_(" ", (ftnlen)1);
    } else {

/*        We shouldn't arrive here. */

	setmsg_("Informational message type # is not supported.", (ftnlen)46);
	errch_("#", inftyp, (ftnlen)1, inftyp_len);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("PRCINF", (ftnlen)6);
	return 0;
    }
    chkout_("PRCINF", (ftnlen)6);
    return 0;
} /* prcinf_ */

