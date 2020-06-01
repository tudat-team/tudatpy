/* dspdsc.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;
static integer c__37 = 37;
static integer c__1 = 1;
static integer c__38 = 38;
static integer c__2 = 2;

/* $Procedure DSPDSC ( Display DSK segment descriptor ) */
/* Subroutine */ int dspdsc_(doublereal *dskdsc, integer *n, integer *items, 
	integer *nsig)
{
    /* Initialized data */

    static char typlst[80*4] = "<Not implemented>                           "
	    "                                    " "Shape model using triangu"
	    "lar plates                                             " "<Not i"
	    "mplemented>                                                     "
	    "          " "Packed integer DEM                                 "
	    "                             ";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    icilist ici__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rnge(char *, integer, 
	    char *, integer);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    char body[36];
    doublereal minx, miny, maxx, maxy, maxz, minz, f;
    integer i__, j;
    char table[132*3];
    extern /* Subroutine */ int etcal_(doublereal *, char *, ftnlen);
    char frame[32];
    extern /* Subroutine */ int chkin_(char *, ftnlen), repmc_(char *, char *,
	     char *, char *, ftnlen, ftnlen, ftnlen, ftnlen), repmi_(char *, 
	    char *, integer *, char *, ftnlen, ftnlen, ftnlen), bodc2n_(
	    integer *, char *, logical *, ftnlen), srfc2s_(integer *, integer 
	    *, char *, logical *, ftnlen);
    extern logical failed_(void);
    integer frmcde;
    char labels[132*3];
    integer srface;
    char begtim[35];
    integer dclass;
    doublereal maxrad, minrad;
    integer bodyid;
    char endtim[35];
    doublereal valcol[3];
    logical isname;
    doublereal minalt, minlat;
    char srfnam[36];
    doublereal maxalt, maxlat, maxlon, minlon;
    char typnam[80], outlin[132];
    doublereal values[6]	/* was [2][3] */;
    extern logical return_(void);
    integer corsys, starts[3], typcde;
    extern /* Subroutine */ int chkout_(char *, ftnlen), tostdo_(char *, 
	    ftnlen), frmnam_(integer *, char *, ftnlen), cortab_(integer *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    char *, ftnlen, ftnlen), setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), sigerr_(char *, ftnlen), suffix_(char *, 
	    integer *, char *, ftnlen, ftnlen);
    logical fnd;
    extern doublereal dpr_(void);

/* $ Abstract */

/*     Display a specified set of attributes from a DSK segment */
/*     descriptor. */

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

/*     DSKBRIEF */

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

/* $ Abstract */

/*     Declare public surface name/ID mapping parameters. */

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

/*     NAIF_IDS */

/* $ Keywords */

/*     CONVERSION */
/*     NAME */
/*     STRING */
/*     SURFACE */

/* $ Restrictions */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman      (JPL) */

/* $ Literature_References */

/*     None. */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 02-DEC-2015 (NJB) */

/* -& */

/*     Maximum number of surface name/ID mapping entries: */


/*     Maximum length of a surface name string: */


/*     End of file srftrn.inc. */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     DSKDSC     I   DSK segment descriptor. */
/*     N          I   Number of descriptor items to display. */
/*     ITEMS      I   Array of descriptor item specifiers. */
/*     NSIG       I   Number of significant digits in floating point */
/*                    output. */

/* $ Detailed_Input */

/*     DSKDSC     is a DSK segment descriptor for which summary */
/*                information is to be displayed. */

/*     N          is the number of DSK segment attribute specifiers */
/*                in the ITEMS array. */

/*     ITEMS      is an array of DSK segment attribute specifiers. */
/*                These specifiers identify the descriptor information */
/*                to display. Information is displayed in the order */
/*                of the corresponding elements of ITEMS. */

/*                Each specifier is the index in the DSK descriptor */
/*                of an attribute. For example, if ITEMS(1) is */
/*                set to the value */

/*                   FRMIDX */

/*                then the first item displayed is the DSK segment's */
/*                reference frame. */

/*                Pairs of coordinate bounds and time bounds are */
/*                indicated by the index for the lower bound alone. For */
/*                example, both the lower and upper latitude bounds are */
/*                displayed as the Nth item if ITEMS(N) is set to MN2IDX */
/*                and if the coordinate system is latitudinal. */

/*                When the coordinate system name is displayed, any */
/*                associated coordinate system parameters are displayed */
/*                as well. */

/*     NSIG       is the number of significant digits in floating point */
/*                numeric output. The range of NSIG is 6:17. */

/* $ Detailed_Output */

/*     None. This routine operates by side effects. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the coordinate system descriptor in DSKDSC is not */
/*         recognized, the error SPICE(NOTSUPPORTED) is signaled. */

/*     2)  If NSIG is outside of the range 6:17, it is replaced with the */
/*         closest value in this range. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine writes to standard output. */

/*     ID codes of bodies, surfaces, and reference frames are */
/*     displayed along with the corresponding names, whenever */
/*     the required ID-name mappings are available. */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    DSKBRIEF Version 1.0.0, 15-MAR-2017 (NJB) */

/*        Previous version 04-OCT-2016 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Initial values */

    if (return_()) {
	return 0;
    }
    chkin_("DSPDSC", (ftnlen)6);

/*     Display items in the order they're listed. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (items[i__ - 1] == 2) {

/*           For historical reasons, the index of the body ID */
/*           is named CTRIDX. */

	    bodyid = i_dnnt(&dskdsc[1]);

/*           Show body ID. */

	    bodc2n_(&bodyid, body, &fnd, (ftnlen)36);
	    if (failed_()) {
		chkout_("DSPDSC", (ftnlen)6);
		return 0;
	    }
	    if (! fnd) {
		s_copy(body, "Name not available", (ftnlen)36, (ftnlen)18);
	    }
	    s_copy(outlin, "Body:                               # (#)", (
		    ftnlen)132, (ftnlen)41);
	    repmi_(outlin, "#", &bodyid, outlin, (ftnlen)132, (ftnlen)1, (
		    ftnlen)132);
	    repmc_(outlin, "#", body, outlin, (ftnlen)132, (ftnlen)1, (ftnlen)
		    36, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	}
	if (items[i__ - 1] == 1) {

/*           Display the surface. */

	    srface = i_dnnt(dskdsc);
	    bodyid = i_dnnt(&dskdsc[1]);
	    srfc2s_(&srface, &bodyid, srfnam, &isname, (ftnlen)36);
	    if (! isname) {
		s_copy(srfnam, "  Name not available", (ftnlen)36, (ftnlen)20)
			;
	    }
	    s_copy(outlin, "  Surface:                          # (#)", (
		    ftnlen)132, (ftnlen)41);
	    repmi_(outlin, "#", &srface, outlin, (ftnlen)132, (ftnlen)1, (
		    ftnlen)132);
	    repmc_(outlin, "#", srfnam, outlin, (ftnlen)132, (ftnlen)1, (
		    ftnlen)36, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	}
	if (items[i__ - 1] == 5) {

/*           Display the reference frame. */

	    frmcde = i_dnnt(&dskdsc[4]);
	    frmnam_(&frmcde, frame, (ftnlen)32);
	    if (s_cmp(frame, " ", (ftnlen)32, (ftnlen)1) != 0) {
		s_copy(outlin, "  Reference frame:                  #", (
			ftnlen)132, (ftnlen)37);
		repmc_(outlin, "#", frame, outlin, (ftnlen)132, (ftnlen)1, (
			ftnlen)32, (ftnlen)132);
	    } else {
		s_copy(outlin, "  Reference frame name N/A; ID code: #", (
			ftnlen)132, (ftnlen)38);
		repmi_(outlin, "#", &frmcde, outlin, (ftnlen)132, (ftnlen)1, (
			ftnlen)132);
	    }
	    tostdo_(outlin, (ftnlen)132);
	}
	if (items[i__ - 1] == 4) {

/*           Display the data type. */

	    typcde = i_dnnt(&dskdsc[3]);
	    if (typcde > 0 && typcde <= 4) {
		s_copy(typnam, typlst + ((i__2 = typcde - 1) < 4 && 0 <= i__2 
			? i__2 : s_rnge("typlst", i__2, "dspdsc_", (ftnlen)
			340)) * 80, (ftnlen)80, (ftnlen)80);
	    } else {
		s_copy(typnam, "  Data type description not available", (
			ftnlen)80, (ftnlen)37);
	    }
	    s_copy(outlin, "  Data type:                        # (#)", (
		    ftnlen)132, (ftnlen)41);
	    repmi_(outlin, "#", &typcde, outlin, (ftnlen)132, (ftnlen)1, (
		    ftnlen)132);
	    repmc_(outlin, "#", typnam, outlin, (ftnlen)132, (ftnlen)1, (
		    ftnlen)80, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	}
	if (items[i__ - 1] == 3) {
	    dclass = i_dnnt(&dskdsc[2]);

/*           Display the data class. */

	    if (dclass == 1) {
		s_copy(outlin, "  Data class:                       1 (Singl"
			"e-valued surface)", (ftnlen)132, (ftnlen)61);
	    } else if (dclass == 2) {
		s_copy(outlin, "  Data class:                       2 (Gener"
			"al surface)", (ftnlen)132, (ftnlen)55);
	    } else {
		s_copy(outlin, "  Data class:                       unknown", 
			(ftnlen)132, (ftnlen)43);
	    }
	    tostdo_(outlin, (ftnlen)132);
	}
	corsys = i_dnnt(&dskdsc[5]);
	if (items[i__ - 1] == 6) {

/*           Display the coordinate system and coordinate system */
/*           parameters, if applicable. */

/*           Display the coordinate bounds as well. */

	    s_copy(outlin, "  Coordinate system:                #", (ftnlen)
		    132, (ftnlen)37);
	    if (corsys == 1) {

/*              The system is LATITUDINAL. */

		repmc_(outlin, "#", "Planetocentric Latitudinal", outlin, (
			ftnlen)132, (ftnlen)1, (ftnlen)26, (ftnlen)132);
		tostdo_(outlin, (ftnlen)132);
	    } else if (corsys == 4) {

/*              The system is PLANETODETIC. */

		repmc_(outlin, "#", "Planetodetic", outlin, (ftnlen)132, (
			ftnlen)1, (ftnlen)12, (ftnlen)132);
		tostdo_(outlin, (ftnlen)132);
		s_copy(labels, "   Equatorial radius (km):", (ftnlen)132, (
			ftnlen)26);
		s_copy(labels + 132, "   Polar radius      (km):", (ftnlen)
			132, (ftnlen)26);
		s_copy(labels + 264, "   Flattening coefficient:", (ftnlen)
			132, (ftnlen)26);
		valcol[0] = dskdsc[6];
		f = dskdsc[7];
		valcol[1] = (1. - f) * valcol[0];
		valcol[2] = f;
		cortab_(&c__3, labels, &c__37, nsig, &c__1, valcol, starts, 
			table, (ftnlen)132, (ftnlen)132);
		for (j = 1; j <= 3; ++j) {
		    tostdo_(table + ((i__2 = j - 1) < 3 && 0 <= i__2 ? i__2 : 
			    s_rnge("table", i__2, "dspdsc_", (ftnlen)420)) * 
			    132, (ftnlen)132);
		}
	    } else if (corsys == 3) {

/*              The system is RECTANGULAR. */

		repmc_(outlin, "#", "Rectangular", outlin, (ftnlen)132, (
			ftnlen)1, (ftnlen)11, (ftnlen)132);
		tostdo_(outlin, (ftnlen)132);
	    } else {
		setmsg_("The coordinate system code # is not recognized.", (
			ftnlen)47);
		errint_("#", &corsys, (ftnlen)1);
		sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
		chkout_("DSPDSC", (ftnlen)6);
		return 0;
	    }
	}
	if (corsys == 1) {
	    minlon = dskdsc[16];
	    maxlon = dskdsc[17];
	    minlat = dskdsc[18];
	    maxlat = dskdsc[19];
	    minrad = dskdsc[20];
	    maxrad = dskdsc[21];
	    s_copy(labels, "    Min, max longitude  (deg):", (ftnlen)132, (
		    ftnlen)30);
	    s_copy(labels + 132, "    Min, max latitude   (deg):", (ftnlen)
		    132, (ftnlen)30);
	    s_copy(labels + 264, "    Min, max radius      (km):", (ftnlen)
		    132, (ftnlen)30);
	    values[0] = minlon * dpr_();
	    values[2] = minlat * dpr_();
	    values[4] = minrad;
	    values[1] = maxlon * dpr_();
	    values[3] = maxlat * dpr_();
	    values[5] = maxrad;
	    cortab_(&c__3, labels, &c__38, nsig, &c__2, values, starts, table,
		     (ftnlen)132, (ftnlen)132);
	    if (items[i__ - 1] == 17) {

/*              Show longitude bounds. */

		tostdo_(table, (ftnlen)132);
	    }
	    if (items[i__ - 1] == 19) {

/*              Show latitude bounds. */

		tostdo_(table + 132, (ftnlen)132);
	    }
	    if (items[i__ - 1] == 21) {

/*              Show radius bounds. */

		tostdo_(table + 264, (ftnlen)132);
	    }
	} else if (corsys == 4) {

/*           The system is PLANETODETIC. */

	    minlon = dskdsc[16];
	    maxlon = dskdsc[17];
	    minlat = dskdsc[18];
	    maxlat = dskdsc[19];
	    minalt = dskdsc[20];
	    maxalt = dskdsc[21];
	    s_copy(labels, "    Min, max longitude  (deg):", (ftnlen)132, (
		    ftnlen)30);
	    s_copy(labels + 132, "    Min, max latitude   (deg):", (ftnlen)
		    132, (ftnlen)30);
	    s_copy(labels + 264, "    Min, max altitude    (km):", (ftnlen)
		    132, (ftnlen)30);
	    values[0] = minlon * dpr_();
	    values[2] = minlat * dpr_();
	    values[4] = minalt;
	    values[1] = maxlon * dpr_();
	    values[3] = maxlat * dpr_();
	    values[5] = maxalt;
	    cortab_(&c__3, labels, &c__38, nsig, &c__2, values, starts, table,
		     (ftnlen)132, (ftnlen)132);
	    if (items[i__ - 1] == 17) {

/*              Show longitude bounds. */

		tostdo_(table, (ftnlen)132);
	    }
	    if (items[i__ - 1] == 19) {

/*              Show latitude bounds. */

		tostdo_(table + 132, (ftnlen)132);
	    }
	    if (items[i__ - 1] == 21) {

/*              Show altitude bounds. */

		tostdo_(table + 264, (ftnlen)132);
	    }
	} else if (corsys == 3) {

/*           The system is RECTANGULAR. */

	    minx = dskdsc[16];
	    maxx = dskdsc[17];
	    miny = dskdsc[18];
	    maxy = dskdsc[19];
	    minz = dskdsc[20];
	    maxz = dskdsc[21];
	    s_copy(labels, "   Min, max X coordinate (km):", (ftnlen)132, (
		    ftnlen)30);
	    s_copy(labels + 132, "   Min, max Y coordinate (km):", (ftnlen)
		    132, (ftnlen)30);
	    s_copy(labels + 264, "   Min, max Z coordinate (km):", (ftnlen)
		    132, (ftnlen)30);
	    values[0] = minx;
	    values[2] = miny;
	    values[4] = minz;
	    values[1] = maxx;
	    values[3] = maxy;
	    values[5] = maxz;
	    cortab_(&c__3, labels, &c__38, nsig, &c__2, values, starts, table,
		     (ftnlen)132, (ftnlen)132);
	    if (items[i__ - 1] == 17) {

/*              Show X bounds. */

		tostdo_(table, (ftnlen)132);
	    }
	    if (items[i__ - 1] == 19) {
		tostdo_(table + 132, (ftnlen)132);
	    }
	    if (items[i__ - 1] == 21) {
		tostdo_(table + 264, (ftnlen)132);
	    }
	} else {
	    setmsg_("The coordinate system code # is not recognized.", (
		    ftnlen)47);
	    errint_("#", &corsys, (ftnlen)1);
	    sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	    chkout_("DSPDSC", (ftnlen)6);
	    return 0;
	}
	if (items[i__ - 1] == 23) {

/*           Display both the start and stop times. */

	    etcal_(&dskdsc[22], begtim, (ftnlen)35);
	    etcal_(&dskdsc[23], endtim, (ftnlen)35);
	    suffix_("TDB", &c__1, begtim, (ftnlen)3, (ftnlen)35);
	    suffix_("TDB", &c__1, endtim, (ftnlen)3, (ftnlen)35);
/* Writing concatenation */
	    i__3[0] = 36, a__1[0] = "  Start time:                       ";
	    i__3[1] = 35, a__1[1] = begtim;
	    s_cat(outlin, a__1, i__3, &c__2, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 132;
	    ici__1.iciunit = outlin;
	    ici__1.icifmt = "(A,(1PE24.16))";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, "    Seconds past J2000 TDB:             ", (ftnlen)
		    40);
	    do_fio(&c__1, (char *)&dskdsc[22], (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    tostdo_(outlin, (ftnlen)132);
/* Writing concatenation */
	    i__3[0] = 36, a__1[0] = "  Stop time:                        ";
	    i__3[1] = 35, a__1[1] = endtim;
	    s_cat(outlin, a__1, i__3, &c__2, (ftnlen)132);
	    tostdo_(outlin, (ftnlen)132);
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 132;
	    ici__1.iciunit = outlin;
	    ici__1.icifmt = "(A,(1PE24.16))";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, "    Seconds past J2000 TDB:             ", (ftnlen)
		    40);
	    do_fio(&c__1, (char *)&dskdsc[23], (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    tostdo_(outlin, (ftnlen)132);
	}
    }
    chkout_("DSPDSC", (ftnlen)6);
    return 0;
} /* dspdsc_ */

