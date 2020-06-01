/* dspgap.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;

/* $Procedure   DSPGAP ( DSKBRIEF, display spatial coverage gaps ) */
/* Subroutine */ int dspgap_(logical *gap, integer *corsys, integer *nsig, 
	integer *ncomp, doublereal *bds1, doublereal *bds2)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer from, i__;
    char table[132*1000];
    doublereal scale;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    char co1str[132], co2str[132];
    integer start1;
    extern logical failed_(void);
    char labels[1*1000];
    extern /* Subroutine */ int cortab_(integer *, char *, integer *, integer 
	    *, integer *, doublereal *, integer *, char *, ftnlen, ftnlen);
    integer remain, nlines;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen);
    doublereal values[4000]	/* was [4][1000] */;
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen);
    char outlin[132];
    extern /* Subroutine */ int tostdo_(char *, ftnlen);
    extern logical return_(void);
    integer starts[4];
    extern doublereal dpr_(void);

/* $ Abstract */

/*     Display spatial coverage gaps of a collection of DSK segments. */
/*     The gaps are expressed as a union of rectangles in a specified */
/*     coordinate system. */

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

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     GAP        I   Flag commanding display of summary or warning. */
/*     CORSYS     I   Coordinate system code. */
/*     NSIG       I   Number of significant digits. */
/*     NCOMP      I   Number of rectangles in gap region. */
/*     BDS1       I   Bounds of first coordinates of rectangles. */
/*     BDS2       I   Bounds of second coordinates of rectangles. */

/* $ Detailed_Input */

/*     GAP        is a logical flag indicating whether a gap summary or */
/*                warning is to be displayed. The summary is displayed */
/*                if and only if GAP is .TRUE. and coverage gaps exist. */

/*     CORSYS     is an integer code indicating the coordinate */
/*                system in which the coverage gaps are represented. */

/*     NSIG       is the number of significant digits to display */
/*                for floating point values. */

/*     NCOMP      is the number of coordinate rectangles---also */
/*                called "components"---comprising the coverage */
/*                gap region. */

/*     BDS1       is an array of bounds of the first coordinates of the */
/*                gap rectangles. If the coordinate system is */
/*                latitudinal or planetodetic, the first coordinate is */
/*                longitude. If the coordinate system is rectangular, */
/*                the first coordinate is X. */

/*                The first coordinate of the Ith rectangle ranges from */

/*                   BDS1(1,I) to BDS1(2,I) */

/*                Units of longitude are radians; rectangular units are */
/*                km. */

/*     BDS2       is an array of bounds of the second coordinates of the */
/*                gap rectangles. If the coordinate system is */
/*                latitudinal or planetodetic, the second coordinate is */
/*                latitude. If the coordinate system is rectangular, the */
/*                second coordinate is Y. */

/*                The second coordinate of the Ith rectangle ranges from */

/*                   BDS2(1,I) to BDS2(2,I) */

/*                Units of latitude are radians; rectangular units are */
/*                km. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If an invalid coordinate system is encountered, this routine */
/*         signals the error SPICE(NOTSUPPORTED). */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine displays summary information for spatial coverage */
/*     gaps of a set of DSK segments. The gaps are represented as union */
/*     of rectangles in a specified coordinate system. */

/*     The display is written to standard output. */

/*     The coverage "gap region" of a set of DSK segments is that */
/*     region, within the rectangle defined by the global extrema of the */
/*     extents of the segments' first and second coordinates, where */
/*     there is no spatial coverage. */

/*     For example, if a set of DSK segments uses latitudinal */
/*     coordinates, the global extrema of the sets coordinates are the */
/*     minimum and maximum longitudes and latitudes, where the extrema */
/*     are taken over all segments in the set. Any point having */
/*     longitude between the set's longitude extrema and latitude */
/*     between the set's latitude extrema, such that point's longitude */
/*     and latitude are not within the coverage of any segment of the */
/*     set, lies in a gap. The coverage gap of the set is the union of */
/*     all such points. This region can be represented as a union of */
/*     rectangles. */

/*     In order for the concept of spatial coverage of a set of DSK */
/*     segments to make sense, the segments must have the common */
/*     attributes: */

/*        Body */
/*        Surface */
/*        Frame */

/*        Coordinate system */

/*           If the coordinate system is planetodetic, the */
/*           system parameters must match exactly. */

/*        Data type */
/*        Data class */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -    DSKBRIEF Version 1.0.0, 10-FEB-2017 (NJB) */

/*        Original version 07-OCT-2016 (NJB) */

/* -& */
/* $ Index_Entries */

/*     display spatial coverage gaps of dsk segment set */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("DSPGAP", (ftnlen)6);
    for (i__ = 1; i__ <= 1000; ++i__) {
	*(unsigned char *)&labels[(i__1 = i__ - 1) < 1000 && 0 <= i__1 ? i__1 
		: s_rnge("labels", i__1, "dspgap_", (ftnlen)246)] = ' ';
	s_copy(table + ((i__1 = i__ - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"table", i__1, "dspgap_", (ftnlen)247)) * 132, " ", (ftnlen)
		132, (ftnlen)1);
    }

/*     Indent the first label by 2 spaces relative to the label above */
/*     it. */

    start1 = 7;

/*     Set up coordinate system-specific column titles */
/*     and scale factors. */

    if (*gap && *ncomp > 0) {
	tostdo_("    Coverage gaps:", (ftnlen)18);
	if (*corsys == 1 || *corsys == 4) {
	    scale = dpr_();
	    s_copy(co1str, "Longitude range (deg)", (ftnlen)132, (ftnlen)21);
	    s_copy(co2str, "Latitude range (deg)", (ftnlen)132, (ftnlen)20);
	} else if (*corsys == 3) {
	    scale = 1.;
	    s_copy(co1str, "X range (km)", (ftnlen)132, (ftnlen)12);
	    s_copy(co2str, "Y range (km)", (ftnlen)132, (ftnlen)12);
	} else {
	    setmsg_("Coordinate system # is not currently supported.", (
		    ftnlen)47);
	    errint_("#", corsys, (ftnlen)1);
	    sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	    chkout_("DSPGAP", (ftnlen)6);
	    return 0;
	}

/*        Create the gap table for the first batch */
/*        of data. (We need to get the column positions */
/*        before we output the header, hence the order */
/*        of operations conducted here.) */

	remain = *ncomp;
	nlines = min(1000,remain);
	i__1 = nlines;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    values[(i__2 = (i__ << 2) - 4) < 4000 && 0 <= i__2 ? i__2 : 
		    s_rnge("values", i__2, "dspgap_", (ftnlen)301)] = bds1[(
		    i__ << 1) - 2] * scale;
	    values[(i__2 = (i__ << 2) - 3) < 4000 && 0 <= i__2 ? i__2 : 
		    s_rnge("values", i__2, "dspgap_", (ftnlen)302)] = bds1[(
		    i__ << 1) - 1] * scale;
	    values[(i__2 = (i__ << 2) - 2) < 4000 && 0 <= i__2 ? i__2 : 
		    s_rnge("values", i__2, "dspgap_", (ftnlen)303)] = bds2[(
		    i__ << 1) - 2] * scale;
	    values[(i__2 = (i__ << 2) - 1) < 4000 && 0 <= i__2 ? i__2 : 
		    s_rnge("values", i__2, "dspgap_", (ftnlen)304)] = bds2[(
		    i__ << 1) - 1] * scale;
	}
	cortab_(&nlines, labels, &start1, nsig, &c__4, values, starts, table, 
		(ftnlen)1, (ftnlen)132);
	if (failed_()) {
	    chkout_("DSPGAP", (ftnlen)6);
	    return 0;
	}
	remain -= nlines;

/*        Display the header. */

	s_copy(outlin, " ", (ftnlen)132, (ftnlen)1);
	s_copy(outlin + (start1 - 1), co1str, 132 - (start1 - 1), (ftnlen)132)
		;
	i__1 = starts[2] - 1;
	s_copy(outlin + i__1, co2str, 132 - i__1, (ftnlen)132);
	tostdo_(outlin, (ftnlen)132);

/*        Display the  first batch of gap information. */

	i__1 = nlines;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tostdo_(table + ((i__2 = i__ - 1) < 1000 && 0 <= i__2 ? i__2 : 
		    s_rnge("table", i__2, "dspgap_", (ftnlen)332)) * 132, (
		    ftnlen)132);
	}

/*        Display the remaining gap information. */

	from = nlines + 1;
	while(remain > 0) {
	    nlines = min(1000,remain);
	    i__1 = nlines;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		values[(i__2 = (i__ << 2) - 4) < 4000 && 0 <= i__2 ? i__2 : 
			s_rnge("values", i__2, "dspgap_", (ftnlen)347)] = 
			bds1[(from << 1) - 2] * scale;
		values[(i__2 = (i__ << 2) - 3) < 4000 && 0 <= i__2 ? i__2 : 
			s_rnge("values", i__2, "dspgap_", (ftnlen)348)] = 
			bds1[(from << 1) - 1] * scale;
		values[(i__2 = (i__ << 2) - 2) < 4000 && 0 <= i__2 ? i__2 : 
			s_rnge("values", i__2, "dspgap_", (ftnlen)349)] = 
			bds2[(from << 1) - 2] * scale;
		values[(i__2 = (i__ << 2) - 1) < 4000 && 0 <= i__2 ? i__2 : 
			s_rnge("values", i__2, "dspgap_", (ftnlen)350)] = 
			bds2[(from << 1) - 1] * scale;
	    }
	    from += nlines;
	    cortab_(&nlines, labels, &start1, nsig, &c__4, values, starts, 
		    table, (ftnlen)1, (ftnlen)132);
	    if (failed_()) {
		chkout_("DSPGAP", (ftnlen)6);
		return 0;
	    }
	    i__1 = nlines;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		tostdo_(table + ((i__2 = i__ - 1) < 1000 && 0 <= i__2 ? i__2 :
			 s_rnge("table", i__2, "dspgap_", (ftnlen)366)) * 132,
			 (ftnlen)132);
	    }
	    remain -= nlines;
	}
    } else if (*ncomp > 0) {

/*        Gap display is not enabled, but there are gaps. */

	tostdo_("    ***Coverage has gaps. Use the -gaps option to display t"
		"hem.***", (ftnlen)66);
    }
    chkout_("DSPGAP", (ftnlen)6);
    return 0;
} /* dspgap_ */

