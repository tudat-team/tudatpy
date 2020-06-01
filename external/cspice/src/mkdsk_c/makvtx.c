/* makvtx.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure MAKVTX ( MKDSK: create vertex from coordinates ) */
/* Subroutine */ int makvtx_(integer *corsys, doublereal *corpar, doublereal *
	coords, doublereal *refval, doublereal *height, doublereal *vertex)
{
    doublereal f;
    extern /* Subroutine */ int chkin_(char *, ftnlen), errdp_(char *, 
	    doublereal *, ftnlen);
    doublereal re;
    extern /* Subroutine */ int georec_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), latrec_(
	    doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal radius;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), setmsg_(char *, 
	    ftnlen), chkout_(char *, ftnlen), errint_(char *, integer *, 
	    ftnlen);
    extern logical return_(void);
    doublereal alt;

/* $ Abstract */

/*     Create a vertex from domain coordinates and height. */

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

/*     MKDSK */

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

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     CORSYS     I   Coordinate system. */
/*     CORPAR     I   Coordinate parameters. */
/*     COORDS     I   Domain coordinates. */
/*     REFVAL     I   Height reference value. */
/*     HEIGHT     I   Height. */
/*     VERTEX     O   Vertex. */

/* $ Detailed_Input */

/*     CORSYS     is a DSK subsystem code designating the coordinate */
/*                system of the input coordinates. */

/*     CORPAR     is an array containing parameters associated with */
/*                the input coordinate system. The contents of the */
/*                array are as described in the DSK include file */
/*                dskdsc.inc. */

/*     COORDS     is a pair of domain coordinates: these may be, */

/*                   - planetocentric longitude and latitude */

/*                   - planetodetic longitude and latitude */

/*                   - X and Y */

/*                For a given coordinate system, the order of the */
/*                elements of COORDS is that of the coordinate names in */
/*                the list above. */

/*     REFVAL     is a reference value to be added to the input height. */
/*                REFVAL is used only for latitudinal and rectangular */
/*                coordinates. */

/*                REFVAL must be non-negative. */

/*                Units are km. */


/*     HEIGHT     is a height datum. Units are km. */


/* $ Detailed_Output */

/*     VERTEX     is a 3-vector corresponding to the input coordinates, */
/*                height, and if applicable, height reference. */

/*                Units are always km. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the coordinate system code is not recognized, the error */
/*         SPICE(NOTSUPPORTED) is signaled. */

/*     2)  If an error occurs while converting planetodetic coordinates */
/*         to rectangular coordinates, the error will be diagnosed by a */
/*         routine in the call tree of this routine. */

/*     3)  If REFVAL is negative, the error SPICE(VALUEOUTOFRANGE) is */
/*         signaled. REFVAL is checked whether or not it is applicable. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     None. */

/* $ Examples */

/*     See usage in the MKDSK routine MKVARR. */

/* $ Restrictions */

/*     1) For use only within program MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 1.0.0, 25-FEB-2017 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("MAKVTX", (ftnlen)6);
    if (*refval < 0.) {
	setmsg_("Reference value # must be non-negative.", (ftnlen)39);
	errdp_("#", refval, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("MAKVTX", (ftnlen)6);
	return 0;
    }
    if (*corsys == 1) {

/*        REFVAL has the same units as HEIGHT. */
/*        HSCALE converts these units to km. */

	radius = *refval + *height;
	latrec_(&radius, coords, &coords[1], vertex);
    } else if (*corsys == 4) {

/*        Height is relative to the system's reference spheroid. */

	re = corpar[0];
	f = corpar[1];
	alt = *height;
	georec_(coords, &coords[1], &alt, &re, &f, vertex);
    } else if (*corsys == 3) {

/*        Height is relative to the reference Z-value. */

/*        REFVAL has the same units as HEIGHT. */
/*        HSCALE converts these units to km. */

	vertex[0] = coords[0];
	vertex[1] = coords[1];
	vertex[2] = *refval + *height;
    } else {
	setmsg_("Coordinate system code # is not recognized.", (ftnlen)43);
	errint_("#", corsys, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("MAKVTX", (ftnlen)6);
	return 0;
    }
    chkout_("MAKVTX", (ftnlen)6);
    return 0;
} /* makvtx_ */

