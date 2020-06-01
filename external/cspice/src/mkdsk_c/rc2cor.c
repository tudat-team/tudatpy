/* rc2cor.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure RC2COR ( MKDSK: map grid row and column to coordinates ) */
/* Subroutine */ int rc2cor_(doublereal *lftcor, doublereal *topcor, 
	doublereal *colstp, doublereal *rowstp, integer *col, integer *row, 
	doublereal *coords)
{
    extern /* Subroutine */ int chkin_(char *, ftnlen), errdp_(char *, 
	    doublereal *, ftnlen), sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen);
    extern logical return_(void);

/* $ Abstract */

/*     Map grid row and column indices to standard coordinates. */

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
/*     LFTCOR     I   Name of input file. */
/*     TOPCOR     I   Maximum number of values to return. */
/*     COLSTP     I   Column step size. */
/*     ROWSTP     I   Row step size. */
/*     COL        I   Column index. */
/*     ROW        I   Row index. */
/*     COORDS     O   Coordinates of input grid point. */

/* $ Detailed_Input */

/*     LFTCOR     is the coordinate of the leftmost column of a */
/*                rectangular data grid. If the coordinate system is */
/*                latitudinal, this is the minimum longitude on the */
/*                grid. If the system is rectangular, this is the */
/*                minimum X value on the grid. */

/*                Units are caller-defined but must be consistent */
/*                with those of COLSTP. */

/*     TOPCOR     is the coordinate of the top row of a rectangular data */
/*                grid. If the coordinate system is latitudinal, this is */
/*                the maximum latitude on the grid. If the system is */
/*                rectangular, this is the maximum Y value on the grid. */

/*                Units are caller-defined but must be consistent */
/*                with those of ROWSTP. */

/*     COLSTP     is the uniform separation between columns of the grid. */

/*                Units are caller-defined but must be consistent */
/*                with those of LFTCOR. */
/*     ROWSTP     is the uniform separation between rows of the grid. */

/*                Units are caller-defined but must be consistent */
/*                with those of TOPCOR. */


/*     COL        is the index of a column in the grid. Indices are */
/*                1-based: the first coordinate of the first column */
/*                is LFTCOR. */

/*     ROW        is the index of a row in the grid. Indices are */
/*                1-based: the second coordinate of the top row */
/*                is TOPCOR. */


/* $ Detailed_Output */

/*     COORDS     is an array containing the coordinates corresponding */
/*                to ROW and COL. The first element of COORDS is the */
/*                coordinate corresponding to COL; the second element */
/*                corresponds to ROW. */

/*                For example, if the coordinate system is latitudinal, */
/*                COORDS(1) is the longitude corresponding to COL and */
/*                COORDS(2) is the latitude corresponding to ROW. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If either the row or column step is not strictly positive, */
/*         the error SPICE(INVALIDSTEP) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     The computation performed by this routine refers to a rectangular */
/*     grid in a latitudinal, planetodetic, or rectangular coordinate */
/*     system. The coordinate system is not specified but must be */
/*     the same for all inputs. */

/*     All rows have equal separation in latitude or Y-coordinate, */
/*     depending on the coordinate system used. All columns have equal */
/*     separation in longitude or X-coordinate. */

/*     Column 1 of the grid corresponds to the coordinate value specified */
/*     by LFTCOR. */

/*     Row 1 of the grid corresponds to the coordinate value specified */
/*     by TOPCOR. Thus the second coordinate of the rows decreases with */
/*     increasing row index. */

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
    chkin_("RC2COR", (ftnlen)6);

/*     Make sure the step sizes are valid. */

    if (*colstp <= 0.) {
	setmsg_("Column step was #; must be strictly positive.", (ftnlen)45);
	errdp_("#", colstp, (ftnlen)1);
	sigerr_("SPICE(INVALIDSTEP)", (ftnlen)18);
	chkout_("RC2COR", (ftnlen)6);
	return 0;
    }
    if (*rowstp <= 0.) {
	setmsg_("Row step was #; must be strictly positive.", (ftnlen)42);
	errdp_("#", rowstp, (ftnlen)1);
	sigerr_("SPICE(INVALIDSTEP)", (ftnlen)18);
	chkout_("RC2COR", (ftnlen)6);
	return 0;
    }

/*     In a lon/lat system, the leftmost column is at minimum longitude; */
/*     the top row is at maximum latitude. */

/*     In the rectangular system, the leftmost column is at minimum */
/*     X; the top row is at maximum Y. */

/*     For all systems, the computation is the same. */

    coords[0] = *lftcor + (*col - 1) * *colstp;
    coords[1] = *topcor - (*row - 1) * *rowstp;
    chkout_("RC2COR", (ftnlen)6);
    return 0;
} /* rc2cor_ */

