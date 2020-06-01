/* mkvarr.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b3 = 1.;
static integer c__1000 = 1000;

/* $Procedure MKVARR ( MKDSK: make vertex array from grid data ) */
/* Subroutine */ int mkvarr_(char *infile__, char *aunits, char *dunits, 
	logical *rowmaj, logical *topdwn, logical *leftrt, integer *corsys, 
	doublereal *corpar, doublereal *refval, doublereal *hscale, integer *
	ncols, integer *nrows, doublereal *lftcor, doublereal *topcor, 
	doublereal *colstp, doublereal *rowstp, integer *maxnv, doublereal *
	verts, ftnlen infile_len, ftnlen aunits_len, ftnlen dunits_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    doublereal href;
    logical done;
    integer i__, j, k, n;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    doublereal hstep, topco;
    integer total;
    doublereal vstep;
    extern /* Subroutine */ int rdgrd5_(char *, integer *, integer *, 
	    doublereal *, logical *, ftnlen), rc2cor_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *);
    integer nc;
    extern logical failed_(void);
    doublereal ascale, dscale;
    integer nr, nv;
    doublereal height, leftco, corscl, coords[2];
    static doublereal values[1000];
    extern /* Subroutine */ int convrt_(doublereal *, char *, char *, 
	    doublereal *, ftnlen, ftnlen);
    extern logical return_(void);
    extern /* Subroutine */ int chkout_(char *, ftnlen), setmsg_(char *, 
	    ftnlen), errint_(char *, integer *, ftnlen), sigerr_(char *, 
	    ftnlen), makvtx_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    integer col, row;

/* $ Abstract */

/*     Make normalized vertex array from a regular height grid. Elements */
/*     of the output array are in row-major, top-down, left-right order. */

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
/*     INFILE     I   Name of input file. */
/*     AUNITS     I   Angular units. */
/*     DUNITS     I   Distance units. */
/*     ROWMAJ     I   Flag indicting whether grid is row-major. */
/*     TOPDWN     I   Flag indicting whether grid is top-down. */
/*     LEFTRT     I   Flag indicting whether grid is left-right. */
/*     CORSYS     I   Coordinate system. */
/*     CORPAR     I   Coordinate parameters. */
/*     REFVAL     I   Height reference value. */
/*     HSCALE     I   Height scale. */
/*     NROWS      I   Number of rows in grid. */
/*     NCOLS      I   Number of columns in grid. */
/*     LFTCOR     I   Coordinate of leftmost column of grid. */
/*     TOPCOR     I   Coordinate of top row of grid. */
/*     COLSTP     I   Column step size. */
/*     ROWSTP     I   Row step size. */
/*     MAXNV      I   Maximum number of vertices to return. */
/*     VERTS      O   Vertex array. */

/* $ Detailed_Input */

/*     INFILE     is the name of an input data file containing height */
/*                grid data. */


/*     AUNITS     is the name of the angular unit associated with the */
/*                grid coordinates, if the grid coordinate system is */
/*                latitudinal or planetodetic. AUNITS must be supported */
/*                by the SPICELIB routine CONVRT. */

/*     DUNITS     is the name of the distance unit associated with the */
/*                grid coordinates. DUNITS must be supported by the */
/*                SPICELIB routine CONVRT. */

/*     ROWMAJ     is a logical flag that is set by the caller to .TRUE. */
/*                when the input grid data are organized in row-major */
/*                order, and that is set to .FALSE. when the data are in */
/*                column-major order. */

/*                "Row-major" means that the Nth consecutive sequence of */
/*                NCOLS tokens in the file represents the Nth row of */
/*                data in the output grid. */

/*                "Column-major" means that the Nth consecutive sequence */
/*                of NROWS tokens in the file represents the Nth column */
/*                of data in the output grid. */

/*                Note that the mapping from a token's position in the */
/*                input file to its position in the output grid is not */
/*                defined by ROWMAJ alone: the values of TOPDWN and */
/*                LEFTRT are needed as well to fully define the mapping. */


/*     TOPDWN     is a logical flag that is set by the caller to .TRUE. */
/*                when the input file contains data in top-down order, */
/*                and is set to .FALSE. when the input data are in */
/*                bottom-up order. Here "top" means "having the highest */
/*                value of the second coordinate." */

/*                TOPDWN is true if and only if the datum of the top */
/*                row and Nth column in the output grid precedes the */
/*                datum for any other row and Nth column in the input */
/*                file. In other words, the data from the input file */
/*                fill in the rows of the output grid in "top-down" */
/*                order. */

/*                When the input file is row-major and TOPDWN is .TRUE., */
/*                the top row of the output grid contains the data from */
/*                the first NCOLS tokens in the input file. When the */
/*                input file is column-major and TOPDWN is .TRUE., the */
/*                columns themselves are in top-down order: the first */
/*                element of each sequence of NROWS tokens belongs to */
/*                the top row of the output grid. */

/*     LEFTRT     is a logical flag that is set by the caller to .TRUE. */
/*                when the input file contains data in left-right order, */
/*                and is set to .FALSE. when the input data are in */
/*                right-left order. Here "left" means "having the lowest */
/*                value of the first coordinate." */

/*                LEFTRT is true if and only if the datum of the left */
/*                column and Nth row in the output grid precedes the */
/*                datum for any other column and Nth row in the input */
/*                file. In other words, the data from the input file */
/*                fill in the columns of the output grid in "left-right" */
/*                order. */

/*                When the input file is column-major and LEFTRT is */
/*                .TRUE., the leftmost column of the output grid */
/*                contains the data from the first NROWS tokens in the */
/*                input file. When the input file is row-major and LEFT */
/*                is .TRUE., the rows themselves are in left-right */
/*                order: the first element of each sequence of NCOLS */
/*                tokens belongs to the left column of the output grid. */

/*     CORSYS     is a DSK subsystem code designating the coordinate */
/*                system of the input coordinates. */

/*     CORPAR     is an array containing parameters associated with */
/*                the input coordinate system. The contents of the */
/*                array are as described in the DSK include file */
/*                dskdsc.inc. */

/*     REFVAL     is a reference value to be added to the input height. */
/*                REFVAL is used only for latitudinal and rectangular */
/*                coordinates. */

/*                REFVAL must be non-negative. */

/*                Units are km. */

/*     HSCALE     is a conversion factor to be applied to height data. */
/*                Multiplying a datum by HSCALE converts the datum to */
/*                kilometers. HSCALE need not correspond to a standard */
/*                unit. */

/*     NROWS      is the number of rows in the output grid. */

/*     NCOLS      is the number of columns in the output grid. */

/*     LFTCOR     is the coordinate of the leftmost column of a */
/*                rectangular data grid. If the coordinate system is */
/*                latitudinal, this is the minimum longitude on the */
/*                grid. If the system is rectangular, this is the */
/*                minimum X value on the grid. */

/*                Units are given by AUNITS or DUNITS, as applicable. */

/*     TOPCOR     is the coordinate of the top row of a rectangular data */
/*                grid. If the coordinate system is latitudinal, this is */
/*                the maximum latitude on the grid. If the system is */
/*                rectangular, this is the maximum Y value on the grid. */

/*                Units are given by AUNITS or DUNITS, as applicable. */

/*     COLSTP     is the uniform separation between columns of the grid. */

/*                Units are given by AUNITS or DUNITS, as applicable. */

/*     ROWSTP     is the uniform separation between rows of the grid. */

/*                Units are given by AUNITS or DUNITS, as applicable. */

/*     MAXNV      is the maximum number of vertices that can be placed */
/*                in the output array. */


/* $ Detailed_Output */

/*     VERTS      is an array of 3-vectors corresponding to the height */
/*                data in the input file. The data in VERTS always */
/*                represent a grid organized in */

/*                   row-major */
/*                   top-down */
/*                   left-right */

/*                order, regardless of the organization of the input */
/*                file. */

/*                VERTS should be declared by the caller as: */

/*                   DOUBLE PRECISION VERTS ( 3, MAXNV ) */

/*                Units are always km. */

/* $ Parameters */

/*     The input file must have a maximum line length of LNSIZE. */

/*     See mkdsk.inc. */

/* $ Exceptions */

/*     1)  If the number of data values read does not match the value */

/*             NROWS * NCOLS */

/*         the error SPICE(INVALIDDATACOUNT) is signaled. */

/*     2)  If an error occurs while converting angular units to radians */
/*         or distance units to km, the error will be diagnosed by a */
/*         routine in the call tree of this routine. */

/*     3)  If MAXNV indicates the output buffer is too small to hold the */
/*         output array, the error SPICE(BUFFERTOOSMALL) is signaled. */

/*     4)  If REFVAL is negative, the error is diagnosed by routines */
/*         in the call tree of this routine. */

/*     5)  If the input coordinate parameters in CORPAR are invalid, the */
/*         error is diagnosed by routines in the call tree of this */
/*         routine. */

/*     6)  If an error occurs while reading the input file, the */
/*         error is diagnosed by routines in the call tree of this */
/*         routine. */

/* $ Files */

/*     The file specified by INFILE can have any of the attributes (one */
/*     choice from each row below): */

/*        row-major  or column-major */
/*        top-down   or bottom-up */
/*        left-right or right-left */

/*     The number of tokens per line may vary. The number need have no */
/*     particular relationship to the row or column dimensions of the */
/*     output grid. */

/*     The file must contain only tokens that can be read as double */
/*     precision values. No non-printing characters can be present in */
/*     the file. */

/*     Tokens can be delimited by blanks or commas. Tokens must not be */
/*     split across lines. */

/*     Blank lines are allowed; however, their use is discouraged */
/*     because they'll cause line numbers in diagnostic messages to */
/*     be out of sync with actual line numbers in the file. */

/*     The file must end with a line terminator. */

/* $ Particulars */

/*     The "output grid" is the two-dimensional array of 3-dimensional */
/*     vertices; the row count is NROWS and the column count is NCOLS. */

/*     The output grid is "normalized" in the sense that it always is */

/*        row-major */
/*        top-down */
/*        left-right */

/* $ Examples */

/*     See usage in the MKDSK routine MKGRID. */

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


/*     Local parameters */


/*     Local variables */


/*     Saved variables */

    if (return_()) {
	return 0;
    }
    chkin_("MKVARR", (ftnlen)6);

/*     Compute unit scales. */

    convrt_(&c_b3, aunits, "RADIANS", &ascale, aunits_len, (ftnlen)7);
    convrt_(&c_b3, dunits, "KM", &dscale, dunits_len, (ftnlen)2);
    if (failed_()) {
	chkout_("MKVARR", (ftnlen)6);
	return 0;
    }

/*     Convert grid boundary and step values to radians and km, */
/*     as needed. Convert the reference value as well. */

    if (*corsys == 3) {
	corscl = dscale;
    } else {
	corscl = ascale;
    }
    leftco = corscl * *lftcor;
    topco = corscl * *topcor;
    hstep = corscl * *colstp;
    vstep = corscl * *rowstp;
    href = dscale * *refval;

/*     Presume the row and column counts are correct. */

    nv = *nrows * *ncols;

/*     Make sure we have room for the output array. */

    if (nv > *maxnv) {
	setmsg_("Room for # vertices is needed; amount available is #.", (
		ftnlen)53);
	errint_("#", &nv, (ftnlen)1);
	errint_("#", maxnv, (ftnlen)1);
	sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
	chkout_("MKVARR", (ftnlen)6);
	return 0;
    }

/*     Fetch data from the input file; distribute it to the */
/*     vertex array. Grab the first buffer of data: */

    rdgrd5_(infile__, &c__1000, &n, values, &done, infile_len);
    if (failed_()) {
	chkout_("MKVARR", (ftnlen)6);
	return 0;
    }
    total = 0;
    while(total < nv && n > 0) {

/*        Distribute the data we've buffered. */

	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++total;

/*           Compute the row and column of the current */
/*           vertex, then covert those indices to a one- */
/*           dimensional index. */

	    if (*rowmaj) {

/*              We fill in the vertex grid one row at a time. */

		if (*topdwn) {

/*                 We fill in rows from the top down. */

		    row = (total - 1) / *ncols + 1;
		    if (*leftrt) {

/*                    We fill in rows from left to right. */

			col = total - (row - 1) * *ncols;
		    } else {

/*                    We fill in rows from right to left. */

			j = total - (row - 1) * *ncols;
			col = *ncols + 1 - j;
		    }
		} else {

/*                 We fill in rows from the bottom up. */

		    nr = (total - 1) / *ncols + 1;
		    row = *nrows + 1 - nr;
		    if (*leftrt) {

/*                    We fill in rows from left to right. */

			col = total - (nr - 1) * *ncols;
		    } else {

/*                    We fill in rows from right to left. */

			j = total - (nr - 1) * *ncols;
			col = *ncols + 1 - j;
		    }
		}
	    } else {

/*              We fill in the vertex grid one column at a time. */

		if (*leftrt) {

/*                 We fill in columns from left to right. */

		    col = (total - 1) / *nrows + 1;
		    if (*topdwn) {

/*                    We fill in columns from the top down. */

			row = total - (col - 1) * *nrows;
		    } else {

/*                    We fill in columns from the bottom up. */

			j = total - (col - 1) * *nrows;
			row = *nrows + 1 - j;
		    }
		} else {

/*                 We fill in columns from right to left. */

		    nc = (total - 1) / *nrows + 1;
		    col = *ncols + 1 - nc;
		    if (*topdwn) {

/*                    We fill in columns from the top down. */

			row = total - (nc - 1) * *nrows;
		    } else {

/*                    We fill in columns from the bottom up. */

			j = total - (nc - 1) * *nrows;
			row = *nrows + 1 - j;
		    }
		}
	    }

/*           Compute the domain coordinates for this vertex. */

	    rc2cor_(&leftco, &topco, &hstep, &vstep, &col, &row, coords);
	    if (failed_()) {
		chkout_("MKVARR", (ftnlen)6);
		return 0;
	    }

/*           Compute the 1-D vertex index. */

	    k = (row - 1) * *ncols + col;

/*           Compute the vertex and insert it into the vertex array. */

	    height = values[(i__2 = i__ - 1) < 1000 && 0 <= i__2 ? i__2 : 
		    s_rnge("values", i__2, "mkvarr_", (ftnlen)588)] * *hscale;
	    makvtx_(corsys, corpar, coords, &href, &height, &verts[k * 3 - 3])
		    ;
	    if (failed_()) {
		chkout_("MKVARR", (ftnlen)6);
		return 0;
	    }
	}

/*        Try to read more data, but only if the EOF hasn't been */
/*        reached. Note that reading after EOF has been reached */
/*        will initiate a new read of the file, starting with the */
/*        first line. */

	if (! done) {
	    rdgrd5_(infile__, &c__1000, &n, values, &done, infile_len);
	} else {
	    n = 0;
	}
	if (failed_()) {
	    chkout_("MKVARR", (ftnlen)6);
	    return 0;
	}
    }

/*     We've read the data. Make sure we got the number of vertices */
/*     we were expecting. */

    if (total != nv) {
	setmsg_("Column count = #; row count = #. Expected height values for"
		" # vertices but found data for #. Setup file and data file  "
		"are inconsistent.", (ftnlen)136);
	errint_("#", ncols, (ftnlen)1);
	errint_("#", nrows, (ftnlen)1);
	errint_("#", &nv, (ftnlen)1);
	errint_("#", &total, (ftnlen)1);
	sigerr_("SPICE(INVALIDDATACOUNT)", (ftnlen)23);
	chkout_("MKVARR", (ftnlen)6);
	return 0;
    }
    chkout_("MKVARR", (ftnlen)6);
    return 0;
} /* mkvarr_ */

