/* mkgrid.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static logical c_true = TRUE_;
static doublereal c_b64 = 0.;
static logical c_false = FALSE_;
static integer c__5 = 5;

/* $Procedure MKGRID ( MKDSK: create plate set from height grid ) */
/* Subroutine */ int mkgrid_(char *infile__, integer *plttyp, char *aunits, 
	char *dunits, integer *corsys, doublereal *corpar, integer *maxnv, 
	integer *maxnp, integer *nv, doublereal *verts, integer *np, integer *
	plates, ftnlen infile_len, ftnlen aunits_len, ftnlen dunits_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer nmid;
    logical wrap;
    extern /* Subroutine */ int zzcapplt_(integer *, logical *, logical *, 
	    integer *, integer *, integer *, integer *), zzgrdplt_(integer *, 
	    integer *, logical *, integer *, integer *);
    integer b, i__, j;
    extern /* Subroutine */ int getg05_(integer *, logical *, logical *, 
	    logical *, logical *, logical *, logical *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), chkin_(char *, ftnlen), errch_(char *
	    , char *, ftnlen, ftnlen), vpack_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), errdp_(char *, doublereal *, ftnlen);
    integer ncols, reqnp, reqnv;
    extern doublereal vnorm_(doublereal *);
    integer nrows;
    extern logical failed_(void);
    doublereal hscale;
    logical mkncap, mkscap;
    doublereal refval;
    integer pltbas;
    doublereal lftcor, colstp, topcor;
    integer nnorth, polidx;
    logical leftrt, rowmaj;
    extern logical return_(void);
    integer nsouth;
    logical topdwn;
    extern /* Subroutine */ int chkout_(char *, ftnlen), setmsg_(char *, 
	    ftnlen), errint_(char *, integer *, ftnlen), sigerr_(char *, 
	    ftnlen), mkvarr_(char *, char *, char *, logical *, logical *, 
	    logical *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    doublereal rowstp, sum;

/* $ Abstract */

/*     Create a DSK type 2 plate set from a height grid provided */
/*     in a file. */

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
/*     INFILE     I   Name of input file. */
/*     PLTTYP     I   MKDSK input file format code. */
/*     AUNITS     I   Angular units. */
/*     DUNITS     I   Distance units. */
/*     CORSYS     I   Coordinate system. */
/*     CORPAR     I   Coordinate parameters. */
/*     MAXNV      I   Maximum number of vertices. */
/*     MAXNP      I   Maximum number of plates. */
/*     NV         O   Number of vertices. */
/*     VERTS      O   Vertex array. */
/*     NP         O   Number of plates. */
/*     PLATES     O   Plate array. */

/* $ Detailed_Input */

/*     INFILE     is the name of an input data file containing height */
/*                grid data. */

/*     PLTTYP     is the MKDSK code indicating the format of the input */
/*                data file. */

/*     AUNITS     is the name of the angular unit associated with the */
/*                grid coordinates, if the grid coordinate system is */
/*                latitudinal or planetodetic. AUNITS must be supported */
/*                by the SPICELIB routine CONVRT. */

/*     DUNITS     is the name of the distance unit associated with the */
/*                grid coordinates. DUNITS must be supported by the */
/*                SPICELIB routine CONVRT. */

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

/*     MAXNV      I   Maximum number of vertices to return. */

/*     MAXNP      I   Maximum number of plates to return. */


/* $ Detailed_Output */

/*     NV         is the number of vertices in VERTS. */

/*     VERTS      is an array of 3-dimensional vertices corresponding */
/*                to the height grid. */

/*                Units are km. */

/*     NP         is the number of plates in PLATES. */

/*     PLATES     is an array of plates representing a tessellation of */
/*                the height grid. */

/* $ Parameters */

/*     See the MKDSK include files */

/*        mkdsk.inc */
/*        mkdsk02.inc */

/*     and the DSK include file */

/*        dskdsc.inc */


/* $ Exceptions */

/*     1)  If either the row or column count is insufficient to define a */
/*         surface, the error SPICE(INVALIDCOUNT) is signaled. */

/*     2)  If longitude wrap is specified for a rectangular coordinate */
/*         system, the error SPICE(SPURIOUSFLAG) is signaled. */

/*     3)  If either the north or south polar cap flag is .TRUE., and */
/*         the coordinate system is rectangular, the error */
/*         SPICE(SPURIOUSFLAG) is signaled. */

/*     4)  If either the row or column step is not strictly positive, */
/*         the error SPICE(INVALIDSTEP) is signaled. */

/*     5)  If the height scale is not is not strictly positive, */
/*         the error SPICE(INVALIDSCALE) is signaled. */

/*     6)  If the coordinate system is latitudinal and the height scale */
/*         is negative, the error SPICE(INVALIDREFVAL) is signaled. */

/*     7)  If the number of vertices that must be created exceeds */
/*         MAXNV, the error SPICE(TOOMANYVERTICES) is signaled. */

/*     8)  If the number of plates that must be created exceeds */
/*         MAXNP, the error SPICE(TOOMANYPLATES) is signaled. */

/*     9)  If the input file format code is not recognized, the */
/*         error SPICE(NOTSUPPORTED) is signaled. */

/*    10)  If an error occurs while reading the input file, the */
/*         error will be diagnosed by routines in the call tree */
/*         of this routine. */

/*    11)  If an error occurs while processing the setup file, the error */
/*         will be diagnosed by routines in the call tree of this */
/*         routine. */

/*    12)  If an error occurs while converting the input data to a */
/*         vertex array, the error will be diagnosed by routines in the */
/*         call tree of this routine. */

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

/*     None. */

/* $ Examples */

/*     See usage in the MKDSK routine ZZWSEG02. */

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
    chkin_("MKGRID", (ftnlen)6);

/*     If we support the requested grid type, process the data. */
/*     The type is contained in the PLTTYP argument. */

    if (*plttyp == 5) {

/*        Fetch grid parameters from the kernel pool. */

	getg05_(corsys, &wrap, &mkncap, &mkscap, &rowmaj, &topdwn, &leftrt, &
		refval, &hscale, &ncols, &nrows, &lftcor, &topcor, &colstp, &
		rowstp);
	if (failed_()) {
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}

/*        For safety, check the parameters that could get us into real */
/*        trouble. */

	if (*corsys == 3) {
	    if (nrows < 2) {
		setmsg_("Number of rows was #; must have at least two rows t"
			"o create a grid using the rectangular coordinate sys"
			"tem.", (ftnlen)107);
		errint_("#", &nrows, (ftnlen)1);
		sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
		chkout_("MKGRID", (ftnlen)6);
		return 0;
	    }
	    if (wrap) {
		setmsg_("Longitude wrap is not applicable to the rectangular"
			" coordinate system.", (ftnlen)70);
		sigerr_("SPICE(SPURIOUSFLAG)", (ftnlen)19);
		chkout_("MKGRID", (ftnlen)6);
		return 0;
	    }
	    if (mkncap || mkscap) {
		setmsg_("Polar cap creation is not applicable to the rectang"
			"ular coordinate system.", (ftnlen)74);
		sigerr_("SPICE(SPURIOUSFLAG)", (ftnlen)19);
		chkout_("MKGRID", (ftnlen)6);
		return 0;
	    }
	} else {
	    if (mkncap || mkscap) {
		if (nrows < 1) {
		    setmsg_("Number of rows was #; must have at  least one r"
			    "ow to create a grid using the # coordinate syste"
			    "m when at least one polar cap is created.", (
			    ftnlen)136);
		    errint_("#", &nrows, (ftnlen)1);
		    if (*corsys == 1) {
			errch_("#", "latitudinal", (ftnlen)1, (ftnlen)11);
		    } else if (*corsys == 4) {
			errch_("#", "planetodetic", (ftnlen)1, (ftnlen)12);
		    } else {
			errint_("#", corsys, (ftnlen)1);
		    }
		    sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
		    chkout_("MKGRID", (ftnlen)6);
		    return 0;
		}
	    } else {
		if (nrows < 2) {
		    setmsg_("Number of rows was #; must have at least two ro"
			    "ws to create a grid using the # coordinate syste"
			    "m when at no polar caps are created.", (ftnlen)
			    131);
		    errint_("#", &nrows, (ftnlen)1);
		    if (*corsys == 1) {
			errch_("#", "latitudinal", (ftnlen)1, (ftnlen)11);
		    } else if (*corsys == 4) {
			errch_("#", "planetodetic", (ftnlen)1, (ftnlen)12);
		    } else {
			errint_("#", corsys, (ftnlen)1);
		    }
		    sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
		    chkout_("MKGRID", (ftnlen)6);
		    return 0;
		}
	    }
	}
	if (ncols < 2) {
	    setmsg_("Number of columns was #; must have at least two columns"
		    " to create a grid.", (ftnlen)73);
	    errint_("#", &ncols, (ftnlen)1);
	    sigerr_("SPICE(INVALIDCOUNT)", (ftnlen)19);
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}
	if (colstp <= 0.) {
	    setmsg_("Column step must be strictly positive but was #.", (
		    ftnlen)48);
	    errdp_("#", &colstp, (ftnlen)1);
	    sigerr_("SPICE(INVALIDSTEP)", (ftnlen)18);
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}
	if (rowstp <= 0.) {
	    setmsg_("Row step must be strictly positive but was #.", (ftnlen)
		    45);
	    errdp_("#", &rowstp, (ftnlen)1);
	    sigerr_("SPICE(INVALIDSTEP)", (ftnlen)18);
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}
	if (hscale <= 0.) {
	    setmsg_("Height scale must be strictly positive but was #.", (
		    ftnlen)49);
	    errdp_("#", &hscale, (ftnlen)1);
	    sigerr_("SPICE(INVALIDSCALE)", (ftnlen)19);
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}
	if (*corsys == 1) {
	    if (refval < 0.) {
		setmsg_("For latitudinal coordinates, the height reference v"
			"alue must be non-negative. It was #.", (ftnlen)87);
		errdp_("#", &refval, (ftnlen)1);
		sigerr_("SPICE(INVALIDREFVAL)", (ftnlen)20);
		chkout_("MKGRID", (ftnlen)6);
		return 0;
	    }
	}

/*        Let REQNV and REQNP be, respectively, the numbers of */
/*        vertices and plates we need to create. Make sure we can handle */
/*        these number. */

	*nv = nrows * ncols;
	reqnv = *nv;
	reqnp = (nrows - 1 << 1) * (ncols - 1);
	if (wrap) {
	    reqnp += nrows - 1 << 1;
	}
	if (mkncap) {
	    ++reqnv;
	    reqnp += ncols;
	    if (wrap) {
		++reqnp;
	    }
	}
	if (mkscap) {
	    ++reqnv;
	    reqnp += ncols;
	    if (wrap) {
		++reqnp;
	    }
	}
	if (reqnv > *maxnv) {
	    setmsg_("The number of vertices that must be created is #. The m"
		    "aximum allowed number is #.", (ftnlen)82);
	    errint_("#", &reqnv, (ftnlen)1);
	    errint_("#", maxnv, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYVERTICES)", (ftnlen)22);
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}
	if (reqnp > *maxnp) {
	    setmsg_("The number of plates that must be created is #. The max"
		    "imum allowed number is #.", (ftnlen)80);
	    errint_("#", &reqnp, (ftnlen)1);
	    errint_("#", maxnp, (ftnlen)1);
	    sigerr_("SPICE(TOOMANYPLATES)", (ftnlen)20);
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}

/*        Create vertices. If we're making a north polar cap, leave */
/*        room for it at the start of the vertex array. */

	if (mkncap) {
	    b = 2;
	    ++(*nv);
	} else {
	    b = 1;
	}
	mkvarr_(infile__, aunits, dunits, &rowmaj, &topdwn, &leftrt, corsys, 
		corpar, &refval, &hscale, &ncols, &nrows, &lftcor, &topcor, &
		colstp, &rowstp, maxnv, &verts[b * 3 - 3], infile_len, 
		aunits_len, dunits_len);
	if (failed_()) {
	    chkout_("MKGRID", (ftnlen)6);
	    return 0;
	}

/*        The output vertices have units of km. */


/*        Make plates. Fill in the polar vertices, if they're needed. */

/*        We create the plates in top-down order, so the polar caps */
/*        will be adjacent to the nearby non-polar plates. */

	pltbas = 1;
	nnorth = 0;
	nsouth = 0;
	*np = 0;
	if (mkncap) {
	    polidx = 1;
	    zzcapplt_(&ncols, &c_true, &wrap, &pltbas, &polidx, &nnorth, 
		    plates);
	    *np = nnorth;

/*           The north vertex magnitude is the average of the */
/*           magnitudes of the vertices in the top row. */

	    sum = 0.;
	    i__1 = ncols;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		j = b - 1 + i__;
		sum += vnorm_(&verts[j * 3 - 3]);
	    }
	    d__1 = sum / ncols;
	    vpack_(&c_b64, &c_b64, &d__1, verts);
	}

/*        Create the non-polar grid, if we have enough rows for it. */

	if (nrows > 1) {
	    zzgrdplt_(&nrows, &ncols, &wrap, &nmid, &plates[(*np + 1) * 3 - 3]
		    );
	    if (failed_()) {
		chkout_("MKGRID", (ftnlen)6);
		return 0;
	    }
	    *np += nmid;
	    if (mkncap) {

/*              Adjust the vertex indices in the plate set. */

		i__1 = *np;
		for (i__ = nnorth + 1; i__ <= i__1; ++i__) {
		    for (j = 1; j <= 3; ++j) {
			++plates[j + i__ * 3 - 4];
		    }
		}
	    }
	} else {

/*           We need to make at least one polar cap, or we won't have */
/*           any output. */

	    if (! mkncap && ! mkscap) {
		setmsg_("We have only one row of data in the input grid, and"
			" no polar caps were commanded to be constructed. Thi"
			"s gives us an empty output plate set.", (ftnlen)140);
		sigerr_("SPICE(NOPLATES)", (ftnlen)15);
		chkout_("MKGRID", (ftnlen)6);
		return 0;
	    }
	}
	if (mkscap) {
	    polidx = *nv + 1;
	    pltbas = b - 1 + (nrows - 1) * ncols;
	    zzcapplt_(&ncols, &c_false, &wrap, &pltbas, &polidx, &nsouth, &
		    plates[(*np + 1) * 3 - 3]);
	    *np += nsouth;

/*           The south vertex magnitude is the average of the */
/*           magnitudes of the vertices in the bottom row. */

	    sum = 0.;
	    i__1 = ncols;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		j = pltbas + i__;
		sum += vnorm_(&verts[j * 3 - 3]);
	    }
	    d__1 = -sum / ncols;
	    vpack_(&c_b64, &c_b64, &d__1, &verts[polidx * 3 - 3]);
	    ++(*nv);
	}
    } else {
	setmsg_("Input data format type is #; only type # is supported.", (
		ftnlen)54);
	errint_("#", plttyp, (ftnlen)1);
	errint_("#", &c__5, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("MKGRID", (ftnlen)6);
	return 0;
    }
    chkout_("MKGRID", (ftnlen)6);
    return 0;
} /* mkgrid_ */

