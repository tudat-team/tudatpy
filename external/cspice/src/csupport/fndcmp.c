/* fndcmp.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;

/* $Procedure FNDCMP ( DSKBRIEF, find rectangular components ) */
/* Subroutine */ int fndcmp_(integer *nrows, integer *ncols, logical *value, 
	integer *maxn, logical *grid, integer *vset, integer *mrkset, integer 
	*tmpset, integer *ncomp, integer *minpxx, integer *maxpxx, integer *
	minpxy, integer *maxpxy)
{
    /* System generated locals */
    integer grid_dim1, grid_dim2, grid_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    integer i__, j;
    extern /* Subroutine */ int diffi_(integer *, integer *, integer *);
    extern integer cardi_(integer *);
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    logical found;
    extern /* Subroutine */ int copyi_(integer *, integer *);
    integer id;
    extern logical failed_(void);
    extern /* Subroutine */ int scardi_(integer *, integer *);
    integer remain;
    extern /* Subroutine */ int appndi_(integer *, integer *);
    integer mincol, maxcol;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen);
    integer colsiz;
    extern /* Subroutine */ int errint_(char *, integer *, ftnlen);
    integer minrow, maxrow;
    extern logical return_(void);
    integer col, row;

/* $ Abstract */

/*     Find rectangular components in a 2-d logical grid. */

/*     The result is a list of rectangles. For each output rectangle, */
/*     the grid value is equal to the input VALUE at each pixel in that */
/*     rectangle. The rectangles are described by the bounds of their X */
/*     and Y pixel coordinates. */

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
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     NROWS      I   Number of rows in the input grid. */
/*     NCOLS      I   Number of columns in the input grid. */
/*     VALUE      I   Value marking pixels to aggregate. */
/*     MAXN       I   Maximum number of components. */
/*     GRID      I-O  Input grid. */
/*     VSET      I-O  Set of pixels set to VALUE. */
/*     MRKSET    I-O  Set of marked pixels. */
/*     TMPSET    I-O  Workspace set. */
/*     NCOMP      O   Number of rectangular components. */
/*     MINPXX     O   Minimum X coordinates of components. */
/*     MAXPXX     O   Maximum X coordinates of components. */
/*     MINPXY     O   Minimum Y coordinates of components. */
/*     MAXPXY     O   Maximum Y coordinates of components. */

/* $ Detailed_Input */

/*     NROWS      is the number of rows in the input pixel grid. */

/*     NCOLS      is number of columns in the input pixel grid. */

/*     VALUE      is the logical value used in the pixel grid to mark */
/*                pixels that are to be grouped into rectangular */
/*                components. */

/*     MAXN       is the maximum number of rectangular components that */
/*                can be returned. */

/*     GRID       is, on input, a pixel grid. GRID is an array of */
/*                logical values; the dimensions of GRID are NROWS x */
/*                NCOLS. */

/*                Caution: GRID is modified by this routine. */

/*     VSET       is, on input, an initialized integer set. VSET will */
/*                be filled with the 1-dimensional indices of pixels */
/*                in GRID that are set to the value VALUE. */

/*     MRKSET     is, on input, an initialized integer set. MRKSET is a */
/*                workspace variable: it will be filled with the */
/*                1-dimensional indices of pixels in GRID that belong to */
/*                a rectangular component that under construction. */

/*     TMPSET     is, on input, an initialized integer set. TMPSET is a */
/*                temporary variable. */

/* $ Detailed_Output */

/*     GRID       is an overwritten version of the input pixel grid. */
/*                Pixels belonging to components are marked by negating */
/*                those pixels' values. */

/*     VSET       is a modified integer set. The contents of VSET are */
/*                undefined. */

/*     MRKSET     is a modified integer set. The contents of MRKSET are */
/*                undefined. */

/*     TMPSET     is a modified integer set. The contents of TMPSET are */
/*                undefined. */

/*     NCOMP      is the number of rectangular components making up */
/*                the set of marked pixels. */

/*     MINPXX     is an array of lower bounds of the X coordinates of */
/*                the rectangular components. The Ith element of MINPXX */
/*                is the lower X bound of the Ith component. */

/*     MAXPXX     is an array of upper bounds of the X coordinates of */
/*                the rectangular components. The Ith element of MAXPXX */
/*                is the upper X bound of the Ith component. */

/*     MINPXY     is an array of lower bounds of the Y coordinates of */
/*                the rectangular components. The Ith element of MINPXY */
/*                is the lower Y bound of the Ith component. */

/*     MAXPXY     is an array of upper bounds of the Y coordinates of */
/*                the rectangular components. The Ith element of MAXPXY */
/*                is the upper Y bound of the Ith component. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If MAXN is smaller than the number of rectangular */
/*         components found, the error SPICE(ARRAYTOOSMALL) is */
/*         signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine supports determination of DSK spatial coverage gaps. */
/*     It can just as easily be used to determine regions of spatial */
/*     coverage. */

/*     The marked pixels in the input pixel grid are grouped by this */
/*     routine into disjoint, rectangular sets called "components." The */
/*     output coordinate bounds returned by this routine define the */
/*     boundaries of these components. */

/*     This routine is designed to operate efficiently. While */
/*     constructing a component, it accumulates pixels in increasing */
/*     order of the pixels' 1-dimensional indices. The allows the */
/*     routine to avoid sorting the indices in order to make the */
/*     component's pixel indices into a set. However, this algorithm */
/*     constrains the way components are identified: their pixels are */
/*     accumulated in column-major order. */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    DSKBRIEF Version 1.0.0, 30-JAN-2017 (NJB) */

/*        Original version 03-SEP-2016 (NJB) */
/* -& */

/*     SPICELIB functions */


/*     Local variables */

    /* Parameter adjustments */
    grid_dim1 = *nrows;
    grid_dim2 = *ncols;
    grid_offset = grid_dim1 + 1;

    /* Function Body */
    if (return_()) {
	return 0;
    }
    chkin_("FNDCMP", (ftnlen)6);
    scardi_(&c__0, vset);
    scardi_(&c__0, mrkset);
    scardi_(&c__0, tmpset);

/*     First step: make a pass through the grid, and store the ID */
/*     of each pixel matching VALUE. */

/*     Proceed in column-major order. */

    i__1 = *ncols;
    for (col = 1; col <= i__1; ++col) {
	i__2 = *nrows;
	for (row = 1; row <= i__2; ++row) {
	    if (grid[(i__3 = row + col * grid_dim1 - grid_offset) < grid_dim1 
		    * grid_dim2 && 0 <= i__3 ? i__3 : s_rnge("grid", i__3, 
		    "fndcmp_", (ftnlen)256)] == *value) {

/*              It's a match. */

/*              ID is the one-dimensional index of the current pixel. */

		id = (col - 1) * *nrows + row;

/*              Since we're traversing the grid in increasing ID */
/*              order, the elements of VSET will automatically be */
/*              in increasing order. We don't need to sort them. */

		appndi_(&id, vset);
		if (failed_()) {
		    chkout_("FNDCMP", (ftnlen)6);
		    return 0;
		}
	    }
	}
    }

/*     Now find rectangular sets of pixels equal to VALUE. */

    scardi_(&c__0, mrkset);
    remain = cardi_(vset);
    *ncomp = 0;
    while(remain > 0) {

/*        Get the row and column coordinates of the first pixel in VSET. */

	id = vset[6];
	col = (id - 1) / *nrows + 1;
	row = id - (col - 1) * *nrows;
	minrow = row;
	mincol = col;

/*        We'll extend the component in the direction of higher row */
/*        indices as far as possible, then in the direction of higher */
/*        column indices  as far as possible. The reason for this is */
/*        that we want to accumulate pixels in increasing order of ID. */

	maxrow = *nrows;
	maxcol = col;
	found = TRUE_;
	while(col <= *ncols && found) {

/*           COL is a valid column number at the top of the loop. */
/*           We increment COL at the bottom of the loop. */

/*           Initialize ROW for a pass through the current column. */

	    row = minrow - 1;

/*           Caution: the value of MAXROW in the loop termination */
/*           condition below changes during loop execution! The */
/*           value is NROWS on the first pass; then it changes */
/*           to the maximum row number of the first column of the */
/*           component. */

	    while(row < maxrow && found) {

/*              Note the .LT. operator in the loop termination */
/*              condition. We increment ROW at the top of the */
/*              loop, so the value of ROW is correct after */
/*              loop termination. */

		++row;
		found = grid[(i__1 = row + col * grid_dim1 - grid_offset) < 
			grid_dim1 * grid_dim2 && 0 <= i__1 ? i__1 : s_rnge(
			"grid", i__1, "fndcmp_", (ftnlen)338)] == *value;
	    }
	    if (col == mincol) {

/*              The index of the last row that matched becomes the */
/*              maximum row index of this component. */

		if (found) {

/*                 The row index reached NROWS. */

		    maxrow = *nrows;
		} else {

/*                 The last matching row was the one preceding ROW. */

		    maxrow = row - 1;

/*                 Set FOUND to .TRUE. so we'll go on to look at the */
/*                 next column. */

		    found = TRUE_;
		}

/*              Now we know the size of the columns of the component. */

		colsiz = maxrow - minrow + 1;

/*              Always go on to look at the next column. FOUND is */
/*              .TRUE. at this point. */

		maxcol = col;
	    } else {

/*              After we process the first column of the component, */
/*              we don't adjust MAXROW again. It's set to the highest */
/*              row number of the first column of the component. */

		if (! found) {

/*                 The current column fails to match in some row of the */
/*                 current column. This column can't be included in the */
/*                 component. */

		    maxcol = col - 1;
		} else {

/*                 The current column matches from row indices MINROW */
/*                 to MAXROW. This column is part of the component. */

		    maxcol = col;

/*                 Set FOUND to .TRUE. so we'll go on to look at the */
/*                 next column. */

		    found = TRUE_;
		}
	    }
	    if (found) {

/*              We've found a column of matching pixels in the */
/*              current component. */

/*              Add the pixels of the column to the marked set. */

/*              Let ID be the ID of the first pixel of the column. */

		id = (maxcol - 1) * *nrows + minrow;
		i__1 = colsiz;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    j = id - 1 + i__;
		    appndi_(&j, mrkset);
		    if (failed_()) {
			chkout_("FNDCMP", (ftnlen)6);
			return 0;
		    }

/*                 Fill in the matching pixels so they won't */
/*                 match again. */

		    grid[(i__2 = minrow - 1 + i__ + maxcol * grid_dim1 - 
			    grid_offset) < grid_dim1 * grid_dim2 && 0 <= i__2 
			    ? i__2 : s_rnge("grid", i__2, "fndcmp_", (ftnlen)
			    432)] = ! (*value);
		}

/*              Note that we've added the IDs to MRKSET in increasing */
/*              order, so MRKSET remains a set. We don't need to sort */
/*              its contents. */

/*              Prepare to examine the next column. */

		remain -= colsiz;
		++col;
	    }
	}

/*        We've finished building a component. */

	++(*ncomp);

/*        Update VSET: subtract the pixels of the new component. */

/*        Note that subtracting one set from another should be an */
/*        efficient process, if done correctly. We trust DIFFI to manage */
/*        this. */

	diffi_(vset, mrkset, tmpset);
	copyi_(tmpset, vset);
	scardi_(&c__0, mrkset);
	scardi_(&c__0, tmpset);
	if (failed_()) {
	    chkout_("FNDCMP", (ftnlen)6);
	    return 0;
	}
/*        The bounds of the component we just found are given by */

/*           MINROW, MAXROW, MINCOL, MAXCOL */

	if (*ncomp <= *maxn) {
	    minpxx[*ncomp - 1] = mincol;
	    maxpxx[*ncomp - 1] = maxcol;
	    minpxy[*ncomp - 1] = minrow;
	    maxpxy[*ncomp - 1] = maxrow;
	} else {

/*           We're out of room. */

	    setmsg_("There are more output rectangles than can be accommodat"
		    "ed in the output rectangle boundary arrays. So far, # co"
		    "mponents have been found; the maximum supported number i"
		    "s #.", (ftnlen)171);
	    errint_("#", ncomp, (ftnlen)1);
	    errint_("#", maxn, (ftnlen)1);
	    sigerr_("SPICE(ARRAYTOOSMALL)", (ftnlen)20);
	    chkout_("FNDCMP", (ftnlen)6);
	    return 0;
	}
    }
    chkout_("FNDCMP", (ftnlen)6);
    return 0;
} /* fndcmp_ */

