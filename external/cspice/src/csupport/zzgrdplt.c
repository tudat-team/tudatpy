/* zzgrdplt.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure ZZGRDPLT ( Create grid of plates ) */
/* Subroutine */ int zzgrdplt_(integer *nrows, integer *ncols, logical *wrap, 
	integer *np, integer *plates)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer bl, br, ub, tl, tr;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), errint_(char *, integer *, 
	    ftnlen);
    extern logical return_(void);

/* $ Abstract */

/*     Generate a set of plates covering the surface spanned by a */
/*     rectangular grid of vertices. */

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
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     NROWS      I   Number of rows of data. */
/*     NCOLS      I   Number of columns of data. */
/*     WRAP       I   Longitude wrap flag: .TRUE. enables wrapping. */
/*     NP         O   Number of plates. */
/*     PLATES     O   Plate array. */

/* $ Detailed_Input */

/*     NROWS, */
/*     NCOLS          are, respectively, the numbers of rows and columns */
/*                    of data in the vertex grid for which a plate set */
/*                    is to be constructed. */

/*     WRAP           is a logical flag indicating whether longitude */
/*                    wrapping is to be performed. Longitude wrapping */
/*                    entails creating an extra column of plates to */
/*                    connect the east and west edges of the vertex */
/*                    grid. */

/* $ Detailed_Output */

/*     NP             is the output plate count. */

/*     PLATES         is a set of plates covering the region spanned */
/*                    by the input vertex grid. Every region having */
/*                    four adjacent vertices at its corners is covered */
/*                    by two plates. */

/*                    Each plate consists of three vertex indices. */

/*                    The indices of the vertices constituting the grid */
/*                    are presumed to range from */

/*                       1  to  NROWS*NCOLS */

/*                    Regardless of the format of the input data file, */
/*                    the vertex grid is presumed to have row major */
/*                    order, and have top-down, left-to-right indexing. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If either the input column count or input row count */
/*        is not at least two, the error SPICE(INVALIDDIMENSION) is */
/*        signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     To avoid portability problems that can arise from */
/*     basing plate-vertex association on computed quantities, */
/*     all plates are constructed according to the following, */
/*     fixed pattern: */


/*        vertex(I,J)      vertex(I,J+1) */
/*            +----------------+ */
/*            |\               | */
/*            | \              | */
/*            |  \             | */
/*            |   \            | */
/*            |    \           | */
/*            |     \          | */
/*            |      \         | */
/*            |       \        | */
/*            |        \       | */
/*            |         \      | */
/*            |          \     | */
/*            |           \    | */
/*            |            \   | */
/*            |             \  | */
/*            |              \ | */
/*            |               \| */
/*            +----------------+ */
/*        vertex(I+1,J)      vertex(I+1,J+1) */



/* $ Examples */

/*     See usage in the MKDSK routine MKGRID. */

/* $ Restrictions */

/*     1) For use only within program MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 1.0.1, 30-APR-2014 (NJB) */

/*        Corrected some comment typos. */

/* -    MKDSK Version 1.0.0, 25-JUL-2011 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local variables */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    }
    chkin_("ZZGRDPLT", (ftnlen)8);

/*     Check row and column dimensions. */

    if (*nrows < 2) {
	setmsg_("Grid must have at least two rows but NROWS is #.", (ftnlen)
		48);
	errint_("#", nrows, (ftnlen)1);
	sigerr_("SPICE(INVALIDDIMENSION)", (ftnlen)23);
	chkout_("ZZGRDPLT", (ftnlen)8);
	return 0;
    }
    if (*ncols < 2) {
	setmsg_("Grid must have at least two columns but NCOLSS is #.", (
		ftnlen)52);
	errint_("#", ncols, (ftnlen)1);
	sigerr_("SPICE(INVALIDDIMENSION)", (ftnlen)23);
	chkout_("ZZGRDPLT", (ftnlen)8);
	return 0;
    }

/*     Set the upper bound on the column loop. If longitude */
/*     wrapping is turned on, the final column connects to */
/*     the first column. */

    if (*wrap) {
	ub = *ncols;
    } else {
	ub = *ncols - 1;
    }

/*     Connect the vertices to generate plates. */

    *np = 0;
    i__1 = *nrows - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ub;
	for (j = 1; j <= i__2; ++j) {

/*           Since the vertices are stored in a 3xN */
/*           array, it will be convenient to use aliases for */
/*           their indices---this cuts down on use of complicated */
/*           index expressions. For each square defined by */
/*           four neighboring vertices, we'll call the vertices */

/*              TL   "top left" */
/*              TR   "top right" */
/*              BL   "bottom left" */
/*              BR   "bottom right" */

/*           Recall that the input pixel grid has dimensions */

/*              NROWS x NCOLS */

/*           The top row is at the highest latitude. */

/*           The leftmost column corresponds to the west */
/*           boundary of the region. */

/*           The top left vertex corresponds to pixel (I,J). */

	    if (*wrap && j == ub) {

/*              Connect the right edge of the grid to the left edge. */

		tl = i__ * *ncols;
		tr = tl - *ncols + 1;
		bl = tl + *ncols;
		br = tr + *ncols;
	    } else {

/*              This is the normal case: the column at index */
/*              J is connected to the column at index J+1. */

		tl = (i__ - 1) * *ncols + j;
		tr = tl + 1;
		bl = tl + *ncols;
		br = bl + 1;
	    }

/*           For each square defined by neighboring pixel centers, */
/*           we must represent the corresponding surface by a pair */
/*           of plates. We have two choices for the diagonal */
/*           common edge connecting these plates: descending or */
/*           ascending to the right. */

/*           We choose the descending diagonal. */

/*           The vertex assignment must be positively */
/*           oriented about the outward normal direction. */

	    ++(*np);
	    plates[*np * 3 - 3] = bl;
	    plates[*np * 3 - 2] = br;
	    plates[*np * 3 - 1] = tl;
	    ++(*np);
	    plates[*np * 3 - 3] = tl;
	    plates[*np * 3 - 2] = br;
	    plates[*np * 3 - 1] = tr;
	}
    }

/*     The plate and vertex counts and arrays have been */
/*     assigned. */

    chkout_("ZZGRDPLT", (ftnlen)8);
    return 0;
} /* zzgrdplt_ */

