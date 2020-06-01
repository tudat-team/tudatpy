/* zzcapplt.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure ZZCAPPLT ( Make polar cap plates ) */
/* Subroutine */ int zzcapplt_(integer *ncols, logical *north, logical *wrap, 
	integer *basidx, integer *polidx, integer *np, integer *plates)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer bl, br, ub, tl, tr;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), errint_(char *, integer *, 
	    ftnlen);
    extern logical return_(void);

/* $ Abstract */

/*     Generate a set of plates constituting a polar cap. */

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
/*     NCOLS      I   Number of columns of data. */
/*     NORTH      I   Pole selector: .TRUE. if pole is north. */
/*     WRAP       I   Longitude wrap flag: .TRUE. enables wrapping. */
/*     BASIDX     I   Base index of vertices in row next to pole. */
/*     POLIDX     I   Index of pole vertex. */
/*     NP         O   Number of plates. */
/*     PLATES     O   Plate array. */

/* $ Detailed_Input */

/*     NCOLS          is the number of columns of data in the vertex */
/*                    grid for which a cap is to be constructed. */

/*     NORTH          is a logical flag indicating whether the cap */
/*                    is to be constructed for the north or south */
/*                    pole. A value of .TRUE. indicates north; .FALSE. */
/*                    indicates south. */

/*     WRAP           is a logical flag indicating whether longitude */
/*                    wrapping is to be performed. Longitude wrapping */
/*                    entails creating an extra plate to connect the */
/*                    east and west edges of the polar cap. */

/*     BASIDX         is the base index of the set of vertices */
/*                    constituting the vertex row adjacent to the */
/*                    designated pole. The indices of the vertices */
/*                    constituting this row range from */

/*                       BASIDX + 1 */

/*                    to */

/*                       BASIDX + NCOLS */

/*                    Regardless of the format of the input data file, */
/*                    the vertex grid is presumed to have row major */
/*                    order, and have top-down, left-to-right indexing. */

/*     POLIDX         is the index of the pole vertex. */

/* $ Detailed_Output */

/*     NP             is the output plate count. */

/*     PLATES         is a set of plates constituting a "cap" for the */
/*                    specified pole. */

/*                    Each plate consists of three vertex indices. The */
/*                    elements of PLATES are stored in row major order, */
/*                    regardless of the organization of the input data */
/*                    set. */

/*                    The dimension of PLATES is ( 3, NP ). */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the input column count is not at least two, the error */
/*        SPICE(INVALIDDIMENSION) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine is used to extend plate sets to cover */
/*     the latitude range -90 : 90 degrees. Many data sets */
/*     cover latitude ranges that exclude one or both poles. */
/*     It can be convenient to fill in the gaps at the poles */
/*     to ensure that certain operations, such as ray-surface */
/*     intercept computations performed to create latitude- */
/*     longitude grids, succeed for the entire surface. */

/*     Vertices in the row adjacent to the pole are assumed to */
/*     have left-to-right order. Vertex order determines the */
/*     direction of the output plates' outward normal vectors. */

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

    if (return_()) {
	return 0;
    }
    chkin_("ZZCAPPLT", (ftnlen)8);

/*     Check column dimensions. */

    if (*ncols < 2) {
	setmsg_("Grid must have at least two columns but NCOLSS is #.", (
		ftnlen)52);
	errint_("#", ncols, (ftnlen)1);
	sigerr_("SPICE(INVALIDDIMENSION)", (ftnlen)23);
	chkout_("ZZCAPPLT", (ftnlen)8);
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
    i__1 = ub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*north) {

/*           Create plates for a north polar cap. */

/*           Set the vertex index of the north pole. */

	    tl = *polidx;

/*           Longitude increases with increasing column index. */

	    if (*wrap && i__ == ub) {

/*              Form a plate by connecting the right edge */
/*              of the surface to the left edge. */

		bl = *basidx + ub;
		br = *basidx + 1;
	    } else {
		bl = *basidx + i__;
		br = bl + 1;
	    }

/*           Create the current plate. */

	    ++(*np);
	    plates[*np * 3 - 3] = tl;
	    plates[*np * 3 - 2] = bl;
	    plates[*np * 3 - 1] = br;
	} else {

/*           Create plates for a south polar cap. */

/*           Set the vertex index of the south pole. */

	    bl = *polidx;

/*           Longitude increases with increasing column index. */

	    if (*wrap && i__ == ub) {

/*              Form a plate by connecting the right edge */
/*              of the surface to the left edge. */

		tl = *basidx + ub;
		tr = *basidx + 1;
	    } else {
		tl = *basidx + i__;
		tr = tl + 1;
	    }

/*           Create the current plate. */

	    ++(*np);
	    plates[*np * 3 - 3] = bl;
	    plates[*np * 3 - 2] = tr;
	    plates[*np * 3 - 1] = tl;
	}
    }

/*     The plate and vertex counts and arrays have been */
/*     assigned. */

    chkout_("ZZCAPPLT", (ftnlen)8);
    return 0;
} /* zzcapplt_ */

