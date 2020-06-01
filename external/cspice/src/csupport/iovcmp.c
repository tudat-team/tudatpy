/* iovcmp.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure IOVCMP ( Inverse order vector with compressed range ) */
/* Subroutine */ int iovcmp_(doublereal *darray, integer *ndim, integer *
	iorder, integer *invord, integer *rngmax)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int orderd_(doublereal *, integer *, integer *);
    integer nupred;

/* $ Abstract */

/*     Create an inverse order vector having a compressed range. */
/*     In this vector, the Ith element is the index of DARRAY(I) */
/*     in the array that would be produced by sorting DARRAY and */
/*     removing duplicates. */

/*     Produce this result in O(N) time. */

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
/*     DARRAY     I   Double precision array. */
/*     NDIM       I   Size of DARRAY. */
/*     IORDER     O   Order vector for DARRAY. */
/*     INVORD     O   Compressed inverse order vector. */
/*     RNGMAX     O   Maximum value in range of INVORD. */

/* $ Detailed_Input */

/*     DARRAY     is an array of double precision numbers. */

/*     NDIM       is the size of DARRAY. */

/* $ Detailed_Output */

/*     IORDER     is an order vector for DARRAY. */

/*     INVORD     is a compressed inverse order vector for DARRAY. */

/*                INVORD(I) is the position that DARRAY(I) would */
/*                have in an ordered set created by sorting DARRAY */
/*                and removing duplicates. */

/*                INVORD has size NDIM. Its elements belong to the */
/*                set {1, ..., RNGMAX}. */

/*     RNGMAX     is the largest element in INVORD. This value is the */
/*                number of distinct values in DARRAY. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     Error free. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine supports efficient determination of spatial */
/*     coverage gaps by DSKBRIEF. */

/* $ Examples */

/*     See usage in DSKBRIEF. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/*     DSKBRIEF Version 1.0.0, 05-OCT-2016 (NJB) */

/* -& */

/*     Local variables */


/*     First step: create an order vector for DARRAY. */

    orderd_(darray, ndim, iorder);

/*     Produce the corresponding inverse order vector. */

    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	invord[iorder[i__ - 1] - 1] = i__;
    }

/*     Step through the order vector, keeping track of the count of */
/*     unique predecessors, in the array that would be produced by */
/*     sorting DARRAY, of each element pointed to by an element of the */
/*     order vector. */

/*     The element of DARRAY at index IORDER(1) has no predecessors, */
/*     and the element INVORD( IORDER(1) ) is already correct. So */
/*     we start at the second element of IORDER (if it exists). */

/*     Initialize NUPRED to the number of unique predecessors of */
/*     the first value. */

    nupred = 0;
    i__1 = *ndim;
    for (i__ = 2; i__ <= i__1; ++i__) {

/*        At this point, NUPRED is the number of unique predecessors of */
/*        DARRAY(I). I is greater than or equal to 2. */

	if (darray[iorder[i__ - 1] - 1] > darray[iorder[i__ - 2] - 1]) {

/*           DARRAY( IORDER(I) ) is strictly greater than, and hence not */
/*           a copy of, its predecessor. It has NUPRED + 1 unique */
/*           predecessors in the array produced by sorting DARRAY. */

	    ++nupred;

/*           The position of DARRAY( IORDER(I) ) in the sorted, */
/*           compressed set derived from DARRAY is one more than the */
/*           count of unique predecessors. */

	    invord[iorder[i__ - 1] - 1] = nupred + 1;
	} else {

/*           DARRAY( IORDER(I) ) is a duplicate. Its position in the */
/*           sorted, compressed array derived from DARRAY is the same */
/*           as that of DARRAY( IORDER(I-1) ). */

	    invord[iorder[i__ - 1] - 1] = invord[iorder[i__ - 2] - 1];
	}
    }

/*     INVORD has been updated so that its elements belong to the */
/*     set { 1 : NUPRED+1 }. */

/*     Set the maximum range value of INVORD. */

    *rngmax = nupred + 1;
    return 0;
} /* iovcmp_ */

