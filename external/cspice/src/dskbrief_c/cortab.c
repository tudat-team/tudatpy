/* cortab.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* $Procedure CORTAB ( DSKBRIEF, display coordinate table ) */
/* Subroutine */ int cortab_(integer *n, char *labels, integer *start1, 
	integer *nsig, integer *ncols, doublereal *values, integer *starts, 
	char *table, ftnlen labels_len, ftnlen table_len)
{
    /* System generated locals */
    integer values_dim1, values_dim2, values_offset, starts_dim1, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1;
    icilist ici__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rnge(char *, integer, char *, integer);
    double d_lg10(doublereal *);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , i_indx(char *, char *, ftnlen, ftnlen), s_cmp(char *, char *, 
	    ftnlen, ftnlen);

    /* Local variables */
    integer bcol, maxa, dpix, maxp;
    doublereal logv;
    integer maxt, i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen), errch_(char *, char *,
	     ftnlen, ftnlen);
    doublereal dpval;
    extern /* Subroutine */ int repmi_(char *, char *, integer *, char *, 
	    ftnlen, ftnlen, ftnlen);
    integer colix, we, wf, wi, wt;
    extern /* Subroutine */ int shiftc_(char *, char *, integer *, char *, 
	    char *, ftnlen, ftnlen, ftnlen, ftnlen), sigerr_(char *, ftnlen), 
	    prefix_(char *, integer *, char *, ftnlen, ftnlen), chkout_(char *
	    , ftnlen), setmsg_(char *, ftnlen), errint_(char *, integer *, 
	    ftnlen);
    extern logical return_(void);
    integer adj;
    char fmt1[30];
    integer pos1, pos2;

/* $ Abstract */

/*     Display a table of coordinates or other double precision */
/*     values. */

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
/*     N          I   is the number of table rows. */
/*     LABELS     I   is an array of row labels. */
/*     START1     I   is the position of the first data column. */
/*     NSIG       I   is the number of significant digits to display. */
/*     NCOLS      I   is the number of data columns. */
/*     VALUES     I   is an array of data values. */
/*     STARTS     O   is an array of data column positions. */
/*     TABLE      O   is the output table. */

/* $ Detailed_Input */

/*     N          is the number of rows in the table to be created. */

/*     LABELS     is an array of labels for the output rows. The labels */
/*                are located at the left side of the table. */

/*     START1     is the character position of the start of the first */
/*                data column in the table. */

/*     NSIG       is the number of significant digits to display for */
/*                the output values. NSIG ranges from 7 to 16. */

/*     NCOLS      is the number of data columns in the table. */

/*     VALUES     is an array of double precision values to be placed in */
/*                the table. VALUES has dimensions */

/*                   NCOLS x N */

/* $ Detailed_Output */

/*     STARTS     is an array containing the starting character positions */
/*                of the output data columns. The first element of STARTS */
/*                is START1. */

/*     TABLE      is the output table. The Ith row of the table starts */
/*                with */

/*                   LABELS(I) */

/*                The values */

/*                   VALUES(1,I), ..., VALUES(NCOLS,I) */

/*                follow, with the Jth value starting at index STARTS(J). */

/*                Values having magnitudes in the range */

/*                   [1, 1e6) */

/*                as well as zero, are expressed in fixed-point format. */
/*                Non-zero values having magnitudes outside this range */
/*                are expressed in scientific notation. */

/*                In each data column, the decimal points of the values */
/*                are aligned. */

/* $ Parameters */

/*     ANGMRG     See the description in dsktol.inc. */

/* $ Exceptions */

/*     None. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine supports coverage gap display. */

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

/*        Original version 18-NOV-2016 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */

    /* Parameter adjustments */
    starts_dim1 = *ncols;
    values_dim1 = *ncols;
    values_dim2 = *n;
    values_offset = values_dim1 + 1;

    /* Function Body */
    if (return_()) {
	return 0;
    }
    chkin_("CORTAB", (ftnlen)6);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(table + (i__ - 1) * table_len, labels + (i__ - 1) * labels_len,
		 table_len, labels_len);
    }

/*     BCOL is the character index of the first column. */

    starts[(i__1 = 0) < starts_dim1 ? i__1 : s_rnge("starts", i__1, "cortab_",
	     (ftnlen)210)] = *start1;
    bcol = starts[(i__1 = 0) < starts_dim1 ? i__1 : s_rnge("starts", i__1, 
	    "cortab_", (ftnlen)211)];
    i__1 = *ncols;
    for (colix = 1; colix <= i__1; ++colix) {

/*        Find maximum character counts of integer and fractional parts */
/*        of the values in the current column. */

/*           MAXA is the maximum adjustment value for all rows. */

/*           MAXP is the maximum decimal point index for all rows. */
/*           The index is relative to the first character of the row. */

/*           MAXT is the maximum total field width for any value in */
/*           the current column, which has index COLIX. */

	maxa = 0;
	maxp = 0;
	maxt = 0;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dpval = values[(i__3 = colix + i__ * values_dim1 - values_offset) 
		    < values_dim1 * values_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "values", i__3, "cortab_", (ftnlen)233)];
	    if (abs(dpval) >= 1. && abs(dpval) < 1e6) {

/*              We can represent the value in fixed-point format. */

/*              Let WI be the width of the integer part; WF is the width */
/*              of the fractional part. WT is the total width of the */
/*              string representing the value. */

		d__1 = abs(dpval);
		logv = d_lg10(&d__1);

/*              Include room for the sign in WI. */

		wi = (integer) logv + 1;
/* Computing MAX */
		i__3 = 0, i__4 = *nsig - wi;
		wf = max(i__3,i__4);
		wt = wi + wf + 2;
		s_copy(fmt1, "(F@W.@F)", (ftnlen)30, (ftnlen)8);
		repmi_(fmt1, "@W", &wt, fmt1, (ftnlen)30, (ftnlen)2, (ftnlen)
			30);
		repmi_(fmt1, "@F", &wf, fmt1, (ftnlen)30, (ftnlen)2, (ftnlen)
			30);
	    } else if (dpval == 0.) {
		wi = 1;
/* Computing MAX */
		i__3 = 0, i__4 = *nsig - wi;
		wf = max(i__3,i__4);
		wt = wi + wf + 1;
		s_copy(fmt1, "(F@W.@F)", (ftnlen)30, (ftnlen)8);
		repmi_(fmt1, "@W", &wt, fmt1, (ftnlen)30, (ftnlen)2, (ftnlen)
			30);
		repmi_(fmt1, "@F", &wf, fmt1, (ftnlen)30, (ftnlen)2, (ftnlen)
			30);
	    } else {

/*              Use exponential notation. WF includes room for the */
/*              exponent, which includes the 'E' symbol, a sign, and */
/*              three digits. */

		we = 5;
		wi = 1;
/* Computing MAX */
		i__3 = 0, i__4 = *nsig - wi;
		wf = max(i__3,i__4);
		wt = wi + wf + 2 + we;
		s_copy(fmt1, "(1PE@W.@F)", (ftnlen)30, (ftnlen)10);
		repmi_(fmt1, "@W", &wt, fmt1, (ftnlen)30, (ftnlen)2, (ftnlen)
			30);
		repmi_(fmt1, "@F", &wf, fmt1, (ftnlen)30, (ftnlen)2, (ftnlen)
			30);
	    }

/*           Write the value to the output table, starting at */
/*           index BCOL. */

	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = table_len - (bcol - 1);
	    ici__1.iciunit = table + ((i__ - 1) * table_len + (bcol - 1));
	    ici__1.icifmt = fmt1;
	    s_wsfi(&ici__1);
	    do_fio(&c__1, (char *)&dpval, (ftnlen)sizeof(doublereal));
	    e_wsfi();

/*           Find the offset of the decimal point within the string */
/*           value just written. */

	    dpix = i_indx(table + ((i__ - 1) * table_len + (bcol - 1)), ".", 
		    table_len - (bcol - 1), (ftnlen)1);
	    maxp = max(maxp,dpix);
	    maxt = max(maxt,wt);
	    if (dpval == 0.) {

/*              Adjust the output string to have a leading zero if it's */
/*              needed. Fortran 77 says leading zeros are optional for */
/*              fixed-point representations of numbers having magnitude */
/*              less than 1. Zero qualifies. */

		if (dpix == 1) {
		    prefix_("0", &c__0, table + ((i__ - 1) * table_len + (
			    bcol - 1)), (ftnlen)1, bcol - 1 + wt - (bcol - 1))
			    ;
		} else if (dpix > 1) {
		    pos1 = dpix - 1;
		    pos2 = bcol - 1 + pos1;
		    if (s_cmp(table + ((i__ - 1) * table_len + (bcol - 1)), 
			    " ", pos2 - (bcol - 1), (ftnlen)1) == 0) {
			*(unsigned char *)&table[(i__ - 1) * table_len + (
				pos2 - 1)] = '0';
		    }
		} else {

/*                 Backstop case. */

		    setmsg_("No decimal point found in string representing z"
			    "ero. I = #; TABLE(I) = #.", (ftnlen)72);
		    errint_("#", &i__, (ftnlen)1);
		    errch_("#", table + (i__ - 1) * table_len, (ftnlen)1, 
			    table_len);
		    sigerr_("SPICE(BUG)", (ftnlen)10);
		    chkout_("CORTAB", (ftnlen)6);
		    return 0;
		}
	    }
	}

/*        Align the values just written. */

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dpix = i_indx(table + ((i__ - 1) * table_len + (bcol - 1)), ".", 
		    table_len - (bcol - 1), (ftnlen)1);
	    adj = maxp - dpix;
	    shiftc_(table + ((i__ - 1) * table_len + (bcol - 1)), "R", &adj, 
		    " ", table + ((i__ - 1) * table_len + (bcol - 1)), 
		    table_len - (bcol - 1), (ftnlen)1, (ftnlen)1, table_len - 
		    (bcol - 1));
	    maxa = max(maxa,adj);
	}
	bcol = bcol + maxt + maxa + 2;
	if (colix == 2) {
	    bcol += 4;
	}

/*        Set the start value for the next column. */

	if (colix < *ncols) {
	    starts[(i__2 = colix) < starts_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "starts", i__2, "cortab_", (ftnlen)366)] = bcol;
	}
    }
    chkout_("CORTAB", (ftnlen)6);
    return 0;
} /* cortab_ */

