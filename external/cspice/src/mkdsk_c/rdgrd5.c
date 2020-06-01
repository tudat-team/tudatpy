/* rdgrd5.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__100 = 100;

/* $Procedure RDGRD5 ( MKDSK: read from a type 5 height grid file ) */
/* Subroutine */ int rdgrd5_(char *infile__, integer *nmax, integer *n, 
	doublereal *values, logical *done, ftnlen infile_len)
{
    /* Initialized data */

    static logical eof = FALSE_;
    static integer ntk = 0;
    static logical newfil = TRUE_;
    static integer from = 0;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer), s_cmp(char *, char *, 
	    ftnlen, ftnlen);

    /* Local variables */
    char line[255];
    integer room, nxfr, i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen), rdnbl_(char *, char *,
	     logical *, ftnlen, ftnlen), errch_(char *, char *, ftnlen, 
	    ftnlen);
    extern logical failed_(void);
    integer to;
    static integer remain, lineno;
    extern /* Subroutine */ int nparsd_(char *, doublereal *, char *, integer 
	    *, ftnlen, ftnlen), sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), lparsm_(char *, char *, integer *, integer *, char *, 
	    ftnlen, ftnlen, ftnlen);
    char errmsg[320];
    extern /* Subroutine */ int setmsg_(char *, ftnlen);
    static char tokens[35*100];
    extern /* Subroutine */ int errint_(char *, integer *, ftnlen);
    extern logical return_(void);
    integer ptr;

/* $ Abstract */

/*     Read data from a MKDSK type 5 height grid file. */

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
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     INFILE     I   Name of input file. */
/*     NMAX       I   Maximum number of values to return. */
/*     N          O   Number of values returned. */
/*     VALUES     O   Data values. */
/*     DONE       O   Flag indicating EOF was reached. */

/* $ Detailed_Input */

/*     INFILE     is the name of an input data file containing height */
/*                grid data. */

/*     NMAX       is the maximum number of values to place in the */
/*                output array VALUES. */

/* $ Detailed_Output */

/*     N          is the number of values in the array VALUES. */

/*     VALUES     is a set of double precision height values */
/*                read from the file INFILE. */

/*     DONE       is a logical flag that is set to .TRUE. if and */
/*                only if the end of file was reached on the */
/*                current read. */

/*                Once DONE is set to .TRUE., a subsequent call */
/*                to this routine with the same value of INFILE */
/*                will cause the file to be read again from the */
/*                beginning. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If a token cannot be parsed as a double precision number, */
/*         the error SPICE(INVALIDDATA) is signaled. */

/*     2)  If the buffer size MAXN is not at least 1, the error */
/*         SPICE(INVALIDSIZE) is signaled. */

/*     3)  Any error that occurs while opening or reading the input */
/*         file will be diagnosed by a routine in the call tree of */
/*         this routine. */

/* $ Files */

/*     The file specified by INFILE must contain only tokens that can be */
/*     read as double precision values. No non-printing characters can */
/*     be present in the file. */

/*     Tokens can be delimited by blanks or commas. Tokens must not be */
/*     split across lines. */

/*     Blank lines are allowed; however, their use is discouraged */
/*     because they'll cause line numbers in diagnostic messages to */
/*     be out of sync with actual line numbers in the file. */

/*     The file must end with a line terminator. */

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


/*     Local parameters */


/*     Comma and space are the allowed delimiters. */


/*     Local variables */


/*     Saved values */


/*     Initial values */

    if (return_()) {
	return 0;
    }
    chkin_("RDGRD5", (ftnlen)6);

/*     Check NMAX. */

    if (*nmax < 1) {
	setmsg_("Buffer size NMAX = #; must be strictly positive.", (ftnlen)
		48);
	errint_("#", nmax, (ftnlen)1);
	sigerr_("SPICE(INVALIDSIZE)", (ftnlen)18);
	chkout_("RDGRD5", (ftnlen)6);
	return 0;
    }
    if (newfil) {

/*        Initialize our local state variables. */

	eof = FALSE_;
	lineno = 0;
	from = 0;
	ntk = 0;
	newfil = FALSE_;
    }

/*     Transfer as many as NMAX values to the output buffer. */
/*     Stop when the buffer fills up or when we run out of */
/*     data. */

    *n = 0;
    to = 0;
    *done = FALSE_;
    while(*n < *nmax && ! (*done)) {

/*        At this point, there's room in the output buffer, and */
/*        we haven't seen the end of the input file. */

/*        We may have buffered data from the last line read. */

	if (from == 0) {

/*           We don't have any buffered data. Read a new non-blank */
/*           line from the file. */

	    rdnbl_(infile__, line, &eof, infile_len, (ftnlen)255);
	    if (failed_()) {
		chkout_("RDGRD5", (ftnlen)6);
		return 0;
	    }
	    if (! eof) {

/*              We have a new line. Parse the values in the line. */

		++lineno;
		lparsm_(line, " ,", &c__100, &ntk, tokens, (ftnlen)255, (
			ftnlen)2, (ftnlen)35);
		if (failed_()) {
		    chkout_("RDGRD5", (ftnlen)6);
		    return 0;
		}

/*              Let NXFR be the number of values to transfer from */
/*              the current input line. */

		room = *nmax - *n;
		nxfr = min(room,ntk);
		i__1 = nxfr;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ++to;
		    nparsd_(tokens + ((i__2 = i__ - 1) < 100 && 0 <= i__2 ? 
			    i__2 : s_rnge("tokens", i__2, "rdgrd5_", (ftnlen)
			    298)) * 35, &values[to - 1], errmsg, &ptr, (
			    ftnlen)35, (ftnlen)320);
		    if (failed_()) {
			chkout_("RDGRD5", (ftnlen)6);
			return 0;
		    }
		    if (s_cmp(errmsg, " ", (ftnlen)320, (ftnlen)1) != 0) {
			setmsg_("Token number # on non-blank line # of data "
				"file <#>: #", (ftnlen)54);
			errint_("#", &i__, (ftnlen)1);
			errint_("#", &lineno, (ftnlen)1);
			errch_("#", infile__, (ftnlen)1, infile_len);
			errch_("#", errmsg, (ftnlen)1, (ftnlen)320);
			sigerr_("SPICE(INVALIDDATA)", (ftnlen)18);
			chkout_("RDGRD5", (ftnlen)6);
			return 0;
		    }
		}
		*n += nxfr;
		remain = ntk - nxfr;
		if (remain > 0) {

/*                 We didn't transfer all tokens. Let FROM */
/*                 be the index of the next token to transfer. */
/*                 We'll transfer the token on the next call. */

		    from = nxfr + 1;
		    chkout_("RDGRD5", (ftnlen)6);
		    return 0;
		}
	    } else {

/*              There are no more data to be had. RDNBL will close */
/*              INFILE. */

		*done = TRUE_;

/*              Get ready for another file. */

		newfil = TRUE_;
	    }
	} else {

/*           We have buffered tokens. Transfer as many of these */
/*           as we can. */

	    room = *nmax - *n;
	    nxfr = min(room,remain);
	    i__1 = from - 1 + nxfr;
	    for (i__ = from; i__ <= i__1; ++i__) {
		++to;
		nparsd_(tokens + ((i__2 = i__ - 1) < 100 && 0 <= i__2 ? i__2 :
			 s_rnge("tokens", i__2, "rdgrd5_", (ftnlen)365)) * 35,
			 &values[to - 1], errmsg, &ptr, (ftnlen)35, (ftnlen)
			320);
		if (failed_()) {
		    chkout_("RDGRD5", (ftnlen)6);
		    return 0;
		}
		if (s_cmp(errmsg, " ", (ftnlen)320, (ftnlen)1) != 0) {
		    setmsg_("Token number # on non-blank line # of data file"
			    " <#>: #", (ftnlen)54);
		    errint_("#", &i__, (ftnlen)1);
		    errint_("#", &lineno, (ftnlen)1);
		    errch_("#", infile__, (ftnlen)1, infile_len);
		    errch_("#", errmsg, (ftnlen)1, (ftnlen)320);
		    sigerr_("SPICE(INVALIDDATA)", (ftnlen)18);
		    chkout_("RDGRD5", (ftnlen)6);
		    return 0;
		}
	    }
	    *n += nxfr;
	    remain -= nxfr;
	    if (remain > 0) {
		from += nxfr;
	    } else {

/*              Indicate the buffer has been exhausted. */

		from = 0;
	    }
	}
    }
    chkout_("RDGRD5", (ftnlen)6);
    return 0;
} /* rdgrd5_ */

