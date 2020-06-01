/* addcom.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__10000 = 10000;
static integer c__2 = 2;

/* $Procedure ADDCOM ( Add comments to DSK file ) */
/* Subroutine */ int addcom_(integer *handle, char *cmdfil, char *inpfn, char 
	*outfn, char *cmtfil, logical *appflg, ftnlen cmdfil_len, ftnlen 
	inpfn_len, ftnlen outfn_len, ftnlen cmtfil_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rnge(char *, integer, 
	    char *, integer), f_clos(cllist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static char line[80];
    static doublereal tvec[6];
    extern /* Subroutine */ int dasac_(integer *, integer *, char *, ftnlen);
    static integer l, m;
    extern /* Subroutine */ int chkin_(char *, ftnlen), dpfmt_(doublereal *, 
	    char *, char *, ftnlen, ftnlen);
    extern integer rtrim_(char *, ftnlen);
    extern logical failed_(void);
    extern /* Subroutine */ int scardc_(integer *, char *, ftnlen), readln_(
	    integer *, char *, logical *, ftnlen);
    static char cmnbuf[255*10006];
    extern /* Subroutine */ int chkout_(char *, ftnlen), ssizec_(integer *, 
	    char *, ftnlen);
    static char astrln[80];
    static integer cmnunt;
    extern /* Subroutine */ int cputim_(doublereal *);
    static char tstamp[80];
    extern integer frstnp_(char *, ftnlen);
    extern logical return_(void);
    extern /* Subroutine */ int txtopr_(char *, integer *, ftnlen);
    static logical eof;

/* $ Abstract */

/*     Add comments to the output DSK file created by MKDSK. */

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

/*     FILES */
/*     TOPOGRAPHY */

/* $ Declarations */
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
/*     CMDFIL     I   Name of setup file. */
/*     HANDLE     I   Handle of DSK file. */
/*     CMTFIL     I   Name of comment file. */
/*     INPFN      I   Input shape file name. */
/*     OUTFN      I   Output DSK file name. */
/*     APPFLG     I   Append flag. */

/* $ Detailed_Input */

/*     CMDFIL         is the name of the MKDSK setup file. */

/*     HANDLE         is the handle associated with the output */
/*                    DSK file. */

/*     CMTFIL         is the name of a text file containing comments */
/*                    to be added to the DSK file.  If CMTFIL is */
/*                    blank, this is interpreted to mean there is no */
/*                    comment file. */

/*     INPFN          is the name of the input shape file to be */
/*                    converted to DSK format. */

/*     OUTFN          is the name of the output DSK file resulting */
/*                    from conversion of the shape file. */

/*     APPFLG         is a logical flag which is .TRUE. if the output */
/*                    data are to be appended to an existing DSK file */
/*                    and .FALSE. if the output DSK file is new. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If an error occurs while reading the comment file, the */
/*        error will be diagnosed by routines in the call tree of */
/*        this routine. */

/*     2) If an error occurs while attempting to add comments to the */
/*        output DSK file, the error will be diagnosed by routines in */
/*        the call tree of this routine. */

/* $ Files */

/*     Normally, the mapping implemented by this routine is defined */
/*     by a kernel variable introduced into the kernel pool by loading */
/*     a SPICE text kernel. */

/* $ Particulars */

/*     This routine adds comments to the comment area of the output */
/*     DSK file. These are: */

/*        - The contents of a comment file, if any, specified in */
/*          the setup file. */

/*        - The run time and date, the names of the setup, input, */
/*          and output files, and an indication of whether the */
/*          output file was new or appended to. */

/*        - The contents of the setup file. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     1) This routine is intended for use only within the MKDSK */
/*        program. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    MKDSK Version 2.0.0, 07-FEB-2017 (NJB) */

/*        Now writes the MKDSK version string to the comment area */
/*        of the output file. */

/* -    MKDSK Version 1.0.0, 15-APR-2010 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Save large buffers to avoid stack problems in */
/*     some C environments. */

    if (return_()) {
	return 0;
    }
    chkin_("ADDCOM", (ftnlen)6);

/*     Set the maximum size of the comment line buffer. */

    ssizec_(&c__10000, cmnbuf, (ftnlen)255);
    l = 0;

/*     Comment area separator line. */

    s_copy(astrln, "********************************************************"
	    "************************", (ftnlen)80, (ftnlen)80);

/*     If the comment file was provided we copy its content to the */
/*     comment area. */

    if (s_cmp(cmtfil, " ", cmtfil_len, (ftnlen)1) != 0) {

/*        We open the comment file, copy text from it to the comment */
/*        buffer line by line, clean non-printing characters from the */
/*        lines on the fly and dump the buffer to the comment area */
/*        when it's full. We repeat until all comments have been copied. */

	txtopr_(cmtfil, &cmnunt, cmtfil_len);
	if (failed_()) {
	    chkout_("ADDCOM", (ftnlen)6);
	    return 0;
	}

/*        Insert top comment separator line. */

	s_copy(cmnbuf + 1530, astrln, (ftnlen)255, (ftnlen)80);
	s_copy(cmnbuf + 1785, " ", (ftnlen)255, (ftnlen)1);
	l = 2;

/*        Get next comment line. */

	readln_(&cmnunt, line, &eof, (ftnlen)80);
	while(! (eof || failed_())) {

/*           Replace non-printing characters with spaces. */

	    while(frstnp_(line, (ftnlen)80) != 0) {
		m = frstnp_(line, (ftnlen)80);
		*(unsigned char *)&line[m - 1] = ' ';
	    }
	    if (l < 10000) {

/*              Store line in the buffer. */

		++l;
		s_copy(cmnbuf + ((i__1 = l + 5) < 10006 && 0 <= i__1 ? i__1 : 
			s_rnge("cmnbuf", i__1, "addcom_", (ftnlen)265)) * 255,
			 line, (ftnlen)255, rtrim_(line, (ftnlen)80));
	    } else {

/*              Buffer is full. Set the cardinality of the comment */
/*              buffer and write it to DSK comment area. Reset */
/*              counter. */

		scardc_(&l, cmnbuf, (ftnlen)255);
		dasac_(handle, &l, cmnbuf + 1530, (ftnlen)255);
		l = 0;
	    }

/*           Get next comment line. */

	    readln_(&cmnunt, line, &eof, (ftnlen)80);
	}

/*        Dump the rest of the buffer into the comment area. */

	if (l != 0) {
	    scardc_(&l, cmnbuf, (ftnlen)255);
	    dasac_(handle, &l, cmnbuf + 1530, (ftnlen)255);
	    l = 0;
	}

/*        Close comment file. */

	cl__1.cerr = 0;
	cl__1.cunit = cmnunt;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     Add a header preceding contents of the setup file and containing */
/*     setup file name and current CPU time. */

    cputim_(tvec);
    s_copy(tstamp, "YYYY-MM-DDTHR:MN:SC", (ftnlen)80, (ftnlen)19);
    dpfmt_(tvec, "0YYY", tstamp, (ftnlen)4, (ftnlen)4);
    dpfmt_(&tvec[1], "0M", tstamp + 5, (ftnlen)2, (ftnlen)2);
    dpfmt_(&tvec[2], "0D", tstamp + 8, (ftnlen)2, (ftnlen)2);
    dpfmt_(&tvec[3], "0h", tstamp + 11, (ftnlen)2, (ftnlen)2);
    dpfmt_(&tvec[4], "0m", tstamp + 14, (ftnlen)2, (ftnlen)2);
    dpfmt_(&tvec[5], "0s", tstamp + 17, (ftnlen)2, (ftnlen)2);
    s_copy(cmnbuf + 1530, " ", (ftnlen)255, (ftnlen)1);
    s_copy(cmnbuf + 1785, astrln, (ftnlen)255, (ftnlen)80);
    s_copy(cmnbuf + 2040, "MKDSK VERSION:       2.0.0, 28-FEB-2017", (ftnlen)
	    255, (ftnlen)39);
/* Writing concatenation */
    i__2[0] = 21, a__1[0] = "MKDSK RUN DATE/TIME: ";
    i__2[1] = rtrim_(tstamp, (ftnlen)80), a__1[1] = tstamp;
    s_cat(cmnbuf + 2295, a__1, i__2, &c__2, (ftnlen)255);
/* Writing concatenation */
    i__2[0] = 21, a__1[0] = "MKDSK SETUP FILE:    ";
    i__2[1] = rtrim_(cmdfil, cmdfil_len), a__1[1] = cmdfil;
    s_cat(cmnbuf + 2550, a__1, i__2, &c__2, (ftnlen)255);
/* Writing concatenation */
    i__2[0] = 21, a__1[0] = "MKDSK INPUT FILE:    ";
    i__2[1] = rtrim_(inpfn, inpfn_len), a__1[1] = inpfn;
    s_cat(cmnbuf + 2805, a__1, i__2, &c__2, (ftnlen)255);
/* Writing concatenation */
    i__2[0] = 21, a__1[0] = "MKDSK OUTPUT FILE:   ";
    i__2[1] = rtrim_(outfn, outfn_len), a__1[1] = outfn;
    s_cat(cmnbuf + 3060, a__1, i__2, &c__2, (ftnlen)255);
    if (*appflg) {
	s_copy(cmnbuf + 3315, "OUTPUT FILE STATUS:    EXISTING FILE", (ftnlen)
		255, (ftnlen)36);
    } else {
	s_copy(cmnbuf + 3315, "OUTPUT FILE STATUS:    NEW FILE", (ftnlen)255, 
		(ftnlen)31);
    }
    s_copy(cmnbuf + 3570, astrln, (ftnlen)255, (ftnlen)80);
    s_copy(cmnbuf + 3825, " ", (ftnlen)255, (ftnlen)1);
    l = 10;

/*     Now we will copy contents of the setup file to the comment area */
/*     using exactly the same procedure: open the setup file, copy */
/*     text from the file into the buffer line by line, clean */
/*     non-printing characters from the lines on the fly and dump the */
/*     buffer to the comment area when it's full. We repeat until all */
/*     setup lines have been copied. */

    txtopr_(cmdfil, &cmnunt, cmdfil_len);
    eof = FALSE_;
    while(! eof) {

/*        Read next line. */

	readln_(&cmnunt, line, &eof, (ftnlen)80);
	if (! eof) {

/*           Replace non-printing character with spaces. */

	    while(frstnp_(line, (ftnlen)80) != 0) {
		m = frstnp_(line, (ftnlen)80);
		*(unsigned char *)&line[m - 1] = ' ';
	    }
	    if (l < 10000) {

/*              Store line on buffer. */

		++l;
		s_copy(cmnbuf + ((i__1 = l + 5) < 10006 && 0 <= i__1 ? i__1 : 
			s_rnge("cmnbuf", i__1, "addcom_", (ftnlen)371)) * 255,
			 line, (ftnlen)255, rtrim_(line, (ftnlen)80));
	    } else {

/*              Buffer is full. Set the cardinality of the comment */
/*              buffer and write it to DSK comment area. Reset counter */
/*              and store the last line that we have obtained in the */
/*              first line of the buffer. */

		scardc_(&l, cmnbuf, (ftnlen)255);
		dasac_(handle, &l, cmnbuf + 1530, (ftnlen)255);
		l = 0;
	    }
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = cmnunt;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     Add "bottom of the comments" separator line. */

    if (l <= 9998) {
	s_copy(cmnbuf + ((i__1 = l + 6) < 10006 && 0 <= i__1 ? i__1 : s_rnge(
		"cmnbuf", i__1, "addcom_", (ftnlen)398)) * 255, " ", (ftnlen)
		255, (ftnlen)1);
	s_copy(cmnbuf + ((i__1 = l + 7) < 10006 && 0 <= i__1 ? i__1 : s_rnge(
		"cmnbuf", i__1, "addcom_", (ftnlen)399)) * 255, astrln, (
		ftnlen)255, (ftnlen)80);
	l += 2;
    } else {

/*        Dump current contents of the comment buffer, first. After */
/*        that stick separator at the top of the buffer. */

	scardc_(&l, cmnbuf, (ftnlen)255);
	dasac_(handle, &l, cmnbuf + 1530, (ftnlen)255);
	s_copy(cmnbuf + 1530, " ", (ftnlen)255, (ftnlen)1);
	s_copy(cmnbuf + 1785, astrln, (ftnlen)255, (ftnlen)80);
	l = 2;
    }

/*     Dump the buffer one more time, if required. */

    if (l > 0) {
	scardc_(&l, cmnbuf, (ftnlen)255);
	dasac_(handle, &l, cmnbuf + 1530, (ftnlen)255);
    }
    chkout_("ADDCOM", (ftnlen)6);
    return 0;
} /* addcom_ */

