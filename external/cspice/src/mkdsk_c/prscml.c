/* prscml.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__12 = 12;
static integer c__1 = 1;

/* $Procedure PRSCML ( Parse MKDSK command line ) */
/* Subroutine */ int prscml_(char *cmdlin, logical *info, char *inftyp, char *
	setup, char *infil, char *outfil, ftnlen cmdlin_len, ftnlen 
	inftyp_len, ftnlen setup_len, ftnlen infil_len, ftnlen outfil_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int chkin_(char *, ftnlen), errch_(char *, char *,
	     ftnlen, ftnlen);
    static logical found;
    static integer nkeys;
    extern integer rtrim_(char *, ftnlen);
    extern /* Subroutine */ int m2chck_(char *, char *, integer *, char *, 
	    char *, ftnlen, ftnlen, ftnlen, ftnlen), m2getc_(char *, char *, 
	    logical *, char *, ftnlen, ftnlen, ftnlen), m2ints_(integer *, 
	    char *, integer *, char *, ftnlen, ftnlen);
    extern logical m2xist_(char *, ftnlen);
    static char loccmd[2000];
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen);
    static char errmsg[1840*2];
    extern /* Subroutine */ int ssizec_(integer *, char *, ftnlen), ssizei_(
	    integer *, integer *), prefix_(char *, integer *, char *, ftnlen, 
	    ftnlen), setmsg_(char *, ftnlen);
    static char synval[2000*18];
    extern logical return_(void);
    extern /* Subroutine */ int prompt_(char *, char *, ftnlen, ftnlen);
    static char synkey[32*18];
    static integer synptr[18];

/* $ Abstract */

/*     Parse the command line arguments supplied to MKDSK. */

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
/*     PCK */
/*     TIME */

/* $ Keywords */

/*     TOPOGRAPHY */
/*     FILES */

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


/*     Include File:  SPICELIB Error Handling Parameters */

/*        errhnd.inc  Version 2    18-JUN-1997 (WLT) */

/*           The size of the long error message was */
/*           reduced from 25*80 to 23*80 so that it */
/*           will be accepted by the Microsoft Power Station */
/*           FORTRAN compiler which has an upper bound */
/*           of 1900 for the length of a character string. */

/*        errhnd.inc  Version 1    29-JUL-1997 (NJB) */



/*     Maximum length of the long error message: */


/*     Maximum length of the short error message: */


/*     End Include File:  SPICELIB Error Handling Parameters */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     CMDLIN     I   Command line to be parsed. */
/*     INFO       O   Flag indicating whether command is an info request. */
/*     INFTYP     O   Information type if information is requested. */
/*     SETUP      O   Setup file name. */
/*     INFIL      O   Input file name. */
/*     OUTFIL     O   Output file name. */

/* $ Detailed_Input */

/*     CMDLIN         is a string containing the command line arguments */
/*                    supplied when MKDSK was invoked. */

/* $ Detailed_Output */

/*     INFO           is a logical flag indicating whether the command */
/*                    is a request for information:  usage, "help,", or */
/*                    the program's version. */

/*                    When INFO is .TRUE., INFTYP will indicate the type */
/*                    of information requested. */

/*     INFTYP         is a string indicating the type of information */
/*                    requested, if the input command requests such. */
/*                    Values of INFTYP are: */

/*                       'HELP' */
/*                       'VERSION' */
/*                       'USAGE' */
/*                       'TEMPLATE' */

/*                    INFTYP is valid only if INFO is .TRUE.  Otherwise, */
/*                    IFNTYP is returned blank. */

/*     SETUP          is the name of the setup file specified on the */
/*                    command line.  If no file is specified, and */
/*                    if the command is not an information request, the */
/*                    user is prompted for a file name. */

/*                    If INFO is returned .TRUE., SETUP is left blank. */

/*     INFIL          is the name of the input file specified on the */
/*                    command line.  If no input file name is specified, */
/*                    INFIL is left blank. */

/*                    If INFO is returned .TRUE., INFIL is left blank. */

/*     OUTFIL         is the name of the output file specified on the */
/*                    command line.  If no input file name is specified, */
/*                    OUTFIL is left blank. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the input command is syntactically invalid, the error */
/*        SPICE(CMDERROR) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     The expected command syntax is: */

/*        mkdsk   [-setup <setup file name>] */
/*                [-input <input data file name>] */
/*                [-output <output SPK file name>] */
/*                [-h|-help] */
/*                [-t|-template] */
/*                [-u|-usage] */
/*                [-v|-version] */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     1) This routine is intended for use solely within the MKDSK */
/*        program. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -    MKDSK Version 1.0.0, 30-MAR-2010 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Non-SPICELIB  functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Initial values */

    if (return_()) {
	return 0;
    }
    chkin_("PRSCML", (ftnlen)6);

/*     Give initial values to the output arguments. */

    *info = FALSE_;
    s_copy(inftyp, " ", inftyp_len, (ftnlen)1);
    s_copy(setup, " ", setup_len, (ftnlen)1);
    s_copy(infil, " ", infil_len, (ftnlen)1);
    s_copy(outfil, " ", outfil_len, (ftnlen)1);
    if (first) {

/*        Initialize the symbol table representing the command */
/*        language. */

	ssizec_(&c__12, synkey, (ftnlen)32);
	ssizec_(&c__12, synval, (ftnlen)2000);
	ssizei_(&c__12, synptr);
	s_copy(synval + 12000, "HELPKEY (1:1){ -h[help]              |      "
		"         -help[help]           |               -v[version]  "
		"         |               -version[version]     |            "
		"   -t[template]          |               -template[template]"
		"   |               -u[usage]             |               -us"
		"age[usage]            }", (ftnlen)2000, (ftnlen)307);
	s_copy(synval + 14000, "CONVKEY         (1:4){ -setup @word[setup]  "
		" |               -input @word[input]   |               -outp"
		"ut @word[output] |               -append[append]          }", 
		(ftnlen)2000, (ftnlen)163);
	nkeys = 2;
	m2ints_(&nkeys, synkey, synptr, synval, (ftnlen)32, (ftnlen)2000);
	first = FALSE_;
    }

/*     See whether we have a blank command. */

    if (s_cmp(cmdlin, " ", cmdlin_len, (ftnlen)1) == 0) {

/*        Prompt for the setup file name, then return. */

	prompt_("SETUP FILE NAME> ", setup, (ftnlen)17, setup_len);
	chkout_("PRSCML", (ftnlen)6);
	return 0;
    }

/*     See whether we have some type of help command. */

    s_copy(loccmd, cmdlin, (ftnlen)2000, cmdlin_len);
    prefix_("HELPKEY", &c__1, loccmd, (ftnlen)7, (ftnlen)2000);
    m2chck_(loccmd, synkey, synptr, synval, errmsg, (ftnlen)2000, (ftnlen)32, 
	    (ftnlen)2000, (ftnlen)1840);
    if (s_cmp(errmsg, " ", (ftnlen)1840, (ftnlen)1) == 0) {

/*        The command matches the HELP syntax. */

	*info = TRUE_;
	s_copy(inftyp, " ", inftyp_len, (ftnlen)1);

/*        We rely on M2CHCK to make sure one of the expected */
/*        verbs is present, so we don't check the FOUND flag. */

	if (m2xist_("help", (ftnlen)4)) {
	    s_copy(inftyp, "HELP", inftyp_len, (ftnlen)4);
	} else if (m2xist_("version", (ftnlen)7)) {
	    s_copy(inftyp, "VERSION", inftyp_len, (ftnlen)7);
	} else if (m2xist_("usage", (ftnlen)5)) {
	    s_copy(inftyp, "USAGE", inftyp_len, (ftnlen)5);
	} else if (m2xist_("template", (ftnlen)8)) {
	    s_copy(inftyp, "TEMPLATE", inftyp_len, (ftnlen)8);
	}

/*        We're done with this command. */

	chkout_("PRSCML", (ftnlen)6);
	return 0;
    }

/*     See whether a conversion has been requested. */

    s_copy(loccmd, cmdlin, (ftnlen)2000, cmdlin_len);
    prefix_("CONVKEY", &c__1, loccmd, (ftnlen)7, (ftnlen)2000);
    s_copy(errmsg, " ", (ftnlen)1840, (ftnlen)1);
    s_copy(errmsg + 1840, " ", (ftnlen)1840, (ftnlen)1);
    m2chck_(loccmd, synkey, synptr, synval, errmsg, (ftnlen)2000, (ftnlen)32, 
	    (ftnlen)2000, (ftnlen)1840);
    if (s_cmp(errmsg, " ", (ftnlen)1840, (ftnlen)1) == 0) {

/*        We have a syntactically correct conversion command. */

	m2getc_("setup", loccmd, &found, setup, (ftnlen)5, (ftnlen)2000, 
		setup_len);
	if (! found) {

/*           Prompt for the setup file name. */

	    prompt_("SETUP FILE NAME> ", setup, (ftnlen)17, setup_len);
	}

/*        Get the input and output file names if they're present. */

	m2getc_("input", loccmd, &found, infil, (ftnlen)5, (ftnlen)2000, 
		infil_len);
	m2getc_("output", loccmd, &found, outfil, (ftnlen)6, (ftnlen)2000, 
		outfil_len);
    } else {
	setmsg_("The command <#> doesn't match any known command syntax.", (
		ftnlen)55);
	errch_("#", cmdlin, (ftnlen)1, rtrim_(cmdlin, cmdlin_len));
	sigerr_("SPICE(CMDERROR)", (ftnlen)15);
	chkout_("PRSCML", (ftnlen)6);
	return 0;
    }
    chkout_("PRSCML", (ftnlen)6);
    return 0;
} /* prscml_ */

