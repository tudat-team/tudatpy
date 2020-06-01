/* mkdsk.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static logical c_false = FALSE_;

/* Main program */ MAIN__(void)
{
    logical info;
    extern /* Subroutine */ int chkin_(char *, ftnlen), reset_(void);
    extern logical eqstr_(char *, char *, ftnlen, ftnlen);
    char input[255], setup[255];
    extern logical failed_(void);
    integer handle;
    extern /* Subroutine */ int delfil_(char *, ftnlen);
    char cmdlin[2000];
    extern /* Subroutine */ int getcml_(char *, ftnlen);
    char cmtfil[255];
    extern /* Subroutine */ int byebye_(char *, ftnlen), erract_(char *, char 
	    *, ftnlen, ftnlen), prcinf_(char *, ftnlen), dskcls_(integer *, 
	    logical *), getset_(char *, char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), chkout_(char *, ftnlen), prscml_(char *, 
	    logical *, char *, char *, char *, char *, ftnlen, ftnlen, ftnlen,
	     ftnlen, ftnlen);
    char inftyp[32];
    extern /* Subroutine */ int tostdo_(char *, ftnlen), wrtdsk_(char *, char 
	    *, char *, char *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical exists_(char *, ftnlen);
    char output[255];

/* $ Abstract */

/*     Convert an input file containing topography data to a SPICE DSK */
/*     file. */

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

/*     FILES */
/*     TOPOGRAPHY */

/* $ Files */

/*     Inputs:  Shape data file, setup file, leapseconds kernel, optional */
/*              comment file. */

/*     Outputs:  DSK file. */

/* $ Particulars */

/*     None. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -    MKDSK Version 1.0.0, 31-MAR-2010 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Global parameters */


/*     Local parameters */

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


/*     Local variables */

    chkin_("MKDSK", (ftnlen)5);

/*     Turn off the "Oh, by the way..." error handling dreck. */

    erract_("SET", "ABORT", (ftnlen)3, (ftnlen)5);

/*     Show version info. */

    prcinf_("VERSION", (ftnlen)7);

/*     Get and parse the command line. */

    getcml_(cmdlin, (ftnlen)2000);
    prscml_(cmdlin, &info, inftyp, setup, input, output, (ftnlen)2000, (
	    ftnlen)32, (ftnlen)255, (ftnlen)255, (ftnlen)255);
    if (info) {

/*        The command is an information request.  Do nothing if */
/*        the version was requested, since we've already displayed it. */

	if (! eqstr_(inftyp, "VERSION", (ftnlen)32, (ftnlen)7)) {
	    prcinf_(inftyp, (ftnlen)32);
	}
    } else {

/*        Process the setup file.  If the input or output files were */
/*        not named on the command line, use the file names specified */
/*        in the setup file.  GETSET will update INPUT and OUTPUT */
/*        if necessary. */

	getset_(setup, input, output, cmtfil, (ftnlen)255, (ftnlen)255, (
		ftnlen)255, (ftnlen)255);

/*        At this point, we're ready to do the conversion.  Set the error */
/*        handling so we'll return to this point if a conversion error */
/*        occurs. */

	erract_("SET", "RETURN", (ftnlen)3, (ftnlen)6);
	wrtdsk_(setup, input, output, cmtfil, &handle, (ftnlen)255, (ftnlen)
		255, (ftnlen)255, (ftnlen)255);
	if (failed_()) {
	    reset_();
	    tostdo_(" ", (ftnlen)1);
	    tostdo_("Conversion failed.", (ftnlen)18);
	    if (exists_(output, (ftnlen)255)) {

/*              If the output file is open, close it before calling */
/*              DELFIL. */

		dskcls_(&handle, &c_false);
		delfil_(output, (ftnlen)255);
		if (! failed_()) {
		    tostdo_("Output file has been deleted.", (ftnlen)29);
		    tostdo_(" ", (ftnlen)1);
		}
	    }
	    tostdo_(" ", (ftnlen)1);
	    byebye_("FAILURE", (ftnlen)7);
	}
    }

/*     Perform termination functions. */

    if (! info) {
	tostdo_(" ", (ftnlen)1);
	tostdo_("All done.", (ftnlen)9);
	tostdo_(" ", (ftnlen)1);
    }
    chkout_("MKDSK", (ftnlen)5);
    byebye_("SUCCESS", (ftnlen)7);
    return 0;
} /* MAIN__ */

/* Main program alias */ int mkdsk_ () { MAIN__ (); return 0; }
