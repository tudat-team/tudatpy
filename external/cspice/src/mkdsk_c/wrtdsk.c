/* wrtdsk.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__10000 = 10000;
static logical c_false = FALSE_;
static logical c_true = TRUE_;

/* $Procedure    WRTDSK ( MKDSK, write DSK file ) */
/* Subroutine */ int wrtdsk_(char *setup, char *input, char *output, char *
	cmtfil, integer *handle, ftnlen setup_len, ftnlen input_len, ftnlen 
	output_len, ftnlen cmtfil_len)
{
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer dtype;
    extern /* Subroutine */ int ljust_(char *, char *, ftnlen, ftnlen), 
	    rjust_(char *, char *, ftnlen, ftnlen);
    extern logical failed_(void);
    extern /* Subroutine */ int addcom_(integer *, char *, char *, char *, 
	    char *, logical *, ftnlen, ftnlen, ftnlen, ftnlen);
    char ifname[60];
    extern /* Subroutine */ int dskcls_(integer *, logical *), sigerr_(char *,
	     ftnlen), chkout_(char *, ftnlen), dskopn_(char *, char *, 
	    integer *, integer *, ftnlen, ftnlen), setmsg_(char *, ftnlen), 
	    errint_(char *, integer *, ftnlen), gettyp_(integer *), tostdo_(
	    char *, ftnlen);
    extern logical return_(void);
    extern /* Subroutine */ int zzwseg02_(char *, integer *, ftnlen);

/* $ Abstract */

/*     Write a new DSK file. */

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

/*     DSK */
/*     FILES */

/* $ Declarations */
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     SETUP      I   Name of MKDSK setup file. */
/*     INPUT      I   Name of MKDSK input shape file. */
/*     OUTPUT     I   Name of DSK file to create. */
/*     CMTFIL     I   Name of comment file. */
/*     HANDLE     O   DAS file handle of DSK. */

/* $ Detailed_Input */

/*     SETUP          is the name of the MKDSK setup file */
/*                    that describes the file to create. */
/*                    See the MKDSK User's Guide for details. */

/*     INPUT          is the name of the shape file containing */
/*                    data to be converted to DSK format. */
/*                    See the MKDSK User's Guide for details. */

/*     OUTPUT         is the name of the DSK file to create. */

/*     CMTFIL         is the name of a text file containing */
/*                    comments to be inserted into the */
/*                    comment area of the DSK file. */

/*                    The caller can set CMTFIL to blank to */
/*                    indicate that no comment file is provided. */

/* $ Detailed_Output */

/*     HANDLE         is the DAS handle of the output DSK file. */
/*                    HANDLE is returned so the caller can */
/*                    close and delete the DSK file if an error */
/*                    occurs. Normally, the DSK file is closed */
/*                    by this routine and HANDLE is not used */
/*                    by the caller. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     This routine is meant to be operated in RETURN SPICE error */
/*     handling mode. The caller is expected to delete the DSK file if */
/*     an error occurs during file creation. */


/*     1) If the setup file requests creation of a segment having */
/*        an unrecognized data type, the error SPICE(NOTSUPPORTED) */
/*        is signaled. */

/*     2) If a new DSK having the specified name cannot be created, */
/*        the error will be diagnosed by routines in the call tree */
/*        of this routine. */

/*     3) If an error occurs while writing comments to the DSK file, */
/*        the error will be diagnosed by routines in the call tree */
/*        of this routine. */

/*     4) If an error occurs while writing data to the DSK file, */
/*        the error will be diagnosed by routines in the call tree */
/*        of this routine. */

/*     5) If an error is present in the setup file, the error will */
/*        be diagnosed, if possible, by routines in the call tree */
/*        of this routine. */

/* $ Files */

/*     See the Detailed_Input section above. */

/* $ Particulars */

/*     This routine executes the high-level file creation */
/*     operations carried out by MKDSK: */

/*        1) Open a new DSK file */
/*        2) Write comments to the file */
/*        3) Write a single segment to the file */
/*        4) Close the file */

/* $ Examples */

/*     See usage in MKDSK. */

/* $ Restrictions */

/*     This routine should be called only from within MKDSK. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.1.0, 03-MAY-2014 (NJB) */

/*        Now calls ZZWSEG02. */

/* -    SPICELIB Version 1.0.0, 15-APR-2010 (NJB) */

/* -& */
/* $ Index_Entries */

/*     write dsk file */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     DAS internal file name length. */


/*     Local variables */

    if (return_()) {
	return 0;
    }
    chkin_("WRTDSK", (ftnlen)6);

/*     Use the trailing IFNLEN characters of the DSK file */
/*     name as the internal file name. */

    rjust_(output, ifname, output_len, (ftnlen)60);
    ljust_(ifname, ifname, (ftnlen)60, (ftnlen)60);

/*     Open the DSK file for write access. */

    dskopn_(output, ifname, &c__10000, handle, output_len, (ftnlen)60);

/*     Write comments to the file. */

    addcom_(handle, setup, input, output, cmtfil, &c_false, setup_len, 
	    input_len, output_len, cmtfil_len);

/*     Fetch segment data type. */

    gettyp_(&dtype);

/*     Call the segment writer of the appropriate type. */

    if (dtype == 2) {
	zzwseg02_(input, handle, input_len);
    } else {
	setmsg_("Segment type # is not supported.", (ftnlen)32);
	errint_("#", &dtype, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("WRTDSK", (ftnlen)6);
	return 0;
    }
    if (failed_()) {
	chkout_("WRTDSK", (ftnlen)6);
	return 0;
    }

/*     Close DSK file. */

    tostdo_("Segregating and closing DSK file...", (ftnlen)35);
    dskcls_(handle, &c_true);
    if (failed_()) {
	chkout_("WRTDSK", (ftnlen)6);
	return 0;
    }
    tostdo_("DSK file was created.", (ftnlen)21);
    chkout_("WRTDSK", (ftnlen)6);
    return 0;
} /* wrtdsk_ */

