/* extmsi.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure    EXTMSI  ( print message with integer and do error exit ) */
/* Subroutine */ int extmsi_(char *messge, char *marker, integer *icode, 
	ftnlen messge_len, ftnlen marker_len)
{
    extern /* Subroutine */ int chkin_(char *, ftnlen), sigerr_(char *, 
	    ftnlen), chkout_(char *, ftnlen), setmsg_(char *, ftnlen), 
	    errint_(char *, integer *, ftnlen);
    extern logical return_(void);

/* $ Abstract */

/*     This subroutine prints an error message with an integer */
/*     and then exits. */

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

/*     ERROR */

/* $ Declarations */
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     MESSGE     I   Character string with message and an imbedded */
/*                    marker. */
/*     MARKER     I   Character string marker to be replaced. */
/*     ICODE      I   Integer that will be inserted into the message. */

/* $ Detailed_Input */

/*     MESSGE     is a character string with a message, and a marker */
/*                where the integer ICODE will be inserted. */

/*     MARKER     is a character string to be replaced by an integer. */

/*     ICODE      is an integer that will be inserted into the message */
/*                at the location of the marker. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     Error free. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This is a common task in MKPLAT. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     J.A. Bytof      (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 08-OCT-2009 (NJB) */

/*        Re-ordered header sections. */

/* -    SPICELIB Version 1.0.0, 10-APR-1997 (JAB) */

/* -& */
/* $ Index_Entries */

/*     Error exit message with integer. */

/* -& */

/*     SPICELIB functions */


/*     Local variables. */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("EXTMSI", (ftnlen)6);
    }
    setmsg_(messge, messge_len);
    errint_(marker, icode, marker_len);
    sigerr_("SPICE(ERROREXIT)", (ftnlen)16);
    chkout_("EXTMSI", (ftnlen)6);
    return 0;
} /* extmsi_ */

