/* rdffdi.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__6 = 6;

/* $Procedure RDFFDI ( read and parse flatfile records ) */
/* Subroutine */ int rdffdi_(char *infile__, integer *nrec, char *format, 
	integer *nd, integer *ni, integer *nc, doublereal *ard, integer *ari, 
	char *arc, logical *eof, ftnlen infile_len, ftnlen format_len, ftnlen 
	arc_len)
{
    /* Initialized data */

    static char fmtsav[255] = "X                                            "
	    "                                                                "
	    "                                                                "
	    "                                                                "
	    "                  ";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rnge(char *, integer, 
	    char *, integer);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    doublereal datd;
    integer dati, ndat;
    static char line[255];
    static integer nfmt;
    integer i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen), errch_(char *, char *,
	     ftnlen, ftnlen);
    extern integer wdcnt_(char *, ftnlen);
    static char error[160];
    extern /* Subroutine */ int replch_(char *, char *, char *, char *, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
    static char datitm[255*6];
    extern integer lastnb_(char *, ftnlen);
    extern /* Subroutine */ int lparse_(char *, char *, integer *, integer *, 
	    char *, ftnlen, ftnlen, ftnlen), nparsd_(char *, doublereal *, 
	    char *, integer *, ftnlen, ftnlen), sigerr_(char *, ftnlen), 
	    nparsi_(char *, integer *, char *, integer *, ftnlen, ftnlen), 
	    chkout_(char *, ftnlen);
    static char fmtitm[255*6];
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), rdtext_(char *, char *, logical *, ftnlen, 
	    ftnlen);
    extern logical return_(void);
    integer ptr;

/* $ Abstract */

/*     Read a flatfile record and parse it into integer and/or double */
/*     precision and/or string values according to a given format */
/*     specification. */

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

/*     FLATFILE */

/* $ Declarations */
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     INFILE     I   Name of input flatfile. */
/*     NREC       I   Current record number in flatfile. */
/*     FORMAT     I   String containing format picture. */
/*     ND         O   Number of doubles found in data record. */
/*     NI         O   Number of integers found in data record. */
/*     NC         O   Number of strings found in data record. */
/*     ARD        O   Array of doubles read from data record. */
/*     ARI        O   Array of integer read from data record. */
/*     ARI        O   Array of strings read from data record. */
/*     EOF        O   End of file flag. */

/* $ Detailed_Input */

/*     INFILE     is the full pathname of the input flatfile. */

/*     NREC       is the current record number in the flatfile, */
/*                provided by the user. */

/*     FORMAT     is a string containing a format description or */
/*                picture of the data to be parsed.  In the present */
/*                scheme, 'I' signifies integer and 'D' double */
/*                precision values.  A typical format string might */
/*                be 'I D D D'.  Up to MAXITM format items are */
/*                permitted. */

/* $ Detailed_Output */

/*     ND         is the number of double precision values found. */

/*     NI         is the number of integer values found. */

/*     NC         the number of strings extracted from the data record. */

/*     ARD        is an array containing ND consecutively read */
/*                double precision values. */

/*     ARI        is an array containing NI consecutively read */
/*                integer values. */

/*     ARC        is an array containing NC consecutively read */
/*                strings. */

/*     EOF        is .TRUE. if end of file was encountered. */

/* $ Parameters */


/*     MAXWD      Maximum allowed length of a single field. */

/*     MAXSTR     Maximum length of an input data record. */

/*     MAXITM     Maximum number of items in format and input */
/*                data records. */

/* $ Exceptions */

/*     1)  If the format string length exceeds MAXSTR, then */
/*         SPICE(FORMATSTRINGTOOLONG) is signaled. */

/*     2)  If the number of format items exceeds MAXITM, then */
/*         SPICE(FORMATITEMLIMITEXCEEDED) is signaled. */

/*     3)  If the number of data items exceeds MAXITM, then */
/*         SPICE(DATAITEMLIMITEXCEEDED) is signaled. */

/*     4)  If an unknown format type occurs in the format string, */
/*         then SPICE(BADFORMATSPECIFIER) is signaled. */

/*     5)  If a string representing numerical data is exactly */
/*         MAXWD in length then SPICE(DATAWIDTHERROR) is signaled. */

/*     6)  If the format string and data string contain a */
/*         different number of items delimited by spaces, */
/*         then SPICE(FORMATDATAMISMATCH) is signaled. */

/*     7)  If an error is returned from NPARSI, then SPICE(BADINTEGER) */
/*         is signaled. */

/*     8)  If an error is returned from NPARSD, then */
/*         SPICE(BADDOUBLEPRECISION) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     1)  It is the responsibility of the programmer to supply the */
/*         correct record number in the file, NREC. */

/*     2)  Format and data items are delimited by spaces. */

/*     3)  ARC must be declared of suffcient length in the calling */
/*         program to hold the string data items, otherwise the */
/*         string is truncated when copied to ARC. */

/* $ Examples */

/*     C */
/*     C     Read a "self-describing" flatfile: */
/*     C */
/*     C     The first record contains a single integer, */
/*     C     which is the number of data records to be read. */
/*     C */
/*     C     The following line contains a format string for the */
/*     C     data records to follow. */
/*     C */
/*     C     After the last data record is read, a new data group */
/*     C     may begin, signified by a line containing an */
/*     C     integer count, followed by a format string and */
/*     C     data records.  Alternatively, the end of the file may */
/*     C     also be encountered. */
/*     C */


/*           LOGICAL               EOF */
/*           CHARACTER*(*)         FNAME */
/*           CHARACTER*(MAXSTR)    FORMAT */
/*           INTEGER               ND */
/*           INTEGER               NI */
/*           INTEGER               NC */
/*           INTEGER               NREC */
/*           INTEGER               NDAT */
/*           INTEGER               ARI     ( MAXITM ) */
/*           DOUBLE PRECISION      ARD     ( MAXITM ) */
/*           CHARACTER*(MAXSTR)    ARC     ( MAXITM ) */


/*           CALL PROMPT ( 'Enter file name ', FNAME ) */

/*     C */
/*     C     Get record count. */
/*     C */

/*           NREC = 1 */
/*           CALL RDFFDI ( FNAME, NREC, 'I', ND, NI, NC, */
/*          .              ARD, ARI, ARC, EOF ) */

/*           DO WHILE ( .NOT. EOF ) */

/*              NDAT = ARI(1) */

/*     C */
/*     C        Get data record format. */
/*     C */

/*              NREC = NREC + 1 */
/*              CALL RDTEXT ( FNAME, FORMAT, EOF ) */

/*              IF ( EOF ) THEN */
/*                 CALL EXTMSI ( 'End of file after line #.', '#', */
/*           .                    NREC-1 ) */
/*              END IF */

/*     C */
/*     C        Read data records. */
/*     C */

/*              DO I = 1, NDAT */
/*                 NREC = NREC + 1 */
/*                 CALL RDFFDI ( FNAME, NREC, FORMAT, ND, NI, NC */
/*           .                   ARD,   ARI,  ARC, EOF ) */
/*                 IF ( EOF ) THEN */
/*                    CALL EXTMSI ( 'End of file after line #.', '#', */
/*           .                       NREC-1 ) */
/*                 END IF */
/*                 CALL HAVFUN ( ND, NI, ARD, ARI ) */
/*              END DO */

/*     C */
/*     C        Get next data record count, if it's there. */
/*     C */

/*              NREC = NREC + 1 */
/*              CALL RDFFDI ( FNAME, NREC, 'I', ND, NI, NC, */
/*           .                ARD, ARI, ARC, EOF ) */

/*           END DO */


/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     J.A. Bytof      (JPL) */

/* $ Version */

/* -    SPICELIB Version 3.0.0, 03-MAY-2014 (NJB) */

/*        Increased supported input file line length to 255. */
/*        Added SAVE statements for several arrays and strings. */

/* -    SPICELIB Version 2.0.1, 08-OCT-2009 (NJB) */

/*        Re-ordered header sections. */

/* -    SPICELIB Version 2.0.0, 25-OCT-2004 (EDW) */

/*        Added capability to process character */
/*        typed data, 'C'. */

/* -    SPICELIB Version 1.0.0, 09-APRIL-1997 (JAB) */

/* -& */
/* $ Index_Entries */

/*     read and parse flatfile records. */

/* -& */

/*     SPICE functions. */


/*     Local variables. */


/*     Maximum allowed length of a single field. */


/*     Maximum length of an input data record. */


/*     Maximum number of items in format and input data records. */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("RDFFDI", (ftnlen)6);
    }

/*     Initialize array counters and logical flags. */

    *nd = 0;
    *ni = 0;
    *nc = 0;
    *eof = FALSE_;

/*     Check format string length. */

    if (lastnb_(format, format_len) > 255) {
	setmsg_("Format string length exceeds limit in line #.", (ftnlen)45);
	errint_("#", nrec, (ftnlen)1);
	sigerr_("SPICE(FORMATSTRINGTOOLONG)", (ftnlen)26);
    }

/*     Parse format string if it has changed.  Verify that format. */

    if (s_cmp(format, fmtsav, format_len, (ftnlen)255) != 0) {

/*     Does the number of format items exceed MAXITM? */

	if (wdcnt_(format, format_len) > 6) {
	    setmsg_("Too many format items at line #", (ftnlen)31);
	    errint_("#", nrec, (ftnlen)1);
	    sigerr_("SPICE(FMTITEMLIMITEXCEEDED)", (ftnlen)27);
	}
	lparse_(format, " ", &c__6, &nfmt, fmtitm, format_len, (ftnlen)1, (
		ftnlen)255);
	i__1 = nfmt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (lastnb_(fmtitm + ((i__2 = i__ - 1) < 6 && 0 <= i__2 ? i__2 : 
		    s_rnge("fmtitm", i__2, "rdffdi_", (ftnlen)381)) * 255, (
		    ftnlen)255) > 1) {
		setmsg_("Format error detected at line #.", (ftnlen)32);
		errint_("#", nrec, (ftnlen)1);
		sigerr_("SPICE(BADFORMATSPECIFIER)", (ftnlen)25);
	    }
	}
	s_copy(fmtsav, format, (ftnlen)255, format_len);
    }

/*     Read data record, if end of file is encountered, */
/*     check out, then return to calling routine. */

    rdtext_(infile__, line, eof, infile_len, (ftnlen)255);
    if (*eof) {
	chkout_("RDFFDI", (ftnlen)6);
	return 0;
    }

/*     Convert any tabs or carriage returns to blanks. */

    replch_(line, "\t", " ", line, (ftnlen)255, (ftnlen)1, (ftnlen)1, (ftnlen)
	    255);
    replch_(line, "\r", " ", line, (ftnlen)255, (ftnlen)1, (ftnlen)1, (ftnlen)
	    255);

/*     Does the number of data items exceed MAXITM? */

    ndat = wdcnt_(line, (ftnlen)255);
    if (ndat > 6) {
	setmsg_("Too many data items, line #", (ftnlen)27);
	errint_("#", nrec, (ftnlen)1);
	sigerr_("SPICE(DATAITEMLIMITEXCEEDED)", (ftnlen)28);
    }

/*     The format specifier and data record might mistakenly */
/*     have a different number of fields. */

    if (nfmt < ndat) {
	setmsg_("Too many data items in line #. The PLATE_TYPE setting may n"
		"ot  match the data file format.", (ftnlen)90);
	errint_("#", nrec, (ftnlen)1);
	sigerr_("SPICE(FORMATDATAMISMATCH)", (ftnlen)25);
    } else if (nfmt > ndat) {
	setmsg_("Too many format items in line # The PLATE_TYPE setting may "
		"not  match the data file format.", (ftnlen)91);
	errint_("#", nrec, (ftnlen)1);
	sigerr_("SPICE(FORMATDATAMISMATCH)", (ftnlen)25);
    }

/*     Parse data, verify that values fit within MAXWD field width. */

    lparse_(line, " ", &c__6, &ndat, datitm, (ftnlen)255, (ftnlen)1, (ftnlen)
	    255);
    i__1 = ndat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (lastnb_(datitm + ((i__2 = i__ - 1) < 6 && 0 <= i__2 ? i__2 : 
		s_rnge("datitm", i__2, "rdffdi_", (ftnlen)451)) * 255, (
		ftnlen)255) >= 40) {
	    setmsg_("Possible data error in line #.", (ftnlen)30);
	    errint_("#", nrec, (ftnlen)1);
	    sigerr_("SPICE(DATAWIDTHERROR)", (ftnlen)21);
	}
    }

/*     Process each data field according to its format. */

    i__1 = ndat;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     Convert the data field to requested type. */

	if (*(unsigned char *)&fmtitm[((i__2 = i__ - 1) < 6 && 0 <= i__2 ? 
		i__2 : s_rnge("fmtitm", i__2, "rdffdi_", (ftnlen)467)) * 255] 
		== 'I') {

/*           We expect an integer. Use NPARSI. */

	    nparsi_(datitm + ((i__2 = i__ - 1) < 6 && 0 <= i__2 ? i__2 : 
		    s_rnge("datitm", i__2, "rdffdi_", (ftnlen)472)) * 255, &
		    dati, error, &ptr, (ftnlen)255, (ftnlen)160);
	    if (ptr != 0) {
		setmsg_("Integer error (#) in line #.", (ftnlen)28);
		errch_("#", error, (ftnlen)1, (ftnlen)160);
		errint_("#", nrec, (ftnlen)1);
		sigerr_("SPICE(BADINTEGER)", (ftnlen)17);
	    }
	    ++(*ni);
	    ari[*ni - 1] = dati;
	} else if (*(unsigned char *)&fmtitm[((i__2 = i__ - 1) < 6 && 0 <= 
		i__2 ? i__2 : s_rnge("fmtitm", i__2, "rdffdi_", (ftnlen)483)) 
		* 255] == 'D') {

/*           We expect a double. Use NPARSD. */

	    nparsd_(datitm + ((i__2 = i__ - 1) < 6 && 0 <= i__2 ? i__2 : 
		    s_rnge("datitm", i__2, "rdffdi_", (ftnlen)488)) * 255, &
		    datd, error, &ptr, (ftnlen)255, (ftnlen)160);
	    if (ptr != 0) {
		setmsg_("D.P. error (#) in line #.", (ftnlen)25);
		errch_("#", error, (ftnlen)1, (ftnlen)160);
		errint_("#", nrec, (ftnlen)1);
		sigerr_("SPICE(BADDOUBLEPRECISION)", (ftnlen)25);
	    }
	    ++(*nd);
	    ard[*nd - 1] = datd;
	} else if (*(unsigned char *)&fmtitm[((i__2 = i__ - 1) < 6 && 0 <= 
		i__2 ? i__2 : s_rnge("fmtitm", i__2, "rdffdi_", (ftnlen)499)) 
		* 255] == 'C') {

/*           We expect a string. No need to parse, just copy */
/*           then count. Notice that if ARC was not declared */
/*           with a string size suffcient to store the */
/*           DATITM value, the copy op will truncate the value. */

	    ++(*nc);
	    s_copy(arc + (*nc - 1) * arc_len, datitm + ((i__2 = i__ - 1) < 6 
		    && 0 <= i__2 ? i__2 : s_rnge("datitm", i__2, "rdffdi_", (
		    ftnlen)508)) * 255, arc_len, (ftnlen)255);
	} else {

/*           Problem with format specifier. */

	    setmsg_("Bad format specifier # used in line #.", (ftnlen)38);
	    errch_("#", fmtitm + ((i__2 = i__ - 1) < 6 && 0 <= i__2 ? i__2 : 
		    s_rnge("fmtitm", i__2, "rdffdi_", (ftnlen)516)) * 255, (
		    ftnlen)1, (ftnlen)255);
	    errint_("#", nrec, (ftnlen)1);
	    sigerr_("SPICE(BADFORMATSPECIFIER)", (ftnlen)25);
	}
    }
    chkout_("RDFFDI", (ftnlen)6);
    return 0;
} /* rdffdi_ */

