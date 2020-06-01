/* prcinf.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/* $Procedure PRCINF ( Process an information request ) */
/* Subroutine */ int prcinf_(char *inftyp, ftnlen inftyp_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    char ch__1[136];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen), errch_(char *, char *,
	     ftnlen, ftnlen);
    extern integer rtrim_(char *, ftnlen);
    extern logical eqstr_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), tostdo_(char *, ftnlen);
    extern logical return_(void);
    static char hlptxt[80*105], verstr[80];
    extern /* Subroutine */ int tkvrsn_(char *, char *, ftnlen, ftnlen);
    static char usgtxt[80*32];

/* $ Abstract */

/*     Process an information request:  display "help," "usage," */
/*     "template, or program version information. */

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
/*     INFTYP     I   Type of information to display. */

/* $ Detailed_Input */

/*     INFTYP         is a character string indicating the type */
/*                    of information to display.  The options are: */

/*                       'HELP'        Dump the introductory */
/*                                     paragraph of the user's guide. */

/*                       'VERSION'     Display the program version */
/*                                     and creation date. */

/*                       'USAGE'       Display a command syntax and */
/*                                     option summary. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the value of INFTYP is not recognized, the error */
/*        SPICE(NOTSUPPORTED) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine centralizes message display functions for */
/*     DSKBRIEF. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     1) For use only within program DSKBRIEF. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/*     DSKBRIEF Version 2.0.0, 07-MAR-2017 (NJB) */

/*        Adapted from MKDSK Version 4.0.0, 22-AUG-2016 (NJB) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Saved variables */


/*     Initial values */

    if (return_()) {
	return 0;
    }
    chkin_("PRCINF", (ftnlen)6);
    if (first) {

/*        This lovely mess was created using Bill Taber's ImportText */
/*        pipe. */

	s_copy(hlptxt, "DSKBRIEF is a command-line utility program that disp"
		"lays a summary of one", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 80, "or more binary DSK files. The program usage is:",
		 (ftnlen)80, (ftnlen)47);
	s_copy(hlptxt + 160, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 240, "      % dskbrief [options] file [file...]", (
		ftnlen)80, (ftnlen)41);
	s_copy(hlptxt + 320, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 400, "   where [file]s are binary DSK files, meta-ke"
		"rnels, or text kernels needed", (ftnlen)80, (ftnlen)75);
	s_copy(hlptxt + 480, "   to support surface name-ID conversion or co"
		"ntaining frame definitions", (ftnlen)80, (ftnlen)72);
	s_copy(hlptxt + 560, "   (FKs), provided in any order. Meta-kernels "
		"may be used to specify sets", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 640, "   of DSK files to summarize.", (ftnlen)80, (
		ftnlen)29);
	s_copy(hlptxt + 720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 800, "   By default, DSKBRIEF summarizes groups of s"
		"egments from each specified", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 880, "   DSK file. Segments having matching attribut"
		"es are grouped together. (See", (ftnlen)80, (ftnlen)75);
	s_copy(hlptxt + 960, "   the section ``DSK segment matching'' below.)"
		, (ftnlen)80, (ftnlen)47);
	s_copy(hlptxt + 1040, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 1120, "   DSKBRIEF can also be commanded to treat al"
		"l DSK files as a single file,", (ftnlen)80, (ftnlen)74);
	s_copy(hlptxt + 1200, "   in which case segments from any file can b"
		"e grouped together if their", (ftnlen)80, (ftnlen)72);
	s_copy(hlptxt + 1280, "   attributes match.", (ftnlen)80, (ftnlen)20);
	s_copy(hlptxt + 1360, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 1440, "   The user can command DSKBRIEF to display s"
		"egment-by-segment summaries", (ftnlen)80, (ftnlen)72);
	s_copy(hlptxt + 1520, "   rather than grouped summaries.", (ftnlen)80,
		 (ftnlen)33);
	s_copy(hlptxt + 1600, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 1680, "   The available options are shown below. The"
		" order of options is not", (ftnlen)80, (ftnlen)69);
	s_copy(hlptxt + 1760, "   significant. The option keys must be lower"
		"case as shown below.", (ftnlen)80, (ftnlen)65);
	s_copy(hlptxt + 1840, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 1920, "         -a       Treat all DSK files as a si"
		"ngle file.", (ftnlen)80, (ftnlen)55);
	s_copy(hlptxt + 2000, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 2080, "         -gaps    Display coverage gaps.", (
		ftnlen)80, (ftnlen)40);
	s_copy(hlptxt + 2160, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 2240, "         -ext     Display extended summaries:"
		" these include data type, data", (ftnlen)80, (ftnlen)75);
	s_copy(hlptxt + 2320, "                  class, and time bounds. Thi"
		"s option applies to summaries", (ftnlen)80, (ftnlen)74);
	s_copy(hlptxt + 2400, "                  of groups of DSK segments.", 
		(ftnlen)80, (ftnlen)44);
	s_copy(hlptxt + 2480, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 2560, "         -tg      Require segment time bounds"
		" to match when grouping", (ftnlen)80, (ftnlen)68);
	s_copy(hlptxt + 2640, "                  segments.", (ftnlen)80, (
		ftnlen)27);
	s_copy(hlptxt + 2720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 2800, "         -seg     Display a segment-by-segmen"
		"t summary.", (ftnlen)80, (ftnlen)55);
	s_copy(hlptxt + 2880, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 2960, "         -full    Display a detailed summary "
		"for each segment, including", (ftnlen)80, (ftnlen)72);
	s_copy(hlptxt + 3040, "                  data-type-specific paramete"
		"rs. This option implies a", (ftnlen)80, (ftnlen)70);
	s_copy(hlptxt + 3120, "                  segment-by-segment summary.",
		 (ftnlen)80, (ftnlen)45);
	s_copy(hlptxt + 3200, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 3280, "         -d <n>   Display n significant digit"
		"s of floating point values.", (ftnlen)80, (ftnlen)72);
	s_copy(hlptxt + 3360, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 3440, "         -v       Display the version of this"
		" program.", (ftnlen)80, (ftnlen)54);
	s_copy(hlptxt + 3520, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 3600, "         -h       Display help text.", (ftnlen)
		80, (ftnlen)36);
	s_copy(hlptxt + 3680, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 3760, "         -u       Display usage text.", (
		ftnlen)80, (ftnlen)37);
	s_copy(hlptxt + 3840, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 3920, "   The options can be provided in any order a"
		"nd can appear before, after,", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 4000, "   or intermixed with file names. The case of"
		" option keys is significant:", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 4080, "   they must be lowercase as shown above.", (
		ftnlen)80, (ftnlen)41);
	s_copy(hlptxt + 4160, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 4240, "   All option combinations are valid; however"
		", some options override", (ftnlen)80, (ftnlen)68);
	s_copy(hlptxt + 4320, "   others:", (ftnlen)80, (ftnlen)10);
	s_copy(hlptxt + 4400, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 4480, "       --   The options -full and -seg both o"
		"verride -a.", (ftnlen)80, (ftnlen)56);
	s_copy(hlptxt + 4560, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 4640, "       --   The option -ext has no effect whe"
		"n -full or -seg are present.", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 4720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 4800, "       --   The option -tg invokes the option"
		" -ext.", (ftnlen)80, (ftnlen)51);
	s_copy(hlptxt + 4880, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 4960, "       --   The option -gaps applies to sets "
		"of DSK files only when -a is", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 5040, "            used. It applies to sets of match"
		"ing segments within a given", (ftnlen)80, (ftnlen)72);
	s_copy(hlptxt + 5120, "            DSK file unless -full or -seg are"
		" used.", (ftnlen)80, (ftnlen)51);
	s_copy(hlptxt + 5200, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 5280, "       --   The program terminates after disp"
		"laying the requested", (ftnlen)80, (ftnlen)65);
	s_copy(hlptxt + 5360, "            information when any of -h, -v, o"
		"r -u are present.", (ftnlen)80, (ftnlen)62);
	s_copy(hlptxt + 5440, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 5520, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 5600, "DSK segment matching", (ftnlen)80, (ftnlen)20);
	s_copy(hlptxt + 5680, "---------------------------------------------"
		"-----------", (ftnlen)80, (ftnlen)56);
	s_copy(hlptxt + 5760, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 5840, "   When DSKBRIEF summarizes groups of segment"
		"s, either within a single DSK", (ftnlen)80, (ftnlen)74);
	s_copy(hlptxt + 5920, "   file, or taken over all specified DSK file"
		"s, the set of segments is", (ftnlen)80, (ftnlen)70);
	s_copy(hlptxt + 6000, "   partitioned into subsets having matching a"
		"ttributes. Summaries are", (ftnlen)80, (ftnlen)69);
	s_copy(hlptxt + 6080, "   produced for these matching subsets.", (
		ftnlen)80, (ftnlen)39);
	s_copy(hlptxt + 6160, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 6240, "   DSK segments ``match'' if they have the sa"
		"me", (ftnlen)80, (ftnlen)47);
	s_copy(hlptxt + 6320, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 6400, "       --   Body", (ftnlen)80, (ftnlen)16);
	s_copy(hlptxt + 6480, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 6560, "       --   Surface", (ftnlen)80, (ftnlen)19);
	s_copy(hlptxt + 6640, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 6720, "       --   Reference frame", (ftnlen)80, (
		ftnlen)27);
	s_copy(hlptxt + 6800, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 6880, "       --   Coordinate system", (ftnlen)80, (
		ftnlen)29);
	s_copy(hlptxt + 6960, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 7040, "       --   Coordinate system parameters, if "
		"applicable", (ftnlen)80, (ftnlen)55);
	s_copy(hlptxt + 7120, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 7200, "       --   Data type", (ftnlen)80, (ftnlen)21)
		;
	s_copy(hlptxt + 7280, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 7360, "       --   Data class", (ftnlen)80, (ftnlen)
		22);
	s_copy(hlptxt + 7440, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 7520, "   Optionally segment time bounds can be adde"
		"d to the list of attributes", (ftnlen)80, (ftnlen)72);
	s_copy(hlptxt + 7600, "   that must match in order for segments to b"
		"e grouped. The", (ftnlen)80, (ftnlen)59);
	s_copy(hlptxt + 7680, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 7760, "      -tg", (ftnlen)80, (ftnlen)9);
	s_copy(hlptxt + 7840, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 7920, "   option invokes this behavior.", (ftnlen)80, 
		(ftnlen)32);
	s_copy(hlptxt + 8000, " ", (ftnlen)80, (ftnlen)1);
	s_copy(hlptxt + 8080, "   Coordinate bounds displayed for such a gro"
		"up are minimum and maximum", (ftnlen)80, (ftnlen)71);
	s_copy(hlptxt + 8160, "   values, taken over the group. It is possib"
		"le for coverage gaps to exist", (ftnlen)80, (ftnlen)74);
	s_copy(hlptxt + 8240, "   within these bounds. The gaps are not disp"
		"layed by default; the option", (ftnlen)80, (ftnlen)73);
	s_copy(hlptxt + 8320, "   -gaps causes DSKBRIEF to display them.", (
		ftnlen)80, (ftnlen)41);
	s_copy(usgtxt, "   DSKBRIEF is a command-line utility program that d"
		"isplays a summary of", (ftnlen)80, (ftnlen)72);
	s_copy(usgtxt + 80, "   one or more binary DSK files. The program us"
		"age is:", (ftnlen)80, (ftnlen)54);
	s_copy(usgtxt + 160, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 240, "      % dskbrief [options] file [file...]", (
		ftnlen)80, (ftnlen)41);
	s_copy(usgtxt + 320, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 400, "   The available options are shown below. The "
		"order of options is not", (ftnlen)80, (ftnlen)69);
	s_copy(usgtxt + 480, "   significant. The option keys must be lowerc"
		"ase as shown below.", (ftnlen)80, (ftnlen)65);
	s_copy(usgtxt + 560, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 640, "         -a       Treat all DSK files as a sin"
		"gle file.", (ftnlen)80, (ftnlen)55);
	s_copy(usgtxt + 720, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 800, "         -gaps    Display coverage gaps. Appli"
		"es only when -a is used.", (ftnlen)80, (ftnlen)70);
	s_copy(usgtxt + 880, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 960, "         -ext     Display extended summaries: "
		"these include data type, data", (ftnlen)80, (ftnlen)75);
	s_copy(usgtxt + 1040, "                  class, and time bounds. Thi"
		"s option applies to summaries", (ftnlen)80, (ftnlen)74);
	s_copy(usgtxt + 1120, "                  of groups of DSK segments.", 
		(ftnlen)80, (ftnlen)44);
	s_copy(usgtxt + 1200, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 1280, "         -tg      Require segment time bounds"
		" to match when grouping", (ftnlen)80, (ftnlen)68);
	s_copy(usgtxt + 1360, "                  segments.", (ftnlen)80, (
		ftnlen)27);
	s_copy(usgtxt + 1440, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 1520, "         -seg     Display a segment-by-segmen"
		"t summary.", (ftnlen)80, (ftnlen)55);
	s_copy(usgtxt + 1600, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 1680, "         -full    Display a detailed summary "
		"for each segment, including", (ftnlen)80, (ftnlen)72);
	s_copy(usgtxt + 1760, "                  data-type-specific paramete"
		"rs. This option implies a", (ftnlen)80, (ftnlen)70);
	s_copy(usgtxt + 1840, "                  segment-by-segment summary.",
		 (ftnlen)80, (ftnlen)45);
	s_copy(usgtxt + 1920, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 2000, "         -d <n>   Display n significant digit"
		"s of floating point values.", (ftnlen)80, (ftnlen)72);
	s_copy(usgtxt + 2080, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 2160, "         -v       Display the version of this"
		" program.", (ftnlen)80, (ftnlen)54);
	s_copy(usgtxt + 2240, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 2320, "         -h       Display help text.", (ftnlen)
		80, (ftnlen)36);
	s_copy(usgtxt + 2400, " ", (ftnlen)80, (ftnlen)1);
	s_copy(usgtxt + 2480, "         -u       Display usage text.", (
		ftnlen)80, (ftnlen)37);
	first = FALSE_;
    }
    if (eqstr_(inftyp, "VERSION", inftyp_len, (ftnlen)7)) {

/*        Create and display "version" message. */

	tkvrsn_("TOOLKIT", verstr, (ftnlen)7, (ftnlen)80);
	tostdo_(" ", (ftnlen)1);
/* Writing concatenation */
	i__1[0] = 56, a__1[0] = "DSKBRIEF Program; Ver. 2.0.0, 07-MAR-2017; "
		"Toolkit Ver. ";
	i__1[1] = rtrim_(verstr, (ftnlen)80), a__1[1] = verstr;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)136);
	tostdo_(ch__1, rtrim_(verstr, (ftnlen)80) + 56);
	tostdo_(" ", (ftnlen)1);
    } else if (eqstr_(inftyp, "HELP", inftyp_len, (ftnlen)4)) {
	for (i__ = 1; i__ <= 105; ++i__) {
	    tostdo_(hlptxt + ((i__2 = i__ - 1) < 105 && 0 <= i__2 ? i__2 : 
		    s_rnge("hlptxt", i__2, "prcinf_", (ftnlen)438)) * 80, (
		    ftnlen)80);
	}
	tostdo_(" ", (ftnlen)1);
    } else if (eqstr_(inftyp, "USAGE", inftyp_len, (ftnlen)5)) {
	for (i__ = 1; i__ <= 32; ++i__) {
	    tostdo_(usgtxt + ((i__2 = i__ - 1) < 32 && 0 <= i__2 ? i__2 : 
		    s_rnge("usgtxt", i__2, "prcinf_", (ftnlen)448)) * 80, (
		    ftnlen)80);
	}
	tostdo_(" ", (ftnlen)1);
    } else {

/*        We shouldn't arrive here. */

	setmsg_("Informational message type # is not supported.", (ftnlen)46);
	errch_("#", inftyp, (ftnlen)1, inftyp_len);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("PRCINF", (ftnlen)6);
	return 0;
    }
    chkout_("PRCINF", (ftnlen)6);
    return 0;
} /* prcinf_ */

