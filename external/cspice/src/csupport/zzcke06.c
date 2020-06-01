/* zzcke06.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static doublereal c_b58 = -.5;

/* $Procedure      ZZCKE06 ( C-Kernel, evaluate, type 6 ) */
/* Subroutine */ int zzcke06_(doublereal *record, doublereal *qstate, 
	doublereal *clkout)
{
    /* Initialized data */

    static integer pktszs[4] = { 8,4,14,7 };

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *), s_rnge(char *, integer, char *, integer);

    /* Local variables */
    doublereal mags, qneg[4], rate;
    integer from;
    extern /* Subroutine */ int vequ_(doublereal *, doublereal *);
    doublereal work[1360]	/* was [680][2] */;
    integer i__, j, n;
    doublereal q[4];
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    doublereal vbuff[6];
    extern /* Subroutine */ int moved_(doublereal *, integer *, doublereal *),
	     errdp_(char *, doublereal *, ftnlen), vsclg_(doublereal *, 
	    doublereal *, integer *, doublereal *);
    doublereal state[8];
    extern doublereal vdotg_(doublereal *, doublereal *, integer *);
    extern /* Subroutine */ int vsubg_(doublereal *, doublereal *, integer *, 
	    doublereal *);
    doublereal dq[4], av[3], ds[4];
    integer to;
    doublereal locrec[340], sclddq[4];
    extern /* Subroutine */ int lgrind_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal sclkdp, radtrm[4];
    integer packsz;
    extern /* Subroutine */ int sigerr_(char *, ftnlen);
    extern doublereal lgrint_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), vdistg_(doublereal *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), chkout_(char *, ftnlen), vminug_(doublereal *,
	     integer *, doublereal *);
    extern doublereal vnormg_(doublereal *, integer *);
    extern /* Subroutine */ int xpsgip_(integer *, integer *, doublereal *), 
	    vsclip_(doublereal *, doublereal *), hrmint_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern logical return_(void);
    integer newptr, xstart, subtyp, ystart, prvptr;
    doublereal qav[4];
    extern /* Subroutine */ int qxq_(doublereal *, doublereal *, doublereal *)
	    ;

/* $ Abstract */

/*     SPICE Private routine intended solely for the support of SPICE */
/*     routines.  Users should not call this routine directly due */
/*     to the volatile nature of this routine. */

/*     Evaluate a single data record from a type 6 CK segment. The */
/*     output is expressed as an interpolated unit quaternion and */
/*     quaternion derivative rather than as a C-matrix and angular */
/*     velocity vector. */

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

/*     CK */

/* $ Keywords */

/*     POINTING */

/* $ Declarations */
/* $ Abstract */

/*     Declare parameters specific to CK type 06. */

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

/*     CK */

/* $ Keywords */

/*     CK */

/* $ Restrictions */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman      (JPL) */
/*     B.V. Semenov      (JPL) */

/* $ Literature_References */

/*     None. */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 10-MAR-2014 (NJB) (BVS) */

/* -& */

/*     Maximum polynomial degree supported by the current */
/*     implementation of this CK type. */


/*     Integer code indicating `true': */


/*     Integer code indicating `false': */


/*     CK type 6 subtype codes: */


/*     Subtype 0:  Hermite interpolation, 8-element packets. Quaternion */
/*                 and quaternion derivatives only, no angular velocity */
/*                 vector provided. Quaternion elements are listed */
/*                 first, followed by derivatives. Angular velocity is */
/*                 derived from the quaternions and quaternion */
/*                 derivatives. */


/*     Subtype 1:  Lagrange interpolation, 4-element packets. Quaternion */
/*                 only. Angular velocity is derived by differentiating */
/*                 the interpolating polynomials. */


/*     Subtype 2:  Hermite interpolation, 14-element packets. */
/*                 Quaternion and angular angular velocity vector, as */
/*                 well as derivatives of each, are provided. The */
/*                 quaternion comes first, then quaternion derivatives, */
/*                 then angular velocity and its derivatives. */


/*     Subtype 3:  Lagrange interpolation, 7-element packets. Quaternion */
/*                 and angular velocity vector provided.  The quaternion */
/*                 comes first. */


/*     Number of subtypes: */


/*     Packet sizes associated with the various subtypes: */


/*     Maximum packet size for type 6: */


/*     Minimum packet size for type 6: */


/*     The CKPFS record size declared in ckparam.inc must be at least as */
/*     large as the maximum possible size of a CK type 6 record. */

/*     The largest possible CK type 6 record has subtype 3 (note that */
/*     records of subtype 2 have half as many epochs as those of subtype */
/*     3, for a given polynomial degree). A subtype 3 record contains */

/*        - The evaluation epoch */
/*        - The subtype and packet count */
/*        - MAXDEG+1 packets of size C06PS3 */
/*        - MAXDEG+1 time tags */


/*     End of file ck06.inc. */

/* $ Abstract */

/*     Declarations of the CK data type specific and general CK low */
/*     level routine parameters. */

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

/*     CK.REQ */

/* $ Keywords */

/*     CK */

/* $ Restrictions */

/*     1) If new CK types are added, the size of the record passed */
/*        between CKRxx and CKExx must be registered as separate */
/*        parameter. If this size will be greater than current value */
/*        of the CKMRSZ parameter (which specifies the maximum record */
/*        size for the record buffer used inside CKPFS) then it should */
/*        be assigned to CKMRSZ as a new value. */

/* $ Author_and_Institution */

/*     N.J. Bachman      (JPL) */
/*     B.V. Semenov      (JPL) */

/* $ Literature_References */

/*     CK Required Reading. */

/* $ Version */

/* -    SPICELIB Version 3.0.0, 27-JAN-2014 (NJB) */

/*        Updated to support CK type 6. Maximum degree for */
/*        type 5 was updated to be consistent with the */
/*        maximum degree for type 6. */

/* -    SPICELIB Version 2.0.0, 19-AUG-2002 (NJB) */

/*        Updated to support CK type 5. */

/* -    SPICELIB Version 1.0.0, 05-APR-1999 (BVS) */

/* -& */

/*     Number of quaternion components and number of quaternion and */
/*     angular rate components together. */


/*     CK Type 1 parameters: */

/*     CK1DTP   CK data type 1 ID; */

/*     CK1RSZ   maximum size of a record passed between CKR01 */
/*              and CKE01. */


/*     CK Type 2 parameters: */

/*     CK2DTP   CK data type 2 ID; */

/*     CK2RSZ   maximum size of a record passed between CKR02 */
/*              and CKE02. */


/*     CK Type 3 parameters: */

/*     CK3DTP   CK data type 3 ID; */

/*     CK3RSZ   maximum size of a record passed between CKR03 */
/*              and CKE03. */


/*     CK Type 4 parameters: */

/*     CK4DTP   CK data type 4 ID; */

/*     CK4PCD   parameter defining integer to DP packing schema that */
/*              is applied when seven number integer array containing */
/*              polynomial degrees for quaternion and angular rate */
/*              components packed into a single DP number stored in */
/*              actual CK records in a file; the value of must not be */
/*              changed or compatibility with existing type 4 CK files */
/*              will be lost. */

/*     CK4MXD   maximum Chebychev polynomial degree allowed in type 4 */
/*              records; the value of this parameter must never exceed */
/*              value of the CK4PCD; */

/*     CK4SFT   number of additional DPs, which are not polynomial */
/*              coefficients, located at the beginning of a type 4 */
/*              CK record that passed between routines CKR04 and CKE04; */

/*     CK4RSZ   maximum size of type 4 CK record passed between CKR04 */
/*              and CKE04; CK4RSZ is computed as follows: */

/*                 CK4RSZ = ( CK4MXD + 1 ) * QAVSIZ + CK4SFT */


/*     CK Type 5 parameters: */


/*     CK5DTP   CK data type 5 ID; */

/*     CK5MXD   maximum polynomial degree allowed in type 5 */
/*              records. */

/*     CK5MET   number of additional DPs, which are not polynomial */
/*              coefficients, located at the beginning of a type 5 */
/*              CK record that passed between routines CKR05 and CKE05; */

/*     CK5MXP   maximum packet size for any subtype.  Subtype 2 */
/*              has the greatest packet size, since these packets */
/*              contain a quaternion, its derivative, an angular */
/*              velocity vector, and its derivative.  See ck05.inc */
/*              for a description of the subtypes. */

/*     CK5RSZ   maximum size of type 5 CK record passed between CKR05 */
/*              and CKE05; CK5RSZ is computed as follows: */

/*                 CK5RSZ = ( CK5MXD + 1 ) * CK5MXP + CK5MET */


/*     CK Type 6 parameters: */


/*     CK6DTP   CK data type 6 ID; */

/*     CK6MXD   maximum polynomial degree allowed in type 6 */
/*              records. */

/*     CK6MET   number of additional DPs, which are not polynomial */
/*              coefficients, located at the beginning of a type 6 */
/*              CK record that passed between routines CKR06 and CKE06; */

/*     CK6MXP   maximum packet size for any subtype.  Subtype 2 */
/*              has the greatest packet size, since these packets */
/*              contain a quaternion, its derivative, an angular */
/*              velocity vector, and its derivative.  See ck06.inc */
/*              for a description of the subtypes. */

/*     CK6RSZ   maximum size of type 6 CK record passed between CKR06 */
/*              and CKE06; CK6RSZ is computed as follows: */

/*                 CK6RSZ = CK6MET + ( CK6MXD + 1 ) * ( CK6PS3 + 1 ) */

/*              where CK6PS3 is equal to the parameter CK06PS3 defined */
/*              in ck06.inc. Note that the subtype having the largest */
/*              packet size (subtype 2) does not give rise to the */
/*              largest record size, because that type is Hermite and */
/*              requires half the window size used by subtype 3 for a */
/*              given polynomial degree. */


/*     The parameter CK6PS3 must be in sync with C06PS3 defined in */
/*     ck06.inc. */



/*     Maximum record size that can be handled by CKPFS. This value */
/*     must be set to the maximum of all CKxRSZ parameters (currently */
/*     CK5RSZ.) */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     RECORD    I-O  Data type 6 record. */
/*     QSTATE     O   Interpolated output record. */
/*     CLKOUT     O   SCLK associated with input record. */

/* $ Detailed_Input */

/*     RECORD      is a record from a type 6 CK segment which, when */
/*                 evaluated at the epoch contained in its first */
/*                 element, will give the attitude and angular velocity */
/*                 of a spacecraft structure or instrument relative to a */
/*                 base reference frame. */

/*                 The structure of the record is as follows: */

/*                    +----------------------+ */
/*                    | evaluation epoch     | */
/*                    +----------------------+ */
/*                    | subtype code         | */
/*                    +----------------------+ */
/*                    | number of packets (n)| */
/*                    +----------------------+ */
/*                    | nominal SCLK rate    | */
/*                    +----------------------+ */
/*                    | packet 1             | */
/*                    +----------------------+ */
/*                    | packet 2             | */
/*                    +----------------------+ */
/*                             . */
/*                             . */
/*                             . */
/*                    +----------------------+ */
/*                    | packet n             | */
/*                    +----------------------+ */
/*                    | epochs 1--n          | */
/*                    +----------------------+ */

/*                See the CK Required Reading or the include file */
/*                ck06.inc for details on CK type 6 packet contents. */


/* $ Detailed_Output */

/*     RECORD     has been modified due to its use as a workspace array. */
/*                The contents are undefined. */

/*     QSTATE     is an interpolated output record, represented as a unit */
/*                quaternion and its derivative with respect to time. */

/*     CLKOUT     is the encoded SCLK associated with the returned */
/*                C-matrix and angular velocity vector. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the input record contains an unrecognized subtype code, */
/*         the error SPICE(NOTSUPPORTED) is signaled. */

/*     2)  If the record subtype is one for which quaternion derivatives */
/*         are stored (subtypes 0 and 2), and if the Ith quaternion in */
/*         the input record is farther than its negative from the (I-1)st */
/*         quaternion in the record, the error SPICE(BADQUATSIGN) */
/*         is signaled. */

/*         For subtypes 1 and 3, this condition is not considered an */
/*         error: the closer to the preceding quaternion of the two */
/*         quaternion representations is used for interpolation. */

/*     3)  If a zero-magnitude quaternion is produced as a result */
/*         of interpolating quaternions, the error SPICE(DIVIDEBYZERO) */
/*         is signaled. */

/*     4)  If the input record contains a non-positive SCLK rate value, */
/*         the error SPICE(INVALIDSCLKRATE) is signaled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine returns quaternion states that, presuming validity */
/*     of the input record, are suitable for interpolation, unlike */
/*     quaternion states obtained by calling CKE06 and then M2Q. The */
/*     latter method of obtaining quaternions is subject to branch */
/*     singularities. */

/*     The exact format and structure of CK type 6 (MEX/Rosetta Attitude */
/*     file interpolation) CK segments is described in the CK Required */
/*     Reading. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     1)  This routine performs minimal error checking. The input data */
/*         are assumed to have been checked when the source CK file was */
/*         created. */

/*     2)  With the exception of the check described in item 2 of */
/*         the Exceptions section above, the input data are assumed to */
/*         be suitable for the interpolation method specified by the */
/*         input record's subtype and packet count (which implies an */
/*         interpolating polynomial degree). */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */
/*     B.V. Semenov   (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 10-AUG-2015 (NJB) (BVS) */

/* -& */
/* $ Index_Entries */

/*     evaluate type_6 ck_segment */

/* -& */
/* $ Revisions */

/*     None. */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Index of evaluation epoch in record: */


/*     Index of subtype code in record: */


/*     Index of packet count in record: */


/*     Index at which packets start; packet base: */


/*     Local variables */


/*     Saved variables */


/*     Initial values */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    }
    chkin_("ZZCKE06", (ftnlen)7);

/*     Transfer the input record's epoch to the output epoch. */

    *clkout = record[0];

/*     Capture the subtype from the record and set the packet size */
/*     accordingly. */

    subtyp = i_dnnt(&record[1]);
    if (subtyp < 0 || subtyp >= 4) {
	setmsg_("Unexpected CK type 6 subtype # found in type 6 segment.", (
		ftnlen)55);
	errint_("#", &subtyp, (ftnlen)1);
	sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	chkout_("ZZCKE06", (ftnlen)7);
	return 0;
    } else {
	packsz = pktszs[(i__1 = subtyp) < 4 && 0 <= i__1 ? i__1 : s_rnge(
		"pktszs", i__1, "zzcke06_", (ftnlen)318)];
    }

/*     Get the packet count and epoch. */

    n = i_dnnt(&record[2]);
    sclkdp = record[0];

/*     Get the nominal clock rate. */

    rate = record[3];
    if (rate <= 0.) {
	setmsg_("SCLK rate is #; rate must be positive.", (ftnlen)38);
	errdp_("#", &rate, (ftnlen)1);
	sigerr_("SPICE(INVALIDSCLKRATE)", (ftnlen)22);
	chkout_("ZZCKE06", (ftnlen)7);
	return 0;
    }

/*     Adjust quaternion "signs" as necessary to minimize distance */
/*     between successive quaternions. This adjustment is performed */
/*     only for subtypes that don't store quaternion derivatives */
/*     (these are the Lagrange subtypes). */

    if (subtyp == 1 || subtyp == 3) {

/*        For these subtypes, only the quaternions themselves need be */
/*        adjusted. */

/*        PRVPTR is the index of the "previous" quaternion---the one to */
/*        which the successor and its negative will be compared. */

	prvptr = 5;
	i__1 = n;
	for (i__ = 2; i__ <= i__1; ++i__) {

/*           NEWPTR points to the quaternion ahead of the one */
/*           pointed to by PRVPTR. */

	    newptr = packsz * (i__ - 1) + 5;
	    vminug_(&record[newptr - 1], &c__4, qneg);

/*           Replace the Ith quaternion with QNEG if QNEG is closer */
/*           than the current quaternion to the previous quaternion. */

	    if (vdistg_(&record[prvptr - 1], qneg, &c__4) < vdistg_(&record[
		    prvptr - 1], &record[newptr - 1], &c__4)) {
		moved_(qneg, &c__4, &record[newptr - 1]);
	    }
	    prvptr = newptr;
	}
    } else {

/*        For the Hermite types, if the quaternions need to be adjusted, */
/*        we have an error condition. */

/*        PRVPTR is the index of the "previous" quaternion---the one to */
/*        which the successor and its negative will be compared. */

	prvptr = 5;
	i__1 = n;
	for (i__ = 2; i__ <= i__1; ++i__) {

/*           NEWPTR points to the quaternion ahead of the one */
/*           pointed to by PRVPTR. */

	    newptr = packsz * (i__ - 1) + 5;
	    vminug_(&record[newptr - 1], &c__4, qneg);

/*           For the Hermite subtypes, it's an error for the current */
/*           quaternion to be closer to QNEG than to the previous */
/*           quaternion. */

	    if (vdistg_(&record[prvptr - 1], qneg, &c__4) < vdistg_(&record[
		    prvptr - 1], &record[newptr - 1], &c__4)) {
		setmsg_("Quaternion sign error: quaternion at index # in the"
			" input record is farther than its negative from the "
			"preceding quaternion in the record. Quaternion is (#"
			", #, #, #); predecessor is (#, #, #, #). This makes "
			"the quaternion sequence unsuitable for Hermite inter"
			"polation. The quaternions, and if applicable, their "
			"derivatives, must be adjusted before they are passed"
			" to this routine.", (ftnlen)380);
		errint_("#", &i__, (ftnlen)1);
		errdp_("#", &record[newptr - 1], (ftnlen)1);
		errdp_("#", &record[newptr], (ftnlen)1);
		errdp_("#", &record[newptr + 1], (ftnlen)1);
		errdp_("#", &record[newptr + 2], (ftnlen)1);
		errdp_("#", &record[prvptr - 1], (ftnlen)1);
		errdp_("#", &record[prvptr], (ftnlen)1);
		errdp_("#", &record[prvptr + 1], (ftnlen)1);
		errdp_("#", &record[prvptr + 2], (ftnlen)1);
		sigerr_("SPICE(BADQUATSIGN)", (ftnlen)18);
		chkout_("ZZCKE06", (ftnlen)7);
		return 0;
	    }
	    prvptr = newptr;
	}
    }
    if (subtyp == 1) {

/*        We perform Lagrange interpolation on each quaternion */
/*        component, and obtain quaternion derivatives from the */
/*        interpolating polynomials. */

/*        We'll transpose the pointing information in the input record so */
/*        that contiguous pieces of it can be shoved directly into the */
/*        interpolation routine LGRIND. */

	n = i_dnnt(&record[2]);
	xpsgip_(&packsz, &n, &record[4]);

/*        We interpolate each state component in turn. */

	xstart = n * packsz + 5;
	i__1 = packsz;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ystart = n * (i__ - 1) + 5;
	    lgrind_(&n, &record[xstart - 1], &record[ystart - 1], work, &
		    sclkdp, &state[(i__2 = i__ - 1) < 8 && 0 <= i__2 ? i__2 : 
		    s_rnge("state", i__2, "zzcke06_", (ftnlen)468)], &state[(
		    i__3 = i__ + 3) < 8 && 0 <= i__3 ? i__3 : s_rnge("state", 
		    i__3, "zzcke06_", (ftnlen)468)]);
	}

/*        The output quaternion is a unitized version of the */
/*        interpolated state. */

	mags = vnormg_(state, &c__4);
	if (mags == 0.) {
	    setmsg_("Quaternion magnitude at SCLK # was zero.", (ftnlen)40);
	    errdp_("#", &sclkdp, (ftnlen)1);
	    sigerr_("SPICE(DIVIDEBYZERO)", (ftnlen)19);
	    chkout_("ZZCKE06", (ftnlen)7);
	    return 0;
	}
	d__1 = 1. / mags;
	vsclg_(&d__1, state, &c__4, q);

/*        Find the time derivative of the unit quaternion: */
/*        Letting S represent the quaternion portion of STATE, we */
/*        have */

/*           Q = S/||S|| */


/*        Then letting < , > denote the 4-dimensional inner product */
/*        operator, we have */


/*                      d(S)/dt      < Q, d(S)/dt > */
/*           d(Q)/dt =  -------  -   -------------- * Q */
/*                       ||S||            ||S|| */


	moved_(&state[4], &c__4, ds);
	d__1 = 1. / mags;
	vsclg_(&d__1, ds, &c__4, sclddq);
	d__1 = vdotg_(q, ds, &c__4) / mags;
	vsclg_(&d__1, q, &c__4, radtrm);
	vsubg_(sclddq, radtrm, &c__4, dq);

/*        Scale the derivative from 1/tick to 1/second. */

	d__1 = 1. / rate;
	vsclg_(&d__1, dq, &c__4, sclddq);
	moved_(q, &c__4, qstate);
	moved_(sclddq, &c__4, &qstate[4]);
    } else if (subtyp == 3) {

/*        This is the easiest case:  we perform Lagrange interpolation */
/*        on each quaternion or angular velocity component. */

/*        We'll transpose the pointing information in the input record so */
/*        that contiguous pieces of it can be shoved directly into the */
/*        interpolation routine LGRINT.  We allow LGRINT to overwrite */
/*        the state values in the input record, since this saves local */
/*        storage and does no harm.  (See the header of LGRINT for a */
/*        description of its work space usage.) */

	n = i_dnnt(&record[2]);
	xpsgip_(&packsz, &n, &record[4]);

/*        We interpolate each state component in turn. */

	xstart = n * packsz + 5;
	i__1 = packsz;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ystart = n * (i__ - 1) + 5;
	    state[(i__2 = i__ - 1) < 8 && 0 <= i__2 ? i__2 : s_rnge("state", 
		    i__2, "zzcke06_", (ftnlen)555)] = lgrint_(&n, &record[
		    xstart - 1], &record[ystart - 1], locrec, &sclkdp);
	}
	mags = vnormg_(state, &c__4);
	if (mags == 0.) {
	    setmsg_("Quaternion magnitude at SCLK # was zero.", (ftnlen)40);
	    errdp_("#", &sclkdp, (ftnlen)1);
	    sigerr_("SPICE(DIVIDEBYZERO)", (ftnlen)19);
	    chkout_("ZZCKE06", (ftnlen)7);
	    return 0;
	}
	d__1 = 1. / mags;
	vsclg_(&d__1, state, &c__4, q);

/*        The angular velocity already is in units of radians/second. */

	vequ_(&state[4], av);

/*        Convert AV to a quaternion derivative. We have from */
/*        the header of QDQ2AV */

/*                       * */
/*           AV =  -2 * Q  * DQ */

/*        so */

/*           DQ =  -1/2 * Q * AV */


	vsclip_(&c_b58, av);
	qav[0] = 0.;
	vequ_(av, &qav[1]);
	qxq_(q, qav, dq);
	moved_(q, &c__4, qstate);
	moved_(dq, &c__4, &qstate[4]);
    } else {

/*        We have a Hermite-style subtype.  Whether it's subtype 0 */
/*        or 2, we perform Hermite interpolation on the quaternions. */

/*        We interpolate each quaternion component in turn.  Attitude and */
/*        angular velocity are interpolated separately. */

	xstart = packsz * n + 5;
	for (i__ = 1; i__ <= 4; ++i__) {
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {

/*              For the Jth input packet, copy the Ith position and */
/*              velocity components into the local record buffer RECORD. */

/*              In order to perform Hermite interpolation, the */
/*              quaternions and quaternion derivatives must have a */
/*              common time scale. So prior to interpolation, we scale */
/*              the units of the quaternion derivatives from radians/sec */
/*              to radians/tick. */

		from = packsz * (j - 1) + 4 + i__;
		to = (j << 1) - 1;
		locrec[(i__2 = to - 1) < 340 && 0 <= i__2 ? i__2 : s_rnge(
			"locrec", i__2, "zzcke06_", (ftnlen)631)] = record[
			from - 1];
		locrec[(i__2 = to) < 340 && 0 <= i__2 ? i__2 : s_rnge("locrec"
			, i__2, "zzcke06_", (ftnlen)632)] = record[from + 3] *
			 rate;
	    }

/*           Interpolate the Ith quaternion and quaternion derivative */
/*           components. */

	    hrmint_(&n, &record[xstart - 1], locrec, &sclkdp, work, &state[(
		    i__1 = i__ - 1) < 8 && 0 <= i__1 ? i__1 : s_rnge("state", 
		    i__1, "zzcke06_", (ftnlen)640)], &state[(i__2 = i__ + 3) <
		     8 && 0 <= i__2 ? i__2 : s_rnge("state", i__2, "zzcke06_",
		     (ftnlen)640)]);
	}

/*        The output quaternion is a unitized version of the */
/*        interpolated state. */

	mags = vnormg_(state, &c__4);
	if (mags == 0.) {
	    setmsg_("Quaternion magnitude at SCLK # was zero.", (ftnlen)40);
	    errdp_("#", &sclkdp, (ftnlen)1);
	    sigerr_("SPICE(DIVIDEBYZERO)", (ftnlen)19);
	    chkout_("ZZCKE06", (ftnlen)7);
	    return 0;
	}
	d__1 = 1. / mags;
	vsclg_(&d__1, state, &c__4, q);
	if (subtyp == 0) {

/*           Find the time derivative of the unit quaternion: */
/*           Letting S represent the quaternion portion of STATE, we */
/*           have */

/*              Q = S/||S|| */


/*           Then letting < , > denote the 4-dimensional inner product */
/*           operator, we have */


/*                         d(S)/dt      < Q, d(S)/dt > */
/*              d(Q)/dt =  -------  -   -------------- * Q */
/*                          ||S||            ||S|| */


	    moved_(&state[4], &c__4, ds);
	    d__1 = 1. / mags;
	    vsclg_(&d__1, ds, &c__4, sclddq);
	    d__1 = vdotg_(q, ds, &c__4) / mags;
	    vsclg_(&d__1, q, &c__4, radtrm);
	    vsubg_(sclddq, radtrm, &c__4, dq);

/*           Scale the derivative from radians/tick to */
/*           radians/second. */

	    d__1 = 1. / rate;
	    vsclg_(&d__1, dq, &c__4, sclddq);

/*           Store Q and DQ in QSTATE. In the process, */

	    moved_(q, &c__4, qstate);
	    moved_(sclddq, &c__4, &qstate[4]);
	} else {

/*           This is subtype 2; we perform Hermite interpolation on */
/*           the angular velocity and its derivative. */

/*           Now interpolate angular velocity, using separate angular */
/*           velocity data and angular acceleration. */

	    for (i__ = 1; i__ <= 3; ++i__) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {

/*                 For the Jth input packet, copy the Ith angular */
/*                 velocity and angular acceleration components into the */
/*                 local record buffer LOCREC.  Note that, as with */
/*                 quaternion derivatives, we must scale angular */
/*                 acceleration from radians/sec**2 to */
/*                 radians/(sec*tick) before interpolating. We would */
/*                 need to scale the angular acceleration to */
/*                 radians/sec**2 for output, if we were returning this */
/*                 quantity. However, we're returning only angular */
/*                 velocity, which is already in the correct units of */
/*                 radians/second. */

		    from = packsz * (j - 1) + 12 + i__;
		    to = (j << 1) - 1;
		    locrec[(i__2 = to - 1) < 340 && 0 <= i__2 ? i__2 : s_rnge(
			    "locrec", i__2, "zzcke06_", (ftnlen)734)] = 
			    record[from - 1];
		    locrec[(i__2 = to) < 340 && 0 <= i__2 ? i__2 : s_rnge(
			    "locrec", i__2, "zzcke06_", (ftnlen)735)] = 
			    record[from + 2] * rate;
		}

/*              Interpolate the Ith angular velocity and angular */
/*              acceleration components of the attitude. We'll */
/*              capture the result in a temporary buffer, then */
/*              transfer the velocity to the output argument AV. */

		hrmint_(&n, &record[xstart - 1], locrec, &sclkdp, work, &
			vbuff[(i__1 = i__ - 1) < 6 && 0 <= i__1 ? i__1 : 
			s_rnge("vbuff", i__1, "zzcke06_", (ftnlen)745)], &
			vbuff[(i__2 = i__ + 2) < 6 && 0 <= i__2 ? i__2 : 
			s_rnge("vbuff", i__2, "zzcke06_", (ftnlen)745)]);
	    }

/*           Fill in the angular velocity in the output angular */
/*           velocity vector using the results of interpolating */
/*           velocity and acceleration. */

/*           The angular velocity is already in units of */
/*           radians/second. */

	    vequ_(vbuff, av);

/*           Convert AV to a quaternion derivative. We have from */
/*           the header of QDQ2AV */

/*                          * */
/*              AV =  -2 * Q  * DQ */

/*           so */

/*              DQ =  -1/2 * Q * AV */


	    vsclip_(&c_b58, av);
	    qav[0] = 0.;
	    vequ_(av, &qav[1]);
	    qxq_(q, qav, dq);
	    moved_(q, &c__4, qstate);
	    moved_(dq, &c__4, &qstate[4]);
	}

/*        We've handled the type 0 and type 2 cases. */


/*        We've computed the angular velocity AV for the Hermite */
/*        subtypes, if a.v. was requested. */

    }

/*     We've handled all four subtypes. */

    chkout_("ZZCKE06", (ftnlen)7);
    return 0;
} /* zzcke06_ */

