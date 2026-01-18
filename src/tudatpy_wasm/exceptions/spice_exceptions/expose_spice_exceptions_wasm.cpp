/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../wasm_module.h"

#include <tudat/interface/spice/spiceException.h>

namespace te = tudat::exceptions;

WASM_MODULE_PATH("exceptions_spice_exceptions")

EMSCRIPTEN_BINDINGS(tudatpy_exceptions_spice_exceptions) {
    using namespace emscripten;

    // Register all SPICE exception types for JavaScript interop
    // Note: In WASM, these exceptions are thrown as JavaScript Error objects
    // The registration allows type-checking and proper error messages

    // Base SpiceError class
    class_<te::SpiceError>("exceptions_SpiceError")
        .function("what", optional_override([](const te::SpiceError& self) {
            return std::string(self.what());
        }));

    // ADDRESSOUTOFBOUNDS exception
    class_<te::SpiceADDRESSOUTOFBOUNDS, base<te::SpiceError>>("exceptions_SpiceADDRESSOUTOFBOUNDS");

    // AGENTLISTOVERFLOW exception
    class_<te::SpiceAGENTLISTOVERFLOW, base<te::SpiceError>>("exceptions_SpiceAGENTLISTOVERFLOW");

    // ALLGONE exception
    class_<te::SpiceALLGONE, base<te::SpiceError>>("exceptions_SpiceALLGONE");

    // AMBIGTEMPL exception
    class_<te::SpiceAMBIGTEMPL, base<te::SpiceError>>("exceptions_SpiceAMBIGTEMPL");

    // ARRAYSIZEMISMATCH exception
    class_<te::SpiceARRAYSIZEMISMATCH, base<te::SpiceError>>("exceptions_SpiceARRAYSIZEMISMATCH");

    // ARRAYTOOSMALL exception
    class_<te::SpiceARRAYTOOSMALL, base<te::SpiceError>>("exceptions_SpiceARRAYTOOSMALL");

    // AVALOUTOFRANGE exception
    class_<te::SpiceAVALOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceAVALOUTOFRANGE");

    // AXISUNDERFLOW exception
    class_<te::SpiceAXISUNDERFLOW, base<te::SpiceError>>("exceptions_SpiceAXISUNDERFLOW");

    // BADACTION exception
    class_<te::SpiceBADACTION, base<te::SpiceError>>("exceptions_SpiceBADACTION");

    // BADADDRESS exception
    class_<te::SpiceBADADDRESS, base<te::SpiceError>>("exceptions_SpiceBADADDRESS");

    // BADANGLE exception
    class_<te::SpiceBADANGLE, base<te::SpiceError>>("exceptions_SpiceBADANGLE");

    // BADANGLEUNITS exception
    class_<te::SpiceBADANGLEUNITS, base<te::SpiceError>>("exceptions_SpiceBADANGLEUNITS");

    // BADANGRATEERROR exception
    class_<te::SpiceBADANGRATEERROR, base<te::SpiceError>>("exceptions_SpiceBADANGRATEERROR");

    // BADANGULARRATE exception
    class_<te::SpiceBADANGULARRATE, base<te::SpiceError>>("exceptions_SpiceBADANGULARRATE");

    // BADANGULARRATEFLAG exception
    class_<te::SpiceBADANGULARRATEFLAG, base<te::SpiceError>>("exceptions_SpiceBADANGULARRATEFLAG");

    // BADARCHITECTURE exception
    class_<te::SpiceBADARCHITECTURE, base<te::SpiceError>>("exceptions_SpiceBADARCHITECTURE");

    // BADARRAYSIZE exception
    class_<te::SpiceBADARRAYSIZE, base<te::SpiceError>>("exceptions_SpiceBADARRAYSIZE");

    // BADATTIME exception
    class_<te::SpiceBADATTIME, base<te::SpiceError>>("exceptions_SpiceBADATTIME");

    // BADATTRIBUTE exception
    class_<te::SpiceBADATTRIBUTE, base<te::SpiceError>>("exceptions_SpiceBADATTRIBUTE");

    // BADATTRIBUTES exception
    class_<te::SpiceBADATTRIBUTES, base<te::SpiceError>>("exceptions_SpiceBADATTRIBUTES");

    // BADAUVALUE exception
    class_<te::SpiceBADAUVALUE, base<te::SpiceError>>("exceptions_SpiceBADAUVALUE");

    // BADAVFLAG exception
    class_<te::SpiceBADAVFLAG, base<te::SpiceError>>("exceptions_SpiceBADAVFLAG");

    // BADAVFRAMEFLAG exception
    class_<te::SpiceBADAVFRAMEFLAG, base<te::SpiceError>>("exceptions_SpiceBADAVFRAMEFLAG");

    // BADAXIS exception
    class_<te::SpiceBADAXIS, base<te::SpiceError>>("exceptions_SpiceBADAXIS");

    // BADAXISLENGTH exception
    class_<te::SpiceBADAXISLENGTH, base<te::SpiceError>>("exceptions_SpiceBADAXISLENGTH");

    // BADAXISNUMBERS exception
    class_<te::SpiceBADAXISNUMBERS, base<te::SpiceError>>("exceptions_SpiceBADAXISNUMBERS");

    // BADBLOCKSIZE exception
    class_<te::SpiceBADBLOCKSIZE, base<te::SpiceError>>("exceptions_SpiceBADBLOCKSIZE");

    // BADBODYID exception
    class_<te::SpiceBADBODYID, base<te::SpiceError>>("exceptions_SpiceBADBODYID");

    // BADBORESIGHTSPEC exception
    class_<te::SpiceBADBORESIGHTSPEC, base<te::SpiceError>>("exceptions_SpiceBADBORESIGHTSPEC");

    // BADBOUNDARY exception
    class_<te::SpiceBADBOUNDARY, base<te::SpiceError>>("exceptions_SpiceBADBOUNDARY");

    // BADCATALOGFILE exception
    class_<te::SpiceBADCATALOGFILE, base<te::SpiceError>>("exceptions_SpiceBADCATALOGFILE");

    // BADCENTERNAME exception
    class_<te::SpiceBADCENTERNAME, base<te::SpiceError>>("exceptions_SpiceBADCENTERNAME");

    // BADCHECKFLAG exception
    class_<te::SpiceBADCHECKFLAG, base<te::SpiceError>>("exceptions_SpiceBADCHECKFLAG");

    // BADCKTYPESPEC exception
    class_<te::SpiceBADCKTYPESPEC, base<te::SpiceError>>("exceptions_SpiceBADCKTYPESPEC");

    // BADCOARSEVOXSCALE exception
    class_<te::SpiceBADCOARSEVOXSCALE, base<te::SpiceError>>("exceptions_SpiceBADCOARSEVOXSCALE");

    // BADCOLUMDECL exception
    class_<te::SpiceBADCOLUMDECL, base<te::SpiceError>>("exceptions_SpiceBADCOLUMDECL");

    // BADCOLUMNCOUNT exception
    class_<te::SpiceBADCOLUMNCOUNT, base<te::SpiceError>>("exceptions_SpiceBADCOLUMNCOUNT");

    // BADCOLUMNDECL exception
    class_<te::SpiceBADCOLUMNDECL, base<te::SpiceError>>("exceptions_SpiceBADCOLUMNDECL");

    // BADCOMMENTAREA exception
    class_<te::SpiceBADCOMMENTAREA, base<te::SpiceError>>("exceptions_SpiceBADCOMMENTAREA");

    // BADCOMPNUMBER exception
    class_<te::SpiceBADCOMPNUMBER, base<te::SpiceError>>("exceptions_SpiceBADCOMPNUMBER");

    // BADCOORDBOUNDS exception
    class_<te::SpiceBADCOORDBOUNDS, base<te::SpiceError>>("exceptions_SpiceBADCOORDBOUNDS");

    // BADCOORDSYS exception
    class_<te::SpiceBADCOORDSYS, base<te::SpiceError>>("exceptions_SpiceBADCOORDSYS");

    // BADCURVETYPE exception
    class_<te::SpiceBADCURVETYPE, base<te::SpiceError>>("exceptions_SpiceBADCURVETYPE");

    // BADDAFTRANSFERFILE exception
    class_<te::SpiceBADDAFTRANSFERFILE, base<te::SpiceError>>("exceptions_SpiceBADDAFTRANSFERFILE");

    // BADDASCOMMENTAREA exception
    class_<te::SpiceBADDASCOMMENTAREA, base<te::SpiceError>>("exceptions_SpiceBADDASCOMMENTAREA");

    // BADDASDIRECTORY exception
    class_<te::SpiceBADDASDIRECTORY, base<te::SpiceError>>("exceptions_SpiceBADDASDIRECTORY");

    // BADDASFILE exception
    class_<te::SpiceBADDASFILE, base<te::SpiceError>>("exceptions_SpiceBADDASFILE");

    // BADDASTRANSFERFILE exception
    class_<te::SpiceBADDASTRANSFERFILE, base<te::SpiceError>>("exceptions_SpiceBADDASTRANSFERFILE");

    // BADDATALINE exception
    class_<te::SpiceBADDATALINE, base<te::SpiceError>>("exceptions_SpiceBADDATALINE");

    // BADDATAORDERTOKEN exception
    class_<te::SpiceBADDATAORDERTOKEN, base<te::SpiceError>>("exceptions_SpiceBADDATAORDERTOKEN");

    // BADDATATYPE exception
    class_<te::SpiceBADDATATYPE, base<te::SpiceError>>("exceptions_SpiceBADDATATYPE");

    // BADDATATYPEFLAG exception
    class_<te::SpiceBADDATATYPEFLAG, base<te::SpiceError>>("exceptions_SpiceBADDATATYPEFLAG");

    // BADDEFAULTVALUE exception
    class_<te::SpiceBADDEFAULTVALUE, base<te::SpiceError>>("exceptions_SpiceBADDEFAULTVALUE");

    // BADDESCRTIMES exception
    class_<te::SpiceBADDESCRTIMES, base<te::SpiceError>>("exceptions_SpiceBADDESCRTIMES");

    // BADDIMENSION exception
    class_<te::SpiceBADDIMENSION, base<te::SpiceError>>("exceptions_SpiceBADDIMENSION");

    // BADDIMENSIONS exception
    class_<te::SpiceBADDIMENSIONS, base<te::SpiceError>>("exceptions_SpiceBADDIMENSIONS");

    // BADDIRECTION exception
    class_<te::SpiceBADDIRECTION, base<te::SpiceError>>("exceptions_SpiceBADDIRECTION");

    // BADDOUBLEPRECISION exception
    class_<te::SpiceBADDOUBLEPRECISION, base<te::SpiceError>>("exceptions_SpiceBADDOUBLEPRECISION");

    // BADDOWNSAMPLINGTOL exception
    class_<te::SpiceBADDOWNSAMPLINGTOL, base<te::SpiceError>>("exceptions_SpiceBADDOWNSAMPLINGTOL");

    // BADECCENTRICITY exception
    class_<te::SpiceBADECCENTRICITY, base<te::SpiceError>>("exceptions_SpiceBADECCENTRICITY");

    // BADENDPOINTS exception
    class_<te::SpiceBADENDPOINTS, base<te::SpiceError>>("exceptions_SpiceBADENDPOINTS");

    // BADEULERANGLEUNITS exception
    class_<te::SpiceBADEULERANGLEUNITS, base<te::SpiceError>>("exceptions_SpiceBADEULERANGLEUNITS");

    // BADFILEFORMAT exception
    class_<te::SpiceBADFILEFORMAT, base<te::SpiceError>>("exceptions_SpiceBADFILEFORMAT");

    // BADFILENAME exception
    class_<te::SpiceBADFILENAME, base<te::SpiceError>>("exceptions_SpiceBADFILENAME");

    // BADFINEVOXELSCALE exception
    class_<te::SpiceBADFINEVOXELSCALE, base<te::SpiceError>>("exceptions_SpiceBADFINEVOXELSCALE");

    // BADFORMATSPECIFIER exception
    class_<te::SpiceBADFORMATSPECIFIER, base<te::SpiceError>>("exceptions_SpiceBADFORMATSPECIFIER");

    // BADFRAME exception
    class_<te::SpiceBADFRAME, base<te::SpiceError>>("exceptions_SpiceBADFRAME");

    // BADFRAMECLASS exception
    class_<te::SpiceBADFRAMECLASS, base<te::SpiceError>>("exceptions_SpiceBADFRAMECLASS");

    // BADFRAMECOUNT exception
    class_<te::SpiceBADFRAMECOUNT, base<te::SpiceError>>("exceptions_SpiceBADFRAMECOUNT");

    // BADFRAMESPEC exception
    class_<te::SpiceBADFRAMESPEC, base<te::SpiceError>>("exceptions_SpiceBADFRAMESPEC");

    // BADFROMTIME exception
    class_<te::SpiceBADFROMTIME, base<te::SpiceError>>("exceptions_SpiceBADFROMTIME");

    // BADFROMTIMESYSTEM exception
    class_<te::SpiceBADFROMTIMESYSTEM, base<te::SpiceError>>("exceptions_SpiceBADFROMTIMESYSTEM");

    // BADFROMTIMETYPE exception
    class_<te::SpiceBADFROMTIMETYPE, base<te::SpiceError>>("exceptions_SpiceBADFROMTIMETYPE");

    // BADGEOMETRY exception
    class_<te::SpiceBADGEOMETRY, base<te::SpiceError>>("exceptions_SpiceBADGEOMETRY");

    // BADGM exception
    class_<te::SpiceBADGM, base<te::SpiceError>>("exceptions_SpiceBADGM");

    // BADHARDSPACE exception
    class_<te::SpiceBADHARDSPACE, base<te::SpiceError>>("exceptions_SpiceBADHARDSPACE");

    // BADHERMITDEGREE exception
    class_<te::SpiceBADHERMITDEGREE, base<te::SpiceError>>("exceptions_SpiceBADHERMITDEGREE");

    // BADINDEX exception
    class_<te::SpiceBADINDEX, base<te::SpiceError>>("exceptions_SpiceBADINDEX");

    // BADINITSTATE exception
    class_<te::SpiceBADINITSTATE, base<te::SpiceError>>("exceptions_SpiceBADINITSTATE");

    // BADINPUTDATALINE exception
    class_<te::SpiceBADINPUTDATALINE, base<te::SpiceError>>("exceptions_SpiceBADINPUTDATALINE");

    // BADINPUTETTIME exception
    class_<te::SpiceBADINPUTETTIME, base<te::SpiceError>>("exceptions_SpiceBADINPUTETTIME");

    // BADINPUTTYPE exception
    class_<te::SpiceBADINPUTTYPE, base<te::SpiceError>>("exceptions_SpiceBADINPUTTYPE");

    // BADINPUTUTCTIME exception
    class_<te::SpiceBADINPUTUTCTIME, base<te::SpiceError>>("exceptions_SpiceBADINPUTUTCTIME");

    // BADINSTRUMENTID exception
    class_<te::SpiceBADINSTRUMENTID, base<te::SpiceError>>("exceptions_SpiceBADINSTRUMENTID");

    // BADINTEGER exception
    class_<te::SpiceBADINTEGER, base<te::SpiceError>>("exceptions_SpiceBADINTEGER");

    // BADKERNELTYPE exception
    class_<te::SpiceBADKERNELTYPE, base<te::SpiceError>>("exceptions_SpiceBADKERNELTYPE");

    // BADKERNELVARTYPE exception
    class_<te::SpiceBADKERNELVARTYPE, base<te::SpiceError>>("exceptions_SpiceBADKERNELVARTYPE");

    // BADLAGRANGEDEGREE exception
    class_<te::SpiceBADLAGRANGEDEGREE, base<te::SpiceError>>("exceptions_SpiceBADLAGRANGEDEGREE");

    // BADLATITUDEBOUNDS exception
    class_<te::SpiceBADLATITUDEBOUNDS, base<te::SpiceError>>("exceptions_SpiceBADLATITUDEBOUNDS");

    // BADLATITUDERANGE exception
    class_<te::SpiceBADLATITUDERANGE, base<te::SpiceError>>("exceptions_SpiceBADLATITUDERANGE");

    // BADLATUSRECTUM exception
    class_<te::SpiceBADLATUSRECTUM, base<te::SpiceError>>("exceptions_SpiceBADLATUSRECTUM");

    // BADLEAPSECONDS exception
    class_<te::SpiceBADLEAPSECONDS, base<te::SpiceError>>("exceptions_SpiceBADLEAPSECONDS");

    // BADLIMBLOCUSMIX exception
    class_<te::SpiceBADLIMBLOCUSMIX, base<te::SpiceError>>("exceptions_SpiceBADLIMBLOCUSMIX");

    // BADLINEPERRECCOUNT exception
    class_<te::SpiceBADLINEPERRECCOUNT, base<te::SpiceError>>("exceptions_SpiceBADLINEPERRECCOUNT");

    // BADLISTFILENAME exception
    class_<te::SpiceBADLISTFILENAME, base<te::SpiceError>>("exceptions_SpiceBADLISTFILENAME");

    // BADLONGITUDERANGE exception
    class_<te::SpiceBADLONGITUDERANGE, base<te::SpiceError>>("exceptions_SpiceBADLONGITUDERANGE");

    // BADMATRIX exception
    class_<te::SpiceBADMATRIX, base<te::SpiceError>>("exceptions_SpiceBADMATRIX");

    // BADMEANMOTION exception
    class_<te::SpiceBADMEANMOTION, base<te::SpiceError>>("exceptions_SpiceBADMEANMOTION");

    // BADMECCENTRICITY exception
    class_<te::SpiceBADMECCENTRICITY, base<te::SpiceError>>("exceptions_SpiceBADMECCENTRICITY");

    // BADMETHODSYNTAX exception
    class_<te::SpiceBADMETHODSYNTAX, base<te::SpiceError>>("exceptions_SpiceBADMETHODSYNTAX");

    // BADMIDNIGHTTYPE exception
    class_<te::SpiceBADMIDNIGHTTYPE, base<te::SpiceError>>("exceptions_SpiceBADMIDNIGHTTYPE");

    // BADMSEMIMAJOR exception
    class_<te::SpiceBADMSEMIMAJOR, base<te::SpiceError>>("exceptions_SpiceBADMSEMIMAJOR");

    // BADMSOPQUATERNION exception
    class_<te::SpiceBADMSOPQUATERNION, base<te::SpiceError>>("exceptions_SpiceBADMSOPQUATERNION");

    // BADNOFDIGITS exception
    class_<te::SpiceBADNOFDIGITS, base<te::SpiceError>>("exceptions_SpiceBADNOFDIGITS");

    // BADNOFSTATES exception
    class_<te::SpiceBADNOFSTATES, base<te::SpiceError>>("exceptions_SpiceBADNOFSTATES");

    // BADNUMBEROFPOINTS exception
    class_<te::SpiceBADNUMBEROFPOINTS, base<te::SpiceError>>("exceptions_SpiceBADNUMBEROFPOINTS");

    // BADOBJECTID exception
    class_<te::SpiceBADOBJECTID, base<te::SpiceError>>("exceptions_SpiceBADOBJECTID");

    // BADOBJECTNAME exception
    class_<te::SpiceBADOBJECTNAME, base<te::SpiceError>>("exceptions_SpiceBADOBJECTNAME");

    // BADOFFSETANGLES exception
    class_<te::SpiceBADOFFSETANGLES, base<te::SpiceError>>("exceptions_SpiceBADOFFSETANGLES");

    // BADOFFSETANGUNITS exception
    class_<te::SpiceBADOFFSETANGUNITS, base<te::SpiceError>>("exceptions_SpiceBADOFFSETANGUNITS");

    // BADOFFSETAXESFORMAT exception
    class_<te::SpiceBADOFFSETAXESFORMAT, base<te::SpiceError>>("exceptions_SpiceBADOFFSETAXESFORMAT");

    // BADOFFSETAXISXYZ exception
    class_<te::SpiceBADOFFSETAXISXYZ, base<te::SpiceError>>("exceptions_SpiceBADOFFSETAXISXYZ");

    // BADORBITALPERIOD exception
    class_<te::SpiceBADORBITALPERIOD, base<te::SpiceError>>("exceptions_SpiceBADORBITALPERIOD");

    // BADOUTPUTSPKTYPE exception
    class_<te::SpiceBADOUTPUTSPKTYPE, base<te::SpiceError>>("exceptions_SpiceBADOUTPUTSPKTYPE");

    // BADOUTPUTTYPE exception
    class_<te::SpiceBADOUTPUTTYPE, base<te::SpiceError>>("exceptions_SpiceBADOUTPUTTYPE");

    // BADPARTNUMBER exception
    class_<te::SpiceBADPARTNUMBER, base<te::SpiceError>>("exceptions_SpiceBADPARTNUMBER");

    // BADPCKVALUE exception
    class_<te::SpiceBADPCKVALUE, base<te::SpiceError>>("exceptions_SpiceBADPCKVALUE");

    // BADPECCENTRICITY exception
    class_<te::SpiceBADPECCENTRICITY, base<te::SpiceError>>("exceptions_SpiceBADPECCENTRICITY");

    // BADPERIAPSEVALUE exception
    class_<te::SpiceBADPERIAPSEVALUE, base<te::SpiceError>>("exceptions_SpiceBADPERIAPSEVALUE");

    // BADPICTURE exception
    class_<te::SpiceBADPICTURE, base<te::SpiceError>>("exceptions_SpiceBADPICTURE");

    // BADPLATECOUNT exception
    class_<te::SpiceBADPLATECOUNT, base<te::SpiceError>>("exceptions_SpiceBADPLATECOUNT");

    // BADPODLOCATION exception
    class_<te::SpiceBADPODLOCATION, base<te::SpiceError>>("exceptions_SpiceBADPODLOCATION");

    // BADPRECVALUE exception
    class_<te::SpiceBADPRECVALUE, base<te::SpiceError>>("exceptions_SpiceBADPRECVALUE");

    // BADPRIORITYSPEC exception
    class_<te::SpiceBADPRIORITYSPEC, base<te::SpiceError>>("exceptions_SpiceBADPRIORITYSPEC");

    // BADQUATSIGN exception
    class_<te::SpiceBADQUATSIGN, base<te::SpiceError>>("exceptions_SpiceBADQUATSIGN");

    // BADQUATTHRESHOLD exception
    class_<te::SpiceBADQUATTHRESHOLD, base<te::SpiceError>>("exceptions_SpiceBADQUATTHRESHOLD");

    // BADRADIUS exception
    class_<te::SpiceBADRADIUS, base<te::SpiceError>>("exceptions_SpiceBADRADIUS");

    // BADRADIUSCOUNT exception
    class_<te::SpiceBADRADIUSCOUNT, base<te::SpiceError>>("exceptions_SpiceBADRADIUSCOUNT");

    // BADRATEFRAMEFLAG exception
    class_<te::SpiceBADRATEFRAMEFLAG, base<te::SpiceError>>("exceptions_SpiceBADRATEFRAMEFLAG");

    // BADRATETHRESHOLD exception
    class_<te::SpiceBADRATETHRESHOLD, base<te::SpiceError>>("exceptions_SpiceBADRATETHRESHOLD");

    // BADRECORDCOUNT exception
    class_<te::SpiceBADRECORDCOUNT, base<te::SpiceError>>("exceptions_SpiceBADRECORDCOUNT");

    // BADREFVECTORSPEC exception
    class_<te::SpiceBADREFVECTORSPEC, base<te::SpiceError>>("exceptions_SpiceBADREFVECTORSPEC");

    // BADROTATIONAXISXYZ exception
    class_<te::SpiceBADROTATIONAXISXYZ, base<te::SpiceError>>("exceptions_SpiceBADROTATIONAXISXYZ");

    // BADROTATIONSORDER exception
    class_<te::SpiceBADROTATIONSORDER, base<te::SpiceError>>("exceptions_SpiceBADROTATIONSORDER");

    // BADROTATIONTYPE exception
    class_<te::SpiceBADROTATIONTYPE, base<te::SpiceError>>("exceptions_SpiceBADROTATIONTYPE");

    // BADROTAXESFORMAT exception
    class_<te::SpiceBADROTAXESFORMAT, base<te::SpiceError>>("exceptions_SpiceBADROTAXESFORMAT");

    // BADROWCOUNT exception
    class_<te::SpiceBADROWCOUNT, base<te::SpiceError>>("exceptions_SpiceBADROWCOUNT");

    // BADSCID exception
    class_<te::SpiceBADSCID, base<te::SpiceError>>("exceptions_SpiceBADSCID");

    // BADSEMIAXIS exception
    class_<te::SpiceBADSEMIAXIS, base<te::SpiceError>>("exceptions_SpiceBADSEMIAXIS");

    // BADSEMILATUS exception
    class_<te::SpiceBADSEMILATUS, base<te::SpiceError>>("exceptions_SpiceBADSEMILATUS");

    // BADSHAPE exception
    class_<te::SpiceBADSHAPE, base<te::SpiceError>>("exceptions_SpiceBADSHAPE");

    // BADSOLDAY exception
    class_<te::SpiceBADSOLDAY, base<te::SpiceError>>("exceptions_SpiceBADSOLDAY");

    // BADSOLINDEX exception
    class_<te::SpiceBADSOLINDEX, base<te::SpiceError>>("exceptions_SpiceBADSOLINDEX");

    // BADSOLTIME exception
    class_<te::SpiceBADSOLTIME, base<te::SpiceError>>("exceptions_SpiceBADSOLTIME");

    // BADSOURCERADIUS exception
    class_<te::SpiceBADSOURCERADIUS, base<te::SpiceError>>("exceptions_SpiceBADSOURCERADIUS");

    // BADSPICEQUATERNION exception
    class_<te::SpiceBADSPICEQUATERNION, base<te::SpiceError>>("exceptions_SpiceBADSPICEQUATERNION");

    // BADSTARINDEX exception
    class_<te::SpiceBADSTARINDEX, base<te::SpiceError>>("exceptions_SpiceBADSTARINDEX");

    // BADSTARTTIME exception
    class_<te::SpiceBADSTARTTIME, base<te::SpiceError>>("exceptions_SpiceBADSTARTTIME");

    // BADSTDIONAME exception
    class_<te::SpiceBADSTDIONAME, base<te::SpiceError>>("exceptions_SpiceBADSTDIONAME");

    // BADSTOPTIME exception
    class_<te::SpiceBADSTOPTIME, base<te::SpiceError>>("exceptions_SpiceBADSTOPTIME");

    // BADSUBSTR exception
    class_<te::SpiceBADSUBSTR, base<te::SpiceError>>("exceptions_SpiceBADSUBSTR");

    // BADSUBSTRINGBOUNDS exception
    class_<te::SpiceBADSUBSTRINGBOUNDS, base<te::SpiceError>>("exceptions_SpiceBADSUBSTRINGBOUNDS");

    // BADSURFACEMAP exception
    class_<te::SpiceBADSURFACEMAP, base<te::SpiceError>>("exceptions_SpiceBADSURFACEMAP");

    // BADTABLEFLAG exception
    class_<te::SpiceBADTABLEFLAG, base<te::SpiceError>>("exceptions_SpiceBADTABLEFLAG");

    // BADTERMLOCUSMIX exception
    class_<te::SpiceBADTERMLOCUSMIX, base<te::SpiceError>>("exceptions_SpiceBADTERMLOCUSMIX");

    // BADTIMEBOUNDS exception
    class_<te::SpiceBADTIMEBOUNDS, base<te::SpiceError>>("exceptions_SpiceBADTIMEBOUNDS");

    // BADTIMECASE exception
    class_<te::SpiceBADTIMECASE, base<te::SpiceError>>("exceptions_SpiceBADTIMECASE");

    // BADTIMECOUNT exception
    class_<te::SpiceBADTIMECOUNT, base<te::SpiceError>>("exceptions_SpiceBADTIMECOUNT");

    // BADTIMEFORMAT exception
    class_<te::SpiceBADTIMEFORMAT, base<te::SpiceError>>("exceptions_SpiceBADTIMEFORMAT");

    // BADTIMEITEM exception
    class_<te::SpiceBADTIMEITEM, base<te::SpiceError>>("exceptions_SpiceBADTIMEITEM");

    // BADTIMEOFFSET exception
    class_<te::SpiceBADTIMEOFFSET, base<te::SpiceError>>("exceptions_SpiceBADTIMEOFFSET");

    // BADTIMESPEC exception
    class_<te::SpiceBADTIMESPEC, base<te::SpiceError>>("exceptions_SpiceBADTIMESPEC");

    // BADTIMESTRING exception
    class_<te::SpiceBADTIMESTRING, base<te::SpiceError>>("exceptions_SpiceBADTIMESTRING");

    // BADTIMETYPE exception
    class_<te::SpiceBADTIMETYPE, base<te::SpiceError>>("exceptions_SpiceBADTIMETYPE");

    // BADTIMETYPEFLAG exception
    class_<te::SpiceBADTIMETYPEFLAG, base<te::SpiceError>>("exceptions_SpiceBADTIMETYPEFLAG");

    // BADTLE exception
    class_<te::SpiceBADTLE, base<te::SpiceError>>("exceptions_SpiceBADTLE");

    // BADTLECOVERAGEPAD exception
    class_<te::SpiceBADTLECOVERAGEPAD, base<te::SpiceError>>("exceptions_SpiceBADTLECOVERAGEPAD");

    // BADTLEPADS exception
    class_<te::SpiceBADTLEPADS, base<te::SpiceError>>("exceptions_SpiceBADTLEPADS");

    // BADTOTIME exception
    class_<te::SpiceBADTOTIME, base<te::SpiceError>>("exceptions_SpiceBADTOTIME");

    // BADTOTIMESYSTEM exception
    class_<te::SpiceBADTOTIMESYSTEM, base<te::SpiceError>>("exceptions_SpiceBADTOTIMESYSTEM");

    // BADTOTIMETYPE exception
    class_<te::SpiceBADTOTIMETYPE, base<te::SpiceError>>("exceptions_SpiceBADTOTIMETYPE");

    // BADTYPESHAPECOMBO exception
    class_<te::SpiceBADTYPESHAPECOMBO, base<te::SpiceError>>("exceptions_SpiceBADTYPESHAPECOMBO");

    // BADVARASSIGN exception
    class_<te::SpiceBADVARASSIGN, base<te::SpiceError>>("exceptions_SpiceBADVARASSIGN");

    // BADVARIABLESIZE exception
    class_<te::SpiceBADVARIABLESIZE, base<te::SpiceError>>("exceptions_SpiceBADVARIABLESIZE");

    // BADVARIABLETYPE exception
    class_<te::SpiceBADVARIABLETYPE, base<te::SpiceError>>("exceptions_SpiceBADVARIABLETYPE");

    // BADVARNAME exception
    class_<te::SpiceBADVARNAME, base<te::SpiceError>>("exceptions_SpiceBADVARNAME");

    // BADVECTOR exception
    class_<te::SpiceBADVECTOR, base<te::SpiceError>>("exceptions_SpiceBADVECTOR");

    // BADVERTEXCOUNT exception
    class_<te::SpiceBADVERTEXCOUNT, base<te::SpiceError>>("exceptions_SpiceBADVERTEXCOUNT");

    // BADVERTEXINDEX exception
    class_<te::SpiceBADVERTEXINDEX, base<te::SpiceError>>("exceptions_SpiceBADVERTEXINDEX");

    // BADWINDOWSIZE exception
    class_<te::SpiceBADWINDOWSIZE, base<te::SpiceError>>("exceptions_SpiceBADWINDOWSIZE");

    // BARRAYTOOSMALL exception
    class_<te::SpiceBARRAYTOOSMALL, base<te::SpiceError>>("exceptions_SpiceBARRAYTOOSMALL");

    // BARYCENTEREPHEM exception
    class_<te::SpiceBARYCENTEREPHEM, base<te::SpiceError>>("exceptions_SpiceBARYCENTEREPHEM");

    // BARYCENTERIDCODE exception
    class_<te::SpiceBARYCENTERIDCODE, base<te::SpiceError>>("exceptions_SpiceBARYCENTERIDCODE");

    // BEFOREBEGSTR exception
    class_<te::SpiceBEFOREBEGSTR, base<te::SpiceError>>("exceptions_SpiceBEFOREBEGSTR");

    // BLANKCOMMANDLINE exception
    class_<te::SpiceBLANKCOMMANDLINE, base<te::SpiceError>>("exceptions_SpiceBLANKCOMMANDLINE");

    // BLANKFILENAME exception
    class_<te::SpiceBLANKFILENAME, base<te::SpiceError>>("exceptions_SpiceBLANKFILENAME");

    // BLANKFILETYPE exception
    class_<te::SpiceBLANKFILETYPE, base<te::SpiceError>>("exceptions_SpiceBLANKFILETYPE");

    // BLANKINPUTFILENAME exception
    class_<te::SpiceBLANKINPUTFILENAME, base<te::SpiceError>>("exceptions_SpiceBLANKINPUTFILENAME");

    // BLANKINPUTTIME exception
    class_<te::SpiceBLANKINPUTTIME, base<te::SpiceError>>("exceptions_SpiceBLANKINPUTTIME");

    // BLANKNAMEASSIGNED exception
    class_<te::SpiceBLANKNAMEASSIGNED, base<te::SpiceError>>("exceptions_SpiceBLANKNAMEASSIGNED");

    // BLANKOUTPTFILENAME exception
    class_<te::SpiceBLANKOUTPTFILENAME, base<te::SpiceError>>("exceptions_SpiceBLANKOUTPTFILENAME");

    // BLANKSCLKSTRING exception
    class_<te::SpiceBLANKSCLKSTRING, base<te::SpiceError>>("exceptions_SpiceBLANKSCLKSTRING");

    // BLANKTIMEFORMAT exception
    class_<te::SpiceBLANKTIMEFORMAT, base<te::SpiceError>>("exceptions_SpiceBLANKTIMEFORMAT");

    // BLOCKSNOTEVEN exception
    class_<te::SpiceBLOCKSNOTEVEN, base<te::SpiceError>>("exceptions_SpiceBLOCKSNOTEVEN");

    // BODIESNOTDISTINCT exception
    class_<te::SpiceBODIESNOTDISTINCT, base<te::SpiceError>>("exceptions_SpiceBODIESNOTDISTINCT");

    // BODYANDCENTERSAME exception
    class_<te::SpiceBODYANDCENTERSAME, base<te::SpiceError>>("exceptions_SpiceBODYANDCENTERSAME");

    // BOGUSENTRY exception
    class_<te::SpiceBOGUSENTRY, base<te::SpiceError>>("exceptions_SpiceBOGUSENTRY");

    // BORESIGHTMISSING exception
    class_<te::SpiceBORESIGHTMISSING, base<te::SpiceError>>("exceptions_SpiceBORESIGHTMISSING");

    // BOUNDARYMISSING exception
    class_<te::SpiceBOUNDARYMISSING, base<te::SpiceError>>("exceptions_SpiceBOUNDARYMISSING");

    // BOUNDARYTOOBIG exception
    class_<te::SpiceBOUNDARYTOOBIG, base<te::SpiceError>>("exceptions_SpiceBOUNDARYTOOBIG");

    // BOUNDSDISAGREE exception
    class_<te::SpiceBOUNDSDISAGREE, base<te::SpiceError>>("exceptions_SpiceBOUNDSDISAGREE");

    // BOUNDSOUTOFORDER exception
    class_<te::SpiceBOUNDSOUTOFORDER, base<te::SpiceError>>("exceptions_SpiceBOUNDSOUTOFORDER");

    // BUFFEROVERFLOW exception
    class_<te::SpiceBUFFEROVERFLOW, base<te::SpiceError>>("exceptions_SpiceBUFFEROVERFLOW");

    // BUFFERSIZESMISMATCH exception
    class_<te::SpiceBUFFERSIZESMISMATCH, base<te::SpiceError>>("exceptions_SpiceBUFFERSIZESMISMATCH");

    // BUFFERTOOSMALL exception
    class_<te::SpiceBUFFERTOOSMALL, base<te::SpiceError>>("exceptions_SpiceBUFFERTOOSMALL");

    // BUG exception
    class_<te::SpiceBUG, base<te::SpiceError>>("exceptions_SpiceBUG");

    // BUGWRITEFAILED exception
    class_<te::SpiceBUGWRITEFAILED, base<te::SpiceError>>("exceptions_SpiceBUGWRITEFAILED");

    // CALLCKBSSFIRST exception
    class_<te::SpiceCALLCKBSSFIRST, base<te::SpiceError>>("exceptions_SpiceCALLCKBSSFIRST");

    // CALLEDOUTOFORDER exception
    class_<te::SpiceCALLEDOUTOFORDER, base<te::SpiceError>>("exceptions_SpiceCALLEDOUTOFORDER");

    // CALLZZDSKBSSFIRST exception
    class_<te::SpiceCALLZZDSKBSSFIRST, base<te::SpiceError>>("exceptions_SpiceCALLZZDSKBSSFIRST");

    // CANNOTFINDGRP exception
    class_<te::SpiceCANNOTFINDGRP, base<te::SpiceError>>("exceptions_SpiceCANNOTFINDGRP");

    // CANNOTGETPACKET exception
    class_<te::SpiceCANNOTGETPACKET, base<te::SpiceError>>("exceptions_SpiceCANNOTGETPACKET");

    // CANNOTMAKEFILE exception
    class_<te::SpiceCANNOTMAKEFILE, base<te::SpiceError>>("exceptions_SpiceCANNOTMAKEFILE");

    // CANNOTPICKFRAME exception
    class_<te::SpiceCANNOTPICKFRAME, base<te::SpiceError>>("exceptions_SpiceCANNOTPICKFRAME");

    // CANTFINDFRAME exception
    class_<te::SpiceCANTFINDFRAME, base<te::SpiceError>>("exceptions_SpiceCANTFINDFRAME");

    // CANTGETROTATIONTYPE exception
    class_<te::SpiceCANTGETROTATIONTYPE, base<te::SpiceError>>("exceptions_SpiceCANTGETROTATIONTYPE");

    // CANTUSEPERIAPEPOCH exception
    class_<te::SpiceCANTUSEPERIAPEPOCH, base<te::SpiceError>>("exceptions_SpiceCANTUSEPERIAPEPOCH");

    // CBNOSUCHSTR exception
    class_<te::SpiceCBNOSUCHSTR, base<te::SpiceError>>("exceptions_SpiceCBNOSUCHSTR");

    // CELLARRAYTOOSMALL exception
    class_<te::SpiceCELLARRAYTOOSMALL, base<te::SpiceError>>("exceptions_SpiceCELLARRAYTOOSMALL");

    // CELLTOOSMALL exception
    class_<te::SpiceCELLTOOSMALL, base<te::SpiceError>>("exceptions_SpiceCELLTOOSMALL");

    // CKBOGUSENTRY exception
    class_<te::SpiceCKBOGUSENTRY, base<te::SpiceError>>("exceptions_SpiceCKBOGUSENTRY");

    // CKDOESNTEXIST exception
    class_<te::SpiceCKDOESNTEXIST, base<te::SpiceError>>("exceptions_SpiceCKDOESNTEXIST");

    // CKFILE exception
    class_<te::SpiceCKFILE, base<te::SpiceError>>("exceptions_SpiceCKFILE");

    // CKNONEXISTREC exception
    class_<te::SpiceCKNONEXISTREC, base<te::SpiceError>>("exceptions_SpiceCKNONEXISTREC");

    // CKTOOMANYFILES exception
    class_<te::SpiceCKTOOMANYFILES, base<te::SpiceError>>("exceptions_SpiceCKTOOMANYFILES");

    // CKUNKNOWNDATATYPE exception
    class_<te::SpiceCKUNKNOWNDATATYPE, base<te::SpiceError>>("exceptions_SpiceCKUNKNOWNDATATYPE");

    // CKWRONGDATATYPE exception
    class_<te::SpiceCKWRONGDATATYPE, base<te::SpiceError>>("exceptions_SpiceCKWRONGDATATYPE");

    // CMDERROR exception
    class_<te::SpiceCMDERROR, base<te::SpiceError>>("exceptions_SpiceCMDERROR");

    // CMDPARSEERROR exception
    class_<te::SpiceCMDPARSEERROR, base<te::SpiceError>>("exceptions_SpiceCMDPARSEERROR");

    // COARSEGRIDOVERFLOW exception
    class_<te::SpiceCOARSEGRIDOVERFLOW, base<te::SpiceError>>("exceptions_SpiceCOARSEGRIDOVERFLOW");

    // COLDESCTABLEFULL exception
    class_<te::SpiceCOLDESCTABLEFULL, base<te::SpiceError>>("exceptions_SpiceCOLDESCTABLEFULL");

    // COLUMNTOOSMALL exception
    class_<te::SpiceCOLUMNTOOSMALL, base<te::SpiceError>>("exceptions_SpiceCOLUMNTOOSMALL");

    // COMMANDTOOLONG exception
    class_<te::SpiceCOMMANDTOOLONG, base<te::SpiceError>>("exceptions_SpiceCOMMANDTOOLONG");

    // COMMENTTOOLONG exception
    class_<te::SpiceCOMMENTTOOLONG, base<te::SpiceError>>("exceptions_SpiceCOMMENTTOOLONG");

    // COMMFILENOTEXIST exception
    class_<te::SpiceCOMMFILENOTEXIST, base<te::SpiceError>>("exceptions_SpiceCOMMFILENOTEXIST");

    // COMPETINGEPOCHSPEC exception
    class_<te::SpiceCOMPETINGEPOCHSPEC, base<te::SpiceError>>("exceptions_SpiceCOMPETINGEPOCHSPEC");

    // COMPETINGFRAMESPEC exception
    class_<te::SpiceCOMPETINGFRAMESPEC, base<te::SpiceError>>("exceptions_SpiceCOMPETINGFRAMESPEC");

    // COORDSYSNOTREC exception
    class_<te::SpiceCOORDSYSNOTREC, base<te::SpiceError>>("exceptions_SpiceCOORDSYSNOTREC");

    // COUNTMISMATCH exception
    class_<te::SpiceCOUNTMISMATCH, base<te::SpiceError>>("exceptions_SpiceCOUNTMISMATCH");

    // COUNTTOOLARGE exception
    class_<te::SpiceCOUNTTOOLARGE, base<te::SpiceError>>("exceptions_SpiceCOUNTTOOLARGE");

    // COVERAGEGAP exception
    class_<te::SpiceCOVERAGEGAP, base<te::SpiceError>>("exceptions_SpiceCOVERAGEGAP");

    // CROSSANGLEMISSING exception
    class_<te::SpiceCROSSANGLEMISSING, base<te::SpiceError>>("exceptions_SpiceCROSSANGLEMISSING");

    // DAFBADCRECLEN exception
    class_<te::SpiceDAFBADCRECLEN, base<te::SpiceError>>("exceptions_SpiceDAFBADCRECLEN");

    // DAFBEGGTEND exception
    class_<te::SpiceDAFBEGGTEND, base<te::SpiceError>>("exceptions_SpiceDAFBEGGTEND");

    // DAFCRNOTFOUND exception
    class_<te::SpiceDAFCRNOTFOUND, base<te::SpiceError>>("exceptions_SpiceDAFCRNOTFOUND");

    // DAFDPWRITEFAIL exception
    class_<te::SpiceDAFDPWRITEFAIL, base<te::SpiceError>>("exceptions_SpiceDAFDPWRITEFAIL");

    // DAFFRNOTFOUND exception
    class_<te::SpiceDAFFRNOTFOUND, base<te::SpiceError>>("exceptions_SpiceDAFFRNOTFOUND");

    // DAFFTFULL exception
    class_<te::SpiceDAFFTFULL, base<te::SpiceError>>("exceptions_SpiceDAFFTFULL");

    // DAFILLEGWRITE exception
    class_<te::SpiceDAFILLEGWRITE, base<te::SpiceError>>("exceptions_SpiceDAFILLEGWRITE");

    // DAFINVALIDACCESS exception
    class_<te::SpiceDAFINVALIDACCESS, base<te::SpiceError>>("exceptions_SpiceDAFINVALIDACCESS");

    // DAFINVALIDPARAMS exception
    class_<te::SpiceDAFINVALIDPARAMS, base<te::SpiceError>>("exceptions_SpiceDAFINVALIDPARAMS");

    // DAFNEGADDR exception
    class_<te::SpiceDAFNEGADDR, base<te::SpiceError>>("exceptions_SpiceDAFNEGADDR");

    // DAFNEWCONFLICT exception
    class_<te::SpiceDAFNEWCONFLICT, base<te::SpiceError>>("exceptions_SpiceDAFNEWCONFLICT");

    // DAFNOIFNMATCH exception
    class_<te::SpiceDAFNOIFNMATCH, base<te::SpiceError>>("exceptions_SpiceDAFNOIFNMATCH");

    // DAFNONAMEMATCH exception
    class_<te::SpiceDAFNONAMEMATCH, base<te::SpiceError>>("exceptions_SpiceDAFNONAMEMATCH");

    // DAFNORESV exception
    class_<te::SpiceDAFNORESV, base<te::SpiceError>>("exceptions_SpiceDAFNORESV");

    // DAFNOSEARCH exception
    class_<te::SpiceDAFNOSEARCH, base<te::SpiceError>>("exceptions_SpiceDAFNOSEARCH");

    // DAFNOSUCHADDR exception
    class_<te::SpiceDAFNOSUCHADDR, base<te::SpiceError>>("exceptions_SpiceDAFNOSUCHADDR");

    // DAFNOSUCHFILE exception
    class_<te::SpiceDAFNOSUCHFILE, base<te::SpiceError>>("exceptions_SpiceDAFNOSUCHFILE");

    // DAFNOSUCHHANDLE exception
    class_<te::SpiceDAFNOSUCHHANDLE, base<te::SpiceError>>("exceptions_SpiceDAFNOSUCHHANDLE");

    // DAFNOSUCHUNIT exception
    class_<te::SpiceDAFNOSUCHUNIT, base<te::SpiceError>>("exceptions_SpiceDAFNOSUCHUNIT");

    // DAFNOWRITE exception
    class_<te::SpiceDAFNOWRITE, base<te::SpiceError>>("exceptions_SpiceDAFNOWRITE");

    // DAFOVERFLOW exception
    class_<te::SpiceDAFOVERFLOW, base<te::SpiceError>>("exceptions_SpiceDAFOVERFLOW");

    // DAFREADFAIL exception
    class_<te::SpiceDAFREADFAIL, base<te::SpiceError>>("exceptions_SpiceDAFREADFAIL");

    // DAFWRITEFAIL exception
    class_<te::SpiceDAFWRITEFAIL, base<te::SpiceError>>("exceptions_SpiceDAFWRITEFAIL");

    // DASFILEREADFAILED exception
    class_<te::SpiceDASFILEREADFAILED, base<te::SpiceError>>("exceptions_SpiceDASFILEREADFAILED");

    // DASFILEWRITEFAILED exception
    class_<te::SpiceDASFILEWRITEFAILED, base<te::SpiceError>>("exceptions_SpiceDASFILEWRITEFAILED");

    // DASFTFULL exception
    class_<te::SpiceDASFTFULL, base<te::SpiceError>>("exceptions_SpiceDASFTFULL");

    // DASINVALIDACCESS exception
    class_<te::SpiceDASINVALIDACCESS, base<te::SpiceError>>("exceptions_SpiceDASINVALIDACCESS");

    // DASINVALIDCOUNT exception
    class_<te::SpiceDASINVALIDCOUNT, base<te::SpiceError>>("exceptions_SpiceDASINVALIDCOUNT");

    // DASINVALIDTYPE exception
    class_<te::SpiceDASINVALIDTYPE, base<te::SpiceError>>("exceptions_SpiceDASINVALIDTYPE");

    // DASNOSUCHADDRESS exception
    class_<te::SpiceDASNOSUCHADDRESS, base<te::SpiceError>>("exceptions_SpiceDASNOSUCHADDRESS");

    // DASNOSUCHFILE exception
    class_<te::SpiceDASNOSUCHFILE, base<te::SpiceError>>("exceptions_SpiceDASNOSUCHFILE");

    // DASNOSUCHHANDLE exception
    class_<te::SpiceDASNOSUCHHANDLE, base<te::SpiceError>>("exceptions_SpiceDASNOSUCHHANDLE");

    // DASNOSUCHUNIT exception
    class_<te::SpiceDASNOSUCHUNIT, base<te::SpiceError>>("exceptions_SpiceDASNOSUCHUNIT");

    // DASNOTEMPTY exception
    class_<te::SpiceDASNOTEMPTY, base<te::SpiceError>>("exceptions_SpiceDASNOTEMPTY");

    // DASREADFAIL exception
    class_<te::SpiceDASREADFAIL, base<te::SpiceError>>("exceptions_SpiceDASREADFAIL");

    // DASWRITEFAIL exception
    class_<te::SpiceDASWRITEFAIL, base<te::SpiceError>>("exceptions_SpiceDASWRITEFAIL");

    // DATAITEMLIMITEXCEEDED exception
    class_<te::SpiceDATAITEMLIMITEXCEEDED, base<te::SpiceError>>("exceptions_SpiceDATAITEMLIMITEXCEEDED");

    // DATAREADFAILED exception
    class_<te::SpiceDATAREADFAILED, base<te::SpiceError>>("exceptions_SpiceDATAREADFAILED");

    // DATAWIDTHERROR exception
    class_<te::SpiceDATAWIDTHERROR, base<te::SpiceError>>("exceptions_SpiceDATAWIDTHERROR");

    // DATEEXPECTED exception
    class_<te::SpiceDATEEXPECTED, base<te::SpiceError>>("exceptions_SpiceDATEEXPECTED");

    // DECODINGERROR exception
    class_<te::SpiceDECODINGERROR, base<te::SpiceError>>("exceptions_SpiceDECODINGERROR");

    // DEGENERATECASE exception
    class_<te::SpiceDEGENERATECASE, base<te::SpiceError>>("exceptions_SpiceDEGENERATECASE");

    // DEGENERATEINTERVAL exception
    class_<te::SpiceDEGENERATEINTERVAL, base<te::SpiceError>>("exceptions_SpiceDEGENERATEINTERVAL");

    // DEGENERATESURFACE exception
    class_<te::SpiceDEGENERATESURFACE, base<te::SpiceError>>("exceptions_SpiceDEGENERATESURFACE");

    // DEGREEOUTOFRANGE exception
    class_<te::SpiceDEGREEOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceDEGREEOUTOFRANGE");

    // DEPENDENTVECTORS exception
    class_<te::SpiceDEPENDENTVECTORS, base<te::SpiceError>>("exceptions_SpiceDEPENDENTVECTORS");

    // DEVICENAMETOOLONG exception
    class_<te::SpiceDEVICENAMETOOLONG, base<te::SpiceError>>("exceptions_SpiceDEVICENAMETOOLONG");

    // DIFFLINETOOLARGE exception
    class_<te::SpiceDIFFLINETOOLARGE, base<te::SpiceError>>("exceptions_SpiceDIFFLINETOOLARGE");

    // DIFFLINETOOSMALL exception
    class_<te::SpiceDIFFLINETOOSMALL, base<te::SpiceError>>("exceptions_SpiceDIFFLINETOOSMALL");

    // DIMENSIONTOOSMALL exception
    class_<te::SpiceDIMENSIONTOOSMALL, base<te::SpiceError>>("exceptions_SpiceDIMENSIONTOOSMALL");

    // DISARRAY exception
    class_<te::SpiceDISARRAY, base<te::SpiceError>>("exceptions_SpiceDISARRAY");

    // DISORDER exception
    class_<te::SpiceDISORDER, base<te::SpiceError>>("exceptions_SpiceDISORDER");

    // DIVIDEBYZERO exception
    class_<te::SpiceDIVIDEBYZERO, base<te::SpiceError>>("exceptions_SpiceDIVIDEBYZERO");

    // DSKBOGUSENTRY exception
    class_<te::SpiceDSKBOGUSENTRY, base<te::SpiceError>>("exceptions_SpiceDSKBOGUSENTRY");

    // DSKDATANOTFOUND exception
    class_<te::SpiceDSKDATANOTFOUND, base<te::SpiceError>>("exceptions_SpiceDSKDATANOTFOUND");

    // DSKTOOMANYFILES exception
    class_<te::SpiceDSKTOOMANYFILES, base<te::SpiceError>>("exceptions_SpiceDSKTOOMANYFILES");

    // DTOUTOFRANGE exception
    class_<te::SpiceDTOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceDTOUTOFRANGE");

    // DUBIOUSMETHOD exception
    class_<te::SpiceDUBIOUSMETHOD, base<te::SpiceError>>("exceptions_SpiceDUBIOUSMETHOD");

    // DUPLICATETIMES exception
    class_<te::SpiceDUPLICATETIMES, base<te::SpiceError>>("exceptions_SpiceDUPLICATETIMES");

    // ECCOUTOFBOUNDS exception
    class_<te::SpiceECCOUTOFBOUNDS, base<te::SpiceError>>("exceptions_SpiceECCOUTOFBOUNDS");

    // ECCOUTOFRANGE exception
    class_<te::SpiceECCOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceECCOUTOFRANGE");

    // EKCOLATTRTABLEFULL exception
    class_<te::SpiceEKCOLATTRTABLEFULL, base<te::SpiceError>>("exceptions_SpiceEKCOLATTRTABLEFULL");

    // EKCOLNUMMISMATCH exception
    class_<te::SpiceEKCOLNUMMISMATCH, base<te::SpiceError>>("exceptions_SpiceEKCOLNUMMISMATCH");

    // EKFILE exception
    class_<te::SpiceEKFILE, base<te::SpiceError>>("exceptions_SpiceEKFILE");

    // EKFILETABLEFULL exception
    class_<te::SpiceEKFILETABLEFULL, base<te::SpiceError>>("exceptions_SpiceEKFILETABLEFULL");

    // EKMISSINGCOLUMN exception
    class_<te::SpiceEKMISSINGCOLUMN, base<te::SpiceError>>("exceptions_SpiceEKMISSINGCOLUMN");

    // EKNOSEGMENTS exception
    class_<te::SpiceEKNOSEGMENTS, base<te::SpiceError>>("exceptions_SpiceEKNOSEGMENTS");

    // EKSEGTABLEFULL exception
    class_<te::SpiceEKSEGTABLEFULL, base<te::SpiceError>>("exceptions_SpiceEKSEGTABLEFULL");

    // EKTABLELISTFULL exception
    class_<te::SpiceEKTABLELISTFULL, base<te::SpiceError>>("exceptions_SpiceEKTABLELISTFULL");

    // ELEMENTSTOOSHORT exception
    class_<te::SpiceELEMENTSTOOSHORT, base<te::SpiceError>>("exceptions_SpiceELEMENTSTOOSHORT");

    // EMPTYINPUTFILE exception
    class_<te::SpiceEMPTYINPUTFILE, base<te::SpiceError>>("exceptions_SpiceEMPTYINPUTFILE");

    // EMPTYSEGMENT exception
    class_<te::SpiceEMPTYSEGMENT, base<te::SpiceError>>("exceptions_SpiceEMPTYSEGMENT");

    // ENDOFFILE exception
    class_<te::SpiceENDOFFILE, base<te::SpiceError>>("exceptions_SpiceENDOFFILE");

    // ENDPOINTSMATCH exception
    class_<te::SpiceENDPOINTSMATCH, base<te::SpiceError>>("exceptions_SpiceENDPOINTSMATCH");

    // ERROREXIT exception
    class_<te::SpiceERROREXIT, base<te::SpiceError>>("exceptions_SpiceERROREXIT");

    // EVECOUTOFRANGE exception
    class_<te::SpiceEVECOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceEVECOUTOFRANGE");

    // EVENHERMITDEGREE exception
    class_<te::SpiceEVENHERMITDEGREE, base<te::SpiceError>>("exceptions_SpiceEVENHERMITDEGREE");

    // EVILBOGUSENTRY exception
    class_<te::SpiceEVILBOGUSENTRY, base<te::SpiceError>>("exceptions_SpiceEVILBOGUSENTRY");

    // EXTERNALOPEN exception
    class_<te::SpiceEXTERNALOPEN, base<te::SpiceError>>("exceptions_SpiceEXTERNALOPEN");

    // FACENOTFOUND exception
    class_<te::SpiceFACENOTFOUND, base<te::SpiceError>>("exceptions_SpiceFACENOTFOUND");

    // FAKESCLKEXISTS exception
    class_<te::SpiceFAKESCLKEXISTS, base<te::SpiceError>>("exceptions_SpiceFAKESCLKEXISTS");

    // FILARCHMISMATCH exception
    class_<te::SpiceFILARCHMISMATCH, base<te::SpiceError>>("exceptions_SpiceFILARCHMISMATCH");

    // FILARCMISMATCH exception
    class_<te::SpiceFILARCMISMATCH, base<te::SpiceError>>("exceptions_SpiceFILARCMISMATCH");

    // FILEALREADYEXISTS exception
    class_<te::SpiceFILEALREADYEXISTS, base<te::SpiceError>>("exceptions_SpiceFILEALREADYEXISTS");

    // FILECURRENTLYOPEN exception
    class_<te::SpiceFILECURRENTLYOPEN, base<te::SpiceError>>("exceptions_SpiceFILECURRENTLYOPEN");

    // FILEDELETEFAILED exception
    class_<te::SpiceFILEDELETEFAILED, base<te::SpiceError>>("exceptions_SpiceFILEDELETEFAILED");

    // FILEDOESNOTEXIST exception
    class_<te::SpiceFILEDOESNOTEXIST, base<te::SpiceError>>("exceptions_SpiceFILEDOESNOTEXIST");

    // FILEEXISTS exception
    class_<te::SpiceFILEEXISTS, base<te::SpiceError>>("exceptions_SpiceFILEEXISTS");

    // FILEISNOTSPK exception
    class_<te::SpiceFILEISNOTSPK, base<te::SpiceError>>("exceptions_SpiceFILEISNOTSPK");

    // FILENAMETOOLONG exception
    class_<te::SpiceFILENAMETOOLONG, base<te::SpiceError>>("exceptions_SpiceFILENAMETOOLONG");

    // FILENOTCONNECTED exception
    class_<te::SpiceFILENOTCONNECTED, base<te::SpiceError>>("exceptions_SpiceFILENOTCONNECTED");

    // FILENOTFOUND exception
    class_<te::SpiceFILENOTFOUND, base<te::SpiceError>>("exceptions_SpiceFILENOTFOUND");

    // FILENOTOPEN exception
    class_<te::SpiceFILENOTOPEN, base<te::SpiceError>>("exceptions_SpiceFILENOTOPEN");

    // FILEOPENCONFLICT exception
    class_<te::SpiceFILEOPENCONFLICT, base<te::SpiceError>>("exceptions_SpiceFILEOPENCONFLICT");

    // FILEOPENERROR exception
    class_<te::SpiceFILEOPENERROR, base<te::SpiceError>>("exceptions_SpiceFILEOPENERROR");

    // FILEOPENFAIL exception
    class_<te::SpiceFILEOPENFAIL, base<te::SpiceError>>("exceptions_SpiceFILEOPENFAIL");

    // FILEOPENFAILED exception
    class_<te::SpiceFILEOPENFAILED, base<te::SpiceError>>("exceptions_SpiceFILEOPENFAILED");

    // FILEREADERROR exception
    class_<te::SpiceFILEREADERROR, base<te::SpiceError>>("exceptions_SpiceFILEREADERROR");

    // FILEREADFAILED exception
    class_<te::SpiceFILEREADFAILED, base<te::SpiceError>>("exceptions_SpiceFILEREADFAILED");

    // FILETABLEFULL exception
    class_<te::SpiceFILETABLEFULL, base<te::SpiceError>>("exceptions_SpiceFILETABLEFULL");

    // FILETRUNCATED exception
    class_<te::SpiceFILETRUNCATED, base<te::SpiceError>>("exceptions_SpiceFILETRUNCATED");

    // FILEWRITEFAILED exception
    class_<te::SpiceFILEWRITEFAILED, base<te::SpiceError>>("exceptions_SpiceFILEWRITEFAILED");

    // FIRSTRECORDMISMATCH exception
    class_<te::SpiceFIRSTRECORDMISMATCH, base<te::SpiceError>>("exceptions_SpiceFIRSTRECORDMISMATCH");

    // FKDOESNTEXIST exception
    class_<te::SpiceFKDOESNTEXIST, base<te::SpiceError>>("exceptions_SpiceFKDOESNTEXIST");

    // FMTITEMLIMITEXCEEDED exception
    class_<te::SpiceFMTITEMLIMITEXCEEDED, base<te::SpiceError>>("exceptions_SpiceFMTITEMLIMITEXCEEDED");

    // FORMATDATAMISMATCH exception
    class_<te::SpiceFORMATDATAMISMATCH, base<te::SpiceError>>("exceptions_SpiceFORMATDATAMISMATCH");

    // FORMATDOESNTAPPLY exception
    class_<te::SpiceFORMATDOESNTAPPLY, base<te::SpiceError>>("exceptions_SpiceFORMATDOESNTAPPLY");

    // FORMATERROR exception
    class_<te::SpiceFORMATERROR, base<te::SpiceError>>("exceptions_SpiceFORMATERROR");

    // FORMATNOTAPPLICABLE exception
    class_<te::SpiceFORMATNOTAPPLICABLE, base<te::SpiceError>>("exceptions_SpiceFORMATNOTAPPLICABLE");

    // FORMATSTRINGTOOLONG exception
    class_<te::SpiceFORMATSTRINGTOOLONG, base<te::SpiceError>>("exceptions_SpiceFORMATSTRINGTOOLONG");

    // FOVTOOWIDE exception
    class_<te::SpiceFOVTOOWIDE, base<te::SpiceError>>("exceptions_SpiceFOVTOOWIDE");

    // FRAMEDATANOTFOUND exception
    class_<te::SpiceFRAMEDATANOTFOUND, base<te::SpiceError>>("exceptions_SpiceFRAMEDATANOTFOUND");

    // FRAMEDEFERROR exception
    class_<te::SpiceFRAMEDEFERROR, base<te::SpiceError>>("exceptions_SpiceFRAMEDEFERROR");

    // FRAMEIDNOTFOUND exception
    class_<te::SpiceFRAMEIDNOTFOUND, base<te::SpiceError>>("exceptions_SpiceFRAMEIDNOTFOUND");

    // FRAMEINFONOTFOUND exception
    class_<te::SpiceFRAMEINFONOTFOUND, base<te::SpiceError>>("exceptions_SpiceFRAMEINFONOTFOUND");

    // FRAMEMISSING exception
    class_<te::SpiceFRAMEMISSING, base<te::SpiceError>>("exceptions_SpiceFRAMEMISSING");

    // FRAMENAMENOTFOUND exception
    class_<te::SpiceFRAMENAMENOTFOUND, base<te::SpiceError>>("exceptions_SpiceFRAMENAMENOTFOUND");

    // FRAMENOTFOUND exception
    class_<te::SpiceFRAMENOTFOUND, base<te::SpiceError>>("exceptions_SpiceFRAMENOTFOUND");

    // FRAMENOTRECOGNIZED exception
    class_<te::SpiceFRAMENOTRECOGNIZED, base<te::SpiceError>>("exceptions_SpiceFRAMENOTRECOGNIZED");

    // FTFULL exception
    class_<te::SpiceFTFULL, base<te::SpiceError>>("exceptions_SpiceFTFULL");

    // FTPXFERERROR exception
    class_<te::SpiceFTPXFERERROR, base<te::SpiceError>>("exceptions_SpiceFTPXFERERROR");

    // GRIDTOOLARGE exception
    class_<te::SpiceGRIDTOOLARGE, base<te::SpiceError>>("exceptions_SpiceGRIDTOOLARGE");

    // HANDLENOTFOUND exception
    class_<te::SpiceHANDLENOTFOUND, base<te::SpiceError>>("exceptions_SpiceHANDLENOTFOUND");

    // HASHISFULL exception
    class_<te::SpiceHASHISFULL, base<te::SpiceError>>("exceptions_SpiceHASHISFULL");

    // HLULOCKFAILED exception
    class_<te::SpiceHLULOCKFAILED, base<te::SpiceError>>("exceptions_SpiceHLULOCKFAILED");

    // IDCODENOTFOUND exception
    class_<te::SpiceIDCODENOTFOUND, base<te::SpiceError>>("exceptions_SpiceIDCODENOTFOUND");

    // IDSTRINGTOOLONG exception
    class_<te::SpiceIDSTRINGTOOLONG, base<te::SpiceError>>("exceptions_SpiceIDSTRINGTOOLONG");

    // ILLEGALCHARACTER exception
    class_<te::SpiceILLEGALCHARACTER, base<te::SpiceError>>("exceptions_SpiceILLEGALCHARACTER");

    // ILLEGALOPTIONNAME exception
    class_<te::SpiceILLEGALOPTIONNAME, base<te::SpiceError>>("exceptions_SpiceILLEGALOPTIONNAME");

    // ILLEGSHIFTDIR exception
    class_<te::SpiceILLEGSHIFTDIR, base<te::SpiceError>>("exceptions_SpiceILLEGSHIFTDIR");

    // ILLEGTEMPL exception
    class_<te::SpiceILLEGTEMPL, base<te::SpiceError>>("exceptions_SpiceILLEGTEMPL");

    // IMMUTABLEVALUE exception
    class_<te::SpiceIMMUTABLEVALUE, base<te::SpiceError>>("exceptions_SpiceIMMUTABLEVALUE");

    // IMPROPERFILE exception
    class_<te::SpiceIMPROPERFILE, base<te::SpiceError>>("exceptions_SpiceIMPROPERFILE");

    // IMPROPEROPEN exception
    class_<te::SpiceIMPROPEROPEN, base<te::SpiceError>>("exceptions_SpiceIMPROPEROPEN");

    // INACTIVEOBJECT exception
    class_<te::SpiceINACTIVEOBJECT, base<te::SpiceError>>("exceptions_SpiceINACTIVEOBJECT");

    // INCOMPATIBLEEOL exception
    class_<te::SpiceINCOMPATIBLEEOL, base<te::SpiceError>>("exceptions_SpiceINCOMPATIBLEEOL");

    // INCOMPATIBLENUMREF exception
    class_<te::SpiceINCOMPATIBLENUMREF, base<te::SpiceError>>("exceptions_SpiceINCOMPATIBLENUMREF");

    // INCOMPATIBLESCALE exception
    class_<te::SpiceINCOMPATIBLESCALE, base<te::SpiceError>>("exceptions_SpiceINCOMPATIBLESCALE");

    // INCOMPATIBLEUNITS exception
    class_<te::SpiceINCOMPATIBLEUNITS, base<te::SpiceError>>("exceptions_SpiceINCOMPATIBLEUNITS");

    // INCOMPLETEFRAME exception
    class_<te::SpiceINCOMPLETEFRAME, base<te::SpiceError>>("exceptions_SpiceINCOMPLETEFRAME");

    // INCONSISTCENTERID exception
    class_<te::SpiceINCONSISTCENTERID, base<te::SpiceError>>("exceptions_SpiceINCONSISTCENTERID");

    // INCONSISTENTTIMES exception
    class_<te::SpiceINCONSISTENTTIMES, base<te::SpiceError>>("exceptions_SpiceINCONSISTENTTIMES");

    // INCONSISTFRAME exception
    class_<te::SpiceINCONSISTFRAME, base<te::SpiceError>>("exceptions_SpiceINCONSISTFRAME");

    // INCONSISTSTARTTIME exception
    class_<te::SpiceINCONSISTSTARTTIME, base<te::SpiceError>>("exceptions_SpiceINCONSISTSTARTTIME");

    // INCONSISTSTOPTIME exception
    class_<te::SpiceINCONSISTSTOPTIME, base<te::SpiceError>>("exceptions_SpiceINCONSISTSTOPTIME");

    // INCORRECTUSAGE exception
    class_<te::SpiceINCORRECTUSAGE, base<te::SpiceError>>("exceptions_SpiceINCORRECTUSAGE");

    // INDEFINITELOCALSECOND exception
    class_<te::SpiceINDEFINITELOCALSECOND, base<te::SpiceError>>("exceptions_SpiceINDEFINITELOCALSECOND");

    // INDEXOUTOFRANGE exception
    class_<te::SpiceINDEXOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceINDEXOUTOFRANGE");

    // INDEXTOOLARGE exception
    class_<te::SpiceINDEXTOOLARGE, base<te::SpiceError>>("exceptions_SpiceINDEXTOOLARGE");

    // INDICESOUTOFORDER exception
    class_<te::SpiceINDICESOUTOFORDER, base<te::SpiceError>>("exceptions_SpiceINDICESOUTOFORDER");

    // INPUTDOESNOTEXIST exception
    class_<te::SpiceINPUTDOESNOTEXIST, base<te::SpiceError>>("exceptions_SpiceINPUTDOESNOTEXIST");

    // INPUTFILENOTEXIST exception
    class_<te::SpiceINPUTFILENOTEXIST, base<te::SpiceError>>("exceptions_SpiceINPUTFILENOTEXIST");

    // INPUTOUTOFBOUNDS exception
    class_<te::SpiceINPUTOUTOFBOUNDS, base<te::SpiceError>>("exceptions_SpiceINPUTOUTOFBOUNDS");

    // INPUTSTOOLARGE exception
    class_<te::SpiceINPUTSTOOLARGE, base<te::SpiceError>>("exceptions_SpiceINPUTSTOOLARGE");

    // INQUIREERROR exception
    class_<te::SpiceINQUIREERROR, base<te::SpiceError>>("exceptions_SpiceINQUIREERROR");

    // INQUIREFAILED exception
    class_<te::SpiceINQUIREFAILED, base<te::SpiceError>>("exceptions_SpiceINQUIREFAILED");

    // INSIDEBODY exception
    class_<te::SpiceINSIDEBODY, base<te::SpiceError>>("exceptions_SpiceINSIDEBODY");

    // INSUFFICIENTANGLES exception
    class_<te::SpiceINSUFFICIENTANGLES, base<te::SpiceError>>("exceptions_SpiceINSUFFICIENTANGLES");

    // INSUFFICIENTDATA exception
    class_<te::SpiceINSUFFICIENTDATA, base<te::SpiceError>>("exceptions_SpiceINSUFFICIENTDATA");

    // INSUFFLEN exception
    class_<te::SpiceINSUFFLEN, base<te::SpiceError>>("exceptions_SpiceINSUFFLEN");

    // INSUFPTRSIZE exception
    class_<te::SpiceINSUFPTRSIZE, base<te::SpiceError>>("exceptions_SpiceINSUFPTRSIZE");

    // INTERVALSTARTNOTFOUND exception
    class_<te::SpiceINTERVALSTARTNOTFOUND, base<te::SpiceError>>("exceptions_SpiceINTERVALSTARTNOTFOUND");

    // INTINDEXTOOSMALL exception
    class_<te::SpiceINTINDEXTOOSMALL, base<te::SpiceError>>("exceptions_SpiceINTINDEXTOOSMALL");

    // INTLENNOTPOS exception
    class_<te::SpiceINTLENNOTPOS, base<te::SpiceError>>("exceptions_SpiceINTLENNOTPOS");

    // INTOUTOFRANGE exception
    class_<te::SpiceINTOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceINTOUTOFRANGE");

    // INVALIDACCESS exception
    class_<te::SpiceINVALIDACCESS, base<te::SpiceError>>("exceptions_SpiceINVALIDACCESS");

    // INVALIDACTION exception
    class_<te::SpiceINVALIDACTION, base<te::SpiceError>>("exceptions_SpiceINVALIDACTION");

    // INVALIDADD exception
    class_<te::SpiceINVALIDADD, base<te::SpiceError>>("exceptions_SpiceINVALIDADD");

    // INVALIDADDRESS exception
    class_<te::SpiceINVALIDADDRESS, base<te::SpiceError>>("exceptions_SpiceINVALIDADDRESS");

    // INVALIDANGLE exception
    class_<te::SpiceINVALIDANGLE, base<te::SpiceError>>("exceptions_SpiceINVALIDANGLE");

    // INVALIDARCHTYPE exception
    class_<te::SpiceINVALIDARCHTYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDARCHTYPE");

    // INVALIDARGUMENT exception
    class_<te::SpiceINVALIDARGUMENT, base<te::SpiceError>>("exceptions_SpiceINVALIDARGUMENT");

    // INVALIDAXIS exception
    class_<te::SpiceINVALIDAXIS, base<te::SpiceError>>("exceptions_SpiceINVALIDAXIS");

    // INVALIDAXISLENGTH exception
    class_<te::SpiceINVALIDAXISLENGTH, base<te::SpiceError>>("exceptions_SpiceINVALIDAXISLENGTH");

    // INVALIDBOUNDS exception
    class_<te::SpiceINVALIDBOUNDS, base<te::SpiceError>>("exceptions_SpiceINVALIDBOUNDS");

    // INVALIDCARDINALITY exception
    class_<te::SpiceINVALIDCARDINALITY, base<te::SpiceError>>("exceptions_SpiceINVALIDCARDINALITY");

    // INVALIDCASE exception
    class_<te::SpiceINVALIDCASE, base<te::SpiceError>>("exceptions_SpiceINVALIDCASE");

    // INVALIDCOLUMN exception
    class_<te::SpiceINVALIDCOLUMN, base<te::SpiceError>>("exceptions_SpiceINVALIDCOLUMN");

    // INVALIDCONSTSTEP exception
    class_<te::SpiceINVALIDCONSTSTEP, base<te::SpiceError>>("exceptions_SpiceINVALIDCONSTSTEP");

    // INVALIDCOUNT exception
    class_<te::SpiceINVALIDCOUNT, base<te::SpiceError>>("exceptions_SpiceINVALIDCOUNT");

    // INVALIDDATA exception
    class_<te::SpiceINVALIDDATA, base<te::SpiceError>>("exceptions_SpiceINVALIDDATA");

    // INVALIDDATACOUNT exception
    class_<te::SpiceINVALIDDATACOUNT, base<te::SpiceError>>("exceptions_SpiceINVALIDDATACOUNT");

    // INVALIDDATATYPE exception
    class_<te::SpiceINVALIDDATATYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDDATATYPE");

    // INVALIDDEGREE exception
    class_<te::SpiceINVALIDDEGREE, base<te::SpiceError>>("exceptions_SpiceINVALIDDEGREE");

    // INVALIDDESCRTIME exception
    class_<te::SpiceINVALIDDESCRTIME, base<te::SpiceError>>("exceptions_SpiceINVALIDDESCRTIME");

    // INVALIDDIMENSION exception
    class_<te::SpiceINVALIDDIMENSION, base<te::SpiceError>>("exceptions_SpiceINVALIDDIMENSION");

    // INVALIDDIRECTION exception
    class_<te::SpiceINVALIDDIRECTION, base<te::SpiceError>>("exceptions_SpiceINVALIDDIRECTION");

    // INVALIDDIVISOR exception
    class_<te::SpiceINVALIDDIVISOR, base<te::SpiceError>>("exceptions_SpiceINVALIDDIVISOR");

    // INVALIDELLIPSE exception
    class_<te::SpiceINVALIDELLIPSE, base<te::SpiceError>>("exceptions_SpiceINVALIDELLIPSE");

    // INVALIDENDPNTSPEC exception
    class_<te::SpiceINVALIDENDPNTSPEC, base<te::SpiceError>>("exceptions_SpiceINVALIDENDPNTSPEC");

    // INVALIDENDPTS exception
    class_<te::SpiceINVALIDENDPTS, base<te::SpiceError>>("exceptions_SpiceINVALIDENDPTS");

    // INVALIDEPOCH exception
    class_<te::SpiceINVALIDEPOCH, base<te::SpiceError>>("exceptions_SpiceINVALIDEPOCH");

    // INVALIDFILETYPE exception
    class_<te::SpiceINVALIDFILETYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDFILETYPE");

    // INVALIDFIXREF exception
    class_<te::SpiceINVALIDFIXREF, base<te::SpiceError>>("exceptions_SpiceINVALIDFIXREF");

    // INVALIDFLAG exception
    class_<te::SpiceINVALIDFLAG, base<te::SpiceError>>("exceptions_SpiceINVALIDFLAG");

    // INVALIDFORMAT exception
    class_<te::SpiceINVALIDFORMAT, base<te::SpiceError>>("exceptions_SpiceINVALIDFORMAT");

    // INVALIDFOV exception
    class_<te::SpiceINVALIDFOV, base<te::SpiceError>>("exceptions_SpiceINVALIDFOV");

    // INVALIDFRAME exception
    class_<te::SpiceINVALIDFRAME, base<te::SpiceError>>("exceptions_SpiceINVALIDFRAME");

    // INVALIDFRAMEDEF exception
    class_<te::SpiceINVALIDFRAMEDEF, base<te::SpiceError>>("exceptions_SpiceINVALIDFRAMEDEF");

    // INVALIDGEOMETRY exception
    class_<te::SpiceINVALIDGEOMETRY, base<te::SpiceError>>("exceptions_SpiceINVALIDGEOMETRY");

    // INVALIDHANDLE exception
    class_<te::SpiceINVALIDHANDLE, base<te::SpiceError>>("exceptions_SpiceINVALIDHANDLE");

    // INVALIDINDEX exception
    class_<te::SpiceINVALIDINDEX, base<te::SpiceError>>("exceptions_SpiceINVALIDINDEX");

    // INVALIDINTEGER exception
    class_<te::SpiceINVALIDINTEGER, base<te::SpiceError>>("exceptions_SpiceINVALIDINTEGER");

    // INVALIDLIMBTYPE exception
    class_<te::SpiceINVALIDLIMBTYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDLIMBTYPE");

    // INVALIDLISTITEM exception
    class_<te::SpiceINVALIDLISTITEM, base<te::SpiceError>>("exceptions_SpiceINVALIDLISTITEM");

    // INVALIDLOCUS exception
    class_<te::SpiceINVALIDLOCUS, base<te::SpiceError>>("exceptions_SpiceINVALIDLOCUS");

    // INVALIDLONEXTENT exception
    class_<te::SpiceINVALIDLONEXTENT, base<te::SpiceError>>("exceptions_SpiceINVALIDLONEXTENT");

    // INVALIDMETADATA exception
    class_<te::SpiceINVALIDMETADATA, base<te::SpiceError>>("exceptions_SpiceINVALIDMETADATA");

    // INVALIDMETHOD exception
    class_<te::SpiceINVALIDMETHOD, base<te::SpiceError>>("exceptions_SpiceINVALIDMETHOD");

    // INVALIDMSGTYPE exception
    class_<te::SpiceINVALIDMSGTYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDMSGTYPE");

    // INVALIDNAME exception
    class_<te::SpiceINVALIDNAME, base<te::SpiceError>>("exceptions_SpiceINVALIDNAME");

    // INVALIDNODE exception
    class_<te::SpiceINVALIDNODE, base<te::SpiceError>>("exceptions_SpiceINVALIDNODE");

    // INVALIDNUMBEROFINTERVALS exception
    class_<te::SpiceINVALIDNUMBEROFINTERVALS, base<te::SpiceError>>("exceptions_SpiceINVALIDNUMBEROFINTERVALS");

    // INVALIDNUMBEROFRECORDS exception
    class_<te::SpiceINVALIDNUMBEROFRECORDS, base<te::SpiceError>>("exceptions_SpiceINVALIDNUMBEROFRECORDS");

    // INVALIDNUMINT exception
    class_<te::SpiceINVALIDNUMINT, base<te::SpiceError>>("exceptions_SpiceINVALIDNUMINT");

    // INVALIDNUMINTS exception
    class_<te::SpiceINVALIDNUMINTS, base<te::SpiceError>>("exceptions_SpiceINVALIDNUMINTS");

    // INVALIDNUMREC exception
    class_<te::SpiceINVALIDNUMREC, base<te::SpiceError>>("exceptions_SpiceINVALIDNUMREC");

    // INVALIDOCCTYPE exception
    class_<te::SpiceINVALIDOCCTYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDOCCTYPE");

    // INVALIDOPERATION exception
    class_<te::SpiceINVALIDOPERATION, base<te::SpiceError>>("exceptions_SpiceINVALIDOPERATION");

    // INVALIDOPTION exception
    class_<te::SpiceINVALIDOPTION, base<te::SpiceError>>("exceptions_SpiceINVALIDOPTION");

    // INVALIDPLANE exception
    class_<te::SpiceINVALIDPLANE, base<te::SpiceError>>("exceptions_SpiceINVALIDPLANE");

    // INVALIDRADII exception
    class_<te::SpiceINVALIDRADII, base<te::SpiceError>>("exceptions_SpiceINVALIDRADII");

    // INVALIDRADIUS exception
    class_<te::SpiceINVALIDRADIUS, base<te::SpiceError>>("exceptions_SpiceINVALIDRADIUS");

    // INVALIDREFFRAME exception
    class_<te::SpiceINVALIDREFFRAME, base<te::SpiceError>>("exceptions_SpiceINVALIDREFFRAME");

    // INVALIDREFVAL exception
    class_<te::SpiceINVALIDREFVAL, base<te::SpiceError>>("exceptions_SpiceINVALIDREFVAL");

    // INVALIDROLLSTEP exception
    class_<te::SpiceINVALIDROLLSTEP, base<te::SpiceError>>("exceptions_SpiceINVALIDROLLSTEP");

    // INVALIDSCALE exception
    class_<te::SpiceINVALIDSCALE, base<te::SpiceError>>("exceptions_SpiceINVALIDSCALE");

    // INVALIDSCLKRATE exception
    class_<te::SpiceINVALIDSCLKRATE, base<te::SpiceError>>("exceptions_SpiceINVALIDSCLKRATE");

    // INVALIDSCLKSTRING exception
    class_<te::SpiceINVALIDSCLKSTRING, base<te::SpiceError>>("exceptions_SpiceINVALIDSCLKSTRING");

    // INVALIDSCLKTIME exception
    class_<te::SpiceINVALIDSCLKTIME, base<te::SpiceError>>("exceptions_SpiceINVALIDSCLKTIME");

    // INVALIDSEARCHSTEP exception
    class_<te::SpiceINVALIDSEARCHSTEP, base<te::SpiceError>>("exceptions_SpiceINVALIDSEARCHSTEP");

    // INVALIDSELECTION exception
    class_<te::SpiceINVALIDSELECTION, base<te::SpiceError>>("exceptions_SpiceINVALIDSELECTION");

    // INVALIDSHADOW exception
    class_<te::SpiceINVALIDSHADOW, base<te::SpiceError>>("exceptions_SpiceINVALIDSHADOW");

    // INVALIDSHAPE exception
    class_<te::SpiceINVALIDSHAPE, base<te::SpiceError>>("exceptions_SpiceINVALIDSHAPE");

    // INVALIDSHAPECOMBO exception
    class_<te::SpiceINVALIDSHAPECOMBO, base<te::SpiceError>>("exceptions_SpiceINVALIDSHAPECOMBO");

    // INVALIDSIZE exception
    class_<te::SpiceINVALIDSIZE, base<te::SpiceError>>("exceptions_SpiceINVALIDSIZE");

    // INVALIDSTARTTIME exception
    class_<te::SpiceINVALIDSTARTTIME, base<te::SpiceError>>("exceptions_SpiceINVALIDSTARTTIME");

    // INVALIDSTATE exception
    class_<te::SpiceINVALIDSTATE, base<te::SpiceError>>("exceptions_SpiceINVALIDSTATE");

    // INVALIDSTEP exception
    class_<te::SpiceINVALIDSTEP, base<te::SpiceError>>("exceptions_SpiceINVALIDSTEP");

    // INVALIDSTEPSIZE exception
    class_<te::SpiceINVALIDSTEPSIZE, base<te::SpiceError>>("exceptions_SpiceINVALIDSTEPSIZE");

    // INVALIDSUBLIST exception
    class_<te::SpiceINVALIDSUBLIST, base<te::SpiceError>>("exceptions_SpiceINVALIDSUBLIST");

    // INVALIDSUBTYPE exception
    class_<te::SpiceINVALIDSUBTYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDSUBTYPE");

    // INVALIDTABLENAME exception
    class_<te::SpiceINVALIDTABLENAME, base<te::SpiceError>>("exceptions_SpiceINVALIDTABLENAME");

    // INVALIDTABLESIZE exception
    class_<te::SpiceINVALIDTABLESIZE, base<te::SpiceError>>("exceptions_SpiceINVALIDTABLESIZE");

    // INVALIDTARGET exception
    class_<te::SpiceINVALIDTARGET, base<te::SpiceError>>("exceptions_SpiceINVALIDTARGET");

    // INVALIDTERMTYPE exception
    class_<te::SpiceINVALIDTERMTYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDTERMTYPE");

    // INVALIDTEXT exception
    class_<te::SpiceINVALIDTEXT, base<te::SpiceError>>("exceptions_SpiceINVALIDTEXT");

    // INVALIDTIMEFORMAT exception
    class_<te::SpiceINVALIDTIMEFORMAT, base<te::SpiceError>>("exceptions_SpiceINVALIDTIMEFORMAT");

    // INVALIDTIMESTRING exception
    class_<te::SpiceINVALIDTIMESTRING, base<te::SpiceError>>("exceptions_SpiceINVALIDTIMESTRING");

    // INVALIDTLEORDER exception
    class_<te::SpiceINVALIDTLEORDER, base<te::SpiceError>>("exceptions_SpiceINVALIDTLEORDER");

    // INVALIDTOLERANCE exception
    class_<te::SpiceINVALIDTOLERANCE, base<te::SpiceError>>("exceptions_SpiceINVALIDTOLERANCE");

    // INVALIDTYPE exception
    class_<te::SpiceINVALIDTYPE, base<te::SpiceError>>("exceptions_SpiceINVALIDTYPE");

    // INVALIDVALUE exception
    class_<te::SpiceINVALIDVALUE, base<te::SpiceError>>("exceptions_SpiceINVALIDVALUE");

    // INVALIDVERTEX exception
    class_<te::SpiceINVALIDVERTEX, base<te::SpiceError>>("exceptions_SpiceINVALIDVERTEX");

    // INVERSTARTSTOPTIME exception
    class_<te::SpiceINVERSTARTSTOPTIME, base<te::SpiceError>>("exceptions_SpiceINVERSTARTSTOPTIME");

    // IRFNOTREC exception
    class_<te::SpiceIRFNOTREC, base<te::SpiceError>>("exceptions_SpiceIRFNOTREC");

    // ITEMNOTFOUND exception
    class_<te::SpiceITEMNOTFOUND, base<te::SpiceError>>("exceptions_SpiceITEMNOTFOUND");

    // ITEMNOTRECOGNIZED exception
    class_<te::SpiceITEMNOTRECOGNIZED, base<te::SpiceError>>("exceptions_SpiceITEMNOTRECOGNIZED");

    // ITERATIONEXCEEDED exception
    class_<te::SpiceITERATIONEXCEEDED, base<te::SpiceError>>("exceptions_SpiceITERATIONEXCEEDED");

    // KERNELNOTLOADED exception
    class_<te::SpiceKERNELNOTLOADED, base<te::SpiceError>>("exceptions_SpiceKERNELNOTLOADED");

    // KERNELPOOLFULL exception
    class_<te::SpiceKERNELPOOLFULL, base<te::SpiceError>>("exceptions_SpiceKERNELPOOLFULL");

    // KERNELVARNOTFOUND exception
    class_<te::SpiceKERNELVARNOTFOUND, base<te::SpiceError>>("exceptions_SpiceKERNELVARNOTFOUND");

    // KERVARSETOVERFLOW exception
    class_<te::SpiceKERVARSETOVERFLOW, base<te::SpiceError>>("exceptions_SpiceKERVARSETOVERFLOW");

    // KERVARTOOBIG exception
    class_<te::SpiceKERVARTOOBIG, base<te::SpiceError>>("exceptions_SpiceKERVARTOOBIG");

    // KEYWORDNOTFOUND exception
    class_<te::SpiceKEYWORDNOTFOUND, base<te::SpiceError>>("exceptions_SpiceKEYWORDNOTFOUND");

    // LBCORRUPTED exception
    class_<te::SpiceLBCORRUPTED, base<te::SpiceError>>("exceptions_SpiceLBCORRUPTED");

    // LBLINETOOLONG exception
    class_<te::SpiceLBLINETOOLONG, base<te::SpiceError>>("exceptions_SpiceLBLINETOOLONG");

    // LBNOSUCHLINE exception
    class_<te::SpiceLBNOSUCHLINE, base<te::SpiceError>>("exceptions_SpiceLBNOSUCHLINE");

    // LBTOOMANYLINES exception
    class_<te::SpiceLBTOOMANYLINES, base<te::SpiceError>>("exceptions_SpiceLBTOOMANYLINES");

    // LOWERBOUNDTOOLOW exception
    class_<te::SpiceLOWERBOUNDTOOLOW, base<te::SpiceError>>("exceptions_SpiceLOWERBOUNDTOOLOW");

    // LSKDOESNTEXIST exception
    class_<te::SpiceLSKDOESNTEXIST, base<te::SpiceError>>("exceptions_SpiceLSKDOESNTEXIST");

    // MALFORMEDSEGMENT exception
    class_<te::SpiceMALFORMEDSEGMENT, base<te::SpiceError>>("exceptions_SpiceMALFORMEDSEGMENT");

    // MARKERNOTFOUND exception
    class_<te::SpiceMARKERNOTFOUND, base<te::SpiceError>>("exceptions_SpiceMARKERNOTFOUND");

    // MESSAGETOOLONG exception
    class_<te::SpiceMESSAGETOOLONG, base<te::SpiceError>>("exceptions_SpiceMESSAGETOOLONG");

    // MISMATCHFROMTIMETYPE exception
    class_<te::SpiceMISMATCHFROMTIMETYPE, base<te::SpiceError>>("exceptions_SpiceMISMATCHFROMTIMETYPE");

    // MISMATCHOUTPUTFORMAT exception
    class_<te::SpiceMISMATCHOUTPUTFORMAT, base<te::SpiceError>>("exceptions_SpiceMISMATCHOUTPUTFORMAT");

    // MISMATCHTOTIMETYPE exception
    class_<te::SpiceMISMATCHTOTIMETYPE, base<te::SpiceError>>("exceptions_SpiceMISMATCHTOTIMETYPE");

    // MISSINGARGUMENTS exception
    class_<te::SpiceMISSINGARGUMENTS, base<te::SpiceError>>("exceptions_SpiceMISSINGARGUMENTS");

    // MISSINGCENTER exception
    class_<te::SpiceMISSINGCENTER, base<te::SpiceError>>("exceptions_SpiceMISSINGCENTER");

    // MISSINGCOLSTEP exception
    class_<te::SpiceMISSINGCOLSTEP, base<te::SpiceError>>("exceptions_SpiceMISSINGCOLSTEP");

    // MISSINGCOORDBOUND exception
    class_<te::SpiceMISSINGCOORDBOUND, base<te::SpiceError>>("exceptions_SpiceMISSINGCOORDBOUND");

    // MISSINGCOORDSYS exception
    class_<te::SpiceMISSINGCOORDSYS, base<te::SpiceError>>("exceptions_SpiceMISSINGCOORDSYS");

    // MISSINGDATA exception
    class_<te::SpiceMISSINGDATA, base<te::SpiceError>>("exceptions_SpiceMISSINGDATA");

    // MISSINGDATACLASS exception
    class_<te::SpiceMISSINGDATACLASS, base<te::SpiceError>>("exceptions_SpiceMISSINGDATACLASS");

    // MISSINGDATAORDERTK exception
    class_<te::SpiceMISSINGDATAORDERTK, base<te::SpiceError>>("exceptions_SpiceMISSINGDATAORDERTK");

    // MISSINGDATATYPE exception
    class_<te::SpiceMISSINGDATATYPE, base<te::SpiceError>>("exceptions_SpiceMISSINGDATATYPE");

    // MISSINGEOT exception
    class_<te::SpiceMISSINGEOT, base<te::SpiceError>>("exceptions_SpiceMISSINGEOT");

    // MISSINGEPOCHTOKEN exception
    class_<te::SpiceMISSINGEPOCHTOKEN, base<te::SpiceError>>("exceptions_SpiceMISSINGEPOCHTOKEN");

    // MISSINGFRAME exception
    class_<te::SpiceMISSINGFRAME, base<te::SpiceError>>("exceptions_SpiceMISSINGFRAME");

    // MISSINGFRAMEVAR exception
    class_<te::SpiceMISSINGFRAMEVAR, base<te::SpiceError>>("exceptions_SpiceMISSINGFRAMEVAR");

    // MISSINGGEOCONSTS exception
    class_<te::SpiceMISSINGGEOCONSTS, base<te::SpiceError>>("exceptions_SpiceMISSINGGEOCONSTS");

    // MISSINGHEIGHTREF exception
    class_<te::SpiceMISSINGHEIGHTREF, base<te::SpiceError>>("exceptions_SpiceMISSINGHEIGHTREF");

    // MISSINGHSCALE exception
    class_<te::SpiceMISSINGHSCALE, base<te::SpiceError>>("exceptions_SpiceMISSINGHSCALE");

    // MISSINGKPV exception
    class_<te::SpiceMISSINGKPV, base<te::SpiceError>>("exceptions_SpiceMISSINGKPV");

    // MISSINGLEFTCOR exception
    class_<te::SpiceMISSINGLEFTCOR, base<te::SpiceError>>("exceptions_SpiceMISSINGLEFTCOR");

    // MISSINGLEFTRTFLAG exception
    class_<te::SpiceMISSINGLEFTRTFLAG, base<te::SpiceError>>("exceptions_SpiceMISSINGLEFTRTFLAG");

    // MISSINGNCAPFLAG exception
    class_<te::SpiceMISSINGNCAPFLAG, base<te::SpiceError>>("exceptions_SpiceMISSINGNCAPFLAG");

    // MISSINGNCOLS exception
    class_<te::SpiceMISSINGNCOLS, base<te::SpiceError>>("exceptions_SpiceMISSINGNCOLS");

    // MISSINGNROWS exception
    class_<te::SpiceMISSINGNROWS, base<te::SpiceError>>("exceptions_SpiceMISSINGNROWS");

    // MISSINGPLATETYPE exception
    class_<te::SpiceMISSINGPLATETYPE, base<te::SpiceError>>("exceptions_SpiceMISSINGPLATETYPE");

    // MISSINGROWMAJFLAG exception
    class_<te::SpiceMISSINGROWMAJFLAG, base<te::SpiceError>>("exceptions_SpiceMISSINGROWMAJFLAG");

    // MISSINGROWSTEP exception
    class_<te::SpiceMISSINGROWSTEP, base<te::SpiceError>>("exceptions_SpiceMISSINGROWSTEP");

    // MISSINGSCAPFLAG exception
    class_<te::SpiceMISSINGSCAPFLAG, base<te::SpiceError>>("exceptions_SpiceMISSINGSCAPFLAG");

    // MISSINGSURFACE exception
    class_<te::SpiceMISSINGSURFACE, base<te::SpiceError>>("exceptions_SpiceMISSINGSURFACE");

    // MISSINGTIMEINFO exception
    class_<te::SpiceMISSINGTIMEINFO, base<te::SpiceError>>("exceptions_SpiceMISSINGTIMEINFO");

    // MISSINGTLEKEYWORD exception
    class_<te::SpiceMISSINGTLEKEYWORD, base<te::SpiceError>>("exceptions_SpiceMISSINGTLEKEYWORD");

    // MISSINGTOPCOR exception
    class_<te::SpiceMISSINGTOPCOR, base<te::SpiceError>>("exceptions_SpiceMISSINGTOPCOR");

    // MISSINGTOPDOWNFLAG exception
    class_<te::SpiceMISSINGTOPDOWNFLAG, base<te::SpiceError>>("exceptions_SpiceMISSINGTOPDOWNFLAG");

    // MISSINGVALUE exception
    class_<te::SpiceMISSINGVALUE, base<te::SpiceError>>("exceptions_SpiceMISSINGVALUE");

    // MISSINGVOXELSCALE exception
    class_<te::SpiceMISSINGVOXELSCALE, base<te::SpiceError>>("exceptions_SpiceMISSINGVOXELSCALE");

    // MISSINGWRAPFLAG exception
    class_<te::SpiceMISSINGWRAPFLAG, base<te::SpiceError>>("exceptions_SpiceMISSINGWRAPFLAG");

    // NAMENOTUNIQUE exception
    class_<te::SpiceNAMENOTUNIQUE, base<te::SpiceError>>("exceptions_SpiceNAMENOTUNIQUE");

    // NAMESNOTRESOLVED exception
    class_<te::SpiceNAMESNOTRESOLVED, base<te::SpiceError>>("exceptions_SpiceNAMESNOTRESOLVED");

    // NAMETABLEFULL exception
    class_<te::SpiceNAMETABLEFULL, base<te::SpiceError>>("exceptions_SpiceNAMETABLEFULL");

    // NARATESFLAG exception
    class_<te::SpiceNARATESFLAG, base<te::SpiceError>>("exceptions_SpiceNARATESFLAG");

    // NEGATIVETOL exception
    class_<te::SpiceNEGATIVETOL, base<te::SpiceError>>("exceptions_SpiceNEGATIVETOL");

    // NOACCEPTABLEDATA exception
    class_<te::SpiceNOACCEPTABLEDATA, base<te::SpiceError>>("exceptions_SpiceNOACCEPTABLEDATA");

    // NOANGULARRATEFLAG exception
    class_<te::SpiceNOANGULARRATEFLAG, base<te::SpiceError>>("exceptions_SpiceNOANGULARRATEFLAG");

    // NOARRAYSTARTED exception
    class_<te::SpiceNOARRAYSTARTED, base<te::SpiceError>>("exceptions_SpiceNOARRAYSTARTED");

    // NOATTIME exception
    class_<te::SpiceNOATTIME, base<te::SpiceError>>("exceptions_SpiceNOATTIME");

    // NOAVDATA exception
    class_<te::SpiceNOAVDATA, base<te::SpiceError>>("exceptions_SpiceNOAVDATA");

    // NOBODYID exception
    class_<te::SpiceNOBODYID, base<te::SpiceError>>("exceptions_SpiceNOBODYID");

    // NOCANDOSPKSPCKS exception
    class_<te::SpiceNOCANDOSPKSPCKS, base<te::SpiceError>>("exceptions_SpiceNOCANDOSPKSPCKS");

    // NOCENTERIDORNAME exception
    class_<te::SpiceNOCENTERIDORNAME, base<te::SpiceError>>("exceptions_SpiceNOCENTERIDORNAME");

    // NOCKSEGMENTTYPE exception
    class_<te::SpiceNOCKSEGMENTTYPE, base<te::SpiceError>>("exceptions_SpiceNOCKSEGMENTTYPE");

    // NOCLASS exception
    class_<te::SpiceNOCLASS, base<te::SpiceError>>("exceptions_SpiceNOCLASS");

    // NOCOMMENTSFILE exception
    class_<te::SpiceNOCOMMENTSFILE, base<te::SpiceError>>("exceptions_SpiceNOCOMMENTSFILE");

    // NOCONVERG exception
    class_<te::SpiceNOCONVERG, base<te::SpiceError>>("exceptions_SpiceNOCONVERG");

    // NOCONVERGENCE exception
    class_<te::SpiceNOCONVERGENCE, base<te::SpiceError>>("exceptions_SpiceNOCONVERGENCE");

    // NOCURRENTARRAY exception
    class_<te::SpiceNOCURRENTARRAY, base<te::SpiceError>>("exceptions_SpiceNOCURRENTARRAY");

    // NODATAORDER exception
    class_<te::SpiceNODATAORDER, base<te::SpiceError>>("exceptions_SpiceNODATAORDER");

    // NODATATYPEFLAG exception
    class_<te::SpiceNODATATYPEFLAG, base<te::SpiceError>>("exceptions_SpiceNODATATYPEFLAG");

    // NODELIMCHARACTER exception
    class_<te::SpiceNODELIMCHARACTER, base<te::SpiceError>>("exceptions_SpiceNODELIMCHARACTER");

    // NODETOOFULL exception
    class_<te::SpiceNODETOOFULL, base<te::SpiceError>>("exceptions_SpiceNODETOOFULL");

    // NODSKSEGMENT exception
    class_<te::SpiceNODSKSEGMENT, base<te::SpiceError>>("exceptions_SpiceNODSKSEGMENT");

    // NODSKSEGMENTS exception
    class_<te::SpiceNODSKSEGMENTS, base<te::SpiceError>>("exceptions_SpiceNODSKSEGMENTS");

    // NOENVVARIABLE exception
    class_<te::SpiceNOENVVARIABLE, base<te::SpiceError>>("exceptions_SpiceNOENVVARIABLE");

    // NOEULERANGLEUNITS exception
    class_<te::SpiceNOEULERANGLEUNITS, base<te::SpiceError>>("exceptions_SpiceNOEULERANGLEUNITS");

    // NOFILENAMES exception
    class_<te::SpiceNOFILENAMES, base<te::SpiceError>>("exceptions_SpiceNOFILENAMES");

    // NOFILES exception
    class_<te::SpiceNOFILES, base<te::SpiceError>>("exceptions_SpiceNOFILES");

    // NOFILESPEC exception
    class_<te::SpiceNOFILESPEC, base<te::SpiceError>>("exceptions_SpiceNOFILESPEC");

    // NOFRAME exception
    class_<te::SpiceNOFRAME, base<te::SpiceError>>("exceptions_SpiceNOFRAME");

    // NOFRAMECONNECT exception
    class_<te::SpiceNOFRAMECONNECT, base<te::SpiceError>>("exceptions_SpiceNOFRAMECONNECT");

    // NOFRAMEDATA exception
    class_<te::SpiceNOFRAMEDATA, base<te::SpiceError>>("exceptions_SpiceNOFRAMEDATA");

    // NOFRAMEINFO exception
    class_<te::SpiceNOFRAMEINFO, base<te::SpiceError>>("exceptions_SpiceNOFRAMEINFO");

    // NOFRAMENAME exception
    class_<te::SpiceNOFRAMENAME, base<te::SpiceError>>("exceptions_SpiceNOFRAMENAME");

    // NOFRAMESKERNELNAME exception
    class_<te::SpiceNOFRAMESKERNELNAME, base<te::SpiceError>>("exceptions_SpiceNOFRAMESKERNELNAME");

    // NOFREELOGICALUNIT exception
    class_<te::SpiceNOFREELOGICALUNIT, base<te::SpiceError>>("exceptions_SpiceNOFREELOGICALUNIT");

    // NOFREENODES exception
    class_<te::SpiceNOFREENODES, base<te::SpiceError>>("exceptions_SpiceNOFREENODES");

    // NOFROMTIME exception
    class_<te::SpiceNOFROMTIME, base<te::SpiceError>>("exceptions_SpiceNOFROMTIME");

    // NOFROMTIMESYSTEM exception
    class_<te::SpiceNOFROMTIMESYSTEM, base<te::SpiceError>>("exceptions_SpiceNOFROMTIMESYSTEM");

    // NOHEADNODE exception
    class_<te::SpiceNOHEADNODE, base<te::SpiceError>>("exceptions_SpiceNOHEADNODE");

    // NOINPUTDATATYPE exception
    class_<te::SpiceNOINPUTDATATYPE, base<te::SpiceError>>("exceptions_SpiceNOINPUTDATATYPE");

    // NOINPUTFILENAME exception
    class_<te::SpiceNOINPUTFILENAME, base<te::SpiceError>>("exceptions_SpiceNOINPUTFILENAME");

    // NOINSTRUMENTID exception
    class_<te::SpiceNOINSTRUMENTID, base<te::SpiceError>>("exceptions_SpiceNOINSTRUMENTID");

    // NOINTERVAL exception
    class_<te::SpiceNOINTERVAL, base<te::SpiceError>>("exceptions_SpiceNOINTERVAL");

    // NOKERNELLOADED exception
    class_<te::SpiceNOKERNELLOADED, base<te::SpiceError>>("exceptions_SpiceNOKERNELLOADED");

    // NOLANDINGTIME exception
    class_<te::SpiceNOLANDINGTIME, base<te::SpiceError>>("exceptions_SpiceNOLANDINGTIME");

    // NOLEAPSECONDS exception
    class_<te::SpiceNOLEAPSECONDS, base<te::SpiceError>>("exceptions_SpiceNOLEAPSECONDS");

    // NOLINESPERRECCOUNT exception
    class_<te::SpiceNOLINESPERRECCOUNT, base<te::SpiceError>>("exceptions_SpiceNOLINESPERRECCOUNT");

    // NOLISTFILENAME exception
    class_<te::SpiceNOLISTFILENAME, base<te::SpiceError>>("exceptions_SpiceNOLISTFILENAME");

    // NOLOADEDDSKFILES exception
    class_<te::SpiceNOLOADEDDSKFILES, base<te::SpiceError>>("exceptions_SpiceNOLOADEDDSKFILES");

    // NOLOADEDFILES exception
    class_<te::SpiceNOLOADEDFILES, base<te::SpiceError>>("exceptions_SpiceNOLOADEDFILES");

    // NOLSKFILENAME exception
    class_<te::SpiceNOLSKFILENAME, base<te::SpiceError>>("exceptions_SpiceNOLSKFILENAME");

    // NOMOREROOM exception
    class_<te::SpiceNOMOREROOM, base<te::SpiceError>>("exceptions_SpiceNOMOREROOM");

    // NONCONICMOTION exception
    class_<te::SpiceNONCONICMOTION, base<te::SpiceError>>("exceptions_SpiceNONCONICMOTION");

    // NONDISTINCTPAIR exception
    class_<te::SpiceNONDISTINCTPAIR, base<te::SpiceError>>("exceptions_SpiceNONDISTINCTPAIR");

    // NONEMPTYENTRY exception
    class_<te::SpiceNONEMPTYENTRY, base<te::SpiceError>>("exceptions_SpiceNONEMPTYENTRY");

    // NONEMPTYTREE exception
    class_<te::SpiceNONEMPTYTREE, base<te::SpiceError>>("exceptions_SpiceNONEMPTYTREE");

    // NONEXISTELEMENTS exception
    class_<te::SpiceNONEXISTELEMENTS, base<te::SpiceError>>("exceptions_SpiceNONEXISTELEMENTS");

    // NONINTEGERFIELD exception
    class_<te::SpiceNONINTEGERFIELD, base<te::SpiceError>>("exceptions_SpiceNONINTEGERFIELD");

    // NONNUMERICSTRING exception
    class_<te::SpiceNONNUMERICSTRING, base<te::SpiceError>>("exceptions_SpiceNONNUMERICSTRING");

    // NONPOSBUFLENGTH exception
    class_<te::SpiceNONPOSBUFLENGTH, base<te::SpiceError>>("exceptions_SpiceNONPOSBUFLENGTH");

    // NONPOSITIVEAXIS exception
    class_<te::SpiceNONPOSITIVEAXIS, base<te::SpiceError>>("exceptions_SpiceNONPOSITIVEAXIS");

    // NONPOSITIVEMASS exception
    class_<te::SpiceNONPOSITIVEMASS, base<te::SpiceError>>("exceptions_SpiceNONPOSITIVEMASS");

    // NONPOSITIVERADIUS exception
    class_<te::SpiceNONPOSITIVERADIUS, base<te::SpiceError>>("exceptions_SpiceNONPOSITIVERADIUS");

    // NONPOSITIVESCALE exception
    class_<te::SpiceNONPOSITIVESCALE, base<te::SpiceError>>("exceptions_SpiceNONPOSITIVESCALE");

    // NONPOSITIVEVALUE exception
    class_<te::SpiceNONPOSITIVEVALUE, base<te::SpiceError>>("exceptions_SpiceNONPOSITIVEVALUE");

    // NONPOSPACKETSIZE exception
    class_<te::SpiceNONPOSPACKETSIZE, base<te::SpiceError>>("exceptions_SpiceNONPOSPACKETSIZE");

    // NONPRINTABLECHARS exception
    class_<te::SpiceNONPRINTABLECHARS, base<te::SpiceError>>("exceptions_SpiceNONPRINTABLECHARS");

    // NONPRINTINGCHAR exception
    class_<te::SpiceNONPRINTINGCHAR, base<te::SpiceError>>("exceptions_SpiceNONPRINTINGCHAR");

    // NONPRINTINGCHARS exception
    class_<te::SpiceNONPRINTINGCHARS, base<te::SpiceError>>("exceptions_SpiceNONPRINTINGCHARS");

    // NONUNITNORMAL exception
    class_<te::SpiceNONUNITNORMAL, base<te::SpiceError>>("exceptions_SpiceNONUNITNORMAL");

    // NOOBJECTIDORNAME exception
    class_<te::SpiceNOOBJECTIDORNAME, base<te::SpiceError>>("exceptions_SpiceNOOBJECTIDORNAME");

    // NOOFFSETANGLEAXES exception
    class_<te::SpiceNOOFFSETANGLEAXES, base<te::SpiceError>>("exceptions_SpiceNOOFFSETANGLEAXES");

    // NOOFFSETANGLEUNITS exception
    class_<te::SpiceNOOFFSETANGLEUNITS, base<te::SpiceError>>("exceptions_SpiceNOOFFSETANGLEUNITS");

    // NOOUTPUTFILENAME exception
    class_<te::SpiceNOOUTPUTFILENAME, base<te::SpiceError>>("exceptions_SpiceNOOUTPUTFILENAME");

    // NOOUTPUTSPKTYPE exception
    class_<te::SpiceNOOUTPUTSPKTYPE, base<te::SpiceError>>("exceptions_SpiceNOOUTPUTSPKTYPE");

    // NOPARTITION exception
    class_<te::SpiceNOPARTITION, base<te::SpiceError>>("exceptions_SpiceNOPARTITION");

    // NOPICTURE exception
    class_<te::SpiceNOPICTURE, base<te::SpiceError>>("exceptions_SpiceNOPICTURE");

    // NOPOLYNOMIALDEGREE exception
    class_<te::SpiceNOPOLYNOMIALDEGREE, base<te::SpiceError>>("exceptions_SpiceNOPOLYNOMIALDEGREE");

    // NOPRECESSIONTYPE exception
    class_<te::SpiceNOPRECESSIONTYPE, base<te::SpiceError>>("exceptions_SpiceNOPRECESSIONTYPE");

    // NOPRODUCERID exception
    class_<te::SpiceNOPRODUCERID, base<te::SpiceError>>("exceptions_SpiceNOPRODUCERID");

    // NOROTATIONORDER exception
    class_<te::SpiceNOROTATIONORDER, base<te::SpiceError>>("exceptions_SpiceNOROTATIONORDER");

    // NOSCID exception
    class_<te::SpiceNOSCID, base<te::SpiceError>>("exceptions_SpiceNOSCID");

    // NOSCLKFILENAMES exception
    class_<te::SpiceNOSCLKFILENAMES, base<te::SpiceError>>("exceptions_SpiceNOSCLKFILENAMES");

    // NOSECONDLINE exception
    class_<te::SpiceNOSECONDLINE, base<te::SpiceError>>("exceptions_SpiceNOSECONDLINE");

    // NOSEGMENTSFOUND exception
    class_<te::SpiceNOSEGMENTSFOUND, base<te::SpiceError>>("exceptions_SpiceNOSEGMENTSFOUND");

    // NOSEPARATION exception
    class_<te::SpiceNOSEPARATION, base<te::SpiceError>>("exceptions_SpiceNOSEPARATION");

    // NOSLKFILENAME exception
    class_<te::SpiceNOSLKFILENAME, base<te::SpiceError>>("exceptions_SpiceNOSLKFILENAME");

    // NOSOLMARKER exception
    class_<te::SpiceNOSOLMARKER, base<te::SpiceError>>("exceptions_SpiceNOSOLMARKER");

    // NOSPACECRAFTID exception
    class_<te::SpiceNOSPACECRAFTID, base<te::SpiceError>>("exceptions_SpiceNOSPACECRAFTID");

    // NOSTARTTIME exception
    class_<te::SpiceNOSTARTTIME, base<te::SpiceError>>("exceptions_SpiceNOSTARTTIME");

    // NOSTOPTIME exception
    class_<te::SpiceNOSTOPTIME, base<te::SpiceError>>("exceptions_SpiceNOSTOPTIME");

    // NOSUCHFILE exception
    class_<te::SpiceNOSUCHFILE, base<te::SpiceError>>("exceptions_SpiceNOSUCHFILE");

    // NOSUCHHANDLE exception
    class_<te::SpiceNOSUCHHANDLE, base<te::SpiceError>>("exceptions_SpiceNOSUCHHANDLE");

    // NOSUCHSYMBOL exception
    class_<te::SpiceNOSUCHSYMBOL, base<te::SpiceError>>("exceptions_SpiceNOSUCHSYMBOL");

    // NOSUNGM exception
    class_<te::SpiceNOSUNGM, base<te::SpiceError>>("exceptions_SpiceNOSUNGM");

    // NOTABINARYKERNEL exception
    class_<te::SpiceNOTABINARYKERNEL, base<te::SpiceError>>("exceptions_SpiceNOTABINARYKERNEL");

    // NOTACKFILE exception
    class_<te::SpiceNOTACKFILE, base<te::SpiceError>>("exceptions_SpiceNOTACKFILE");

    // NOTADAFFILE exception
    class_<te::SpiceNOTADAFFILE, base<te::SpiceError>>("exceptions_SpiceNOTADAFFILE");

    // NOTADASFILE exception
    class_<te::SpiceNOTADASFILE, base<te::SpiceError>>("exceptions_SpiceNOTADASFILE");

    // NOTADPNUMBER exception
    class_<te::SpiceNOTADPNUMBER, base<te::SpiceError>>("exceptions_SpiceNOTADPNUMBER");

    // NOTANDPNUMBER exception
    class_<te::SpiceNOTANDPNUMBER, base<te::SpiceError>>("exceptions_SpiceNOTANDPNUMBER");

    // NOTANINTEGER exception
    class_<te::SpiceNOTANINTEGER, base<te::SpiceError>>("exceptions_SpiceNOTANINTEGER");

    // NOTANINTEGERNUMBER exception
    class_<te::SpiceNOTANINTEGERNUMBER, base<te::SpiceError>>("exceptions_SpiceNOTANINTEGERNUMBER");

    // NOTANINTNUMBER exception
    class_<te::SpiceNOTANINTNUMBER, base<te::SpiceError>>("exceptions_SpiceNOTANINTNUMBER");

    // NOTAPCKFILE exception
    class_<te::SpiceNOTAPCKFILE, base<te::SpiceError>>("exceptions_SpiceNOTAPCKFILE");

    // NOTAROTATION exception
    class_<te::SpiceNOTAROTATION, base<te::SpiceError>>("exceptions_SpiceNOTAROTATION");

    // NOTATEXTFILE exception
    class_<te::SpiceNOTATEXTFILE, base<te::SpiceError>>("exceptions_SpiceNOTATEXTFILE");

    // NOTATRANSFERFILE exception
    class_<te::SpiceNOTATRANSFERFILE, base<te::SpiceError>>("exceptions_SpiceNOTATRANSFERFILE");

    // NOTCOMPUTABLE exception
    class_<te::SpiceNOTCOMPUTABLE, base<te::SpiceError>>("exceptions_SpiceNOTCOMPUTABLE");

    // NOTDIMENSIONALLYEQUIV exception
    class_<te::SpiceNOTDIMENSIONALLYEQUIV, base<te::SpiceError>>("exceptions_SpiceNOTDIMENSIONALLYEQUIV");

    // NOTDISJOINT exception
    class_<te::SpiceNOTDISJOINT, base<te::SpiceError>>("exceptions_SpiceNOTDISJOINT");

    // NOTDISTINCT exception
    class_<te::SpiceNOTDISTINCT, base<te::SpiceError>>("exceptions_SpiceNOTDISTINCT");

    // NOTENOUGHPEAS exception
    class_<te::SpiceNOTENOUGHPEAS, base<te::SpiceError>>("exceptions_SpiceNOTENOUGHPEAS");

    // NOTIMETYPEFLAG exception
    class_<te::SpiceNOTIMETYPEFLAG, base<te::SpiceError>>("exceptions_SpiceNOTIMETYPEFLAG");

    // NOTINDEXED exception
    class_<te::SpiceNOTINDEXED, base<te::SpiceError>>("exceptions_SpiceNOTINDEXED");

    // NOTINITIALIZED exception
    class_<te::SpiceNOTINITIALIZED, base<te::SpiceError>>("exceptions_SpiceNOTINITIALIZED");

    // NOTINPART exception
    class_<te::SpiceNOTINPART, base<te::SpiceError>>("exceptions_SpiceNOTINPART");

    // NOTLEDATAFOROBJECT exception
    class_<te::SpiceNOTLEDATAFOROBJECT, base<te::SpiceError>>("exceptions_SpiceNOTLEDATAFOROBJECT");

    // NOTLEGALCB exception
    class_<te::SpiceNOTLEGALCB, base<te::SpiceError>>("exceptions_SpiceNOTLEGALCB");

    // NOTOTIME exception
    class_<te::SpiceNOTOTIME, base<te::SpiceError>>("exceptions_SpiceNOTOTIME");

    // NOTOTIMESYSTEM exception
    class_<te::SpiceNOTOTIMESYSTEM, base<te::SpiceError>>("exceptions_SpiceNOTOTIMESYSTEM");

    // NOTRANSLATION exception
    class_<te::SpiceNOTRANSLATION, base<te::SpiceError>>("exceptions_SpiceNOTRANSLATION");

    // NOTRECOGNIZED exception
    class_<te::SpiceNOTRECOGNIZED, base<te::SpiceError>>("exceptions_SpiceNOTRECOGNIZED");

    // NOTSEMCHECKED exception
    class_<te::SpiceNOTSEMCHECKED, base<te::SpiceError>>("exceptions_SpiceNOTSEMCHECKED");

    // NOTSUPPORTED exception
    class_<te::SpiceNOTSUPPORTED, base<te::SpiceError>>("exceptions_SpiceNOTSUPPORTED");

    // NOTTWOFIELDSCLK exception
    class_<te::SpiceNOTTWOFIELDSCLK, base<te::SpiceError>>("exceptions_SpiceNOTTWOFIELDSCLK");

    // NOTTWOMODULI exception
    class_<te::SpiceNOTTWOMODULI, base<te::SpiceError>>("exceptions_SpiceNOTTWOMODULI");

    // NOTTWOOFFSETS exception
    class_<te::SpiceNOTTWOOFFSETS, base<te::SpiceError>>("exceptions_SpiceNOTTWOOFFSETS");

    // NOUNITSPEC exception
    class_<te::SpiceNOUNITSPEC, base<te::SpiceError>>("exceptions_SpiceNOUNITSPEC");

    // NUMBEREXPECTED exception
    class_<te::SpiceNUMBEREXPECTED, base<te::SpiceError>>("exceptions_SpiceNUMBEREXPECTED");

    // NUMCOEFFSNOTPOS exception
    class_<te::SpiceNUMCOEFFSNOTPOS, base<te::SpiceError>>("exceptions_SpiceNUMCOEFFSNOTPOS");

    // NUMERICOVERFLOW exception
    class_<te::SpiceNUMERICOVERFLOW, base<te::SpiceError>>("exceptions_SpiceNUMERICOVERFLOW");

    // NUMPACKETSNOTPOS exception
    class_<te::SpiceNUMPACKETSNOTPOS, base<te::SpiceError>>("exceptions_SpiceNUMPACKETSNOTPOS");

    // NUMPARTSUNEQUAL exception
    class_<te::SpiceNUMPARTSUNEQUAL, base<te::SpiceError>>("exceptions_SpiceNUMPARTSUNEQUAL");

    // NUMSTATESNOTPOS exception
    class_<te::SpiceNUMSTATESNOTPOS, base<te::SpiceError>>("exceptions_SpiceNUMSTATESNOTPOS");

    // OBJECTLISTFULL exception
    class_<te::SpiceOBJECTLISTFULL, base<te::SpiceError>>("exceptions_SpiceOBJECTLISTFULL");

    // OBJECTSTOOCLOSE exception
    class_<te::SpiceOBJECTSTOOCLOSE, base<te::SpiceError>>("exceptions_SpiceOBJECTSTOOCLOSE");

    // ORBITDECAY exception
    class_<te::SpiceORBITDECAY, base<te::SpiceError>>("exceptions_SpiceORBITDECAY");

    // OUTOFPLACEDELIMITER exception
    class_<te::SpiceOUTOFPLACEDELIMITER, base<te::SpiceError>>("exceptions_SpiceOUTOFPLACEDELIMITER");

    // OUTOFRANGE exception
    class_<te::SpiceOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceOUTOFRANGE");

    // OUTOFROOM exception
    class_<te::SpiceOUTOFROOM, base<te::SpiceError>>("exceptions_SpiceOUTOFROOM");

    // OUTPUTFILEEXISTS exception
    class_<te::SpiceOUTPUTFILEEXISTS, base<te::SpiceError>>("exceptions_SpiceOUTPUTFILEEXISTS");

    // OUTPUTISNOTSPK exception
    class_<te::SpiceOUTPUTISNOTSPK, base<te::SpiceError>>("exceptions_SpiceOUTPUTISNOTSPK");

    // OUTPUTTOOLONG exception
    class_<te::SpiceOUTPUTTOOLONG, base<te::SpiceError>>("exceptions_SpiceOUTPUTTOOLONG");

    // OUTPUTTOOSHORT exception
    class_<te::SpiceOUTPUTTOOSHORT, base<te::SpiceError>>("exceptions_SpiceOUTPUTTOOSHORT");

    // PARSERNOTREADY exception
    class_<te::SpicePARSERNOTREADY, base<te::SpiceError>>("exceptions_SpicePARSERNOTREADY");

    // PARTIALFRAMESPEC exception
    class_<te::SpicePARTIALFRAMESPEC, base<te::SpiceError>>("exceptions_SpicePARTIALFRAMESPEC");

    // PASTENDSTR exception
    class_<te::SpicePASTENDSTR, base<te::SpiceError>>("exceptions_SpicePASTENDSTR");

    // PATHMISMATCH exception
    class_<te::SpicePATHMISMATCH, base<te::SpiceError>>("exceptions_SpicePATHMISMATCH");

    // PATHTOOLONG exception
    class_<te::SpicePATHTOOLONG, base<te::SpiceError>>("exceptions_SpicePATHTOOLONG");

    // PCKDOESNTEXIST exception
    class_<te::SpicePCKDOESNTEXIST, base<te::SpiceError>>("exceptions_SpicePCKDOESNTEXIST");

    // PCKFILE exception
    class_<te::SpicePCKFILE, base<te::SpiceError>>("exceptions_SpicePCKFILE");

    // PCKFILETABLEFULL exception
    class_<te::SpicePCKFILETABLEFULL, base<te::SpiceError>>("exceptions_SpicePCKFILETABLEFULL");

    // PCKKRECTOOLARGE exception
    class_<te::SpicePCKKRECTOOLARGE, base<te::SpiceError>>("exceptions_SpicePCKKRECTOOLARGE");

    // PLATELISTTOOSMALL exception
    class_<te::SpicePLATELISTTOOSMALL, base<te::SpiceError>>("exceptions_SpicePLATELISTTOOSMALL");

    // POINTEROUTOFRANGE exception
    class_<te::SpicePOINTEROUTOFRANGE, base<te::SpiceError>>("exceptions_SpicePOINTEROUTOFRANGE");

    // POINTERSETTOOBIG exception
    class_<te::SpicePOINTERSETTOOBIG, base<te::SpiceError>>("exceptions_SpicePOINTERSETTOOBIG");

    // POINTERTABLEFULL exception
    class_<te::SpicePOINTERTABLEFULL, base<te::SpiceError>>("exceptions_SpicePOINTERTABLEFULL");

    // POINTNOTFOUND exception
    class_<te::SpicePOINTNOTFOUND, base<te::SpiceError>>("exceptions_SpicePOINTNOTFOUND");

    // POINTNOTINSEGMENT exception
    class_<te::SpicePOINTNOTINSEGMENT, base<te::SpiceError>>("exceptions_SpicePOINTNOTINSEGMENT");

    // POINTNOTONSURFACE exception
    class_<te::SpicePOINTNOTONSURFACE, base<te::SpiceError>>("exceptions_SpicePOINTNOTONSURFACE");

    // POINTOFFSURFACE exception
    class_<te::SpicePOINTOFFSURFACE, base<te::SpiceError>>("exceptions_SpicePOINTOFFSURFACE");

    // POINTONZAXIS exception
    class_<te::SpicePOINTONZAXIS, base<te::SpiceError>>("exceptions_SpicePOINTONZAXIS");

    // POINTTOOSMALL exception
    class_<te::SpicePOINTTOOSMALL, base<te::SpiceError>>("exceptions_SpicePOINTTOOSMALL");

    // PTRARRAYTOOSMALL exception
    class_<te::SpicePTRARRAYTOOSMALL, base<te::SpiceError>>("exceptions_SpicePTRARRAYTOOSMALL");

    // QPARAMOUTOFRANGE exception
    class_<te::SpiceQPARAMOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceQPARAMOUTOFRANGE");

    // QUERYFAILURE exception
    class_<te::SpiceQUERYFAILURE, base<te::SpiceError>>("exceptions_SpiceQUERYFAILURE");

    // QUERYNOTPARSED exception
    class_<te::SpiceQUERYNOTPARSED, base<te::SpiceError>>("exceptions_SpiceQUERYNOTPARSED");

    // RADIIOUTOFORDER exception
    class_<te::SpiceRADIIOUTOFORDER, base<te::SpiceError>>("exceptions_SpiceRADIIOUTOFORDER");

    // RAYISZEROVECTOR exception
    class_<te::SpiceRAYISZEROVECTOR, base<te::SpiceError>>("exceptions_SpiceRAYISZEROVECTOR");

    // READFAILED exception
    class_<te::SpiceREADFAILED, base<te::SpiceError>>("exceptions_SpiceREADFAILED");

    // RECORDNOTFOUND exception
    class_<te::SpiceRECORDNOTFOUND, base<te::SpiceError>>("exceptions_SpiceRECORDNOTFOUND");

    // RECURSIONTOODEEP exception
    class_<te::SpiceRECURSIONTOODEEP, base<te::SpiceError>>("exceptions_SpiceRECURSIONTOODEEP");

    // RECURSIVELOADING exception
    class_<te::SpiceRECURSIVELOADING, base<te::SpiceError>>("exceptions_SpiceRECURSIVELOADING");

    // REFANGLEMISSING exception
    class_<te::SpiceREFANGLEMISSING, base<te::SpiceError>>("exceptions_SpiceREFANGLEMISSING");

    // REFVALNOTINTEGER exception
    class_<te::SpiceREFVALNOTINTEGER, base<te::SpiceError>>("exceptions_SpiceREFVALNOTINTEGER");

    // REFVECTORMISSING exception
    class_<te::SpiceREFVECTORMISSING, base<te::SpiceError>>("exceptions_SpiceREFVECTORMISSING");

    // REQUESTOUTOFBOUNDS exception
    class_<te::SpiceREQUESTOUTOFBOUNDS, base<te::SpiceError>>("exceptions_SpiceREQUESTOUTOFBOUNDS");

    // REQUESTOUTOFORDER exception
    class_<te::SpiceREQUESTOUTOFORDER, base<te::SpiceError>>("exceptions_SpiceREQUESTOUTOFORDER");

    // RWCONFLICT exception
    class_<te::SpiceRWCONFLICT, base<te::SpiceError>>("exceptions_SpiceRWCONFLICT");

    // SBINSUFPTRSIZE exception
    class_<te::SpiceSBINSUFPTRSIZE, base<te::SpiceError>>("exceptions_SpiceSBINSUFPTRSIZE");

    // SBTOOMANYSTRS exception
    class_<te::SpiceSBTOOMANYSTRS, base<te::SpiceError>>("exceptions_SpiceSBTOOMANYSTRS");

    // SCLKDOESNTEXIST exception
    class_<te::SpiceSCLKDOESNTEXIST, base<te::SpiceError>>("exceptions_SpiceSCLKDOESNTEXIST");

    // SCLKTRUNCATED exception
    class_<te::SpiceSCLKTRUNCATED, base<te::SpiceError>>("exceptions_SpiceSCLKTRUNCATED");

    // SEGIDTOOLONG exception
    class_<te::SpiceSEGIDTOOLONG, base<te::SpiceError>>("exceptions_SpiceSEGIDTOOLONG");

    // SEGMENTNOTFOUND exception
    class_<te::SpiceSEGMENTNOTFOUND, base<te::SpiceError>>("exceptions_SpiceSEGMENTNOTFOUND");

    // SEGMENTTABLEFULL exception
    class_<te::SpiceSEGMENTTABLEFULL, base<te::SpiceError>>("exceptions_SpiceSEGMENTTABLEFULL");

    // SEGTABLETOOSMALL exception
    class_<te::SpiceSEGTABLETOOSMALL, base<te::SpiceError>>("exceptions_SpiceSEGTABLETOOSMALL");

    // SEGTYPECONFLICT exception
    class_<te::SpiceSEGTYPECONFLICT, base<te::SpiceError>>("exceptions_SpiceSEGTYPECONFLICT");

    // SETEXCESS exception
    class_<te::SpiceSETEXCESS, base<te::SpiceError>>("exceptions_SpiceSETEXCESS");

    // SETTOOSMALL exception
    class_<te::SpiceSETTOOSMALL, base<te::SpiceError>>("exceptions_SpiceSETTOOSMALL");

    // SETUPDOESNOTEXIST exception
    class_<te::SpiceSETUPDOESNOTEXIST, base<te::SpiceError>>("exceptions_SpiceSETUPDOESNOTEXIST");

    // SHAPEMISSING exception
    class_<te::SpiceSHAPEMISSING, base<te::SpiceError>>("exceptions_SpiceSHAPEMISSING");

    // SHAPENOTSUPPORTED exception
    class_<te::SpiceSHAPENOTSUPPORTED, base<te::SpiceError>>("exceptions_SpiceSHAPENOTSUPPORTED");

    // SIZEMISMATCH exception
    class_<te::SpiceSIZEMISMATCH, base<te::SpiceError>>("exceptions_SpiceSIZEMISMATCH");

    // SIZEOUTOFRANGE exception
    class_<te::SpiceSIZEOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceSIZEOUTOFRANGE");

    // SPACETOONARROW exception
    class_<te::SpiceSPACETOONARROW, base<te::SpiceError>>("exceptions_SpiceSPACETOONARROW");

    // SPCRFLNOTCALLED exception
    class_<te::SpiceSPCRFLNOTCALLED, base<te::SpiceError>>("exceptions_SpiceSPCRFLNOTCALLED");

    // SPICEISTIRED exception
    class_<te::SpiceSPICEISTIRED, base<te::SpiceError>>("exceptions_SpiceSPICEISTIRED");

    // SPKDOESNTEXIST exception
    class_<te::SpiceSPKDOESNTEXIST, base<te::SpiceError>>("exceptions_SpiceSPKDOESNTEXIST");

    // SPKFILE exception
    class_<te::SpiceSPKFILE, base<te::SpiceError>>("exceptions_SpiceSPKFILE");

    // SPKFILETABLEFULL exception
    class_<te::SpiceSPKFILETABLEFULL, base<te::SpiceError>>("exceptions_SpiceSPKFILETABLEFULL");

    // SPKINSUFFDATA exception
    class_<te::SpiceSPKINSUFFDATA, base<te::SpiceError>>("exceptions_SpiceSPKINSUFFDATA");

    // SPKINVALIDOPTION exception
    class_<te::SpiceSPKINVALIDOPTION, base<te::SpiceError>>("exceptions_SpiceSPKINVALIDOPTION");

    // SPKNOTASUBSET exception
    class_<te::SpiceSPKNOTASUBSET, base<te::SpiceError>>("exceptions_SpiceSPKNOTASUBSET");

    // SPKRECTOOLARGE exception
    class_<te::SpiceSPKRECTOOLARGE, base<te::SpiceError>>("exceptions_SpiceSPKRECTOOLARGE");

    // SPKREFNOTSUPP exception
    class_<te::SpiceSPKREFNOTSUPP, base<te::SpiceError>>("exceptions_SpiceSPKREFNOTSUPP");

    // SPKSTRUCTUREERROR exception
    class_<te::SpiceSPKSTRUCTUREERROR, base<te::SpiceError>>("exceptions_SpiceSPKSTRUCTUREERROR");

    // SPKTYPENOTSUPP exception
    class_<te::SpiceSPKTYPENOTSUPP, base<te::SpiceError>>("exceptions_SpiceSPKTYPENOTSUPP");

    // SPKTYPENOTSUPPORTD exception
    class_<te::SpiceSPKTYPENOTSUPPORTD, base<te::SpiceError>>("exceptions_SpiceSPKTYPENOTSUPPORTD");

    // SPURIOUSKEYWORD exception
    class_<te::SpiceSPURIOUSKEYWORD, base<te::SpiceError>>("exceptions_SpiceSPURIOUSKEYWORD");

    // STFULL exception
    class_<te::SpiceSTFULL, base<te::SpiceError>>("exceptions_SpiceSTFULL");

    // STRINGTOOSHORT exception
    class_<te::SpiceSTRINGTOOSHORT, base<te::SpiceError>>("exceptions_SpiceSTRINGTOOSHORT");

    // STRINGTOOSMALL exception
    class_<te::SpiceSTRINGTOOSMALL, base<te::SpiceError>>("exceptions_SpiceSTRINGTOOSMALL");

    // STRINGTRUNCATED exception
    class_<te::SpiceSTRINGTRUNCATED, base<te::SpiceError>>("exceptions_SpiceSTRINGTRUNCATED");

    // SUBORBITAL exception
    class_<te::SpiceSUBORBITAL, base<te::SpiceError>>("exceptions_SpiceSUBORBITAL");

    // SUBPOINTNOTFOUND exception
    class_<te::SpiceSUBPOINTNOTFOUND, base<te::SpiceError>>("exceptions_SpiceSUBPOINTNOTFOUND");

    // SYNTAXERROR exception
    class_<te::SpiceSYNTAXERROR, base<te::SpiceError>>("exceptions_SpiceSYNTAXERROR");

    // SYSTEMCALLFAILED exception
    class_<te::SpiceSYSTEMCALLFAILED, base<te::SpiceError>>("exceptions_SpiceSYSTEMCALLFAILED");

    // TABLENOTLOADED exception
    class_<te::SpiceTABLENOTLOADED, base<te::SpiceError>>("exceptions_SpiceTABLENOTLOADED");

    // TIMECONFLICT exception
    class_<te::SpiceTIMECONFLICT, base<te::SpiceError>>("exceptions_SpiceTIMECONFLICT");

    // TIMEOUTOFBOUNDS exception
    class_<te::SpiceTIMEOUTOFBOUNDS, base<te::SpiceError>>("exceptions_SpiceTIMEOUTOFBOUNDS");

    // TIMESDONTMATCH exception
    class_<te::SpiceTIMESDONTMATCH, base<te::SpiceError>>("exceptions_SpiceTIMESDONTMATCH");

    // TIMESOUTOFORDER exception
    class_<te::SpiceTIMESOUTOFORDER, base<te::SpiceError>>("exceptions_SpiceTIMESOUTOFORDER");

    // TIMEZONEERROR exception
    class_<te::SpiceTIMEZONEERROR, base<te::SpiceError>>("exceptions_SpiceTIMEZONEERROR");

    // TOOFEWINPUTLINES exception
    class_<te::SpiceTOOFEWINPUTLINES, base<te::SpiceError>>("exceptions_SpiceTOOFEWINPUTLINES");

    // TOOFEWPACKETS exception
    class_<te::SpiceTOOFEWPACKETS, base<te::SpiceError>>("exceptions_SpiceTOOFEWPACKETS");

    // TOOFEWPLATES exception
    class_<te::SpiceTOOFEWPLATES, base<te::SpiceError>>("exceptions_SpiceTOOFEWPLATES");

    // TOOFEWSTATES exception
    class_<te::SpiceTOOFEWSTATES, base<te::SpiceError>>("exceptions_SpiceTOOFEWSTATES");

    // TOOFEWVERTICES exception
    class_<te::SpiceTOOFEWVERTICES, base<te::SpiceError>>("exceptions_SpiceTOOFEWVERTICES");

    // TOOFEWWINDOWS exception
    class_<te::SpiceTOOFEWWINDOWS, base<te::SpiceError>>("exceptions_SpiceTOOFEWWINDOWS");

    // TOOMANYBASEFRAMES exception
    class_<te::SpiceTOOMANYBASEFRAMES, base<te::SpiceError>>("exceptions_SpiceTOOMANYBASEFRAMES");

    // TOOMANYFIELDS exception
    class_<te::SpiceTOOMANYFIELDS, base<te::SpiceError>>("exceptions_SpiceTOOMANYFIELDS");

    // TOOMANYFILESOPEN exception
    class_<te::SpiceTOOMANYFILESOPEN, base<te::SpiceError>>("exceptions_SpiceTOOMANYFILESOPEN");

    // TOOMANYHITS exception
    class_<te::SpiceTOOMANYHITS, base<te::SpiceError>>("exceptions_SpiceTOOMANYHITS");

    // TOOMANYITERATIONS exception
    class_<te::SpiceTOOMANYITERATIONS, base<te::SpiceError>>("exceptions_SpiceTOOMANYITERATIONS");

    // TOOMANYKEYWORDS exception
    class_<te::SpiceTOOMANYKEYWORDS, base<te::SpiceError>>("exceptions_SpiceTOOMANYKEYWORDS");

    // TOOMANYPAIRS exception
    class_<te::SpiceTOOMANYPAIRS, base<te::SpiceError>>("exceptions_SpiceTOOMANYPAIRS");

    // TOOMANYPARTS exception
    class_<te::SpiceTOOMANYPARTS, base<te::SpiceError>>("exceptions_SpiceTOOMANYPARTS");

    // TOOMANYPEAS exception
    class_<te::SpiceTOOMANYPEAS, base<te::SpiceError>>("exceptions_SpiceTOOMANYPEAS");

    // TOOMANYPLATES exception
    class_<te::SpiceTOOMANYPLATES, base<te::SpiceError>>("exceptions_SpiceTOOMANYPLATES");

    // TOOMANYSURFACES exception
    class_<te::SpiceTOOMANYSURFACES, base<te::SpiceError>>("exceptions_SpiceTOOMANYSURFACES");

    // TOOMANYVERTICES exception
    class_<te::SpiceTOOMANYVERTICES, base<te::SpiceError>>("exceptions_SpiceTOOMANYVERTICES");

    // TOOMANYWATCHES exception
    class_<te::SpiceTOOMANYWATCHES, base<te::SpiceError>>("exceptions_SpiceTOOMANYWATCHES");

    // TRANSFERFILE exception
    class_<te::SpiceTRANSFERFILE, base<te::SpiceError>>("exceptions_SpiceTRANSFERFILE");

    // TRANSFERFORMAT exception
    class_<te::SpiceTRANSFERFORMAT, base<te::SpiceError>>("exceptions_SpiceTRANSFERFORMAT");

    // TWOSCLKFILENAMES exception
    class_<te::SpiceTWOSCLKFILENAMES, base<te::SpiceError>>("exceptions_SpiceTWOSCLKFILENAMES");

    // TYPEMISMATCH exception
    class_<te::SpiceTYPEMISMATCH, base<te::SpiceError>>("exceptions_SpiceTYPEMISMATCH");

    // TYPENOTSUPPORTED exception
    class_<te::SpiceTYPENOTSUPPORTED, base<te::SpiceError>>("exceptions_SpiceTYPENOTSUPPORTED");

    // TYPESMISMATCH exception
    class_<te::SpiceTYPESMISMATCH, base<te::SpiceError>>("exceptions_SpiceTYPESMISMATCH");

    // UNALLOCATEDNODE exception
    class_<te::SpiceUNALLOCATEDNODE, base<te::SpiceError>>("exceptions_SpiceUNALLOCATEDNODE");

    // UNBALANCEDGROUP exception
    class_<te::SpiceUNBALANCEDGROUP, base<te::SpiceError>>("exceptions_SpiceUNBALANCEDGROUP");

    // UNBALANCEDPAIR exception
    class_<te::SpiceUNBALANCEDPAIR, base<te::SpiceError>>("exceptions_SpiceUNBALANCEDPAIR");

    // UNDEFINEDFRAME exception
    class_<te::SpiceUNDEFINEDFRAME, base<te::SpiceError>>("exceptions_SpiceUNDEFINEDFRAME");

    // UNEQUALTIMESTEP exception
    class_<te::SpiceUNEQUALTIMESTEP, base<te::SpiceError>>("exceptions_SpiceUNEQUALTIMESTEP");

    // UNINITIALIZED exception
    class_<te::SpiceUNINITIALIZED, base<te::SpiceError>>("exceptions_SpiceUNINITIALIZED");

    // UNINITIALIZEDHASH exception
    class_<te::SpiceUNINITIALIZEDHASH, base<te::SpiceError>>("exceptions_SpiceUNINITIALIZEDHASH");

    // UNINITIALIZEDVALUE exception
    class_<te::SpiceUNINITIALIZEDVALUE, base<te::SpiceError>>("exceptions_SpiceUNINITIALIZEDVALUE");

    // UNITSMISSING exception
    class_<te::SpiceUNITSMISSING, base<te::SpiceError>>("exceptions_SpiceUNITSMISSING");

    // UNITSNOTREC exception
    class_<te::SpiceUNITSNOTREC, base<te::SpiceError>>("exceptions_SpiceUNITSNOTREC");

    // UNKNONWNTIMESYSTEM exception
    class_<te::SpiceUNKNONWNTIMESYSTEM, base<te::SpiceError>>("exceptions_SpiceUNKNONWNTIMESYSTEM");

    // UNKNOWNBFF exception
    class_<te::SpiceUNKNOWNBFF, base<te::SpiceError>>("exceptions_SpiceUNKNOWNBFF");

    // UNKNOWNCKMETA exception
    class_<te::SpiceUNKNOWNCKMETA, base<te::SpiceError>>("exceptions_SpiceUNKNOWNCKMETA");

    // UNKNOWNCOMPARE exception
    class_<te::SpiceUNKNOWNCOMPARE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNCOMPARE");

    // UNKNOWNDATATYPE exception
    class_<te::SpiceUNKNOWNDATATYPE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNDATATYPE");

    // UNKNOWNFILARC exception
    class_<te::SpiceUNKNOWNFILARC, base<te::SpiceError>>("exceptions_SpiceUNKNOWNFILARC");

    // UNKNOWNFRAME exception
    class_<te::SpiceUNKNOWNFRAME, base<te::SpiceError>>("exceptions_SpiceUNKNOWNFRAME");

    // UNKNOWNFRAMESPEC exception
    class_<te::SpiceUNKNOWNFRAMESPEC, base<te::SpiceError>>("exceptions_SpiceUNKNOWNFRAMESPEC");

    // UNKNOWNFRAMETYPE exception
    class_<te::SpiceUNKNOWNFRAMETYPE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNFRAMETYPE");

    // UNKNOWNID exception
    class_<te::SpiceUNKNOWNID, base<te::SpiceError>>("exceptions_SpiceUNKNOWNID");

    // UNKNOWNINCLUSION exception
    class_<te::SpiceUNKNOWNINCLUSION, base<te::SpiceError>>("exceptions_SpiceUNKNOWNINCLUSION");

    // UNKNOWNINDEXTYPE exception
    class_<te::SpiceUNKNOWNINDEXTYPE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNINDEXTYPE");

    // UNKNOWNKERNELTYPE exception
    class_<te::SpiceUNKNOWNKERNELTYPE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNKERNELTYPE");

    // UNKNOWNKEY exception
    class_<te::SpiceUNKNOWNKEY, base<te::SpiceError>>("exceptions_SpiceUNKNOWNKEY");

    // UNKNOWNMETAITEM exception
    class_<te::SpiceUNKNOWNMETAITEM, base<te::SpiceError>>("exceptions_SpiceUNKNOWNMETAITEM");

    // UNKNOWNMODE exception
    class_<te::SpiceUNKNOWNMODE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNMODE");

    // UNKNOWNOP exception
    class_<te::SpiceUNKNOWNOP, base<te::SpiceError>>("exceptions_SpiceUNKNOWNOP");

    // UNKNOWNPCKTYPE exception
    class_<te::SpiceUNKNOWNPCKTYPE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNPCKTYPE");

    // UNKNOWNREFDIR exception
    class_<te::SpiceUNKNOWNREFDIR, base<te::SpiceError>>("exceptions_SpiceUNKNOWNREFDIR");

    // UNKNOWNSPKTYPE exception
    class_<te::SpiceUNKNOWNSPKTYPE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNSPKTYPE");

    // UNKNOWNSYSTEM exception
    class_<te::SpiceUNKNOWNSYSTEM, base<te::SpiceError>>("exceptions_SpiceUNKNOWNSYSTEM");

    // UNKNOWNTYPE exception
    class_<te::SpiceUNKNOWNTYPE, base<te::SpiceError>>("exceptions_SpiceUNKNOWNTYPE");

    // UNKNOWNUNITS exception
    class_<te::SpiceUNKNOWNUNITS, base<te::SpiceError>>("exceptions_SpiceUNKNOWNUNITS");

    // UNMATCHENDPTS exception
    class_<te::SpiceUNMATCHENDPTS, base<te::SpiceError>>("exceptions_SpiceUNMATCHENDPTS");

    // UNNATURALACT exception
    class_<te::SpiceUNNATURALACT, base<te::SpiceError>>("exceptions_SpiceUNNATURALACT");

    // UNNATURALRELATION exception
    class_<te::SpiceUNNATURALRELATION, base<te::SpiceError>>("exceptions_SpiceUNNATURALRELATION");

    // UNORDEREDREFS exception
    class_<te::SpiceUNORDEREDREFS, base<te::SpiceError>>("exceptions_SpiceUNORDEREDREFS");

    // UNORDEREDTIMES exception
    class_<te::SpiceUNORDEREDTIMES, base<te::SpiceError>>("exceptions_SpiceUNORDEREDTIMES");

    // UNPARSEDQUERY exception
    class_<te::SpiceUNPARSEDQUERY, base<te::SpiceError>>("exceptions_SpiceUNPARSEDQUERY");

    // UNPARSEDTIME exception
    class_<te::SpiceUNPARSEDTIME, base<te::SpiceError>>("exceptions_SpiceUNPARSEDTIME");

    // UNRECOGNAPPFLAG exception
    class_<te::SpiceUNRECOGNAPPFLAG, base<te::SpiceError>>("exceptions_SpiceUNRECOGNAPPFLAG");

    // UNRECOGNDATATYPE exception
    class_<te::SpiceUNRECOGNDATATYPE, base<te::SpiceError>>("exceptions_SpiceUNRECOGNDATATYPE");

    // UNRECOGNDELIMITER exception
    class_<te::SpiceUNRECOGNDELIMITER, base<te::SpiceError>>("exceptions_SpiceUNRECOGNDELIMITER");

    // UNRECOGNIZABLEFILE exception
    class_<te::SpiceUNRECOGNIZABLEFILE, base<te::SpiceError>>("exceptions_SpiceUNRECOGNIZABLEFILE");

    // UNRECOGNIZEDACTION exception
    class_<te::SpiceUNRECOGNIZEDACTION, base<te::SpiceError>>("exceptions_SpiceUNRECOGNIZEDACTION");

    // UNRECOGNIZEDFORMAT exception
    class_<te::SpiceUNRECOGNIZEDFORMAT, base<te::SpiceError>>("exceptions_SpiceUNRECOGNIZEDFORMAT");

    // UNRECOGNIZEDFRAME exception
    class_<te::SpiceUNRECOGNIZEDFRAME, base<te::SpiceError>>("exceptions_SpiceUNRECOGNIZEDFRAME");

    // UNRECOGNIZEDTYPE exception
    class_<te::SpiceUNRECOGNIZEDTYPE, base<te::SpiceError>>("exceptions_SpiceUNRECOGNIZEDTYPE");

    // UNRECOGNPRECTYPE exception
    class_<te::SpiceUNRECOGNPRECTYPE, base<te::SpiceError>>("exceptions_SpiceUNRECOGNPRECTYPE");

    // UNRESOLVEDNAMES exception
    class_<te::SpiceUNRESOLVEDNAMES, base<te::SpiceError>>("exceptions_SpiceUNRESOLVEDNAMES");

    // UNRESOLVEDTIMES exception
    class_<te::SpiceUNRESOLVEDTIMES, base<te::SpiceError>>("exceptions_SpiceUNRESOLVEDTIMES");

    // UNSUPPORTEDARCH exception
    class_<te::SpiceUNSUPPORTEDARCH, base<te::SpiceError>>("exceptions_SpiceUNSUPPORTEDARCH");

    // UNSUPPORTEDBFF exception
    class_<te::SpiceUNSUPPORTEDBFF, base<te::SpiceError>>("exceptions_SpiceUNSUPPORTEDBFF");

    // UNSUPPORTEDMETHOD exception
    class_<te::SpiceUNSUPPORTEDMETHOD, base<te::SpiceError>>("exceptions_SpiceUNSUPPORTEDMETHOD");

    // UNSUPPORTEDSPEC exception
    class_<te::SpiceUNSUPPORTEDSPEC, base<te::SpiceError>>("exceptions_SpiceUNSUPPORTEDSPEC");

    // UNTITLEDHELP exception
    class_<te::SpiceUNTITLEDHELP, base<te::SpiceError>>("exceptions_SpiceUNTITLEDHELP");

    // UPDATEPENDING exception
    class_<te::SpiceUPDATEPENDING, base<te::SpiceError>>("exceptions_SpiceUPDATEPENDING");

    // USAGEERROR exception
    class_<te::SpiceUSAGEERROR, base<te::SpiceError>>("exceptions_SpiceUSAGEERROR");

    // UTFULL exception
    class_<te::SpiceUTFULL, base<te::SpiceError>>("exceptions_SpiceUTFULL");

    // VALUEOUTOFRANGE exception
    class_<te::SpiceVALUEOUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceVALUEOUTOFRANGE");

    // VALUETABLEFULL exception
    class_<te::SpiceVALUETABLEFULL, base<te::SpiceError>>("exceptions_SpiceVALUETABLEFULL");

    // VARIABLENOTFOUND exception
    class_<te::SpiceVARIABLENOTFOUND, base<te::SpiceError>>("exceptions_SpiceVARIABLENOTFOUND");

    // VARNAMETOOLONG exception
    class_<te::SpiceVARNAMETOOLONG, base<te::SpiceError>>("exceptions_SpiceVARNAMETOOLONG");

    // VECTORTOOBIG exception
    class_<te::SpiceVECTORTOOBIG, base<te::SpiceError>>("exceptions_SpiceVECTORTOOBIG");

    // VERSIONMISMATCH exception
    class_<te::SpiceVERSIONMISMATCH, base<te::SpiceError>>("exceptions_SpiceVERSIONMISMATCH");

    // VERTEXNOTINGRID exception
    class_<te::SpiceVERTEXNOTINGRID, base<te::SpiceError>>("exceptions_SpiceVERTEXNOTINGRID");

    // VOXELGRIDTOOBIG exception
    class_<te::SpiceVOXELGRIDTOOBIG, base<te::SpiceError>>("exceptions_SpiceVOXELGRIDTOOBIG");

    // WIDTHTOOSMALL exception
    class_<te::SpiceWIDTHTOOSMALL, base<te::SpiceError>>("exceptions_SpiceWIDTHTOOSMALL");

    // WINDOWEXCESS exception
    class_<te::SpiceWINDOWEXCESS, base<te::SpiceError>>("exceptions_SpiceWINDOWEXCESS");

    // WINDOWSTOOSMALL exception
    class_<te::SpiceWINDOWSTOOSMALL, base<te::SpiceError>>("exceptions_SpiceWINDOWSTOOSMALL");

    // WINDOWTOOSMALL exception
    class_<te::SpiceWINDOWTOOSMALL, base<te::SpiceError>>("exceptions_SpiceWINDOWTOOSMALL");

    // WORKSPACETOOSMALL exception
    class_<te::SpiceWORKSPACETOOSMALL, base<te::SpiceError>>("exceptions_SpiceWORKSPACETOOSMALL");

    // WRITEERROR exception
    class_<te::SpiceWRITEERROR, base<te::SpiceError>>("exceptions_SpiceWRITEERROR");

    // WRITEFAILED exception
    class_<te::SpiceWRITEFAILED, base<te::SpiceError>>("exceptions_SpiceWRITEFAILED");

    // WRONGARCHITECTURE exception
    class_<te::SpiceWRONGARCHITECTURE, base<te::SpiceError>>("exceptions_SpiceWRONGARCHITECTURE");

    // WRONGCKTYPE exception
    class_<te::SpiceWRONGCKTYPE, base<te::SpiceError>>("exceptions_SpiceWRONGCKTYPE");

    // WRONGCONIC exception
    class_<te::SpiceWRONGCONIC, base<te::SpiceError>>("exceptions_SpiceWRONGCONIC");

    // WRONGDATATYPE exception
    class_<te::SpiceWRONGDATATYPE, base<te::SpiceError>>("exceptions_SpiceWRONGDATATYPE");

    // WRONGSEGMENT exception
    class_<te::SpiceWRONGSEGMENT, base<te::SpiceError>>("exceptions_SpiceWRONGSEGMENT");

    // WRONGSPKTYPE exception
    class_<te::SpiceWRONGSPKTYPE, base<te::SpiceError>>("exceptions_SpiceWRONGSPKTYPE");

    // YEAROUTOFRANGE exception
    class_<te::SpiceYEAROUTOFRANGE, base<te::SpiceError>>("exceptions_SpiceYEAROUTOFRANGE");

    // ZEROBORESIGHT exception
    class_<te::SpiceZEROBORESIGHT, base<te::SpiceError>>("exceptions_SpiceZEROBORESIGHT");

    // ZEROBOUNDSEXTENT exception
    class_<te::SpiceZEROBOUNDSEXTENT, base<te::SpiceError>>("exceptions_SpiceZEROBOUNDSEXTENT");

    // ZEROFRAMEID exception
    class_<te::SpiceZEROFRAMEID, base<te::SpiceError>>("exceptions_SpiceZEROFRAMEID");

    // ZEROLENGTHCOLUMN exception
    class_<te::SpiceZEROLENGTHCOLUMN, base<te::SpiceError>>("exceptions_SpiceZEROLENGTHCOLUMN");

    // ZEROPOSITION exception
    class_<te::SpiceZEROPOSITION, base<te::SpiceError>>("exceptions_SpiceZEROPOSITION");

    // ZEROQUATERNION exception
    class_<te::SpiceZEROQUATERNION, base<te::SpiceError>>("exceptions_SpiceZEROQUATERNION");

    // ZEROSTEP exception
    class_<te::SpiceZEROSTEP, base<te::SpiceError>>("exceptions_SpiceZEROSTEP");

    // ZEROVECTOR exception
    class_<te::SpiceZEROVECTOR, base<te::SpiceError>>("exceptions_SpiceZEROVECTOR");

    // ZEROVELOCITY exception
    class_<te::SpiceZEROVELOCITY, base<te::SpiceError>>("exceptions_SpiceZEROVELOCITY");

    // ZZHOLDDGETFAILED exception
    class_<te::SpiceZZHOLDDGETFAILED, base<te::SpiceError>>("exceptions_SpiceZZHOLDDGETFAILED");
}

#endif

