/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPICE_EXCEPTIONS_H
#define TUDAT_SPICE_EXCEPTIONS_H

#include <iostream>
#include <string>

namespace tudat
{

namespace exceptions
{

/// @brief Generic exception class for SPICE.
///
/// This class is a generic base class for all SPICE exceptions.
/// It inherits from std::runtime_error.
class SpiceError : public std::runtime_error
{
private:
public:
    /// @brief Default constructor for SpiceError.
    /// @param errorMessage Error message to be displayed.
    SpiceError( const std::string& shortMessage,
                const std::string& explanation,
                const std::string& longMessage,
                const std::string& traceback ):
        std::runtime_error( shortMessage + "\n" + explanation + "\n" + longMessage + "\n" + "Traceback: " + traceback )
    { }
    ~SpiceError( ) { }
};

class SpiceADDRESSOUTOFBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceAGENTLISTOVERFLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceALLGONE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceAMBIGTEMPL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceARRAYSIZEMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceARRAYTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceAVALOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceAXISUNDERFLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADACTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADADDRESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADANGLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADANGLEUNITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADANGRATEERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADANGULARRATE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADANGULARRATEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADARCHITECTURE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADARRAYSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADATTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADATTRIBUTE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADATTRIBUTES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADAUVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADAVFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADAVFRAMEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADAXIS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADAXISLENGTH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADAXISNUMBERS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADBLOCKSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADBODYID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADBORESIGHTSPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADBOUNDARY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCATALOGFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCENTERNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCHECKFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCKTYPESPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOARSEVOXSCALE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOLUMDECL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOLUMNCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOLUMNDECL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOMMENTAREA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOMPNUMBER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOORDBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCOORDSYS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADCURVETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDAFTRANSFERFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDASCOMMENTAREA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDASDIRECTORY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDASFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDASTRANSFERFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDATALINE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDATAORDERTOKEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDATATYPEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDEFAULTVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDESCRTIMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDIMENSION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDIMENSIONS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDIRECTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDOUBLEPRECISION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADDOWNSAMPLINGTOL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADECCENTRICITY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADENDPOINTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADEULERANGLEUNITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFILEFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFINEVOXELSCALE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFORMATSPECIFIER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFRAMECLASS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFRAMECOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFRAMESPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFROMTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFROMTIMESYSTEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADFROMTIMETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADGEOMETRY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADGM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADHARDSPACE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADHERMITDEGREE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINDEX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINITSTATE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINPUTDATALINE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINPUTETTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINPUTTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINPUTUTCTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINSTRUMENTID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADINTEGER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADKERNELTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADKERNELVARTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLAGRANGEDEGREE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLATITUDEBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLATITUDERANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLATUSRECTUM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLEAPSECONDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLIMBLOCUSMIX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLINEPERRECCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLISTFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADLONGITUDERANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADMATRIX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADMEANMOTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADMECCENTRICITY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADMETHODSYNTAX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADMIDNIGHTTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADMSEMIMAJOR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADMSOPQUATERNION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADNOFDIGITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADNOFSTATES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADNUMBEROFPOINTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOBJECTID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOBJECTNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOFFSETANGLES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOFFSETANGUNITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOFFSETAXESFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOFFSETAXISXYZ : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADORBITALPERIOD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOUTPUTSPKTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADOUTPUTTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPARTNUMBER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPCKVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPECCENTRICITY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPERIAPSEVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPICTURE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPLATECOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPODLOCATION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPRECVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADPRIORITYSPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADQUATSIGN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADQUATTHRESHOLD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADRADIUS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADRADIUSCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADRATEFRAMEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADRATETHRESHOLD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADRECORDCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADREFVECTORSPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADROTATIONAXISXYZ : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADROTATIONSORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADROTATIONTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADROTAXESFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADROWCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSCID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSEMIAXIS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSEMILATUS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSHAPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSOLDAY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSOLINDEX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSOLTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSOURCERADIUS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSPICEQUATERNION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSTARINDEX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSTARTTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSTDIONAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSTOPTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSUBSTR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSUBSTRINGBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADSURFACEMAP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTABLEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTERMLOCUSMIX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMEBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMECASE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMECOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMEFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMEITEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMEOFFSET : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMESPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMESTRING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTIMETYPEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTLECOVERAGEPAD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTLEPADS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTOTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTOTIMESYSTEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTOTIMETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADTYPESHAPECOMBO : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADVARASSIGN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADVARIABLESIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADVARIABLETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADVARNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADVECTOR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADVERTEXCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADVERTEXINDEX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBADWINDOWSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBARRAYTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBARYCENTEREPHEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBARYCENTERIDCODE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBEFOREBEGSTR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKCOMMANDLINE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKFILETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKINPUTFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKINPUTTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKNAMEASSIGNED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKOUTPTFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKSCLKSTRING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLANKTIMEFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBLOCKSNOTEVEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBODIESNOTDISTINCT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBODYANDCENTERSAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBOGUSENTRY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBORESIGHTMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBOUNDARYMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBOUNDARYTOOBIG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBOUNDSDISAGREE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBOUNDSOUTOFORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBUFFEROVERFLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBUFFERSIZESMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBUFFERTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBUG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceBUGWRITEFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCALLCKBSSFIRST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCALLEDOUTOFORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCALLZZDSKBSSFIRST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCANNOTFINDGRP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCANNOTGETPACKET : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCANNOTMAKEFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCANNOTPICKFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCANTFINDFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCANTGETROTATIONTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCANTUSEPERIAPEPOCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCBNOSUCHSTR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCELLARRAYTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCELLTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCKBOGUSENTRY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCKDOESNTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCKFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCKNONEXISTREC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCKTOOMANYFILES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCKUNKNOWNDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCKWRONGDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCMDERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCMDPARSEERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOARSEGRIDOVERFLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOLDESCTABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOLUMNTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOMMANDTOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOMMENTTOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOMMFILENOTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOMPETINGEPOCHSPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOMPETINGFRAMESPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOORDSYSNOTREC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOUNTMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOUNTTOOLARGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCOVERAGEGAP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceCROSSANGLEMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFBADCRECLEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFBEGGTEND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFCRNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFDPWRITEFAIL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFFRNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFFTFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFILLEGWRITE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFINVALIDACCESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFINVALIDPARAMS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNEGADDR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNEWCONFLICT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNOIFNMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNONAMEMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNORESV : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNOSEARCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNOSUCHADDR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNOSUCHFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNOSUCHHANDLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNOSUCHUNIT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFNOWRITE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFOVERFLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFREADFAIL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDAFWRITEFAIL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASFILEREADFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASFILEWRITEFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASFTFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASINVALIDACCESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASINVALIDCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASINVALIDTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASNOSUCHADDRESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASNOSUCHFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASNOSUCHHANDLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASNOSUCHUNIT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASNOTEMPTY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASREADFAIL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDASWRITEFAIL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDATAITEMLIMITEXCEEDED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDATAREADFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDATAWIDTHERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDATEEXPECTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDECODINGERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDEGENERATECASE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDEGENERATEINTERVAL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDEGENERATESURFACE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDEGREEOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDEPENDENTVECTORS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDEVICENAMETOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDIFFLINETOOLARGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDIFFLINETOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDIMENSIONTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDISARRAY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDISORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDIVIDEBYZERO : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDSKBOGUSENTRY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDSKDATANOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDSKTOOMANYFILES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDTOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDUBIOUSMETHOD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceDUPLICATETIMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceECCOUTOFBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceECCOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKCOLATTRTABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKCOLNUMMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKFILETABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKMISSINGCOLUMN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKNOSEGMENTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKSEGTABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEKTABLELISTFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceELEMENTSTOOSHORT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEMPTYINPUTFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEMPTYSEGMENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceENDOFFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceENDPOINTSMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceERROREXIT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEVECOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEVENHERMITDEGREE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEVILBOGUSENTRY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceEXTERNALOPEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFACENOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFAKESCLKEXISTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILARCHMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILARCMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEALREADYEXISTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILECURRENTLYOPEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEDELETEFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEDOESNOTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEEXISTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEISNOTSPK : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILENAMETOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILENOTCONNECTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILENOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILENOTOPEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEOPENCONFLICT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEOPENERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEOPENFAIL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEOPENFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEREADERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEREADFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILETABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILETRUNCATED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFILEWRITEFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFIRSTRECORDMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFKDOESNTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFMTITEMLIMITEXCEEDED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFORMATDATAMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFORMATDOESNTAPPLY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFORMATERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFORMATNOTAPPLICABLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFORMATSTRINGTOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFOVTOOWIDE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMEDATANOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMEDEFERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMEIDNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMEINFONOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMEMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMENAMENOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMENOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFRAMENOTRECOGNIZED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFTFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceFTPXFERERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceGRIDTOOLARGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceHANDLENOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceHASHISFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceHLULOCKFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceIDCODENOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceIDSTRINGTOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceILLEGALCHARACTER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceILLEGALOPTIONNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceILLEGSHIFTDIR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceILLEGTEMPL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceIMMUTABLEVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceIMPROPERFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceIMPROPEROPEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINACTIVEOBJECT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCOMPATIBLEEOL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCOMPATIBLENUMREF : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCOMPATIBLESCALE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCOMPATIBLEUNITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCOMPLETEFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCONSISTCENTERID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCONSISTENTTIMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCONSISTFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCONSISTSTARTTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCONSISTSTOPTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINCORRECTUSAGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINDEFINITELOCALSECOND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINDEXOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINDEXTOOLARGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINDICESOUTOFORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINPUTDOESNOTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINPUTFILENOTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINPUTOUTOFBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINPUTSTOOLARGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINQUIREERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINQUIREFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINSIDEBODY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINSUFFICIENTANGLES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINSUFFICIENTDATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINSUFFLEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINSUFPTRSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINTERVALSTARTNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINTINDEXTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINTLENNOTPOS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINTOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDACCESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDACTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDADD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDADDRESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDANGLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDARCHTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDARGUMENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDAXIS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDAXISLENGTH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDCARDINALITY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDCASE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDCOLUMN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDCONSTSTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDATACOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDEGREE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDESCRTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDIMENSION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDIRECTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDDIVISOR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDELLIPSE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDENDPNTSPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDENDPTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDEPOCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDFILETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDFIXREF : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDFOV : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDFRAMEDEF : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDGEOMETRY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDHANDLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDINDEX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDINTEGER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDLIMBTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDLISTITEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDLOCUS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDLONEXTENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDMETADATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDMETHOD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDMSGTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDNODE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDNUMBEROFINTERVALS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDNUMBEROFRECORDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDNUMINT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDNUMINTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDNUMREC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDOCCTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDOPERATION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDOPTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDPLANE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDRADII : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDRADIUS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDREFFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDREFVAL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDROLLSTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSCALE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSCLKRATE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSCLKSTRING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSCLKTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSEARCHSTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSELECTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSHADOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSHAPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSHAPECOMBO : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSTARTTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSTATE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSTEPSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSUBLIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDSUBTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTABLENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTABLESIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTARGET : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTERMTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTEXT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTIMEFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTIMESTRING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTLEORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTOLERANCE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVALIDVERTEX : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceINVERSTARTSTOPTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceIRFNOTREC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceITEMNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceITEMNOTRECOGNIZED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceITERATIONEXCEEDED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceKERNELNOTLOADED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceKERNELPOOLFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceKERNELVARNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceKERVARSETOVERFLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceKERVARTOOBIG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceKEYWORDNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceLBCORRUPTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceLBLINETOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceLBNOSUCHLINE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceLBTOOMANYLINES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceLOWERBOUNDTOOLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceLSKDOESNTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMALFORMEDSEGMENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMARKERNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMESSAGETOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISMATCHFROMTIMETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISMATCHOUTPUTFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISMATCHTOTIMETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGARGUMENTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGCENTER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGCOLSTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGCOORDBOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGCOORDSYS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGDATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGDATACLASS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGDATAORDERTK : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGEOT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGEPOCHTOKEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGFRAMEVAR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGGEOCONSTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGHEIGHTREF : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGHSCALE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGKPV : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGLEFTCOR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGLEFTRTFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGNCAPFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGNCOLS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGNROWS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGPLATETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGROWMAJFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGROWSTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGSCAPFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGSURFACE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGTIMEINFO : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGTLEKEYWORD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGTOPCOR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGTOPDOWNFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGVOXELSCALE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceMISSINGWRAPFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNAMENOTUNIQUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNAMESNOTRESOLVED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNAMETABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNARATESFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNEGATIVETOL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOACCEPTABLEDATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOANGULARRATEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOARRAYSTARTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOATTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOAVDATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOBODYID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCANDOSPKSPCKS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCENTERIDORNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCKSEGMENTTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCLASS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCOMMENTSFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCONVERG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCONVERGENCE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOCURRENTARRAY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNODATAORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNODATATYPEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNODELIMCHARACTER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNODETOOFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNODSKSEGMENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNODSKSEGMENTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOENVVARIABLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOEULERANGLEUNITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFILENAMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFILES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFILESPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFRAMECONNECT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFRAMEDATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFRAMEINFO : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFRAMENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFRAMESKERNELNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFREELOGICALUNIT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFREENODES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFROMTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOFROMTIMESYSTEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOHEADNODE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOINPUTDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOINPUTFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOINSTRUMENTID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOINTERVAL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOKERNELLOADED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOLANDINGTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOLEAPSECONDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOLINESPERRECCOUNT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOLISTFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOLOADEDDSKFILES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOLOADEDFILES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOLSKFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOMOREROOM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONCONICMOTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONDISTINCTPAIR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONEMPTYENTRY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONEMPTYTREE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONEXISTELEMENTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONINTEGERFIELD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONNUMERICSTRING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPOSBUFLENGTH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPOSITIVEAXIS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPOSITIVEMASS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPOSITIVERADIUS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPOSITIVESCALE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPOSITIVEVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPOSPACKETSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPRINTABLECHARS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPRINTINGCHAR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONPRINTINGCHARS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNONUNITNORMAL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOOBJECTIDORNAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOOFFSETANGLEAXES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOOFFSETANGLEUNITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOOUTPUTFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOOUTPUTSPKTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOPARTITION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOPICTURE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOPOLYNOMIALDEGREE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOPRECESSIONTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOPRODUCERID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOROTATIONORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSCID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSCLKFILENAMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSECONDLINE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSEGMENTSFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSEPARATION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSLKFILENAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSOLMARKER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSPACECRAFTID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSTARTTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSTOPTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSUCHFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSUCHHANDLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSUCHSYMBOL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOSUNGM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTABINARYKERNEL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTACKFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTADAFFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTADASFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTADPNUMBER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTANDPNUMBER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTANINTEGER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTANINTEGERNUMBER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTANINTNUMBER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTAPCKFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTAROTATION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTATEXTFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTATRANSFERFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTCOMPUTABLE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTDIMENSIONALLYEQUIV : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTDISJOINT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTDISTINCT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTENOUGHPEAS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTIMETYPEFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTINDEXED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTINITIALIZED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTINPART : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTLEDATAFOROBJECT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTLEGALCB : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTOTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTOTIMESYSTEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTRANSLATION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTRECOGNIZED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTSEMCHECKED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTSUPPORTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTTWOFIELDSCLK : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTTWOMODULI : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOTTWOOFFSETS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNOUNITSPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNUMBEREXPECTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNUMCOEFFSNOTPOS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNUMERICOVERFLOW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNUMPACKETSNOTPOS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNUMPARTSUNEQUAL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceNUMSTATESNOTPOS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOBJECTLISTFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOBJECTSTOOCLOSE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceORBITDECAY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOUTOFPLACEDELIMITER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOUTOFROOM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOUTPUTFILEEXISTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOUTPUTISNOTSPK : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOUTPUTTOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceOUTPUTTOOSHORT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePARSERNOTREADY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePARTIALFRAMESPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePASTENDSTR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePATHMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePATHTOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePCKDOESNTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePCKFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePCKFILETABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePCKKRECTOOLARGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePLATELISTTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTEROUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTERSETTOOBIG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTERTABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTNOTINSEGMENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTNOTONSURFACE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTOFFSURFACE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTONZAXIS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePOINTTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpicePTRARRAYTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceQPARAMOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceQUERYFAILURE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceQUERYNOTPARSED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceRADIIOUTOFORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceRAYISZEROVECTOR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceREADFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceRECORDNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceRECURSIONTOODEEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceRECURSIVELOADING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceREFANGLEMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceREFVALNOTINTEGER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceREFVECTORMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceREQUESTOUTOFBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceREQUESTOUTOFORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceRWCONFLICT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSBINSUFPTRSIZE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSBTOOMANYSTRS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSCLKDOESNTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSCLKTRUNCATED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSEGIDTOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSEGMENTNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSEGMENTTABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSEGTABLETOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSEGTYPECONFLICT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSETEXCESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSETTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSETUPDOESNOTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSHAPEMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSHAPENOTSUPPORTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSIZEMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSIZEOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPACETOONARROW : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPCRFLNOTCALLED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPICEISTIRED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKDOESNTEXIST : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKFILETABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKINSUFFDATA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKINVALIDOPTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKNOTASUBSET : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKRECTOOLARGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKREFNOTSUPP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKSTRUCTUREERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKTYPENOTSUPP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPKTYPENOTSUPPORTD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSPURIOUSKEYWORD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSTFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSTRINGTOOSHORT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSTRINGTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSTRINGTRUNCATED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSUBORBITAL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSUBPOINTNOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSYNTAXERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceSYSTEMCALLFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTABLENOTLOADED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTIMECONFLICT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTIMEOUTOFBOUNDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTIMESDONTMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTIMESOUTOFORDER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTIMEZONEERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOFEWINPUTLINES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOFEWPACKETS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOFEWPLATES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOFEWSTATES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOFEWVERTICES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOFEWWINDOWS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYBASEFRAMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYFIELDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYFILESOPEN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYHITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYITERATIONS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYKEYWORDS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYPAIRS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYPARTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYPEAS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYPLATES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYSURFACES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYVERTICES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTOOMANYWATCHES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTRANSFERFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTRANSFERFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTWOSCLKFILENAMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTYPEMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTYPENOTSUPPORTED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceTYPESMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNALLOCATEDNODE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNBALANCEDGROUP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNBALANCEDPAIR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNDEFINEDFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNEQUALTIMESTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNINITIALIZED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNINITIALIZEDHASH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNINITIALIZEDVALUE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNITSMISSING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNITSNOTREC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNONWNTIMESYSTEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNBFF : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNCKMETA : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNCOMPARE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNFILARC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNFRAMESPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNFRAMETYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNINCLUSION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNINDEXTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNKERNELTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNKEY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNMETAITEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNMODE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNOP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNPCKTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNREFDIR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNSPKTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNSYSTEM : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNKNOWNUNITS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNMATCHENDPTS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNNATURALACT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNNATURALRELATION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNORDEREDREFS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNORDEREDTIMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNPARSEDQUERY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNPARSEDTIME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNAPPFLAG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNDELIMITER : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNIZABLEFILE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNIZEDACTION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNIZEDFORMAT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNIZEDFRAME : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNIZEDTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRECOGNPRECTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRESOLVEDNAMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNRESOLVEDTIMES : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNSUPPORTEDARCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNSUPPORTEDBFF : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNSUPPORTEDMETHOD : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNSUPPORTEDSPEC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUNTITLEDHELP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUPDATEPENDING : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUSAGEERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceUTFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVALUEOUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVALUETABLEFULL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVARIABLENOTFOUND : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVARNAMETOOLONG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVECTORTOOBIG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVERSIONMISMATCH : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVERTEXNOTINGRID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceVOXELGRIDTOOBIG : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWIDTHTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWINDOWEXCESS : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWINDOWSTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWINDOWTOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWORKSPACETOOSMALL : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRITEERROR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRITEFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRONGARCHITECTURE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRONGCKTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRONGCONIC : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRONGDATATYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRONGSEGMENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceWRONGSPKTYPE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceYEAROUTOFRANGE : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROBORESIGHT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROBOUNDSEXTENT : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROFRAMEID : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROLENGTHCOLUMN : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROPOSITION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROQUATERNION : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROSTEP : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROVECTOR : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZEROVELOCITY : public SpiceError
{
public:
    using SpiceError::SpiceError;
};
class SpiceZZHOLDDGETFAILED : public SpiceError
{
public:
    using SpiceError::SpiceError;
};

}  // namespace exceptions

}  // namespace tudat

#endif  // TUDAT_SPICE_EXCEPTIONS_H