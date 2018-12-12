#include "Interpret.h"

Interpret::Interpret(void)
{
  setSourceFileName("Interpret()");
  setStandardSettings();
  allocateHitArray();
  allocateHitBufferArray();
  allocateTriggerStatusCounterArray();
  allocateEventStatusCounterArray();
  allocateTdcValueArray();
  allocateTdcTriggerDistanceArray();
  allocateServiceRecordCounterArray();
  reset();
}

Interpret::~Interpret(void)
{
  debug("~Interpret()");
  deleteHitArray();
  deleteHitBufferArray();
  deleteTriggerStatusCounterArray();
  deleteEventStatusCounterArray();
  deleteTdcValueArray();
  deleteTdcTriggerDistanceArray();
  deleteServiceRecordCounterArray();
}

void Interpret::setStandardSettings()
{
  info("setStandardSettings()");
  _hitInfoSize = 1000000;
  _hitInfo = 0;
  _hitIndex = 0;
  _startDebugEvent = 0;
  _stopDebugEvent = 0;
  _NbCID = 16;
  _maxTot = 13;
  _fEI4B = true;
  _metaDataSet = false;
  _debugEvents = false;
  _lastMetaIndexNotSet = 0;
  _lastWordIndexSet = 0;
  _metaEventIndexLength = 0;
  _metaEventIndex = 0;
  _startWordIndex = 0;
  _createMetaDataWordIndex = false;
  _createEmptyEventHits = false;
  _isMetaTableV2 = true;
  _alignAtTriggerNumber = false;
  _TriggerDataFormat = TRIGGER_FROMAT_TRIGGER_NUMBER;
  _haveTdcTriggerTimeStamp = false;
  _haveTdcTriggerDistance = false;
  _maxTdcDelay = 255;
  _alignAtTdcWord = false;
  _dataWordIndex = 0;
  _maxTriggerNumber = (2 ^ 31) - 1;
}

bool Interpret::interpretRawData(unsigned int* pDataWords, const unsigned int& pNdataWords)
{
  if (Basis::debugSet()) {
    std::stringstream tDebug;
    tDebug << "interpretRawData with " << pNdataWords << " words at total word " << _nDataWords;
    debug(tDebug.str());
  }
  _hitIndex = 0;
  _actualMetaWordIndex = 0;

  int tActualCol1 = 0;  // column position of the first hit in the actual data record
  int tActualRow1 = 0;  // row position of the first hit in the actual data record
  int tActualTot1 = -1;  // tot value of the first hit in the actual data record
  int tActualCol2 = 0;  // column position of the second hit in the actual data record
  int tActualRow2 = 0;  // row position of the second hit in the actual data record
  int tActualTot2 = -1;  // tot value of the second hit in the actual data record

  for (unsigned int iWord = 0; iWord < pNdataWords; ++iWord) {  // loop over the SRAM words
    if (_debugEvents) {
      if (_nEvents >= _startDebugEvent && _nEvents <= _stopDebugEvent)
        setDebugOutput();
      else
        setDebugOutput(false);
      setInfoOutput(false);
      setWarningOutput(false);  // TODO: do not always set to false
    }

    _nDataWords++;
    unsigned int tActualWord = pDataWords[iWord];  // take the actual SRAM word
    tActualTot1 = -1;  // TOT1 value stays negative if it can not be set properly in getHitsfromDataRecord()
    tActualTot2 = -1;  // TOT2 value stays negative if it can not be set properly in getHitsfromDataRecord()
    if (getTimefromDataHeader(tActualWord, tActualLVL1ID, tActualBCID)) {  // data word is data header if true is returned
      _nDataHeaders++;  // increase global data header counter
      if (tNdataHeader >= _NbCID) {  // maximum event window is reached (tNdataHeader > BCIDs, mostly tNdataHeader > 15)
        if (_alignAtTriggerNumber) {  // do not create new event
          addEventStatus(__TRUNC_EVENT);
          if (Basis::warningSet())
            warning("interpretRawData: " + IntToStr(_nDataWords) + " DH WORD " + IntToStr(tActualWord) + " - " + IntToStr(tNdataHeader) + ">" + IntToStr(_NbCID - 1) + " at event " + LongIntToStr(_nEvents) + " aligning at trigger number, too many data headers (set __TRUNC_EVENT)");
        }
        else {  // create new event
          addEvent();
        }
      }
      if (tNdataHeader == 0) {  // set the BCID of the first data header
        tStartBCID = tActualBCID;
        tStartLVL1ID = tActualLVL1ID;
      }
      else {
        tDbCID++;  // increase relative BCID counter [0:15]
        if (_fEI4B) {
          if (tStartBCID + tDbCID > __BCIDCOUNTERSIZE_FEI4B - 1)  // BCID counter overflow for FEI4B (10 bit BCID counter)
            tStartBCID = tStartBCID - __BCIDCOUNTERSIZE_FEI4B;
        }
        else {
          if (tStartBCID + tDbCID > __BCIDCOUNTERSIZE_FEI4A - 1)  // BCID counter overflow for FEI4A (8 bit BCID counter)
            tStartBCID = tStartBCID - __BCIDCOUNTERSIZE_FEI4A;
        }

        if (tStartBCID + tDbCID != tActualBCID) {  // check if BCID is increasing by 1 in the event window, if not close actual event and create new event with actual data header
          if (tActualLVL1ID == tStartLVL1ID) {  // happens sometimes, non inc. BCID, FE feature, only abort if the LVL1ID is not constant (if no external trigger is used or)
            addEventStatus(__BCID_JUMP);
            if (Basis::infoSet())
              info("interpretRawData: " + IntToStr(_nDataWords) + " DH WORD " + IntToStr(tActualWord) + " - " + IntToStr(tStartBCID + tDbCID) + "!=" + IntToStr(tActualBCID) + " at event " + LongIntToStr(_nEvents) + " BCID jumping");
          } else if (_alignAtTriggerNumber || _alignAtTdcWord) {  // rely here on the trigger number or TDC word and do not start a new event
            addEventStatus(__BCID_JUMP);
            if (Basis::infoSet())
              info("interpretRawData: " + IntToStr(_nDataWords) + " DH WORD " + IntToStr(tActualWord) + " - " + IntToStr(tStartBCID + tDbCID) + "!=" + IntToStr(tActualBCID) + " at event " + LongIntToStr(_nEvents) + " BCID jumping");
          } else {
            tBCIDerror = true;  // BCID number wrong, abort event and take actual data header for the first hit of the new event
            addEventStatus(__EVENT_INCOMPLETE);
            if (Basis::infoSet())
              info("interpretRawData: " + IntToStr(_nDataWords) + " DH WORD " + IntToStr(tActualWord) + " - " + IntToStr(tStartBCID + tDbCID) + "!=" + IntToStr(tActualBCID) + " at event " + LongIntToStr(_nEvents) + " event incomplete");
          }
        }
        if (!tBCIDerror && tActualLVL1ID != tStartLVL1ID) {  // LVL1ID not constant, is expected for CMOS pulse trigger/HitOR self-trigger, but not for trigger word triggering
          addEventStatus(__NON_CONST_LVL1ID);
          if (Basis::infoSet())
            info("interpretRawData: " + IntToStr(_nDataWords) + " DH WORD " + IntToStr(tActualWord) + " - " + IntToStr(tActualLVL1ID) + "!=" + IntToStr(tStartLVL1ID) + " at event " + LongIntToStr(_nEvents) + " LVL1 is not constant");
        }
      }
      tNdataHeader++;  // increase event data header counter
      if (Basis::debugSet())
        debug(std::string(" ") + IntToStr(_nDataWords) + " DH WORD " + IntToStr(tActualWord) + " - " + "LVL1ID/BCID " + IntToStr(tActualLVL1ID) + "/" + IntToStr(tActualBCID) + " at event " + LongIntToStr(_nEvents));
    }
    else if (isTriggerWord(tActualWord)) {  // data word is trigger word, is first word of the event data if external trigger is present
      _nTriggers++;  // increase global trigger word counter
      if (_alignAtTriggerNumber) {  // use trigger number for event building, first word is trigger word in event data stream
        // check for _firstTriggerNrSet, prevent building new event for the very first trigger word
        if (!_firstTriggerNrSet && tNdataHeader >= _NbCID) {  // for old data where trigger word (first raw data word) might be missing
          if (Basis::infoSet())
            info("interpretRawData: " + IntToStr(_nDataWords) + " TW WORD " + IntToStr(tActualWord) + " - " + IntToStr(tNdataHeader) + ">" + IntToStr(_NbCID) + " at event " + LongIntToStr(_nEvents) +  " missing trigger (adding new event)");
          addEventStatus(__NO_TRG_WORD);
          addEvent();
        }
        else if (_firstTriggerNrSet && tNdataHeader < _NbCID) {  // when data headers are missing
          if (Basis::infoSet())
            info("interpretRawData: " + IntToStr(_nDataWords) + " TW WORD " + IntToStr(tActualWord) + " - " + IntToStr(tNdataHeader) + "<" + IntToStr(_NbCID) + " at event " + LongIntToStr(_nEvents) + " event incomplete (adding new event)");
          addEventStatus(__EVENT_INCOMPLETE);
          addEvent();
        }
        else if (_firstTriggerNrSet) {  // usually the case
          addEvent();
        }
      } else if (tNdataHeader >= _NbCID) {  // use trigger word as indicator for a new event, otherwise trigger number gets assigned to wrong event
          addEvent();
      }

      tTriggerWord++;  // increase event trigger word counter
      if (_TriggerDataFormat == TRIGGER_FROMAT_TRIGGER_NUMBER) {  // trigger number
        tTriggerNumber = TRIGGER_DATA_MACRO(tActualWord);  // 31bit trigger number
        tTriggerTimeStamp = 0;  // 31bit time stamp
      }
      else if (_TriggerDataFormat == TRIGGER_FROMAT_TIME_STAMP) {  // time stamp
        // for compatibility assign to tTriggerNumber
        tTriggerNumber = 0;  // 31bit trigger number
        tTriggerTimeStamp = TRIGGER_DATA_MACRO(tActualWord);  // 31bit time stamp
      }
      else if (_TriggerDataFormat == TRIGGER_FROMAT_COMBINED) {  // combined
        tTriggerNumber = TRIGGER_NUMBER_COMBINED_MACRO(tActualWord);  // 16bit trigger number
        tTriggerTimeStamp = TRIGGER_TIME_STAMP_COMBINED_MACRO(tActualWord);  // 15bit time stamp
      }
      else {
        throw std::out_of_range("Invalid mode for trigger data format.");
      }
      if (Basis::debugSet()) {
        if (_TriggerDataFormat == TRIGGER_FROMAT_TRIGGER_NUMBER) {  // trigger number
          debug(std::string(" ") + IntToStr(_nDataWords) + " TW NUMBER " + IntToStr(tTriggerNumber) + " WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
        }
        else if (_TriggerDataFormat == TRIGGER_FROMAT_TIME_STAMP) {  // time stamp
          debug(std::string(" ") + IntToStr(_nDataWords) + " TW TIME STAMP " + IntToStr(tTriggerTimeStamp) + " WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
        }
        else if (_TriggerDataFormat == TRIGGER_FROMAT_COMBINED) {  // combined
          debug(std::string(" ") + IntToStr(_nDataWords) + " TW TIME STAMP " + IntToStr(tTriggerTimeStamp) + " TW NUMBER " + IntToStr(tTriggerNumber) + " WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
          }
      }
      // TLU error handling
      if (!_firstTriggerNrSet) {
        _firstTriggerNrSet = true;
      } else if ((_TriggerDataFormat != TRIGGER_FROMAT_TIME_STAMP) && (_lastTriggerNumber + 1 != tTriggerNumber) && !(_lastTriggerNumber == _maxTriggerNumber && tTriggerNumber == 0)) {
        addTriggerStatus(__TRG_NUMBER_INC_ERROR);
        if (Basis::warningSet())
          warning("interpretRawData: Trigger Number not increasing by 1 (old/new): " + IntToStr(_lastTriggerNumber) + "/" + IntToStr(tTriggerNumber) + " at event " + LongIntToStr(_nEvents));
      }

      if (tTriggerWord == 1) {  // event trigger number is trigger number of first trigger word within the event
        tEventTriggerNumber = tTriggerNumber;
        tEventTriggerTimeStamp = tTriggerTimeStamp;
      }
      // store for next event in case of missing trigger word
      _lastTriggerNumber = tTriggerNumber;
      _lastTriggerTimeStamp = tTriggerTimeStamp;
    } else if (getInfoFromServiceRecord(tActualWord, tActualSRcode, tActualSRcounter)) {  // data word is service record
      if (Basis::debugSet())
        debug(std::string(" ") + IntToStr(_nDataWords) + " SR " + IntToStr(tActualSRcode) + " (" + IntToStr(tActualSRcounter) + ") at event " + LongIntToStr(_nEvents));
      addServiceRecord(tActualSRcode, tActualSRcounter);
      addEventStatus(__HAS_SR);
      _nServiceRecords++;
    } else if (isTdcWord(tActualWord)) {  // data word is a TDC word
      addTdcValue(TDC_VALUE_MACRO(tActualWord));
      // TDC trigger distance
      if (_haveTdcTriggerDistance) {
        addTdcTriggerDistanceValue(TDC_TRIG_DIST_MACRO(tActualWord));
      }
      _nTDCWords++;
      if (_haveTdcTriggerDistance && (TDC_TRIG_DIST_MACRO(tActualWord) > _maxTdcDelay)) {  // if TDC trigger distance > _maxTdcDelay, ignore TDC word
        if (Basis::debugSet())
          debug(std::string(" ") + IntToStr(_nDataWords) + " TDC WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents) + " TDC TRIGGER DISTANCE " + IntToStr(TDC_TRIG_DIST_MACRO(tActualWord)) + " MAX DELAY REJECTED");
        continue;
      }

      // create new event if the option to align at TDC words is active AND the previous event has seen already all data headers OR the previous event had no TDC word
      if (_alignAtTdcWord && _firstTdcSet && ((tNdataHeader > _NbCID - 1) || ((tEventStatus & __TDC_WORD) != __TDC_WORD))) {
        addEvent();
      }
      _firstTdcSet = true;
      // if the event has already a valid TDC word
      if (((tEventStatus & __TDC_WORD) == __TDC_WORD) && (tTdcValue != 0) && (_haveTdcTriggerDistance && (tTdcTriggerDistance != 255))) {
        // already got valid TDC word, set correponding flag if current TDC word is also valid; set proper event status code and keep values of first TDC word
        if ((TDC_VALUE_MACRO(tActualWord) != 0) && (_haveTdcTriggerDistance && (TDC_TRIG_DIST_MACRO(tActualWord) != 255))) {
          addEventStatus(__MORE_THAN_ONE_TDC_WORD);
        }
      } else {  // first TDC word in event or the event has already an invalid TDC word, update values
        addEventStatus(__TDC_WORD);
        tTdcValue = TDC_VALUE_MACRO(tActualWord);
        if (_haveTdcTriggerTimeStamp && _haveTdcTriggerDistance) {
          tTdcTimeStamp = TDC_SHORT_TIME_STAMP_MACRO(tActualWord);
          tTdcTriggerDistance = TDC_TRIG_DIST_MACRO(tActualWord);
        } else if (_haveTdcTriggerTimeStamp && !_haveTdcTriggerDistance) {
          tTdcTimeStamp = TDC_TIME_STAMP_MACRO(tActualWord);
        } else if (!_haveTdcTriggerTimeStamp && _haveTdcTriggerDistance) {
          tTdcTriggerDistance = TDC_TRIG_DIST_MACRO(tActualWord);
        }
        if (Basis::debugSet()) {
          if (_haveTdcTriggerTimeStamp && _haveTdcTriggerDistance) {
            debug(std::string(" ") + IntToStr(_nDataWords) + " TDC WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents) + " TDC VALUE " + IntToStr(tTdcValue) + " TDC TRIGGER DISTANCE " + IntToStr(tTdcTriggerDistance) + " TDC TIME STAMP " + IntToStr(tTdcTimeStamp));
          } else if (_haveTdcTriggerTimeStamp && !_haveTdcTriggerDistance) {
            debug(std::string(" ") + IntToStr(_nDataWords) + " TDC WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents) + " TDC VALUE " + IntToStr(tTdcValue) + " TDC TIME STAMP " + IntToStr(TDC_TIME_STAMP_MACRO(tActualWord)));
          } else if (!_haveTdcTriggerTimeStamp && _haveTdcTriggerDistance) {
            debug(std::string(" ") + IntToStr(_nDataWords) + " TDC WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents) + " TDC VALUE " + IntToStr(tTdcValue) + " TDC TRIGGER DISTANCE " + IntToStr(tTdcTriggerDistance) + " TDC WORD COUNTER " + IntToStr(TDC_SHORT_TIME_STAMP_MACRO(tActualWord)));
          } else {
            debug(std::string(" ") + IntToStr(_nDataWords) + " TDC WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents) + " TDC VALUE " + IntToStr(tTdcValue) + " TDC WORD COUNTER " + IntToStr(TDC_TIME_STAMP_MACRO(tActualWord)));
          }
        }
      }
    } else if (isDataRecord(tActualWord)) {  // data word is data record if true is returned
      if (getHitsfromDataRecord(tActualWord, tActualCol1, tActualRow1, tActualTot1, tActualCol2, tActualRow2, tActualTot2)) {
        tNdataRecord++;  // increase data record counter for this event
        _nDataRecords++;  // increase total data record counter
        if (tActualTot1 >= 0)               // add hit if hit info is reasonable (TOT1 >= 0)
          if (!(addHit(tDbCID, tActualLVL1ID, tActualCol1, tActualRow1, tActualTot1, tActualBCID)))
            if (Basis::warningSet())
              warning("interpretRawData: " + IntToStr(_nDataWords) + " DR " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents) + " too many data records");
        if (tActualTot2 >= 0)               // add hit if hit info is reasonable and set (TOT2 >= 0)
          if (!(addHit(tDbCID, tActualLVL1ID, tActualCol2, tActualRow2, tActualTot2, tActualBCID)))
            if (Basis::warningSet())
              warning("interpretRawData: " + IntToStr(_nDataWords) + " DR " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents) + " too many data records");
        if (Basis::debugSet()) {
          std::stringstream tDebug;
          tDebug << " " << _nDataWords << " DR COL1/ROW1/TOT1  COL2/ROW2/TOT2 " << tActualCol1 << "/" << tActualRow1 << "/" << tActualTot1 << "  " << tActualCol2 << "/" << tActualRow2 << "/" << tActualTot2 << " rBCID " << tDbCID << " at event " << _nEvents;
          debug(tDebug.str());
        }
      } else {
        if (Basis::warningSet())
          warning("interpretRawData: " + IntToStr(_nDataWords) + " UNKNOWN WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
        if (Basis::debugSet())
          debug(std::string(" ") + IntToStr(_nDataWords) + " UNKNOWN WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
      }
    } else if (isAddressRecord(tActualWord)) {  // data word is address record if true is returned
      _nAddressRecords++;
      if (Basis::debugSet()) {
        unsigned int tAddress = 0;
        bool isShiftRegister = false;
        if (isAddressRecord(tActualWord, tAddress, isShiftRegister)) {
          if (isShiftRegister)
            debug(std::string(" ") + IntToStr(_nDataWords) + " ADDRESS RECORD SHIFT REG. " + IntToStr(tAddress) + " WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
          else
            debug(std::string(" ") + IntToStr(_nDataWords) + " ADDRESS RECORD GLOBAL REG. " + IntToStr(tAddress) + " WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
        }
      }
    } else if (isValueRecord(tActualWord)) {  // data word is value record if true is returned
      _nValueRecords++;
      if (Basis::debugSet()) {
        unsigned int tValue = 0;
        if (isValueRecord(tActualWord, tValue)) {
          debug(std::string(" ") + IntToStr(_nDataWords) + " VALUE RECORD " + IntToStr(tValue) + " at event " + LongIntToStr(_nEvents));
        }
      }
    } else {
      if (isOtherWord(tActualWord)) {  // other data words
        addEventStatus(__OTHER_WORD);
        _nOtherWords++;
        if (Basis::debugSet()) {
          debug(std::string(" ") + IntToStr(_nDataWords) + " OTHER WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
        }
      } else {  // remaining data words, unknown words
        addEventStatus(__UNKNOWN_WORD);
        _nUnknownWords++;
        if (Basis::warningSet())
          warning("interpretRawData: " + IntToStr(_nDataWords) + " UNKNOWN WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
        if (Basis::debugSet())
          debug(std::string(" ") + IntToStr(_nDataWords) + " UNKNOWN WORD " + IntToStr(tActualWord) + " at event " + LongIntToStr(_nEvents));
      }
    }
    if (tBCIDerror) {  // tBCIDerror is raised if BCID is not increasing by 1, most likely due to incomplete data transmission, so start new event, actual word is data header here
      if (Basis::warningSet())
        warning("interpretRawData " + IntToStr(_nDataWords) + " BCID ERROR at event " + LongIntToStr(_nEvents));
      addEvent();
      _nIncompleteEvents++;
      getTimefromDataHeader(tActualWord, tActualLVL1ID, tStartBCID);
      tNdataHeader = 1;  // tNdataHeader is already 1, because actual word is first data of new event
      tStartBCID = tActualBCID;
      tStartLVL1ID = tActualLVL1ID;
    }
    correlateMetaWordIndex(_nEvents, _dataWordIndex);
    _dataWordIndex++;
    tNdataWords++;
  }
  return true;
}

bool Interpret::setMetaData(MetaInfo* &rMetaInfo, const unsigned int& tLength)
{
  info("setMetaData with " + IntToStr(tLength) + " entries");
  _isMetaTableV2 = false;
  _metaInfo = rMetaInfo;
  if (tLength == 0) {
    warning("setMetaWordIndex: data is empty");
    return false;
  }
  // sanity check
  for (unsigned int i = 0; i < tLength - 1; ++i) {
    if (_metaInfo[i].startIndex + _metaInfo[i].length != _metaInfo[i].stopIndex)
      throw std::out_of_range("Meta word index out of range.");
    if (_metaInfo[i].stopIndex != _metaInfo[i + 1].startIndex && _metaInfoV2[i + 1].startIndex != 0)
      throw std::out_of_range("Meta word index out of range.");
    if (_metaInfo[i].timeStamp > _metaInfo[i + 1].timeStamp)
      throw std::out_of_range("Time stamp not increasing.");
  }
  if (_metaInfo[tLength - 1].startIndex + _metaInfo[tLength - 1].length != _metaInfo[tLength - 1].stopIndex)
    throw std::out_of_range("Meta word index out of range.");

  _metaEventIndexLength = tLength;
  _metaDataSet = true;

  return true;
}

bool Interpret::setMetaDataV2(MetaInfoV2* &rMetaInfo, const unsigned int& tLength)
{
  info("setMetaDataV2 with " + IntToStr(tLength) + " entries");
  _isMetaTableV2 = true;
  _metaInfoV2 = rMetaInfo;
  if (tLength == 0) {
    warning(std::string("setMetaWordIndex: data is empty"));
    return false;
  }
  // sanity check
  for (unsigned int i = 0; i < tLength - 1; ++i) {
    if (_metaInfoV2[i].startIndex + _metaInfoV2[i].length != _metaInfoV2[i].stopIndex)
      throw std::out_of_range("Meta word index out of range.");
    if (_metaInfoV2[i].stopIndex != _metaInfoV2[i + 1].startIndex && _metaInfoV2[i + 1].startIndex != 0)
      throw std::out_of_range("Meta word index out of range.");
    if (_metaInfoV2[i].startTimeStamp > _metaInfoV2[i].stopTimeStamp || _metaInfoV2[i].stopTimeStamp > _metaInfoV2[i + 1].startTimeStamp)
      throw std::out_of_range("Time stamp not increasing.");
  }
  if (_metaInfoV2[tLength - 1].startIndex + _metaInfoV2[tLength - 1].length != _metaInfoV2[tLength - 1].stopIndex)
    throw std::out_of_range("Meta word index out of range.");

  _metaEventIndexLength = tLength;
  _metaDataSet = true;

  return true;
}

void Interpret::getHits(HitInfo*& rHitInfo, unsigned int& rSize, bool copy)
{
  debug("getHits(...)");
  if (copy)
    std::copy(_hitInfo, _hitInfo + _hitInfoSize, rHitInfo);
  else
    rHitInfo = _hitInfo;
  rSize = _hitIndex;
}

void Interpret::setHitsArraySize(const unsigned int &rSize)
{
  info("setHitsArraySize(...) with size " + IntToStr(rSize));
  deleteHitArray();
  _hitInfoSize = rSize;
  allocateHitArray();
}

void Interpret::setMetaDataEventIndex(uint64_t*& rEventNumber, const unsigned int& rSize)
{
  info("setMetaDataEventIndex(...) with length " + IntToStr(rSize));
  _metaEventIndex = rEventNumber;
  _metaEventIndexLength = rSize;
}

void Interpret::setMetaDataWordIndex(MetaWordInfoOut*& rWordNumber, const unsigned int& rSize)
{
  info("setMetaDataWordIndex(...) with length " + IntToStr(rSize));
  _metaWordIndex = rWordNumber;
  _metaWordIndexLength = rSize;
}

void Interpret::resetCounters()
{
  info("resetCounters()");
  _nDataWords = 0;
  _nTriggers = 0;
  _nEvents = 0;
  _nIncompleteEvents = 0;
  _nDataRecords = 0;
  _nDataHeaders = 0;
  _nAddressRecords = 0;
  _nValueRecords = 0;
  _nServiceRecords = 0;
  _nUnknownWords = 0;
  _nTDCWords = 0;
  _nOtherWords = 0;
  _nHits = 0;
  _nSmallHits = 0;
  _nEmptyEvents = 0;
  _nMaxHitsPerEvent = 0;
  _firstTriggerNrSet = false;
  _firstTdcSet = false;
  _lastTriggerNumber = 0;
  _lastTriggerTimeStamp = 0;
  _dataWordIndex = 0;
  resetTriggerStatusCounterArray();
  resetEventStatusCounterArray();
  resetTdcValueArray();
  resetTdcTriggerDistanceArray();
  resetServiceRecordCounterArray();
}

void Interpret::resetEventVariables()
{
  tNdataWords = 0;
  tNdataHeader = 0;
  tNdataRecord = 0;
  tDbCID = 0;
  tTriggerStatus = 0;
  tEventStatus = 0;
  tServiceRecord = 0;
  tBCIDerror = false;
  tTriggerWord = 0;
  tTdcValue = 0;
  tTdcTimeStamp = 0;
  tTdcTriggerDistance = 0;
  tTriggerNumber = 0;
  tTriggerTimeStamp = 0;
  tEventTriggerNumber = 0;
  tEventTriggerTimeStamp = 0;
  tStartBCID = 0;
  tStartLVL1ID = 0;
  tHitBufferIndex = 0;
  tTotalHits = 0;
}

void Interpret::resetHistograms()
{
  resetTriggerStatusCounterArray();
  resetEventStatusCounterArray();
  resetTdcValueArray();
  resetTdcTriggerDistanceArray();
  resetServiceRecordCounterArray();
}

void Interpret::createMetaDataWordIndex(bool CreateMetaDataWordIndex)
{
  debug("createMetaDataWordIndex");
  _createMetaDataWordIndex = CreateMetaDataWordIndex;
}

void Interpret::createEmptyEventHits(bool CreateEmptyEventHits)
{
  debug("createEmptyEventHits");
  _createEmptyEventHits = CreateEmptyEventHits;
}

void Interpret::setNbCIDs(const unsigned int& NbCIDs)
{
  _NbCID = NbCIDs;
}

void Interpret::setMaxTot(const unsigned int& rMaxTot)
{
  _maxTot = rMaxTot;
}

void Interpret::setMaxTdcDelay(const unsigned int& rtMaxTdcDelay)
{
  _maxTdcDelay = rtMaxTdcDelay;
}

void Interpret::alignAtTriggerNumber(bool alignAtTriggerNumber)
{
  info("alignAtTriggerNumber()");
  _alignAtTriggerNumber = alignAtTriggerNumber;
}

void Interpret::setMaxTriggerNumber(const unsigned int& rMaxTriggerNumber)
{
  _maxTriggerNumber = rMaxTriggerNumber;
}

void Interpret::alignAtTdcWord(bool alignAtTdcWord)
{
  info("alignAtTdcWord()");
  _alignAtTdcWord = alignAtTdcWord;
}

void Interpret::setTriggerDataFormat(const unsigned int& rTriggerDataFormat)
{
  info("setTriggerDataFormat()");
  _TriggerDataFormat = rTriggerDataFormat;
}

void Interpret::setTdcTriggerTimeStamp(bool haveTdcTriggerTimeStamp)
{
  info("setTdcTriggerTimeStamp()");
  _haveTdcTriggerTimeStamp = haveTdcTriggerTimeStamp;
}

void Interpret::setTdcTriggerDistance(bool haveTdcTriggerDistance)
{
  info("setTdcTriggerDistance()");
  _haveTdcTriggerDistance = haveTdcTriggerDistance;
}

void Interpret::getServiceRecordsCounters(unsigned int*& rServiceRecordsCounter, unsigned int& rNserviceRecords, bool copy)
{
  debug("getServiceRecordsCounters(...)");
  if (copy)
    std::copy(_serviceRecordCounter, _serviceRecordCounter + __NSERVICERECORDS, rServiceRecordsCounter);
  else
    rServiceRecordsCounter = _serviceRecordCounter;

  rNserviceRecords = __NSERVICERECORDS;
}

void Interpret::getEventStatusCounters(unsigned int*& rEventStatusCounter, unsigned int& rNeventStatusCounters, bool copy)
{
  debug("getEventStatusCounters(...)");
  if (copy)
    std::copy(_eventStatusCounter, _eventStatusCounter + __N_ERROR_CODES, rEventStatusCounter);
  else
    rEventStatusCounter = _eventStatusCounter;

  rNeventStatusCounters = __N_ERROR_CODES;
}

void Interpret::getTriggerStatusCounters(unsigned int*& rTriggerStatusCounter, unsigned int& rNTriggerStatusCounters, bool copy)
{
  debug(std::string("getTriggerStatusCounters(...)"));
  if (copy)
    std::copy(_triggerStatusCounter, _triggerStatusCounter + __TRG_N_ERROR_CODES, rTriggerStatusCounter);
  else
    rTriggerStatusCounter = _triggerStatusCounter;

  rNTriggerStatusCounters = __TRG_N_ERROR_CODES;
}

void Interpret::getTdcValues(unsigned int*& rTdcValue, unsigned int& rNtdcValues, bool copy)
{
  debug("getEventStatusCounters(...)");
  if (copy)
    std::copy(_tdcValue, _tdcValue + __N_TDC_VALUES, rTdcValue);
  else
    rTdcValue = _tdcValue;

  rNtdcValues = __N_TDC_VALUES;
}

void Interpret::getTdcTriggerDistances(unsigned int*& rTdcTriggerDistance, unsigned int& rNtdcTriggerDistance, bool copy)
{
  debug("getEventStatusCounters(...)");
  if (copy)
    std::copy(_tdcTriggerDistance, _tdcTriggerDistance + __N_TDC_TRG_DIST_VALUES, rTdcTriggerDistance);
  else
    rTdcTriggerDistance = _tdcTriggerDistance;

  rNtdcTriggerDistance = __N_TDC_TRG_DIST_VALUES;
}

unsigned int Interpret::getNwords()
{
  return _nDataWords;
}

void Interpret::printSummary()
{
  std::cout << "# Data Words        " << std::right << std::setw(15) << _nDataWords << "\n";
  std::cout << "# Data Headers      " << std::right << std::setw(15) << _nDataHeaders << "\n";
  std::cout << "# Data Records      " << std::right << std::setw(15) << _nDataRecords << "\n";
  std::cout << "# Address Records   " << std::right << std::setw(15) << _nAddressRecords << "\n";
  std::cout << "# Value Records     " << std::right << std::setw(15) << _nValueRecords << "\n";
  std::cout << "# Service Records   " << std::right << std::setw(15) << _nServiceRecords << "\n";
  std::cout << "# TDC Words         " << std::right << std::setw(15) << _nTDCWords << "\n";
  std::cout << "# Trigger Words     " << std::right << std::setw(15) << _nTriggers << "\n";
  std::cout << "# Other Words       " << std::right << std::setw(15) << _nOtherWords << "\n";
  std::cout << "# Unknown Words     " << std::right << std::setw(15) << _nUnknownWords << "\n\n";

  std::cout << "# Events            " << std::right << std::setw(15) << _nEvents << "\n";
  std::cout << "# Empty Events      " << std::right << std::setw(15) << _nEmptyEvents << "\n";
  std::cout << "# Incomplete Events " << std::right << std::setw(15) << _nIncompleteEvents << "\n\n";

  std::cout << "# Hits              " << std::right << std::setw(15) << _nHits << "\n";
  std::cout << "# Small/Late Hits   " << std::right << std::setw(15) << _nSmallHits << "\n";
  std::cout << "# MaxHitsPerEvent   " << std::right << std::setw(15) << _nMaxHitsPerEvent << "\n\n";

  std::cout << "# Event Status\n";
  std::cout << "0 " << std::right << std::setw(15) << _eventStatusCounter[0] << " Events with SR\n";
  std::cout << "1 " << std::right << std::setw(15) << _eventStatusCounter[1] << " Events without trigger word\n";
  std::cout << "2 " << std::right << std::setw(15) << _eventStatusCounter[2] << " Events without constant LVL1ID\n";
  std::cout << "3 " << std::right << std::setw(15) << _eventStatusCounter[3] << " Events that were incomplete (# BCIDs wrong)\n";
  std::cout << "4 " << std::right << std::setw(15) << _eventStatusCounter[4] << " Events with unknown words\n";
  std::cout << "5 " << std::right << std::setw(15) << _eventStatusCounter[5] << " Events with jumping BCIDs\n";
  std::cout << "6 " << std::right << std::setw(15) << _eventStatusCounter[6] << " Events with trigger status != 0\n";
  std::cout << "7 " << std::right << std::setw(15) << _eventStatusCounter[7] << " Events that were truncated due to too many data headers or data records\n";
  std::cout << "8 " << std::right << std::setw(15) << _eventStatusCounter[8] << " Events with TDC words\n";
  std::cout << "9 " << std::right << std::setw(15) << _eventStatusCounter[9] << " Events with more than one TDC word\n";
  std::cout << "10" << std::right << std::setw(15) << _eventStatusCounter[10] << " Events with TDC overflow\n";
  std::cout << "11" << std::right << std::setw(15) << _eventStatusCounter[11] << " Events without hits\n";
  std::cout << "12" << std::right << std::setw(15) << _eventStatusCounter[12] << " Events with other data words\n";
  std::cout << "13" << std::right << std::setw(15) << _eventStatusCounter[13] << " Events with more than one hit\n";
  std::cout << "14" << std::right << std::setw(15) << _eventStatusCounter[14] << " Events with invalid TDC values\n\n";

  std::cout << "# Trigger Status\n";
  std::cout << "0 " << std::right << std::setw(15) << _triggerStatusCounter[0] << " Trigger number not increasing by 1\n";
  std::cout << "1 " << std::right << std::setw(15) << _triggerStatusCounter[1] << " # Trigger per event > 1\n\n";

  std::cout << "# Service Records\n";
  bool no_service_record = true;
  for (unsigned int i = 0; i < __NSERVICERECORDS; ++i) {
    if (_serviceRecordCounter[i] > 0) {
      no_service_record = false;
      std::cout << std::left << std::setw(2) << i << std::right << std::setw(15) << _serviceRecordCounter[i] << "\n";
    }
  }
  if (no_service_record) {
    std::cout << "NO SERVICE RECORD" << "\n\n";
  } else {
    std::cout << "\n";
  }
}

void Interpret::printStatus()
{
  std::cout << "configuration parameters\n";
  std::cout << "_NbCID " << _NbCID << "\n";
  std::cout << "_maxTot " << _maxTot << "\n";
  std::cout << "_fEI4B " << _fEI4B << "\n";
  std::cout << "_debugEvents " << _debugEvents << "\n";
  std::cout << "_startDebugEvent " << _startDebugEvent << "\n";
  std::cout << "_stopDebugEvent " << _stopDebugEvent << "\n";
  std::cout << "_alignAtTriggerNumber " << _alignAtTriggerNumber << "\n";
  std::cout << "_alignAtTdcWord " << _alignAtTdcWord << "\n";
  std::cout << "_TriggerDataFormat " << _TriggerDataFormat << "\n";
  std::cout << "_haveTdcTriggerTimeStamp " << _haveTdcTriggerTimeStamp << "\n";
  std::cout << "_haveTdcTriggerDistance " << _haveTdcTriggerDistance << "\n";
  std::cout << "_maxTdcDelay " << _maxTdcDelay << "\n";

  std::cout << "\ncurrent event variables\n";
  std::cout << "tNdataWords " << tNdataWords << "\n";
  std::cout << "tNdataHeader " << tNdataHeader << "\n";
  std::cout << "tNdataRecord " << tNdataRecord << "\n";
  std::cout << "tStartBCID " << tStartBCID << "\n";
  std::cout << "tStartLVL1ID " << tStartLVL1ID << "\n";
  std::cout << "tDbCID " << tDbCID << "\n";
  std::cout << "tTriggerStatus " << tTriggerStatus << "\n";
  std::cout << "tEventStatus " << tEventStatus << "\n";
  std::cout << "tServiceRecord " << tServiceRecord << "\n";
  std::cout << "tTriggerNumber " << tTriggerNumber << "\n";
  std::cout << "tTriggerTimeStamp " << tTriggerTimeStamp << "\n";
  std::cout << "tTotalHits " << tTotalHits << "\n";
  std::cout << "tBCIDerror " << tBCIDerror << "\n";
  std::cout << "tTriggerWord " << tTriggerWord << "\n";
  std::cout << "tTdcValue " << tTdcValue << "\n";
  std::cout << "tTdcTimeStamp" << tTdcTimeStamp << "\n";
  std::cout << "tTdcTriggerDistance" << tTdcTriggerDistance << "\n";
  std::cout << "_lastTriggerNumber " << _lastTriggerNumber << "\n";
  std::cout << "_lastTriggerTimeStamp " << _lastTriggerTimeStamp << "\n";

  std::cout << "\ncounters/flags for the raw data processing\n";
  std::cout << "_nTriggers " << _nTriggers << "\n";
  std::cout << "_nEvents " << _nEvents << "\n";
  std::cout << "_nMaxHitsPerEvent " << _nMaxHitsPerEvent << "\n";
  std::cout << "_nEmptyEvents " << _nEmptyEvents << "\n";
  std::cout << "_nIncompleteEvents " << _nIncompleteEvents << "\n";
  std::cout << "_nOtherWords " << _nOtherWords << "\n";
  std::cout << "_nUnknownWords " << _nUnknownWords << "\n";
  std::cout << "_nTDCWords " << _nTDCWords << "\n\n";
  std::cout << "_nServiceRecords " << _nServiceRecords << "\n";
  std::cout << "_nDataRecords " << _nDataRecords << "\n";
  std::cout << "_nDataHeaders " << _nDataHeaders << "\n";
  std::cout << "_nAddressRecords " << _nAddressRecords << "\n";
  std::cout << "_nValueRecords " << _nValueRecords << "\n";
  std::cout << "_nHits " << _nHits << "\n";
  std::cout << "_nSmallHits " << _nSmallHits << "\n";
  std::cout << "_nDataWords " << _nDataWords << "\n";
  std::cout << "_firstTriggerNrSet " << _firstTriggerNrSet << "\n";
  std::cout << "_firstTdcSet " << _firstTdcSet << "\n";
}

void Interpret::debugEvents(const unsigned int& rStartEvent, const unsigned int& rStopEvent, const bool& debugEvents)
{
  _debugEvents = debugEvents;
  _startDebugEvent = rStartEvent;
  _stopDebugEvent = rStopEvent;
}

unsigned int Interpret::getHitSize()
{
  return sizeof(HitInfo);
}

void Interpret::reset()
{
  info("reset()");
  resetCounters();
  resetEventVariables();
  _lastMetaIndexNotSet = 0;
  _lastWordIndexSet = 0;
  _metaEventIndexLength = 0;
  _metaEventIndex = 0;
  _startWordIndex = 0;
  // initialize SRAM variables to 0
  tTriggerNumber = 0;
  tTriggerTimeStamp = 0;
  tActualLVL1ID = 0;
  tActualBCID = 0;
  tActualSRcode= 0;
  tActualSRcounter = 0;
}

void Interpret::resetMetaDataCounter()
{
  _lastWordIndexSet = 0;
  _dataWordIndex = 0;
}

// private

bool Interpret::addHit(const unsigned char& pRelBCID, const unsigned short int& pLVL1ID, const unsigned char& pColumn, const unsigned short int& pRow, const unsigned char& pTot, const unsigned short int& pBCID)  // add hit with event number, column, row, relative BCID [0:15], tot, trigger ID
{
  if (tHitBufferIndex < __MAXHITBUFFERSIZE) {
    _hitBuffer[tHitBufferIndex].event_number = _nEvents;
    _hitBuffer[tHitBufferIndex].trigger_number = tEventTriggerNumber;
    _hitBuffer[tHitBufferIndex].trigger_time_stamp = tEventTriggerTimeStamp;
    _hitBuffer[tHitBufferIndex].relative_BCID = pRelBCID;
    _hitBuffer[tHitBufferIndex].LVL1ID = pLVL1ID;
    _hitBuffer[tHitBufferIndex].column = pColumn;
    _hitBuffer[tHitBufferIndex].row = pRow;
    _hitBuffer[tHitBufferIndex].tot = pTot;
    _hitBuffer[tHitBufferIndex].BCID = pBCID;
    _hitBuffer[tHitBufferIndex].TDC = tTdcValue;
    _hitBuffer[tHitBufferIndex].TDC_time_stamp = tTdcTimeStamp;
    _hitBuffer[tHitBufferIndex].TDC_trigger_distance = tTdcTriggerDistance;
    _hitBuffer[tHitBufferIndex].service_record = tServiceRecord;
    _hitBuffer[tHitBufferIndex].trigger_status = tTriggerStatus;
    _hitBuffer[tHitBufferIndex].event_status = tEventStatus;
    if ((tEventStatus & __NO_HIT) != __NO_HIT)  // only count non-virtual hits
      tTotalHits++;
    tHitBufferIndex++;
    return true;
  } else {
    addEventStatus(__TRUNC_EVENT);  // too many hits in the event, abort this event, add truncated flag
    // addEvent();
    if (Basis::warningSet())
      warning(std::string("addHit: Hit buffer overflow prevented by ignoring hits at event " + LongIntToStr(_nEvents)), __LINE__);
  }
  return false;
}

void Interpret::storeHit(HitInfo& rHit)
{
  _nHits++;
  if (_hitIndex < _hitInfoSize) {
    if (_hitInfo != 0) {
      _hitInfo[_hitIndex] = rHit;
      _hitIndex++;
    } else {
      throw std::runtime_error("Output hit array not set.");
    }
  } else {
    if (Basis::errorSet())
      error("storeHit: _hitIndex = " + IntToStr(_hitIndex), __LINE__);
    throw std::out_of_range("Hit index out of range.");
  }
}

void Interpret::addEvent()
{
  if (Basis::debugSet()) {
    std::stringstream tDebug;
    tDebug << "addEvent() " << _nEvents;
    debug(tDebug.str());
  }
  if (tTotalHits == 0) {
    _nEmptyEvents++;
    if (_createEmptyEventHits) {
      addEventStatus(__NO_HIT);
      addHit(0, 0, 0, 0, 0, 0);
    }
  } else if (tTotalHits > 1) {
    addEventStatus(__MORE_THAN_ONE_HIT);
  }
  if (tTriggerWord == 0) {
    addEventStatus(__NO_TRG_WORD);
    if (_firstTriggerNrSet) {  // set the last existing trigger number for events without trigger number if trigger numbers exist
      tEventTriggerNumber = _lastTriggerNumber;
      tEventTriggerTimeStamp = _lastTriggerTimeStamp;
    }
  }
  if (tTriggerWord > 1) {
    addTriggerStatus(__TRG_NUMBER_MORE_ONE);
    if (Basis::warningSet())
      warning(std::string("addEvent: # trigger words > 1 at event " + LongIntToStr(_nEvents)));
  }
  // 4095 is TDC value overflow, 254 is TDC trigger distance overflow
  if ((tTdcValue == (__N_TDC_VALUES - 1)) || (_haveTdcTriggerDistance && (tTdcTriggerDistance == 254))) {
    addEventStatus(__TDC_OVERFLOW);
  }
  // 0 is invalid TDC value, 255 is invalied TDC trigger distance
  if ((tTdcValue == 0) || (_haveTdcTriggerDistance && (tTdcTriggerDistance == 255))) {
    addEventStatus(__TDC_INVALID);
  }

  storeEventHits();
  if (tTotalHits > _nMaxHitsPerEvent)
    _nMaxHitsPerEvent = tTotalHits;
  histogramTriggerStatusCode();
  histogramErrorCode();
  if (_createMetaDataWordIndex) {
    if (_actualMetaWordIndex < _metaWordIndexLength) {
      _metaWordIndex[_actualMetaWordIndex].eventIndex = _nEvents;
      _metaWordIndex[_actualMetaWordIndex].startWordIdex = _startWordIndex;
      _metaWordIndex[_actualMetaWordIndex].stopWordIdex = _startWordIndex + tNdataWords;  // excluding stop word index
      _startWordIndex = _nDataWords - 1;
      _actualMetaWordIndex++;
    } else {
      std::stringstream tInfo;
      tInfo << "Interpret::addEvent(): meta word index array is too small " << _actualMetaWordIndex << ">=" << _metaWordIndexLength;
      throw std::out_of_range(tInfo.str());
    }
  }
  _nEvents++;
  resetEventVariables();
}

void Interpret::storeEventHits()
{
  for (unsigned int i = 0; i < tHitBufferIndex; ++i) {
    // duplicate certain values for all hits in an event
    _hitBuffer[i].trigger_number = tEventTriggerNumber;
    _hitBuffer[i].trigger_time_stamp = tEventTriggerTimeStamp;
    _hitBuffer[i].TDC = tTdcValue;
    _hitBuffer[i].TDC_time_stamp = tTdcTimeStamp;
    _hitBuffer[i].TDC_trigger_distance = tTdcTriggerDistance;
    _hitBuffer[i].service_record = tServiceRecord;
    // set status bits at the very end
    _hitBuffer[i].trigger_status = tTriggerStatus;
    _hitBuffer[i].event_status = tEventStatus;
    storeHit(_hitBuffer[i]);
  }
}

void Interpret::correlateMetaWordIndex(const uint64_t& pEventNumber, const unsigned int& pDataWordIndex)
{
  if (_metaDataSet && pDataWordIndex == _lastWordIndexSet) {  // this check is to speed up the _metaEventIndex access by using the fact that the index has to increase for consecutive events
//    std::cout<<"_lastMetaIndexNotSet "<<_lastMetaIndexNotSet<<"\n";
    _metaEventIndex[_lastMetaIndexNotSet] = pEventNumber;
    if (_isMetaTableV2 == true) {
      _lastWordIndexSet = _metaInfoV2[_lastMetaIndexNotSet].stopIndex;
      _lastMetaIndexNotSet++;
      while (_metaInfoV2[_lastMetaIndexNotSet - 1].length == 0 && _lastMetaIndexNotSet < _metaEventIndexLength) {
        info("correlateMetaWordIndex: more than one readout during one event, correcting meta info");
//        std::cout<<"correlateMetaWordIndex: pEventNumber "<<pEventNumber<<" _lastWordIndexSet "<<_lastWordIndexSet<<" _lastMetaIndexNotSet "<<_lastMetaIndexNotSet<<"\n";
        _metaEventIndex[_lastMetaIndexNotSet] = pEventNumber;
        _lastWordIndexSet = _metaInfoV2[_lastMetaIndexNotSet].stopIndex;
        _lastMetaIndexNotSet++;
//        std::cout<<"correlateMetaWordIndex: pEventNumber "<<pEventNumber<<" _lastWordIndexSet "<<_lastWordIndexSet<<" _lastMetaIndexNotSet "<<_lastMetaIndexNotSet<<"\n";
//        std::cout<<" finished\n";
      }
    }
    else {
      _lastWordIndexSet = _metaInfo[_lastMetaIndexNotSet].stopIndex;
      _lastMetaIndexNotSet++;
      while (_metaInfo[_lastMetaIndexNotSet - 1].length == 0 && _lastMetaIndexNotSet < _metaEventIndexLength) {
        info("correlateMetaWordIndex: more than one readout during one event, correcting meta info");
//        std::cout<<"correlateMetaWordIndex: pEventNumber "<<pEventNumber<<" _lastWordIndexSet "<<_lastWordIndexSet<<" _lastMetaIndexNotSet "<<_lastMetaIndexNotSet<<"\n";
        _metaEventIndex[_lastMetaIndexNotSet] = pEventNumber;
        _lastWordIndexSet = _metaInfo[_lastMetaIndexNotSet].stopIndex;
        _lastMetaIndexNotSet++;
//        std::cout<<"correlateMetaWordIndex: pEventNumber "<<pEventNumber<<" _lastWordIndexSet "<<_lastWordIndexSet<<" _lastMetaIndexNotSet "<<_lastMetaIndexNotSet<<"\n";
//        std::cout<<" finished\n";
      }
    }
  }
}

bool Interpret::getTimefromDataHeader(const unsigned int& pSRAMWORD, unsigned int& pLVL1ID, unsigned int& pBCID)
{
  if (DATA_HEADER_MACRO(pSRAMWORD)) {
    if (_fEI4B) {
      pLVL1ID = DATA_HEADER_LV1ID_MACRO_FEI4B(pSRAMWORD);
      pBCID = DATA_HEADER_BCID_MACRO_FEI4B(pSRAMWORD);
    }
    else {
      pLVL1ID = DATA_HEADER_LV1ID_MACRO(pSRAMWORD);
      pBCID = DATA_HEADER_BCID_MACRO(pSRAMWORD);
    }
    return true;
  }
  return false;
}

bool Interpret::isDataRecord(const unsigned int& pSRAMWORD)
{
  if (DATA_RECORD_MACRO(pSRAMWORD)) {
    return true;
  }
  return false;
}

bool Interpret::isTdcWord(const unsigned int& pSRAMWORD)
{
  if (TDC_WORD_MACRO(pSRAMWORD))
    return true;
  return false;
}

bool Interpret::getHitsfromDataRecord(const unsigned int& pSRAMWORD, int& pColHit1, int& pRowHit1, int& pTotHit1, int& pColHit2, int& pRowHit2, int& pTotHit2)
{
  // check if the hit values are reasonable
  if ((DATA_RECORD_TOT1_MACRO(pSRAMWORD) == 0xF) || (DATA_RECORD_COLUMN1_MACRO(pSRAMWORD) < RAW_DATA_MIN_COLUMN) || (DATA_RECORD_COLUMN1_MACRO(pSRAMWORD) > RAW_DATA_MAX_COLUMN) || (DATA_RECORD_ROW1_MACRO(pSRAMWORD) < RAW_DATA_MIN_ROW) || (DATA_RECORD_ROW1_MACRO(pSRAMWORD) > RAW_DATA_MAX_ROW)) {
    warning(std::string("getHitsfromDataRecord: data record values (1. Hit) out of bounds at event " + LongIntToStr(_nEvents)));
    return false;
  }
  if ((DATA_RECORD_TOT2_MACRO(pSRAMWORD) != 0xF) && ((DATA_RECORD_COLUMN2_MACRO(pSRAMWORD) < RAW_DATA_MIN_COLUMN) || (DATA_RECORD_COLUMN2_MACRO(pSRAMWORD) > RAW_DATA_MAX_COLUMN) || (DATA_RECORD_ROW2_MACRO(pSRAMWORD) < RAW_DATA_MIN_ROW) || (DATA_RECORD_ROW2_MACRO(pSRAMWORD) > RAW_DATA_MAX_ROW))) {
    warning(std::string("getHitsfromDataRecord: data record values (2. Hit) out of bounds at event " + LongIntToStr(_nEvents)));
    return false;
  }

  // set first hit values
  if (DATA_RECORD_TOT1_MACRO(pSRAMWORD) <= _maxTot) {  // ommit late/small hit and no hit TOT values for the TOT(1) hit
    pColHit1 = DATA_RECORD_COLUMN1_MACRO(pSRAMWORD);
    pRowHit1 = DATA_RECORD_ROW1_MACRO(pSRAMWORD);
    pTotHit1 = DATA_RECORD_TOT1_MACRO(pSRAMWORD);
  }
  if (DATA_RECORD_TOT1_MACRO(pSRAMWORD) == 14) {
    _nSmallHits++;
  }

  // set second hit values
  if (DATA_RECORD_TOT2_MACRO(pSRAMWORD) <= _maxTot) {  // ommit late/small hit and no hit (15) tot values for the TOT(2) hit
    pColHit2 = DATA_RECORD_COLUMN2_MACRO(pSRAMWORD);
    pRowHit2 = DATA_RECORD_ROW2_MACRO(pSRAMWORD);
    pTotHit2 = DATA_RECORD_TOT2_MACRO(pSRAMWORD);
  }
  if (DATA_RECORD_TOT2_MACRO(pSRAMWORD) == 14) {
    _nSmallHits++;
  }
  return true;
}

bool Interpret::getInfoFromServiceRecord(const unsigned int& pSRAMWORD, unsigned int& pSRcode, unsigned int& pSRcount)
{
  if (SERVICE_RECORD_MACRO(pSRAMWORD)) {
    pSRcode = SERVICE_RECORD_CODE_MACRO(pSRAMWORD);
    if (_fEI4B) {
      if (pSRcode == 14)
        pSRcount = 1;
      else if (pSRcode == 16)
        pSRcount = SERVICE_RECORD_ETC_MACRO_FEI4B(pSRAMWORD);
      else
        pSRcount = SERVICE_RECORD_COUNTER_MACRO(pSRAMWORD);
    }
    else {
      pSRcount = SERVICE_RECORD_COUNTER_MACRO(pSRAMWORD);
    }
    return true;
  }
  return false;
}

bool Interpret::isTriggerWord(const unsigned int& pSRAMWORD)
{
  if (TRIGGER_WORD_MACRO(pSRAMWORD))  // data word is trigger word
    return true;
  return false;
}

bool Interpret::isAddressRecord(const unsigned int& pSRAMWORD, unsigned int& rAddress, bool& isShiftRegister)
{
  if (ADDRESS_RECORD_MACRO(pSRAMWORD)) {
    if (ADDRESS_RECORD_TYPE_SET_MACRO(pSRAMWORD))
      isShiftRegister = true;
    rAddress = ADDRESS_RECORD_ADDRESS_MACRO(pSRAMWORD);
    return true;
  }
  return false;
}

bool Interpret::isAddressRecord(const unsigned int& pSRAMWORD)
{
  if (ADDRESS_RECORD_MACRO(pSRAMWORD)) {
    return true;
  }
  return false;
}

bool Interpret::isValueRecord(const unsigned int& pSRAMWORD, unsigned int& rValue)
{
  if (VALUE_RECORD_MACRO(pSRAMWORD)) {
    rValue = VALUE_RECORD_VALUE_MACRO(pSRAMWORD);
    return true;
  }
  return false;
}

bool Interpret::isValueRecord(const unsigned int& pSRAMWORD)
{
  if (VALUE_RECORD_MACRO(pSRAMWORD)) {
    return true;
  }
  return false;
}

bool Interpret::isOtherWord(const unsigned int& pSRAMWORD)
{
  if (OTHER_WORD_MACRO(pSRAMWORD))
    return true;
  return false;
}

void Interpret::addTriggerStatus(const unsigned char& pStatus)
{
  if (Basis::debugSet()) {
    debug(std::string(" ") + "addTriggerStatus: " + IntToBin(pStatus) + " at event " + LongIntToStr(_nEvents));
  }
  addEventStatus(__TRG_ERROR);
  tTriggerStatus |= pStatus;
}

void Interpret::addEventStatus(const unsigned short& pStatus)
{
  if ((tEventStatus & pStatus) != pStatus) {  // only add event error code if it hasn't been set
    if (Basis::debugSet()) {
      debug(std::string(" ") + "addEventStatus: " + IntToBin(pStatus) + " at event " + LongIntToStr(_nEvents));
    }
    tEventStatus |= pStatus;
  }
}

void Interpret::histogramTriggerStatusCode()
{
  unsigned int tBitPosition = 0;
  for (unsigned char iErrorCode = tTriggerStatus; iErrorCode != 0; iErrorCode = iErrorCode >> 1) {
    if (iErrorCode & 0x1)
      _triggerStatusCounter[tBitPosition] += 1;
    tBitPosition++;
  }
}

void Interpret::histogramErrorCode()
{
  unsigned int tBitPosition = 0;
  for (unsigned short int iErrorCode = tEventStatus; iErrorCode != 0; iErrorCode = iErrorCode >> 1) {
    if (iErrorCode & 0x1)
      _eventStatusCounter[tBitPosition] += 1;
    tBitPosition++;
  }
}

void Interpret::addServiceRecord(const unsigned char& pSRcode, const unsigned int& pSRcounter)
{
  tServiceRecord |= pSRcode;
  if (pSRcode < __NSERVICERECORDS)
    _serviceRecordCounter[pSRcode] += pSRcounter;
}

void Interpret::addTdcValue(const unsigned int& pTdcValue)
{
  if (pTdcValue < __N_TDC_VALUES)
    _tdcValue[pTdcValue] += 1;
}

void Interpret::addTdcTriggerDistanceValue(const unsigned int& pTdcTriggerDistanceValue)
{
  if (pTdcTriggerDistanceValue < __N_TDC_TRG_DIST_VALUES)
    _tdcTriggerDistance[pTdcTriggerDistanceValue] += 1;
}

void Interpret::allocateHitArray()
{
  debug(std::string("allocateHitArray()"));
  try {
    _hitInfo = new HitInfo[_hitInfoSize];
  } catch (std::bad_alloc& exception) {
    error(std::string("allocateHitArray(): ") + std::string(exception.what()));
    throw;
  }
}

void Interpret::deleteHitArray()
{
  debug(std::string("deleteHitArray()"));
  if (_hitInfo == 0)
    return;
  delete[] _hitInfo;
  _hitInfo = 0;
}

void Interpret::allocateHitBufferArray()
{
  debug(std::string("allocateHitBufferArray()"));
  try {
    _hitBuffer = new HitInfo[__MAXHITBUFFERSIZE];
  } catch (std::bad_alloc& exception) {
    error(std::string("allocateHitBufferArray(): ") + std::string(exception.what()));
    throw;
  }
}

void Interpret::deleteHitBufferArray()
{
  debug(std::string("deleteHitBufferArray()"));
  if (_hitBuffer == 0)
    return;
  delete[] _hitBuffer;
  _hitBuffer = 0;
}

void Interpret::allocateTriggerStatusCounterArray()
{
  debug(std::string("allocateTriggerStatusCounterArray()"));
  try {
    _triggerStatusCounter = new unsigned int[__TRG_N_ERROR_CODES];
  } catch (std::bad_alloc& exception) {
    error(std::string("allocateTriggerStatusCounterArray(): ") + std::string(exception.what()));
  }
}

void Interpret::resetTriggerStatusCounterArray()
{
  for (unsigned int i = 0; i < __TRG_N_ERROR_CODES; ++i)
    _triggerStatusCounter[i] = 0;
}

void Interpret::deleteTriggerStatusCounterArray()
{
  debug(std::string("deleteTriggerStatusCounterArray()"));
  if (_triggerStatusCounter == 0)
    return;
  delete[] _triggerStatusCounter;
  _triggerStatusCounter = 0;
}

void Interpret::allocateEventStatusCounterArray()
{
  debug(std::string("allocateEventStatusCounterArray()"));
  try {
    _eventStatusCounter = new unsigned int[__N_ERROR_CODES];
  } catch (std::bad_alloc& exception) {
    error(std::string("allocateEventStatusCounterArray(): ") + std::string(exception.what()));
  }
}

void Interpret::allocateTdcValueArray()
{
  debug(std::string("allocateTdcValueArray()"));
  try {
    _tdcValue = new unsigned int[__N_TDC_VALUES];
  } catch (std::bad_alloc& exception) {
    error(std::string("allocateTdcValueArray(): ") + std::string(exception.what()));
  }
}

void Interpret::allocateTdcTriggerDistanceArray()
{
  debug(std::string("allocateTdcTriggerDistanceArray()"));
  try {
    _tdcTriggerDistance = new unsigned int[__N_TDC_TRG_DIST_VALUES];
  } catch (std::bad_alloc& exception) {
    error(std::string("allocateTdcTriggerDistanceArray(): ") + std::string(exception.what()));
  }
}

void Interpret::resetEventStatusCounterArray()
{
  for (unsigned int i = 0; i < __N_ERROR_CODES; ++i)
    _eventStatusCounter[i] = 0;
}

void Interpret::resetTdcValueArray()
{
  for (unsigned int i = 0; i < __N_TDC_VALUES; ++i)
    _tdcValue[i] = 0;
}

void Interpret::resetTdcTriggerDistanceArray()
{
  for (unsigned int i = 0; i < __N_TDC_TRG_DIST_VALUES; ++i)
    _tdcTriggerDistance[i] = 0;
}

void Interpret::deleteEventStatusCounterArray()
{
  debug(std::string("deleteEventStatusCounterArray()"));
  if (_eventStatusCounter == 0)
    return;
  delete[] _eventStatusCounter;
  _eventStatusCounter = 0;
}

void Interpret::deleteTdcValueArray()
{
  debug(std::string("deleteTdcValueArray()"));
  if (_tdcValue == 0)
    return;
  delete[] _tdcValue;
  _tdcValue = 0;
}

void Interpret::deleteTdcTriggerDistanceArray()
{
  debug(std::string("deleteTdcTriggerDistanceArray()"));
  if (_tdcTriggerDistance == 0)
    return;
  delete[] _tdcTriggerDistance;
  _tdcTriggerDistance = 0;
}

void Interpret::allocateServiceRecordCounterArray()
{
  debug(std::string("allocateServiceRecordCounterArray()"));
  try {
    _serviceRecordCounter = new unsigned int[__NSERVICERECORDS];
  } catch (std::bad_alloc& exception) {
    error(std::string("allocateServiceRecordCounterArray(): ") + std::string(exception.what()));
  }
}

void Interpret::resetServiceRecordCounterArray()
{
  for (unsigned int i = 0; i < __NSERVICERECORDS; ++i)
    _serviceRecordCounter[i] = 0;
}

void Interpret::deleteServiceRecordCounterArray()
{
  debug(std::string("deleteServiceRecordCounterArray()"));
  if (_serviceRecordCounter == 0)
    return;
  delete[] _serviceRecordCounter;
  _serviceRecordCounter = 0;
}
