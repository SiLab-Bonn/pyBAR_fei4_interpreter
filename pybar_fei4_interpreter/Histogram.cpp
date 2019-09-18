#include "Histogram.h"

Histogram::Histogram(void)
{
  setSourceFileName("Histogram()");
//  setDebugOutput(true);
  setStandardSettings();
}

Histogram::~Histogram(void)
{
  debug("~Histogram()");
  deleteOccupancyArray();
  deleteTotArray();
  deleteTdcValueArray();
  deleteTdcTriggerDistanceArray();
  deleteRelBcidArray();
  deleteTotPixelArray();
  deleteTdcPixelArray();
  deleteMeanTotArray();
}

void Histogram::setStandardSettings()
{
  info("setStandardSettings()");
    _metaEventIndex = 0;
  _parInfo = 0;
  _lastMetaEventIndex = 0;
  _metaEventIndex = 0;
  _occupancy = 0;
  _relBcid = 0;
  _tot = 0;
  _meanTot = 0;
  _tdcValue = 0;
  _tdcTriggerDistance = 0;
  _totPixel = 0;
  _tdcPixel = 0;
  _NparameterValues = 1;
  _createOccHist = false;
  _createRelBCIDhist = false;
  _createTotHist = false;
  _createMeanTotHist = false;
  _createTdcValueHist = false;
  _createTdcTriggerDistanceHist = false;
  _createTdcPixelHist = false;
  _createTotPixelHist = false;
  _maxTot = 13;
}

void Histogram::createOccupancyHist(bool createOccHist)
{
  _createOccHist = createOccHist;
  if (!createOccHist)
    deleteOccupancyArray();
}

void Histogram::createRelBCIDHist(bool createRelBCIDHist)
{
  _createRelBCIDhist = createRelBCIDHist;
  if (_createRelBCIDhist) {
    allocateRelBcidArray();
    resetRelBcidArray();
  }
  else
    deleteRelBcidArray();
}

void Histogram::createTotHist(bool createTotHist)
{
  _createTotHist = createTotHist;
  if (_createTotHist) {
    allocateTotArray();
    resetTotArray();
  }
  else
    deleteTotArray();
}

void Histogram::createMeanTotHist(bool createMeanTotHist)
{
  _createMeanTotHist = createMeanTotHist;
  if (_createMeanTotHist)
    createOccupancyHist(true);
  else
    deleteMeanTotArray();
}

void Histogram::createTdcValueHist(bool createTdcValueHist)
{
  _createTdcValueHist = createTdcValueHist;
  if (_createTdcValueHist) {
    allocateTdcValueArray();
    resetTdcValueArray();
  }
  else
    deleteTdcValueArray();
}

void Histogram::createTdcTriggerDistanceHist(bool createTdcTriggerDistanceHist)
{
  _createTdcTriggerDistanceHist = createTdcTriggerDistanceHist;
  if (_createTdcTriggerDistanceHist) {
    allocateTdcTriggerDistanceArray();
    resetTdcTriggerDistanceArray();
  }
  else
    deleteTdcTriggerDistanceArray();
}

void Histogram::createTdcPixelHist(bool createTdcPixelHist)
{
  _createTdcPixelHist = createTdcPixelHist;
  if (_createTdcPixelHist) {
    allocateTdcPixelArray();
    resetTdcPixelArray();
  }
  else
    deleteTdcPixelArray();
}

void Histogram::createTotPixelHist(bool createTotPixelHist)
{
  _createTotPixelHist = createTotPixelHist;
  if (_createTotPixelHist) {
    allocateTotPixelArray();
    resetTotPixelArray();
  }
  else
    deleteTotPixelArray();
}

void Histogram::setMaxTot(const unsigned int& rMaxTot)
{
  _maxTot = rMaxTot;
}

void Histogram::addHits(HitInfo*& rHitInfo, const unsigned int& rNhits)
{
  debug("addHits()");
  for (unsigned int i = 0; i<rNhits; ++i) {
    if (((rHitInfo[i].event_status & __NO_HIT) == __NO_HIT) || (rHitInfo[i].column == 0) || (rHitInfo[i].row == 0))  // ignore virtual hits
      continue;
    unsigned short tColumnIndex = rHitInfo[i].column-1;
    if (tColumnIndex > RAW_DATA_MAX_COLUMN-1)
      throw std::out_of_range("Column index out of range.");
    unsigned int tRowIndex = rHitInfo[i].row-1;
    if (tRowIndex > RAW_DATA_MAX_ROW-1)
      throw std::out_of_range("Row index out of range.");
    unsigned int tTot = rHitInfo[i].tot;
    if (tTot > __MAXHITTOT)
      throw std::out_of_range("ToT index out of range.");
    unsigned int tTdc = rHitInfo[i].TDC;
    if (tTdc >= __N_TDC_VALUES)
      throw std::out_of_range("TDC value " + IntToStr(tTdc) + " index out of range.");
    unsigned int tTdcTriggerDistance = rHitInfo[i].TDC_trigger_distance;
    if (tTdcTriggerDistance >= __N_TDC_TRG_DIST_VALUES)
      throw std::out_of_range("TDC trigger distance " + IntToStr(tTdc) + " index out of range.");
    unsigned int tRelBcid = rHitInfo[i].relative_BCID;
    if (tRelBcid >= __MAXBCID)
      throw std::out_of_range("Relative BCID index out of range.");

    unsigned int tParIndex = getParIndex(rHitInfo[i].event_number);

    if (tParIndex >= getNparameters()) {
      error("addHits: tParIndex "+IntToStr(tParIndex)+"\t> "+IntToStr(_NparameterValues));
      throw std::out_of_range("Parameter index out of range.");
    }
    if (_createOccHist) {
      if (tTot <= _maxTot) {
        if (_occupancy!=0) {
          _occupancy[(size_t)tColumnIndex + (size_t)tRowIndex * (size_t)RAW_DATA_MAX_COLUMN + (size_t)tParIndex * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] += 1;
        } else {
          throw std::runtime_error("Occupancy array not initialized. Set scan parameter first!.");
        }
        if (_createMeanTotHist) {
          if (_meanTot!=0) {
            float tOccupancy = (float)_occupancy[(size_t)tColumnIndex + (size_t)tRowIndex * (size_t)RAW_DATA_MAX_COLUMN + (size_t)tParIndex * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW];
            float tMeanTot = _meanTot[(size_t)tColumnIndex + (size_t)tRowIndex * (size_t)RAW_DATA_MAX_COLUMN + (size_t)tParIndex * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW];
            if (tMeanTot != tMeanTot)  // check for NAN, _meanTot initialized with NAN
              tMeanTot = 0.0;
            _meanTot[(size_t)tColumnIndex + (size_t)tRowIndex * (size_t)RAW_DATA_MAX_COLUMN + (size_t)tParIndex * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] = (tMeanTot * (tOccupancy - 1) + tTot) / tOccupancy;
          } else {
            throw std::runtime_error("Mean ToT array not initialized. Set scan parameter first!.");
          }
        }
      }
    }
    if (_createRelBCIDhist)
      if (tTot <= _maxTot)
        _relBcid[tRelBcid] += 1;
    if (_createTotHist)
      if (tTot <= __MAXHITTOT)
        _tot[tTot] += 1;
    if (((rHitInfo[i].event_status & (__NO_HIT | __MORE_THAN_ONE_HIT | __MORE_THAN_ONE_TDC_WORD | __TDC_INVALID)) == 0) && ((rHitInfo[i].event_status & __TDC_WORD) == __TDC_WORD)) {  // get TDC values from single hit events and unambiguous TDC word
      if (_createTdcValueHist) {  // get TDC values from single hit events
        _tdcValue[tTdc] += 1;
      }
      if (_createTdcTriggerDistanceHist) {  // get TDC values from single hit events
        _tdcTriggerDistance[tTdcTriggerDistance] += 1;
      }
      if (_createTdcPixelHist) {
        if (_tdcPixel != 0) {
          if (tTdc >= __N_TDC_VALUES) {  // get TDC values from single hit events
            info("TDC value out of range:" + IntToStr(tTdc) + ">" + IntToStr(__N_TDC_VALUES));
            tTdc = 0;
          }
          _tdcPixel[(size_t)tColumnIndex + (size_t)tRowIndex * (size_t)RAW_DATA_MAX_COLUMN + (size_t)tTdc * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] += 1;
        } else {
          throw std::runtime_error("Output TDC pixel array array not set.");
        }
      }
    }
    if (_createTotPixelHist) {
      if (tTot <= __MAXHITTOT)
        _totPixel[(size_t)tColumnIndex + (size_t)tRowIndex * (size_t)RAW_DATA_MAX_COLUMN + (size_t)tTot * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] += 1;
    }
  }
}

void Histogram::addClusterSeedHits(ClusterInfo*& rClusterInfo, const unsigned int& rNcluster)
{
  if (Basis::debugSet())
    debug("addClusterSeedHits(...,rNcluster="+IntToStr(rNcluster)+")");
  for (unsigned int i = 0; i<rNcluster; ++i) {
    unsigned short tColumnIndex = rClusterInfo[i].seed_column-1;
    if (tColumnIndex > RAW_DATA_MAX_COLUMN-1)
      throw std::out_of_range("Column index out of range.");
    unsigned int tRowIndex = rClusterInfo[i].seed_row-1;
    if (tRowIndex > RAW_DATA_MAX_ROW-1)
      throw std::out_of_range("Row index out of range.");

    unsigned int tParIndex = getParIndex(rClusterInfo[i].event_number);

    if (tParIndex >= getNparameters()) {
      error("addClusterSeedHits: tParIndex "+IntToStr(tParIndex)+"\t> "+IntToStr(_NparameterValues));
      throw std::out_of_range("Parameter index out of range.");
    }
    if (_createOccHist) {
      if (_occupancy!=0)
        _occupancy[(size_t)tColumnIndex + (size_t)tRowIndex * (size_t)RAW_DATA_MAX_COLUMN + (size_t)tParIndex * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] += 1;
      else
        throw std::runtime_error("Occupancy array not initialized. Set scan parameter first!.");
    }
  }
}

unsigned int Histogram::getParIndex(int64_t& rEventNumber)
{
  if (_parInfo == 0)
    return 0;
  for (uint64_t i=_lastMetaEventIndex; i<_nMetaEventIndexLength-1; ++i) {
    if (_metaEventIndex[i+1] > (uint64_t) rEventNumber || _metaEventIndex[i+1] < _metaEventIndex[i]) {  // second case: meta event data not set yet (std value = 0), event number has to increase
      _lastMetaEventIndex = i;
      if (i < _nParInfoLength) {
        return _parInfo[i];
      } else {
        error("Scan parameter index " + LongIntToStr(i) + " out of range");
        throw std::out_of_range("Scan parameter index out of range.");
      }
    }
  }
  if (_metaEventIndex[_nMetaEventIndexLength-1] <= (uint64_t) rEventNumber)  // last read outs
    return _parInfo[_nMetaEventIndexLength-1];
  error("getScanParameter: Correlation issues at event "+LongIntToStr(rEventNumber)+"\n_metaEventIndex[_nMetaEventIndexLength-1] "+LongIntToStr(_metaEventIndex[_nMetaEventIndexLength-1])+"\n_lastMetaEventIndex "+LongIntToStr(_lastMetaEventIndex));
  throw std::logic_error("Event parameter correlation issues.");
  return 0;
}

void Histogram::addScanParameter(int*& rParInfo, const unsigned int& rNparInfoLength)
{
  debug("addScanParameter");
  _nParInfoLength = rNparInfoLength;
  _parInfo = rParInfo;

  std::vector<unsigned int> tParameterValues;

  for(unsigned int i = 0; i < _nParInfoLength; ++i)
    tParameterValues.push_back(_parInfo[i]);

  std::sort(tParameterValues.begin(), tParameterValues.end());  // sort from lowest to highest value
  std::set<unsigned int> tSet(tParameterValues.begin(), tParameterValues.end());  // delete duplicates
  tParameterValues.assign(tSet.begin(), tSet.end() );

  for(unsigned int i = 0; i < tParameterValues.size(); ++i)
    _parameterValues[tParameterValues[i]] = i;

  _NparameterValues = (unsigned int) tSet.size();

  if (_createOccHist) {
    allocateOccupancyArray();
    resetOccupancyArray();
    if (_createMeanTotHist) {
      allocateMeanTotArray();
      resetMeanTotArray();
    }
  }
  if (Basis::debugSet()) {
    for(unsigned int i=0; i<_nParInfoLength; ++i)
      std::cout<<"index "<<i<<"\t parameter index "<<_parInfo[i]<<"\n";
  }
}

void Histogram::addMetaEventIndex(uint64_t*& rMetaEventIndex, const unsigned int& rNmetaEventIndexLength)
{
  debug("addMetaEventIndex()");
  _nMetaEventIndexLength = rNmetaEventIndexLength;
  _metaEventIndex = rMetaEventIndex;
  if (Basis::debugSet())
    for (unsigned int i=0; i<_nMetaEventIndexLength; ++i)
      std::cout<<"index "<<i<<"\t event number "<<_metaEventIndex[i]<<"\n";
}

void Histogram::allocateOccupancyArray()
{
  debug("allocateOccupancyArray() with "+IntToStr(getNparameters())+" parameters");
  deleteOccupancyArray();
  try {
    _occupancy = new unsigned int[(size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW * (size_t)getNparameters()];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateOccupancyArray: ")+std::string(exception.what()));
  }
}

void Histogram::deleteOccupancyArray()
{
  debug("deleteOccupancyArray()");
  if (_occupancy != 0)
    delete[] _occupancy;
  _occupancy = 0;
}

void Histogram::resetOccupancyArray()
{
  info("resetOccupancyArray()");
  if (_occupancy != 0) {
    for (size_t i = 0; i < RAW_DATA_MAX_COLUMN; i++)
      for (size_t j = 0; j < RAW_DATA_MAX_ROW; j++)
        for (size_t k = 0; k < getNparameters(); k++)
          _occupancy[i + j * (size_t)RAW_DATA_MAX_COLUMN + k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] = 0;
  }
}

void Histogram::allocateMeanTotArray()
{
  debug("allocateMeanTotArray() with "+IntToStr(getNparameters())+" parameters");
  deleteMeanTotArray();
  try {
    _meanTot = new float[(size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW * (size_t)getNparameters()];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateMeanTotArray: ")+std::string(exception.what()));
  }
}

void Histogram::deleteMeanTotArray()
{
  debug("deleteMeanTotArray()");
  if (_meanTot != 0) {
    delete[] _meanTot;
  }
  _meanTot = 0;
}

void Histogram::resetMeanTotArray()
{
  info("resetMeanTotArray()");
  if (_meanTot != 0) {
    for (size_t i = 0; i < RAW_DATA_MAX_COLUMN; i++)
      for (size_t j = 0; j < RAW_DATA_MAX_ROW; j++)
        for (size_t k = 0; k < getNparameters(); k++)
          _meanTot[i + j * (size_t)RAW_DATA_MAX_COLUMN + k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] = NAN;
  }
}

void Histogram::resetTdcPixelArray()
{
  info("resetTdcPixelArray()");
  if (_createTdcPixelHist) {
    if (_tdcPixel != 0) {
      for (size_t i = 0; i < RAW_DATA_MAX_COLUMN; i++)
        for (size_t j = 0; j < RAW_DATA_MAX_ROW; j++)
          for (size_t k = 0; k < __N_TDC_VALUES; k++)
            _tdcPixel[i + j * (size_t)RAW_DATA_MAX_COLUMN + k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] = 0;
    } else {
      throw std::runtime_error("Output TDC pixel array array not set.");
    }
  }
}

void Histogram::resetTotPixelArray()
{
  info("resetTotPixelArray()");
  if (_createTotPixelHist) {
    if (_totPixel != 0) {
      for (size_t i = 0; i < RAW_DATA_MAX_COLUMN; i++) {
        for (size_t j = 0; j < RAW_DATA_MAX_ROW; j++) {
          for (size_t k = 0; k <= __MAXHITTOT; k++) {
            _totPixel[i + j * (size_t)RAW_DATA_MAX_COLUMN + k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW] = 0;
          }
        }
      }
    } else {
      throw std::runtime_error("Output ToT pixel array array not set.");
    }
  }
}

void Histogram::allocateTotArray()
{
  debug("allocateTotArray()");
  deleteTotArray();
  try {
    _tot = new unsigned int[__MAXHITTOT + 1];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateTotArray: ")+std::string(exception.what()));
  }
}

void Histogram::allocateTdcValueArray()
{
  debug("allocateTdcValueArray()");
  deleteTdcValueArray();
  try {
    _tdcValue = new unsigned int[__N_TDC_VALUES];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateTdcValueArray: ")+std::string(exception.what()));
  }
}

void Histogram::allocateTdcTriggerDistanceArray()
{
  debug("allocateTdcTriggerDistanceArray()");
  deleteTdcTriggerDistanceArray();
  try {
    _tdcTriggerDistance = new unsigned int[__N_TDC_TRG_DIST_VALUES];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateTdcTriggerDistanceArray: ")+std::string(exception.what()));
  }
}

void Histogram::resetTotArray()
{
  info("resetTotArray()");
  if (_createTotHist) {
    if (_tot != 0) {
      for (size_t i = 0; i <= __MAXHITTOT; i++) {
        _tot[i] = 0;
      }
    }
  }
}

void Histogram::resetTdcValueArray()
{
  info("resetTdcValueArray()");
  if (_createTdcValueHist) {
    if (_tdcValue != 0) {
      for (size_t i = 0; i < __N_TDC_VALUES; i++) {
        _tdcValue[i] = 0;
      }
    }
  }
}

void Histogram::resetTdcTriggerDistanceArray()
{
  info("resetTdcTriggerDistanceArray()");
  if (_createTdcTriggerDistanceHist) {
    if (_tdcTriggerDistance != 0) {
      for (size_t i = 0; i < __N_TDC_TRG_DIST_VALUES; i++) {
        _tdcTriggerDistance[i] = 0;
      }
    }
  }
}

void Histogram::deleteTotArray()
{
  debug("deleteTotArray()");
  if (_tot != 0) {
    delete[] _tot;
  }
  _tot = 0;
}

void Histogram::deleteTdcValueArray()
{
  debug("deleteTdcValueArray()");
  if (_tdcValue != 0) {
    delete[] _tdcValue;
  }
  _tdcValue = 0;
}

void Histogram::deleteTdcTriggerDistanceArray()
{
  debug("deleteTdcTriggerDistanceArray()");
  if (_tdcTriggerDistance != 0) {
    delete[] _tdcTriggerDistance;
  }
  _tdcTriggerDistance = 0;
}

void Histogram::allocateRelBcidArray()
{
  debug("allocateRelBcidArray");
  deleteRelBcidArray();
  try {
    _relBcid = new unsigned int[__MAXBCID];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateRelBcidArray: ")+std::string(exception.what()));
  }
}

void Histogram::resetRelBcidArray()
{
  info("resetRelBcidArray()");
  if (_createRelBCIDhist) {
    if (_relBcid != 0) {
      for (size_t i = 0; i < __MAXBCID; i++)
        _relBcid[i] = 0;
    }
  }
}

void Histogram::deleteRelBcidArray()
{
  debug("deleteRelBcidArray");
  if (_relBcid != 0)
    delete[] _relBcid;
  _relBcid = 0;
}

void Histogram::allocateTotPixelArray()
{
  debug("allocateTotPixelArray()");
  deleteTotPixelArray();
  try {
    _totPixel = new unsigned short[(size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW * ((size_t)__MAXHITTOT + 1)];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateTotPixelArray: ")+std::string(exception.what()));
  }
}

void Histogram::allocateTdcPixelArray()
{
  debug("allocateTdcPixelArray()");
  deleteTdcPixelArray();
  try {
    _tdcPixel = new unsigned short[(size_t) RAW_DATA_MAX_COLUMN * (size_t) RAW_DATA_MAX_ROW * (size_t)__N_TDC_VALUES];
  } catch(std::bad_alloc& exception) {
    error(std::string("allocateTdcPixelArray: ")+std::string(exception.what()));
  }
}

void Histogram::deleteTotPixelArray()
{
  debug("deleteTotPixelArray");
  if (_totPixel != 0)
    delete[] _totPixel;
  _totPixel = 0;
}

void Histogram::deleteTdcPixelArray()
{
  debug("deleteTdcPixelArray");
  if (_tdcPixel != 0)
    delete[] _tdcPixel;
  _tdcPixel = 0;
}

unsigned int Histogram::getNparameters()
{
  return _NparameterValues;
}

void Histogram::getOccupancy(unsigned int& rNparameterValues, unsigned int*& rOccupancy, bool copy)
{
  debug("getOccupancy(...)");
  if (copy) {
    size_t tArrayLength = (size_t)(RAW_DATA_MAX_COLUMN - 1) + (size_t)(RAW_DATA_MAX_ROW - 1) * (size_t)RAW_DATA_MAX_COLUMN + (size_t)(_NparameterValues - 1) * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW + 1;
    std::copy(_occupancy, _occupancy + tArrayLength, rOccupancy);
  } else {
    rOccupancy = _occupancy;
  }
  rNparameterValues = _NparameterValues;
}

void Histogram::getTotHist(unsigned int*& rTotHist, bool copy)
{
  debug("getTotHist(...)");
  if (copy) {
    std::copy(_tot, _tot + __MAXHITTOT + 1, rTotHist);
  } else {
    rTotHist = _tot;
  }
}

void Histogram::getMeanTot(unsigned int& rNparameterValues, float*& rMeanTot, bool copy)
{
  debug("getMeanTot(...)");
  if (copy) {
    size_t tArrayLength = (size_t)(RAW_DATA_MAX_COLUMN - 1) + (size_t)(RAW_DATA_MAX_ROW - 1) * (size_t)RAW_DATA_MAX_COLUMN + (size_t)(_NparameterValues - 1) * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW + 1;
    std::copy(_meanTot, _meanTot + tArrayLength, rMeanTot);
  } else {
    rMeanTot = _meanTot;
  }
  rNparameterValues = _NparameterValues;
}

void Histogram::getTdcValuesHist(unsigned int*& rTdcValueHist, bool copy)
{
  debug("getTdcValuesHist(...)");
  if (copy) {
    std::copy(_tdcValue, _tdcValue + __N_TDC_VALUES, rTdcValueHist);
  } else {
    rTdcValueHist = _tdcValue;
  }
}

void Histogram::getTdcTriggerDistancesHist(unsigned int*& rTdcTriggerDistanceHist, bool copy)
{
  debug("getTdcTriggerDistancesHist(...)");
  if (copy) {
    std::copy(_tdcTriggerDistance, _tdcTriggerDistance + __N_TDC_TRG_DIST_VALUES, rTdcTriggerDistanceHist);
  } else {
    rTdcTriggerDistanceHist = _tdcTriggerDistance;
  }
}

void Histogram::getRelBcidHist(unsigned int*& rRelBcidHist, bool copy)
{
  debug("getRelBcidHist(...)");
  if (copy) {
    std::copy(_relBcid, _relBcid + __MAXBCID, rRelBcidHist);
  } else {
    rRelBcidHist = _relBcid;
  }
}

void Histogram::getTotPixelHist(unsigned short*& rTotPixelHist, bool copy)
{
  debug("getTotPixelHist(...)");
  if (copy) {
    std::copy(_totPixel, _totPixel + __MAXHITTOT + 1, rTotPixelHist);
  } else {
    rTotPixelHist = _totPixel;
  }
}

void Histogram::getTdcPixelHist(unsigned short*& rTdcPixelHist, bool copy)
{
  debug("getTdcPixelHist(...)");
  if (copy) {
    std::copy(_tdcPixel, _tdcPixel + __N_TDC_VALUES, rTdcPixelHist);
  } else {
    rTdcPixelHist = _tdcPixel;
  }
}

void Histogram::calculateThresholdScanArrays(double rMuArray[], double rSigmaArray[], const unsigned int& rMaxInjections, const unsigned int& min_parameter, const unsigned int& max_parameter)
{
  debug("calculateThresholdScanArrays(...)");
  // fast algorithm from M. Mertens, PhD thesis, FZ Juelich 2010

  if (_occupancy==0)
    throw std::runtime_error("Occupancy array not initialized. Set scan parameter first!.");

  if (_NparameterValues<2)  // a minimum number of different scans is needed
    return;

  unsigned int q_min = min_parameter;
  unsigned int q_max = max_parameter;
  unsigned int A = rMaxInjections;
  double d = ((double)q_max - (double)q_min)/(double)(getNparameters()-1);

  for(unsigned int i=0; i<RAW_DATA_MAX_COLUMN; ++i) {
    for(unsigned int j=0; j<RAW_DATA_MAX_ROW; ++j) {
      unsigned int M = 0;

      for(unsigned int k=0; k<getNparameters(); ++k) {
        M += _occupancy[(size_t)i + (size_t)j * (size_t)RAW_DATA_MAX_COLUMN  + (size_t)k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW];
      }
      double threshold = (double)q_max+d/2 - d*(double)M/(double)A;
      rMuArray[i+j*RAW_DATA_MAX_COLUMN] = threshold;

      double mu1 = 0;
      double mu2 = 0;
//      bool do_math = true;
      for(unsigned int k=0; k<getNparameters(); ++k) {
        if(((double)k*d+(double)q_min) < threshold) {
          mu1 += (double)_occupancy[(size_t)i + (size_t)j * (size_t)RAW_DATA_MAX_COLUMN  + (size_t)k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW];
        }else{
//          if(do_math) {
//          double G = 0.5 - (threshold - ((double)k*d+(double)q_min))/d;
//            mu1 += G*(double)_occupancy[(size_t)i + (size_t)j * (size_t)RAW_DATA_MAX_COLUMN  + (size_t)k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW];
//            mu2 += (1-G)*((double)A-(double)_occupancy[(size_t)i + (size_t)j * (size_t)RAW_DATA_MAX_COLUMN  + (size_t)k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW]);
//            do_math = false;
//          }else{
            mu2 += ((double)A-(double)_occupancy[(size_t)i + (size_t)j * (size_t)RAW_DATA_MAX_COLUMN  + (size_t)k * (size_t)RAW_DATA_MAX_COLUMN * (size_t)RAW_DATA_MAX_ROW]);
//          }
        }
      }
      double noise = d*(mu1+mu2)/(double)A*sqrt(3.14159265358979323846/2);
      rSigmaArray[i+j*RAW_DATA_MAX_COLUMN] = noise;
    }
  }
}

void Histogram::setNoScanParameter()
{
  debug("setNoScanParameter()");
  deleteOccupancyArray();
  _NparameterValues = 1;
  allocateOccupancyArray();
  resetOccupancyArray();
  allocateMeanTotArray();
  resetMeanTotArray();
}

void Histogram::reset()
{
  info("reset()");
  resetOccupancyArray();
  resetTotArray();
  resetMeanTotArray();
  resetTdcValueArray();
  resetTdcTriggerDistanceArray();
  resetTotPixelArray();
  resetTdcPixelArray();
  resetRelBcidArray();
  _parInfo = 0;
}
