
# pyBAR_fei4_interpreter [![Build Status](https://travis-ci.org/SiLab-Bonn/pyBAR_fei4_interpreter.svg?branch=master)](https://travis-ci.org/SiLab-Bonn/pyBAR_fei4_interpreter) [![Build Status](https://ci.appveyor.com/api/projects/status/github/SiLab-Bonn/pyBAR_fei4_interpreter?svg=true)](https://ci.appveyor.com/project/DavidLP/pyBAR_fei4_interpreter-71xwl)

pyBAR_fei4_interpreter - An ATLAS FE-I4 raw data interpreter in Python and C++

This package can be used to interpred raw data from the ATLAS FE-I4 taken with the readout framework pyBAR. It also contains histogramming functions and interpretation tools. The interpretation takes place in fast C++ code to increase the speed.

## Installation

The following packages are required for pyBAR's ATLAS FE-I4 interpreter:
  ```
  numpy cython tables
  ```

## Usage
```
from pybar_fei4_interpreter.data_interpreter import PyDataInterpreter
interpreter = PyDataInterpreter()  # Initialize interpretation module
raw_data = np.array([73175087, 73044495, 73058863, 73194895, 73197919, 73093151], np.uint32)  # Some raw data to interpret
interpreter.interpret_raw_data(raw_data)  # Start the raw data interpretation
print interpreter.get_hits()  # Print the hits in the raw data
```

Als take a look at the example folder.
## Support

To subscribe to the pyBAR mailing list, click [here](https://e-groups.cern.ch/e-groups/EgroupsSubscription.do?egroupName=pybar-devel). Please ask questions on the pyBAR mailing list [pybar-devel@cern.ch](mailto:pybar-devel@cern.ch?subject=bug%20report%20%2F%20feature%20request) (subscription required) or file a new bug report / feature request [here](https://github.com/SiLab-Bonn/pyBAR_fei4_interpreter/issues/new).
