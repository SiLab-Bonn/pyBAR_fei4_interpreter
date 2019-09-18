
# pyBAR FE-I4 Interpreter [![Build Status](https://travis-ci.org/SiLab-Bonn/pyBAR_fei4_interpreter.svg?branch=master)](https://travis-ci.org/SiLab-Bonn/pyBAR_fei4_interpreter) [![Build status](https://ci.appveyor.com/api/projects/status/247ws3a200c8tnq1?svg=true)](https://ci.appveyor.com/project/laborleben/pybar-fei4-interpreter)

pyBAR_fei4_interpreter - An ATLAS FE-I4 raw data interpreter in Python and C++

Interpreter for ATLAS FE-I4A/B raw data for the readout framework pyBAR. It also provides histogramming functions. The interpreter is written in C++ to achieve a high throughput.

## Installation

The following packages are required for pyBAR's ATLAS FE-I4 interpreter:
  ```
  cython numpy pytables
  ```

## Usage
```
from pybar_fei4_interpreter.data_interpreter import PyDataInterpreter
interpreter = PyDataInterpreter()  # Initialize interpretation module
raw_data = np.array([73175087, 73044495, 73058863, 73194895, 73197919, 73093151], np.uint32)  # Some raw data to interpret
interpreter.interpret_raw_data(raw_data)  # Start the raw data interpretation
print interpreter.get_hits()  # Print the hits in the raw data
```

Also take a look at the example folder.

## Support

Please use GitHub's [issue tracker](https://github.com/SiLab-Bonn/pyBAR_fei4_interpreter/issues) for bug reports/feature requests/questions.

*For CERN users*: Feel free to subscribe to the [pyBAR mailing list](https://e-groups.cern.ch/e-groups/EgroupsSubscription.do?egroupName=pybar-devel).
