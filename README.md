# Note

This is a personal **unofficial** repository of Cheetah.
The official Website is http://www.desy.de/~barty/cheetah/,
and the repository is https://github.com/antonbarty/cheetah/.

Here (eiger-zmq branch) functions for Eiger streaming data specific for SPring-8 MX beamlines are being developed.
We are using Cheetah's peakfind functions and expose them to Python using Boost.Python.
We started development by merging eiger branch of [biochem-fan's forked repository](https://github.com/biochem-fan/cheetah).

Everything related to this branch is in eiger-zmq/ directory.
Peakfind functions in source/libcheetah/src/peakfinders.cpp and related header files in original code were moved to eiger-zmq/.
They are exported to python by eiger-zmq/boost_python/cheetah_ext.cpp.

# Building
You need Python and its packages (PyUblas and numpy). Replace XXX and followings to match your environment.
```
cd eiger-zmq/
g++ -fPIC -O3 -c peakfinders.cpp
cd boost_python
g++ -O3 -c cheetah_ext.cpp -shared -fPIC  -I../ -I/XXX/python/python-2.7.6/include/python2.7/ \
 -I/XXX/python/python-2.7.6/lib/python2.7/site-packages/PyUblas-2013.1-py2.7-linux-x86_64.egg/pyublas/include \
 -I/XXX/python/python-2.7.6/lib/python2.7/site-packages/numpy/core/include
g++ -shared cheetah_ext.o  -o cheetah_ext.so -L/oys/xtal/cctbx/build/lib -lboost_python ../peakfinders.o 
```

# Usage
Currently it is designed to use at SPring-8 BL32XU only, but it may easily be customized to other environment.

cheetah_client.py receives frames from Eiger via Zmq and find peaks. The results are sent to other program via Zmq, which will be available in yamtbx.
