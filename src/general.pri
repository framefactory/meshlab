# this is the common include for anything compiled inside the meshlab pro


# This is the main coord type inside meshlab
# it can be double or float according to user needs.
DEFINES += MESHLAB_SCALAR=float

VCGDIR   = ../../../vcglib
CONFIG(system_eigen3): EIGENDIR = /usr/include/eigen3
!CONFIG(system_eigen3):EIGENDIR = $$VCGDIR/eigenlib
!CONFIG(system_glew):  GLEWDIR = ../external/glew-2.1.0

CONFIG += c++11

macx:QMAKE_CXXFLAGS += -Wno-inconsistent-missing-override
macx:CONFIG(release, debug|release):QMAKE_CXXFLAGS += -O3 -DNDEBUG
macx:CONFIG(debug, debug|release):QMAKE_CXXFLAGS += -O0 -g

# uncomment these three lines for using the latest clang compiler on OSX to use openmp
# using brew install clang from llvm and libomp
macx:QMAKE_CXX = /usr/local/opt/llvm/bin/clang++
macx:QMAKE_CXXFLAGS += -fopenmp -I/usr/local/opt/llvm/include
macx:QMAKE_LFLAGS += -L/usr/local/opt/llvm/lib -lomp


MACLIBDIR = ../../external/lib/macx64

# the following line is needed to avoid mismatch between
# the awful min/max macros of windows and the limits max
win32:DEFINES += NOMINMAX

linux-g++:QMAKE_CXXFLAGS+=-Wno-unknown-pragmas

# Set up library search paths
linux:QMAKE_LFLAGS+=-L$$PWD/external/lib/linux
linux-g++:QMAKE_LFLAGS+=-L$$PWD/external/lib/linux-g++
linux-g++-32:QMAKE_LFLAGS += -L$$PWD/external/lib/linux-g++-32
linux-g++-64:QMAKE_LFLAGS += -L$$PWD/external/lib/linux-g++-64
