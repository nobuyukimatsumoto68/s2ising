CXX=g++
CXXFLAGS=-Wall -Wextra -std=c++23 -march=native -O3
INCLUDES=-I"/mnt/hdd_barracuda/opt/eigen/" -I"/mnt/hdd_barracuda/opt/HighFive/include/" -I"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/include/" -I"~/opt/eigen/"
LDFLAGS=-L"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/lib/" -L"/usr/lib/" -L"/usr/local/lib/" -lhdf5 # -lhdf5_cpp
# INCLUDES=-I"/Users/nobuyukimatsumoto/opt/eigen/"
# LDFLAGS=

all:
	$(CXX) geometry.cc $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o geometry.o
