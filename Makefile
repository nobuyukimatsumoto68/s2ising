# CXX=g++
CXX=g++-13
CXXFLAGS=-Wall -Wextra -O3 -std=c++23 # -march=native -O3 -fopenmp
INCLUDES=-I"/mnt/hdd_barracuda/opt/eigen/" -I"/mnt/hdd_barracuda/opt/HighFive/include/" -I"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/include/" -I"~/opt/eigen/"
LDFLAGS=-L"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/lib/" -L"/usr/lib/" -L"/usr/local/lib/" -lhdf5 -lgsl -lgslcblas -lm # -lhdf5_cpp
# INCLUDES=-I"/Users/nobuyukimatsumoto/opt/eigen/" -I"/opt/gsl-2.8/"
# LDFLAGS=-L"/usr/local/lib/" -lgsl -lgslcblas -lm



SRCS := $(wildcard *.cc)
OBJS := $(SRCS:%.cc=%.o)

all: $(OBJS)

%.d: %.cc
	$(CXX) $(INCLUDES) -M $< -o $@

include $(SRCS:.cc=.d)

%.o: %.cc
	$(CXX) $< -o $@ $(CXXFLAGS) $(INCLUDES) $(LDFLAGS)


# SRCS := $(wildcard *.cc)
# OBJS := $(SRCS:%.cc=%.o)

# all: $(OBJS)

# %.d: %.cc
# 	$(CXX) $(INCLUDES_CUDA) -M $< -o $@

# include $(SRCS:.cu=.d)

# all:
# 	# $(CXX) test.cc $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o test.o
# 	$(CXX) angle_opt.cc $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o angle_opt.o
# 	$(CXX) gsl.cc $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o gsl.o
