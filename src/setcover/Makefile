#CXXFLAGS := -I. -O3 -g -std=c++14 -mavx2

#setcover: setcover.cpp
#	$(CXX) -o setcover setcover.cpp $(CXXFLAGS) 

CXXFLAGS := -I. -O3 -g -std=c++14
UNAME_S := $(shell uname -s)
ARCH := $(shell uname -m)

ifeq ($(ARCH),x86_64)
    CXXFLAGS += -mavx2
endif

setcover: setcover.cpp
	$(CXX) -o setcover setcover.cpp $(CXXFLAGS)
