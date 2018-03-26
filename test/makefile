#Completely stolen.
CPPFLAGS=-I $(BOOST_INC) \
         -I $(CANVAS_INC) \
         -I $(CETLIB_INC) \
         -I $(FHICLCPP_INC) \
         -I $(GALLERY_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
         -I $(NUSIMDATA_INC) \
	 -I $(UBOONECODE_INC) \
         -I $(ROOT_INC)

UBOONECODE_LIB="/uboone/app/users/mastbaum/uboonecode-v06_26_01_12/localProducts_larsoft_v06_26_01_10_e10_prof/uboonecode/v06_26_01_12/slf6.x86_64.e10.prof/lib"

CXXFLAGS=-std=c++14 -Wall
CXX=g++
LDFLAGS=$$(root-config --libs) \
        -L $(CANVAS_LIB) -l canvas_Utilities -l canvas_Persistency_Common -l canvas_Persistency_Provenance \
        -L $(CETLIB_LIB) -l cetlib \
        -L $(GALLERY_LIB) -l gallery \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
        -L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase \
        -L $(UBOONECODE_LIB) -l EventWeight 


UNAME := $(shell uname -s)

all: check_weights arborist

check_weights: check_weights.cc
	 @echo Building $(EXEC)
	 @$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

arborist: arborist.cpp
	 @echo Building $(EXEC)
	 @$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

