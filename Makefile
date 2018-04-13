CCC         = g++
OPT         = -O0 -g -DDEBUG
CFLAGS="-I/usr/include/gsl"
LDFLAGS=-lgsl -lgslcblas

bin	= taun mn_test fh_test
lib	= libTAUn.so
base	= taun_vars TauData

src	= $(base:=.cc)
inc	= $(base:=.h)
obj	= $(base:=.o)

DEBUG	  ?= -O3 -DDEBUG
CXXFLAGS  += $(DEBUG) -fPIC
SHFLAGS   ?= -shared -Wl,-soname,

ROOTCINT  = rootcint
ROOTFLAGS ?= $(shell root-config --cflags)
ROOTLIBS  ?= $(shell root-config --glibs)
LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs) 


all: bin lib
bin: $(bin)
lib: $(lib)

clean:
	rm -f *.o *_cc.d *_cc.so *_C.d *_C.so *Dict.*
veryclean: clean
	rm -f $(bin) $(lib)


$(bin): %: %.o $(obj) $(lib:.so=_Dict.o)
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTLIBS) $^

$(lib): %.so: $(obj) %_Dict.o
	$(CXX) -o $@ $(CXXFLAGS) $(SHFLAGS)$@ $^


libTAUn_Dict.cc: $(inc) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $^
%_Dict.cc: %.h %_LinkDef.h
	$(ROOTCINT) -f $@ -c -p $^

$(obj):	%.o: %.h
TauData.o: $(inc)
taun.o:	TauData.h
taun:	libTAUn_Dict.o

detsigave: detsigave.cc TauData.cc
	$(CCC) $(OPT) -o detsigave  $(ROOTFLAGS) $^ -L, -l TreePlayer $(ROOTLIBS)

multi_series: multi_series.cc TauData.cc
	$(CCC) $(OPT) -o multi_series  $(ROOTFLAGS) $^ -L, -l TreePlayer $(ROOTLIBS)

multi_seriesv: multi_seriesv.cc TauData.cc
	$(CCC) $(OPT) -o multi_seriesv  $(ROOTFLAGS) $^ -L, -l TreePlayer $(ROOTLIBS)

monitoring: monitoring.cc TauData.cc
	$(CCC) $(OPT) -o monitoring  $(ROOTFLAGS) $^ -L, -l TreePlayer $(ROOTLIBS)

clean:
	rm -f *.o *~ */*.e4* */*.pe4* */*.po4* detsigave

%.o:%.cc
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $(ROOTFLAGS) $<
%.o:%.C
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $(ROOTFLAGS) $<
