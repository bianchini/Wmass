INC_DIR_CMSSW   = /home/users/bianchini/CMSSW_8_0_25/src/Wmass/src/
INC_DIR         = ./
CXX		= g++
LD		= g++
CXXFLAGS	= -O2 -ggdb -std=gnu++0x -Wall
LDFLAGS		= -lz -lm
SOFLAGS		= -fPIC -shared 
SHELL		= bash
###
SrcSuf        = cc
HeadSuf       = h
ObjSuf        = o
DepSuf        = d
DllSuf        = so
StatSuf       = a

CppTestFiles=$(wildcard test/*.$(SrcSuf))
Packages=$(patsubst test/%.$(SrcSuf),%,$(CppTestFiles) )

CppSrcFiles=$(wildcard src/*.$(SrcSuf))
Objects=$(patsubst src/%.$(SrcSuf),%,$(CppSrcFiles))

LibName = Wmass

### ----- OPTIONS ABOVE ----- ####

InfoLine = compiling $(1)

BASEDIR=$(shell pwd)
BINDIR=$(BASEDIR)/bin
LIBDIR=$(BASEDIR)/libs
SRCDIR = $(BASEDIR)/src
HDIR = $(BASEDIR)/interface
#MEM  = -lMathMore -lLHAPDF -lTreePlayer -lRooFit -lRooFitCore -lHtml -lMinuit

### Main Target, first
.PHONY: all
all: info $(Packages) | $(BINDIR)

CXXFLAGS	+=`root-config --cflags`
LDFLAGS 	+=`root-config --libs`

#CXXFLAGS        += $(MEM)

BINOBJ	=$(patsubst %,$(BINDIR)/%.$(ObjSuf),$(Objects) )
SRCFILES=$(patsubst %,$(SRCDIR)/%.$(SrcSuf),$(Objects) )
HFILES	=$(patsubst %,$(HDIR)/%.$(HeadSuf),$(Objects) )
StatLib		=$(BINDIR)/Code.$(StatSuf)
SoLib		=$(BINDIR)/Code.$(DllSuf)

.PRECIOUS:*.ObjSuf *.DepSuf *.DllSuf

Deps=$(patsubst %,$(BINDIR)/%.$(DepSuf),$(Objects) $(Packages) )

############### EXPLICIT RULES ###############

$(BINDIR):
	mkdir -p $(BINDIR)

info:
	@echo "--------------------------"
	@echo "Compile on $(shell hostname)"
	@echo "Packages are: $(Packages)"
	@echo "Objects are: $(Objects)"
	@echo "--------------------------"
	@echo "DEBUG:"

$(StatLib): $(BINOBJ)
	ar rcs $@ $(BINOBJ)
.PHONY: soLib
soLib: $(SoLib)

$(SoLib): $(StatLib)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^

.PHONY: $(Packages) 
$(Packages): % : $(BINDIR)/% | $(BINDIR)
	@echo $(call InfoLine , $@ )

#$(BINDIR)/$(Packages): $(BINDIR)/% : $(BASEDIR)/test/%.$(SrcSuf) $(StatLib) | $(BINDIR)
$(addprefix $(BINDIR)/,$(Packages)): $(BINDIR)/% : $(BASEDIR)/test/%.$(SrcSuf) $(StatLib) | $(BINDIR)
	@echo $(call InfoLine , $@ )
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(StatLib) -I$(INC_DIR) -I$(HDIR) -I$(INC_DIR_CMSSW) -L$(LIBDIR) -L$(ROOTSYS)/lib

.PHONY: clean
clean:
	-rm -v bin/*.$(ObjSuf)
	-rm -v bin/*.$(DllSuf)
	-rm -v bin/*.$(StatSuf)


############### IMPLICIT RULES ###############

#.o
%.$(ObjSuf): $(BINDIR)/%.$(ObjSuf)

#.o
$(BINDIR)/%.$(ObjSuf): $(SRCDIR)/%.$(SrcSuf) $(HDIR)/%.$(HeadSuf)
	@echo $(call InfoLine , $@ )
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $(SRCDIR)/$*.$(SrcSuf) -I$(INC_DIR) -I$(HDIR) -I$(INC_DIR_CMSSW) -L$(LIBDIR) -L$(ROOTSYS)/lib

#.d
$(BINDIR)/%.$(DepSuf): $(SRCDIR)/%.$(SrcSuf) $(HDIR)/%.$(HeadSuf)
	@echo $(call InfoLine , $@ )
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -M -o $@ $(SRCDIR)/$*.$(SrcSuf) -I$(INC_DIR) -I$(HDIR) -I$(INC_DIR_CMSSW) -L$(LIBDIR) -L$(ROOTSYS)/lib
	sed -i'' "s|^.*:|& Makefile $(BINDIR)/&|g" $@
