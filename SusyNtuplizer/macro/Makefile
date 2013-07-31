# Original Author:  Dongwook Jang
# $Id: Makefile,v 1.4 2013/04/01 09:53:24 yiiyama Exp $

#--------------------------------------------------
# Change file names for your need
#--------------------------------------------------
TARGET  = libSusyEvent.so

INC     = -I. -I$(shell root-config --incdir)
LIBS    = $(shell root-config --libs)

OBJDIR  = ../src
skip_files = $(OBJDIR)/SusyNtuplizer.cc $(OBJDIR)/SusyTriggerNtuplizer.cc
OBJECTS    = $(filter-out $(skip_files), $(wildcard $(OBJDIR)/*.cc))
DICTS   = Event_dict

#---------------------------------------------------
# compiler dependent flags
#---------------------------------------------------
ARCH         := $(shell root-config --arch)
#
ifeq ($(ARCH),linuxx8664gcc)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared
endif

ifeq ($(ARCH),macosx64)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = clang++
CXXFLAGS      = -O2 -Wall -fPIC
LD            = clang++
LDFLAGS       = -O2
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxkcc)
# Linux with the KAI compiler
CXX           = KCC --one_instantiation_per_object
CXXFLAGS      = -fPIC +K0 --backend -w --no_exceptions
LD            = KCC
LDFLAGS       = -O $(shell root-config --cflags)
SOFLAGS       =
endif

ifeq ($(ARCH),linuxicc)
# Linux with Intel icc compiler
CXX           = icc
CXXFLAGS      = -O
LD            = icc
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),win32)
# Windows with the VC++ compiler
CXX           = cl
CXXOPT        = -O2
#CXXOPT        = -Z7
CXXFLAGS      = $(CXXOPT) -G5 -GR -GX -MD -DWIN32 -D_WINDOWS -nologo \
                -DVISUAL_CPLUSPLUS -D_X86_=1 -D_DLL
LD            = link
LDOPT         = -opt:ref
#LDOPT         = -debug
LDFLAGS       = $(LDOPT) -pdb:none -nologo
SOFLAGS       = -DLL

LIBS  = $(shell root-config --nonew --libs)
endif

ifeq ($(ARCH),win32gcc)
# Windows with gcc
CXX           = g++
CXXFLAGS      = -O2 -march=pentium4 -mfpmath=sse -Wall -Woverloaded-virtual -I/usr/X11R6/include
LD            = g++
LDFLAGS       = -O2 -Wl,--enable-auto-import
SOFLAGS       = -shared -D_DLL -Wl,--export-all-symbols
endif


#-------------------------------------------
# Rules
#-------------------------------------------

all: $(TARGET) 
	@@echo "Done!"

so: $(TARGET)
	@@echo "Done!"

$(DICTS): 
	@@echo "Generating dictionary ..."
	@@rootcint -f $(DICTS).cc -c $(OBJDIR)/SusyEvent.h $(OBJDIR)/SusyNtuplizer_LinkDef.h


$(TARGET): $(DICTS)
	@@echo "Making $(TARGET) ..."
	@@$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $(TARGET) $(INC) $(LIBS) $(OBJECTS) $(DICTS).cc
	@@/bin/rm -f $(DICTS).*

clean:
	@@/bin/rm -f $(TARGET) $(DICTS).* *_cc.d *_cc.so *~
