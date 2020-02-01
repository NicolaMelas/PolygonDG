#############################################
# Directories set by the user
#############################################

# Enter the directory where the Eigen library headers are:
EIGEN_INC  = /vagrant/eigen

# Enter the directory where the Eigen unsupported library headers are:
EIGEN_INC_UN  = /vagrant/eigen/unsupported

# Enter muparserx include directory
MUPX_INC = /vagrant/pacs-examples/Extras/muparserx/muparserx/parser/
MUPX_LIB = /vagrant/pacs-examples/Extras/muparserx/muparserx

# Enter cgal include directory
CG_INC = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/cgal/4.9/include/
CG_LIB = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/cgal/4.9/lib/

# Enter boost include directory
BOOST_INC = /vagrant/boost

# Enter Qhull include directory
QHULL_INC = /vagrant/qhull-2019.1/src
QHULL_LIB = /vagrant/qhull-lib

# Enter the directory where the GetPot header is (needed only by the executables):
# change if you want to use your own version
GETPOT_INC = ./externalinclude

#enter other directories if needed
OTHER_INC=
OTHER_LIB=

# Enter the directory where you want to install the library
POLYGONDG_PATH = /vagrant/TestLib

##########
# Compiler
##########

CXX = g++

##########################################################################################
#################################	NOT MODIFIABLE PART ###################################
##########################################################################################


#muParserXInterface header :
EXT_INC = ./externalinclude


###################
# Compilation flags
###################
WFLAGS  += -Wall
STDFLAG += -std=c++14

ifdef RELEASE
 OPTFLAGS += -O3 -DNDEBUG
else
 OPTFLAGS += -g
 WFLAGS   += -Wextra -pedantic
endif

CXXFLAGS += $(WFLAGS) $(STDFLAG) $(OPTFLAGS)

ALL_INC_DIR= $(EIGEN_INC) $(MUPX_INC) $(EIGEN_INC_UN) $(QHULL_INC) $(CG_INC) $(BOOST_INC) $(OTHER_INC) $(EXT_INC)
INCLS := $(foreach inc,$(ALL_INC_DIR),-I$(inc) )

ALL_LIB_DIR=$(MUPX_LIB) $(CG_LIB) $(QHULL_LIB) $(QHULL_LIB) $(OTHER_LIB)
LIBD := $(foreach lib,$(ALL_LIB_DIR),-L$(lib) )

CPPFLAGS += $(INCLS)

LD_LIBS  +=  -lPolygonDG -lqhullcpp -lqhull_r  -lmuparserx  -lCGAL -lgmp  -lstdc++ -lm
LD_FLAGS += $(LIBD)


#############
# Directories
#############
POLYGONDG_INC = $(POLYGONDG_PATH)/PolygonDgLib/include
POLYGONDG_LIB = $(POLYGONDG_PATH)/PolygonDgLib/lib
