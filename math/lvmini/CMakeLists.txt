############################################################################
# CMakeLists.txt file for building ROOT math/lvmini package
############################################################################

ROOT_USE_PACKAGE(math/mathcore)
ROOT_USE_PACKAGE(hist/hist)

add_definitions(-DWARNINGMSG -DUSE_ROOT_ERROR)

ROOT_GENERATE_DICTIONARY(G__LVMini *.h LVMini/*.h LINKDEF LinkDef.h)
ROOT_GENERATE_ROOTMAP(LVMini LINKDEF LinkDef.h)
ROOT_LINKER_LIBRARY(LVMini lvmini.F *.cxx G__LVMini.cxx DEPENDENCIES MathCore Hist)
#ROOT_LINKER_LIBRARY(LVMini lvmini.cpp *.cxx G__LVMini.cxx DEPENDENCIES MathCore Hist f2c)
ROOT_INSTALL_HEADERS()
