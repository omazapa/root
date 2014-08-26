find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIRS})

set(Libraries Core RIO Net Hist Graf Graf3d Gpad Tree
              Rint Postscript Matrix Physics MathCore Thread)

AddCompilerFlag("-DUSE_VC")
set(USEAVX)

build_example(genvector vectorOperation.cxx LIBS Core MathCore Matrix GenVector)

