find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIRS})


#ROOT_EXECUTABLE(testKalman testKalman.cxx LIBRARIES Core MathCore Matrix Vc)
#ROOT_ADD_TEST(example_kalman_autovec COMMAND testKalman )

#ROOT_EXECUTABLE(testOperations testOperations.cxx LIBRARIES Core MathCore Matrix Vc)
#ROOT_ADD_TEST(example_smatrix_autovec COMMAND testOperations )

AddCompilerFlag("-DUSE_VC")
set(USEAVX)

build_example(smatrix testOperations.cxx LIBS Core MathCore Matrix)
build_example(kalman testKalman.cxx LIBS Core MathCore Matrix)

