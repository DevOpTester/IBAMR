#include <gtest/gtest.h>
#include "example.cpp"

int ex_argc;
char** ex_argv;
bool ex_runs;
bool run_example(int, char**);

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME IBFE_ex0_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME IBFE_ex0_3d
#endif

TEST(IBFE_ex0, 2d) {
    ex_runs = run_example(ex_argc, ex_argv);
    EXPECT_EQ(ex_runs, true);
}

int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv ); 
    ex_argc = argc;
    ex_argv = argv;
    return RUN_ALL_TESTS( );
}
