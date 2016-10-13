#include <gtest/gtest.h>
#include "example.cpp"
#include <cmath>

int ex_argc;
char** ex_argv;
bool ex_runs;
bool run_example(int, char**);

TEST(navier_stokes_ex4, 2d) { 
    EXPECT_EQ(ex_runs, true);
}

TEST(navier_stokes_ex4, 3d) {
    EXPECT_EQ(ex_runs, true);
}

int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv ); 
    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv);
    return RUN_ALL_TESTS( );
}