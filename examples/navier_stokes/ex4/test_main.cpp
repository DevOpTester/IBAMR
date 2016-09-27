#include <gtest/gtest.h>
#include "example.cpp"
#include <cmath>

int ex_argc;
char** ex_argv;
bool ex_runs;
double run_example(int, char**);
double u_error;
double std_2d_u_error = 2.83046e-05;
double std_3d_u_error = 0.00553941;
double EPSILON = 0.01;

TEST(navier_stokes_ex4, 2d_uMax_norm_error) { 
    EXPECT_LE(std::abs(u_error - std_2d_u_error), EPSILON);
}

TEST(navier_stokes_ex4, 3d_uMax_norm_error) {
    EXPECT_LE(std::abs(u_error - std_3d_u_error), EPSILON);
}

int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv ); 
    ex_argc = argc;
    ex_argv = argv;
    u_error = run_example(ex_argc, ex_argv);
    return RUN_ALL_TESTS( );
}
