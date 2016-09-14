#include <gtest/gtest.h>
#include "example.cpp"
#include <cmath>
#include <vector>

int ex_argc;
char** ex_argv;
bool result;
double FACTOR = 10;
std::vector<double> bench_u_err, bench_p_err;
std::vector<double> u_err, p_err;
double bench;
double actual;

// tests for u error

TEST(navier_stokes_ex0, 2d_u_L1Norm) {
	bench  = bench_u_err[0];
	actual = u_err[0];
    EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 2d_u_L2Norm) {
    bench  = bench_u_err[1];
    actual = u_err[1];
    EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 2d_u_MaxNorm) {
   bench  = bench_u_err[2];
   actual = u_err[2];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 3d_u_L1Norm) {
    bench  = bench_u_err[3];
    actual = u_err[0];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 3d_u_L2Norm) {
    bench  = bench_u_err[4];
    actual = u_err[1];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 3d_u_MaxNorm) {
 	bench  = bench_u_err[5];
    actual = u_err[2];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

// tests for p error

TEST(navier_stokes_ex0, 2d_p_L1Norm) {
	bench  = bench_p_err[0];
	actual = p_err[0];
    EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 2d_p_L2Norm) {
    bench  = bench_p_err[1];
    actual = p_err[1];
    EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 2d_p_MaxNorm) {
   bench  = bench_p_err[2];
   actual = p_err[2];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 3d_p_L1Norm) {
    bench  = bench_p_err[3];
    actual = p_err[0];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 3d_p_L2Norm) {
    bench  = bench_p_err[4];
    actual = p_err[1];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}

TEST(navier_stokes_ex0, 3d_p_MaxNorm) {
 	bench  = bench_p_err[5];
    actual = p_err[2];
   EXPECT_LE(std::abs((actual - bench)), (bench/FACTOR));
}


int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv );
    
    //error recorded from main2d running input2d.test Oct 4, 2016
    // benchmark error in u 
    bench_u_err.resize(6);
    bench_u_err[0] = 0.00357601;  //2d L1Norm
    bench_u_err[1] = 0.00439633;  //2d L2Norm
    bench_u_err[2] = 0.0417876;   //2d maxNorm
    bench_u_err[3] = 0.000510578; //3d L1Norm
    bench_u_err[4] = 0.000477978; //3d L2Norm
    bench_u_err[5] = 0.0055429;   //3d maxNorm
    // benchmark error in p
    bench_p_err.resize(6);
    bench_p_err[0] = 0.0219484;  //2d L1Norm
    bench_p_err[1] = 0.0295763;  //2d L2Norm
    bench_p_err[2] = 0.220349;   //2d maxNorm
    bench_p_err[3] = 0.0021559;  //3d L1Norm
    bench_p_err[4] = 0.00249808; //3d L2Norm
    bench_p_err[5] = 0.00843954; //3d maxNorm

    ex_argc = argc;
    ex_argv = argv;
    run_example(ex_argc, ex_argv, u_err, p_err);
    return RUN_ALL_TESTS( );
}
