#include <gtest/gtest.h>
#include "example.cpp"

int example_argc;
char** example_argv;
bool ExampleDoesRun;
bool runExample(int, char**);

TEST(DummyTest, OneIsOne) {
    ASSERT_EQ(1.0, 1.0);
}

TEST(ExampleDoesRun, Bool) {
    ExampleDoesRun = runExample(example_argc, example_argv);
    ASSERT_EQ(ExampleDoesRun, true);
}

int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv ); 
    example_argc = argc;
    example_argv = argv;
    return RUN_ALL_TESTS( );
}
