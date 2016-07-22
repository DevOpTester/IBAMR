#include <gtest/gtest.h>
#include "example.cpp"

bool runExample(int, char*);

TEST(DummyTest, OneIsOne) {
    ASSERT_EQ(1.0, 1.0);
}

int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS( );
}
