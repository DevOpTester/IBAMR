
##IBAMR automated tests -- a brief guide. ##

## Dependencies ##

To run IBAMR's automated tests using `make gtest` or `make gtest-long` you need to have available the [Google Test Framework](https://github.com/google/googletest) installed using the same compiler and compiler flags as the IBAMR project. 

Discussion of how this is accomplished and why it is necessary is explored in the [Google Test README](googletest/googletest/README.md).

><font size=3>***Note:*** *The Google Test Framwork is often referred to as simply* **`gtest`** </font size>

Once the gtest library is installed, IBAMR must be configured with the LDFLAGS variable pointing to it. 

Example configure invocation:
```bash
./configure \
  CFLAGS="-Wall" \
  CXXFLAGS="-Wall" \
  FCFLAGS="-Wall" \
  CC="ccache $HOME/linux/openmpi/1.10.2/bin/mpicc" \
  CXX="ccache $HOME/linux/openmpi/1.10.2/bin/mpicxx" \
  FC=$HOME/linux/openmpi/1.10.2/bin/mpif90 \
  CPPFLAGS="-DOMPI_SKIP_MPICXX" \
  LDFLAGS="-lgtest" \
  --with-hypre=$PETSC_DIR/$PETSC_ARCH \
  --with-samrai=$SAMRAI_DIR \
  --with-hdf5=$HOME/linux/hdf5/1.8.16 \
  --with-blitz=$HOME/linux/blitz/0.10 \
  --with-silo=$HOME/linux/silo/4.10 \
  --with-boost=$HOME/linux/boost/1.61.0 \
  --enable-libmesh \
  --with-libmesh=$HOME/linux/libmesh/1.0.0/1.0.0-debug \
  --with-libmesh-method=dbg
```

## Invocation ##

Once the gtest library is installed and IBAMR is configured properly, all the tests can be run by invoking `make gtest` or `make gtest-long` in the root build directory. 

#### `make gtest` vs. `make gtest-long`####
* `make gtest` is intended as a smoke test/sanity check that can be run in less than half an hour on an average machine used for development. `make gtest` compiles a selection of 2D and 3D tests, and of these selected only runs the 2D tests.

*  `make gtest-long` compiles and runs every single test (2D and 3D) and takes a significant amount of time (exact runtime depends on your system resources).

## Running Individual Tests ##

Each test currently included in the library is based off of example applications and uses the same source code as the examples, but runs the simulation from within a gtest application and analyzes the results. To run an individul test, navigate to the example of interest and run `make gtest` from the command line.

## Interpreting Results ##

Every gtest application included in the IBAMR library returns a result indicating whether or not the example runs. This 

Below is some annotated output from an individual test application:

```
...( output from simulation )...

[==========] Running 7 tests from 1 test case.
[----------] Global test environment set-up.
[----------] 7 tests from IBFE_ex0_2d
[ RUN      ] IBFE_ex0_2d.example_runs
[       OK ] IBFE_ex0_2d.example_runs (0 ms)
[ RUN      ] IBFE_ex0_2d.u_L1Norm
[       OK ] IBFE_ex0_2d.u_L1Norm (0 ms)
[ RUN      ] IBFE_ex0_2d.u_L2Norm
[       OK ] IBFE_ex0_2d.u_L2Norm (0 ms)
[ RUN      ] IBFE_ex0_2d.u_MaxNorm
[       OK ] IBFE_ex0_2d.u_MaxNorm (0 ms)
[ RUN      ] IBFE_ex0_2d.p_L1Norm
[       OK ] IBFE_ex0_2d.p_L1Norm (0 ms)
[ RUN      ] IBFE_ex0_2d.p_L2Norm
[       OK ] IBFE_ex0_2d.p_L2Norm (0 ms)
[ RUN      ] IBFE_ex0_2d.p_MaxNorm
[       OK ] IBFE_ex0_2d.p_MaxNorm (0 ms)
[----------] 7 tests from IBFE_ex0_2d (0 ms total)

[----------] Global test environment tear-down
[==========] 7 tests from 1 test case ran. (1 ms total)
[  PASSED  ] 7 tests.
```

Notice that each test name follows the convention:
`NAME_OF_TEST_APPLICATION.NAME_OF_QUALITY_TESTED`

**Example interpretation 1:**
```
[ RUN      ] IBFE_ex0_2d.example_runs
[       OK ] IBFE_ex0_2d.example_runs (0 ms)
```
Indicates that the application run was the 2D case of the application found in examples/IBFE/ex0, and the the test was whether or not the application was able to run to completion and did not experience any run time errors. 

This is a base line test to make sure any changes to the API of different IBAMR classes is accurately reflected in the example applications.

**Example interpretation 2:**
```
[ RUN      ] IBFE_ex0_2d.u_MaxNorm
[       OK ] IBFE_ex0_2d.u_MaxNorm (0 ms)
[ RUN      ] IBFE_ex0_2d.p_L1Norm
[       OK ] IBFE_ex0_2d.p_L1Norm (0 ms)
```

Again, each test indicates it was run on the 2D case of the applicaiton found in examples/IBFE/ex0.  Here the first test `IBFE_ex0_2d.u_MaxNorm` is observing whether or not the test application's "Max Norm" (also known as the infinity norm) of the velocity is within an acceptable margin of error from an analytical or laboratory result. The second test is performing a similar calculation except this time with the L1 norm of the pressure.

#### _Notation_ ####
**L1_Norm**: The one-norm (also known as the L1-norm, `1 norm, or mean norm)  is defined as the sum of the absolute values of its components.

**L2_Norm**:  The two-norm (also known as the L2-norm, 2-norm, mean-square norm, or least-squares norm) is defined as the square root of the sum of the squares of the absolute values of its components.

**MaxNorm**:  The infinity norm (also known as the L∞-norm, `∞-norm, max norm, or uniform norm) is defined as the maximum of the absolute values of its
components.

### Writing new tests ###

#### Adding tests to existing applications ####

Adding new tests to examples with existing gtest applications, for example to `$IBAMR_DIR/examples/IBFE/explicit/ex1/`, can be done by adding tests to the `test_main.cpp` file in that directory following the pattern:

```
TEST(TEST_CASE_NAME, your_unique_test_name_here) {
	// do stuff
	// make some assertions
}
```
Many types  of assertions are available and are well documented in the [Google Test Primer](https://github.com/google/googletest/blob/master/googletest/docs/Primer.md) and the [Google Test Advanced Guide](https://github.com/google/googletest/blob/master/googletest/docs/AdvancedGuide.md).

If you desire to collect additional data from the example, it is suggested you pass a reference to the object that you'd like to test to the example. An example of this can be found in `$IBAMR_DIR/examples/IBFE/explicit/ex0/test_main.cpp `.

#### Creating new gtest applications ####

Two considerations must be kept in mind when writing new gtest applications. 

1. Users should not be required to have the Google Test Framework installed in order to use the rest of the library.

2. In order for your new tests to be run by Jenkins, you must update all necessary parent Makefile.am files.

## Output ##
By default, the gtest applications will only output to stdout.

The **`$GTEST_OUTPUT`** environment variable is available to make the gtest applications generate xml reports of test results. Its behavior is detailed in [the advanced guide to gtest](https://github.com/google/googletest/blob/48ee8e98abc950abd8541e15550b18f8f6cfb3a9/googletest/docs/V1_7_AdvancedGuide.md#generating-an-xml-report). 

In short,

> To generate the XML report, set the `GTEST_OUTPUT` environment variable or the `--gtest_output` flag to the string `"xml:_path_to_output_file_"`, which will create the file at the given location. You can also just use the string `"xml"`, in which case the output can be found in the `test_detail.xml` file in the current directory.

This is especially useful when using the [JUnit plugin for Jenkins](https://wiki.jenkins-ci.org/display/JENKINS/JUnit+Plugin) which knows how to parse these xml files. 

The following is a build script for a job running tests and reporting the results using the [JUnit plugin](https://wiki.jenkins-ci.org/display/JENKINS/JUnit+Plugin):
```bash
export HOME=/srv/sfw
export BOOST_ROOT=$HOME/linux/boost/1.61.0
export PETSC_ARCH=linux-debug
export PETSC_DIR=$HOME/petsc/3.7.2
export GTEST_OUTPUT="xml:$WORKSPACE/"
export SAMRAI_DIR="$HOME/samrai/2.4.4/linux-debug"
rm -rf *.xml
./.jenkins_quicktest

```

With this set, all test results will be amalgamated into one report the the plugin can use to determine the build status and generate graphical reports in the job view. Old test results must be removed at the beginning of each build in order to avoid cluttering up the workspace.

> _**Note:**_`$WORKSPACE` *is one of the [environment variable is available from Jenkins](https://wiki.jenkins-ci.org/display/JENKINS/Building+a+software+project#Buildingasoftwareproject-JenkinsSetEnvironmentVariables) and resolves to the absolute path to the workspace*

