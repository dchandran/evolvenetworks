.. _adding_regression_tests:

Adding Regression Tests
=======================

This page describes how to add regression tests for a Boost library in
the CMake-based build system. Before adding regression tests, make
sure you have already followed the directions for
[wiki:CMakeAddingALibrary adding a library to CMake], so that the
CMake system recognizes your Boost library, and (if necessary)
[wiki:CMakeAddingCompiledLibrary adding a compiled library to
CMake]. We also assume that you have already configured your build

.. index:: BUILD_TESTING

.. _BUILD_TESTING:

tree for regression testing, by enabling the ``BUILD_TESTING`` option
described in the section :ref:`testing`.

In this page, we will assume that your library resides in the
subdirectory ``libs/libname``, and that tests for this library are
stored in ``libs/libname/test``. The test directory should be listed
via ``TESTDIRS`` in the use of the [wiki:CMakeLibraryProject
boost_library_project macro], as described in an earlier section,
[wiki:CMakeAddingALibrary "Adding a Library to CMake"]. Follow these
steps to add this new library into Boost's build system. If your
library has multiple testing directories listed after ``TESTDIRS``,
follow these steps for each one.

1. Create a new file ``libs/libname/test/CMakeLists.txt`` file with
your favorite text editor. This file will contain instructions for
building and running each of the regression tests for your library.

1a. If your regression test depends on any other part of boost then
you will need to inform the build system of such with the following
line::

  boost_additional_test_dependencies(libname BOOST_DEPENDS test fusion)

where 'libname' is the name of your library that you are testing.

2. For each test that only needs to be compiled (but not executed),
add a ``compile`` or ``compile_fail`` test using the
[wiki:CMakeTestCompile boost_test_compile] or
[wiki:CMakeTestCompileFail boost_test_compile_fail] macros,
respectively. The most basic usage of these macros provides only the
test name, e.g., ::

  boost_test_compile(compile_test)
  boost_test_compile_fail(compile_fail_test)

This code will create two regression tests. The first,
``compile_test``, will try to compile the source file
``compile_test.cpp`` in the current source directory. If the compile
is successful, the regression test passes. If the compile fails, the
regression test fails. The second regression test works the opposite
way: it will try to compile ``compile_fail_test.cpp``: if the
compilation is successful, the regression test fails. When you run the
regression tests (e.g., by calling ``ctest`` from the build
directory), the regression tests will execute and produce output like
the following::

  
  Running tests...
  Start processing tests
  Test project /Users/dgregor/Projects/boost-darwin
    1/  2 Testing libname::compile_test            Passed
    2/  2 Testing libname::compile_fail_test     ***Failed - supposed to fail
  
  100% tests passed, 0 tests failed out of 2

3. For any tests that need to be built and executed, use the
[wiki:CMakeTestRun boost_test_run] or [wiki:CMakeTestRunFail
boost_test_run_fail] macros. Both tests will build, link and execute a
regression test. The ``boost_test_run`` macro expects that executable
to return an exit code of zero, while the ``boost_test_run_fail``
macro expects that executable to return a non-zero exit code. For
example, we might build a simple test ``simple_test`` from the source
file ``simple_test.cpp``: {{{ boost_test_run(simple_test) }}}

4. Often, we'll want to link against our own Boost library, which we
do using the ``DEPENDS`` argument to ``boost_test_run``::

  boost_test_run(big_test big_test1.cpp big_test2.cpp
    DEPENDS boost_libname-static
    )

Here, we have created a test ``big_test``, built from the source files
``big_test1.cpp`` and ``big_test2.cpp``, which will link against the
static library for ``boost_libname``. We could create a similar test
that links against the shared library for ``boost_libname``, passing
along compilation flags specific to the shared library::

  boost_test_run(big_test_dll big_test1.cpp big_test2.cpp
    DEPENDS boost_libname-shared
    COMPILE_FLAGS "-DBOOST_LIBNAME_DYN_LINK=1"
    )

Some tests require command-line arguments. For example, say we want to
pass ``-loop 1000`` to a randomized test. We can do so using the
``ARGS`` argument to ``boost_test_run`` (or ``boost_test_run_fail``)::

  boost_test_run(random_test ARGS "-loop" "1000" DEPENDS boost_libname-static)

Once you have finished describing your regression tests to the CMake
system, you're done! Your library will now build, test, and install
with CMake and this behavior should be portable across many different
platforms.

