## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "DCA")
set(CTEST_NIGHTLY_START_TIME "00:00:00 GMT")

set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "cdash.cscs.ch")
set(CTEST_DROP_LOCATION "/submit.php?project=DCA")
set(CTEST_DROP_SITE_CDASH TRUE)

