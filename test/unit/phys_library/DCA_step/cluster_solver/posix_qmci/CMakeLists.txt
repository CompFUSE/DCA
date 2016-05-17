set(TEST_INCLUDES ${PROJECT_SOURCE_DIR}/src)
set(TEST_LIBS posix_qmci)

add_gtest(NAME thread_task_handler_test
          GTEST_MAIN
          INCLUDES "${TEST_INCLUDES}"
          LIBS "${TEST_LIBS}")
