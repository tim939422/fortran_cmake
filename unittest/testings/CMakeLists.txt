add_executable(testapp)

target_sources(
    testapp PRIVATE
    test_simple.f90
    testapp.f90
)

target_link_libraries(
    testapp PRIVATE mylib Fortuno::fortuno_serial
)

add_test(NAME test_simple COMMAND testapp)