! This file is part of Fortuno.
! Licensed under the BSD-2-Clause Plus Patent license.
! SPDX-License-Identifier: BSD-2-Clause-Patent

!> Test app with command line interface, collecting and executing the tests.
program testapp
    use fortuno_serial, only : execute_serial_cmd_app, test_list
    use test_simple, only : simple_tests => tests
    implicit none
  
    ! Creating and executing a command line app with the tests to be included.
    ! Note: this function does not return but stops the code with the right exit code.
    ! (0 on success, non-zero otherwise)
    call execute_serial_cmd_app(test_list([ simple_tests()]))
  
  end program testapp