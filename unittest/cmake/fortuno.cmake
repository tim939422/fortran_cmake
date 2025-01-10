# Duosi Fan Jan 10, 2025
# requirement must include(FetchContent) in the main CMakeLists.txt

set(FORTUNO_VERSION "0.1.0")
set(FORTUNO_TARBALL "fortuno-${FORTUNO_VERSION}.tar.gz")

if(DEFINED THIRD_PARTY_TARBALL_DIR)
  message(STATUS "Trying to use local tarball for Fortuno: ${FORTUNO_TARBALL}")
  set(FORTUNO_TARBALL_PATH "${THIRD_PARTY_TARBALL_DIR}/${FORTUNO_TARBALL}")

  if(EXISTS ${FORTUNO_TARBALL_PATH})
    message(STATUS "Local tarball found at ${FORTUNO_TARBALL_PATH}")
    FetchContent_Declare(
      Fortuno
      URL ${FORTUNO_TARBALL_PATH}
    )
  else()
    message(WARNING "Local tarball not found at ${FORTUNO_TARBALL_PATH}. Fetching from Git repository...")
    FetchContent_Declare(
      Fortuno
      GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno"
      GIT_TAG "v${FORTUNO_VERSION}"
    )
  endif()
else()
  message(STATUS "THIRD_PARTY_TARBALL_DIR not set. Fetching Fortuno from Git repository...")
  FetchContent_Declare(
    Fortuno
    GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno"
    GIT_TAG "v${FORTUNO_VERSION}"
  )
endif()
FetchContent_MakeAvailable(Fortuno)