# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.0.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.0.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build

# Include any dependencies generated for this target.
include CMakeFiles/example_doubleIntegrator_qp.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/example_doubleIntegrator_qp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/example_doubleIntegrator_qp.dir/flags.make

CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o: CMakeFiles/example_doubleIntegrator_qp.dir/flags.make
CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o: ../examples/doubleIntegrator_qp.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o   -c /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/doubleIntegrator_qp.c

CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/doubleIntegrator_qp.c > CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.i

CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/doubleIntegrator_qp.c -o CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.s

CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.requires:
.PHONY : CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.requires

CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.provides: CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.requires
	$(MAKE) -f CMakeFiles/example_doubleIntegrator_qp.dir/build.make CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.provides.build
.PHONY : CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.provides

CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.provides.build: CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o

# Object files for target example_doubleIntegrator_qp
example_doubleIntegrator_qp_OBJECTS = \
"CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o"

# External object files for target example_doubleIntegrator_qp
example_doubleIntegrator_qp_EXTERNAL_OBJECTS =

bin/doubleIntegrator_qp: CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o
bin/doubleIntegrator_qp: CMakeFiles/example_doubleIntegrator_qp.dir/build.make
bin/doubleIntegrator_qp: lib/libqpdunes.a
bin/doubleIntegrator_qp: CMakeFiles/example_doubleIntegrator_qp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/doubleIntegrator_qp"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_doubleIntegrator_qp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/example_doubleIntegrator_qp.dir/build: bin/doubleIntegrator_qp
.PHONY : CMakeFiles/example_doubleIntegrator_qp.dir/build

CMakeFiles/example_doubleIntegrator_qp.dir/requires: CMakeFiles/example_doubleIntegrator_qp.dir/examples/doubleIntegrator_qp.c.o.requires
.PHONY : CMakeFiles/example_doubleIntegrator_qp.dir/requires

CMakeFiles/example_doubleIntegrator_qp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/example_doubleIntegrator_qp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/example_doubleIntegrator_qp.dir/clean

CMakeFiles/example_doubleIntegrator_qp.dir/depend:
	cd /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles/example_doubleIntegrator_qp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/example_doubleIntegrator_qp.dir/depend

