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
include CMakeFiles/example_singleIntegrator_b.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/example_singleIntegrator_b.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/example_singleIntegrator_b.dir/flags.make

CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o: CMakeFiles/example_singleIntegrator_b.dir/flags.make
CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o: ../examples/singleIntegrator_b.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o   -c /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/singleIntegrator_b.c

CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/singleIntegrator_b.c > CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.i

CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/singleIntegrator_b.c -o CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.s

CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.requires:
.PHONY : CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.requires

CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.provides: CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.requires
	$(MAKE) -f CMakeFiles/example_singleIntegrator_b.dir/build.make CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.provides.build
.PHONY : CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.provides

CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.provides.build: CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o

# Object files for target example_singleIntegrator_b
example_singleIntegrator_b_OBJECTS = \
"CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o"

# External object files for target example_singleIntegrator_b
example_singleIntegrator_b_EXTERNAL_OBJECTS =

bin/singleIntegrator_b: CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o
bin/singleIntegrator_b: CMakeFiles/example_singleIntegrator_b.dir/build.make
bin/singleIntegrator_b: lib/libqpdunes.a
bin/singleIntegrator_b: CMakeFiles/example_singleIntegrator_b.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/singleIntegrator_b"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_singleIntegrator_b.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/example_singleIntegrator_b.dir/build: bin/singleIntegrator_b
.PHONY : CMakeFiles/example_singleIntegrator_b.dir/build

CMakeFiles/example_singleIntegrator_b.dir/requires: CMakeFiles/example_singleIntegrator_b.dir/examples/singleIntegrator_b.c.o.requires
.PHONY : CMakeFiles/example_singleIntegrator_b.dir/requires

CMakeFiles/example_singleIntegrator_b.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/example_singleIntegrator_b.dir/cmake_clean.cmake
.PHONY : CMakeFiles/example_singleIntegrator_b.dir/clean

CMakeFiles/example_singleIntegrator_b.dir/depend:
	cd /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles/example_singleIntegrator_b.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/example_singleIntegrator_b.dir/depend
