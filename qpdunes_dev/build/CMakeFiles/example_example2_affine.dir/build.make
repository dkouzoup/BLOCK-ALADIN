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
include CMakeFiles/example_example2_affine.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/example_example2_affine.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/example_example2_affine.dir/flags.make

CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o: CMakeFiles/example_example2_affine.dir/flags.make
CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o: ../examples/example2_affine.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o   -c /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/example2_affine.c

CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/example2_affine.c > CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.i

CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/example2_affine.c -o CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.s

CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.requires:
.PHONY : CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.requires

CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.provides: CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.requires
	$(MAKE) -f CMakeFiles/example_example2_affine.dir/build.make CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.provides.build
.PHONY : CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.provides

CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.provides.build: CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o

# Object files for target example_example2_affine
example_example2_affine_OBJECTS = \
"CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o"

# External object files for target example_example2_affine
example_example2_affine_EXTERNAL_OBJECTS =

bin/example2_affine: CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o
bin/example2_affine: CMakeFiles/example_example2_affine.dir/build.make
bin/example2_affine: lib/libqpdunes.a
bin/example2_affine: CMakeFiles/example_example2_affine.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/example2_affine"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_example2_affine.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/example_example2_affine.dir/build: bin/example2_affine
.PHONY : CMakeFiles/example_example2_affine.dir/build

CMakeFiles/example_example2_affine.dir/requires: CMakeFiles/example_example2_affine.dir/examples/example2_affine.c.o.requires
.PHONY : CMakeFiles/example_example2_affine.dir/requires

CMakeFiles/example_example2_affine.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/example_example2_affine.dir/cmake_clean.cmake
.PHONY : CMakeFiles/example_example2_affine.dir/clean

CMakeFiles/example_example2_affine.dir/depend:
	cd /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles/example_example2_affine.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/example_example2_affine.dir/depend
