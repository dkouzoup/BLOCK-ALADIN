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
include CMakeFiles/example_wang2010.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/example_wang2010.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/example_wang2010.dir/flags.make

CMakeFiles/example_wang2010.dir/examples/wang2010.c.o: CMakeFiles/example_wang2010.dir/flags.make
CMakeFiles/example_wang2010.dir/examples/wang2010.c.o: ../examples/wang2010.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/example_wang2010.dir/examples/wang2010.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/example_wang2010.dir/examples/wang2010.c.o   -c /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/wang2010.c

CMakeFiles/example_wang2010.dir/examples/wang2010.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example_wang2010.dir/examples/wang2010.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/wang2010.c > CMakeFiles/example_wang2010.dir/examples/wang2010.c.i

CMakeFiles/example_wang2010.dir/examples/wang2010.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example_wang2010.dir/examples/wang2010.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/examples/wang2010.c -o CMakeFiles/example_wang2010.dir/examples/wang2010.c.s

CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.requires:
.PHONY : CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.requires

CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.provides: CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.requires
	$(MAKE) -f CMakeFiles/example_wang2010.dir/build.make CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.provides.build
.PHONY : CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.provides

CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.provides.build: CMakeFiles/example_wang2010.dir/examples/wang2010.c.o

# Object files for target example_wang2010
example_wang2010_OBJECTS = \
"CMakeFiles/example_wang2010.dir/examples/wang2010.c.o"

# External object files for target example_wang2010
example_wang2010_EXTERNAL_OBJECTS =

bin/wang2010: CMakeFiles/example_wang2010.dir/examples/wang2010.c.o
bin/wang2010: CMakeFiles/example_wang2010.dir/build.make
bin/wang2010: lib/libqpdunes.a
bin/wang2010: CMakeFiles/example_wang2010.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/wang2010"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_wang2010.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/example_wang2010.dir/build: bin/wang2010
.PHONY : CMakeFiles/example_wang2010.dir/build

CMakeFiles/example_wang2010.dir/requires: CMakeFiles/example_wang2010.dir/examples/wang2010.c.o.requires
.PHONY : CMakeFiles/example_wang2010.dir/requires

CMakeFiles/example_wang2010.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/example_wang2010.dir/cmake_clean.cmake
.PHONY : CMakeFiles/example_wang2010.dir/clean

CMakeFiles/example_wang2010.dir/depend:
	cd /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build /Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/CMakeFiles/example_wang2010.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/example_wang2010.dir/depend

