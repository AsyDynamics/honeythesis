# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/asy/Documents/sim/simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/asy/Documents/sim/simulation/build

# Include any dependencies generated for this target.
include CMakeFiles/run.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/run.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/run.dir/flags.make

CMakeFiles/run.dir/src/main.cpp.o: CMakeFiles/run.dir/flags.make
CMakeFiles/run.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/asy/Documents/sim/simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/run.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/src/main.cpp.o -c /home/asy/Documents/sim/simulation/src/main.cpp

CMakeFiles/run.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/asy/Documents/sim/simulation/src/main.cpp > CMakeFiles/run.dir/src/main.cpp.i

CMakeFiles/run.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/asy/Documents/sim/simulation/src/main.cpp -o CMakeFiles/run.dir/src/main.cpp.s

CMakeFiles/run.dir/src/main.cpp.o.requires:

.PHONY : CMakeFiles/run.dir/src/main.cpp.o.requires

CMakeFiles/run.dir/src/main.cpp.o.provides: CMakeFiles/run.dir/src/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/run.dir/build.make CMakeFiles/run.dir/src/main.cpp.o.provides.build
.PHONY : CMakeFiles/run.dir/src/main.cpp.o.provides

CMakeFiles/run.dir/src/main.cpp.o.provides.build: CMakeFiles/run.dir/src/main.cpp.o


# Object files for target run
run_OBJECTS = \
"CMakeFiles/run.dir/src/main.cpp.o"

# External object files for target run
run_EXTERNAL_OBJECTS =

../bin/run: CMakeFiles/run.dir/src/main.cpp.o
../bin/run: CMakeFiles/run.dir/build.make
../bin/run: CMakeFiles/run.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/asy/Documents/sim/simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/run"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/run.dir/build: ../bin/run

.PHONY : CMakeFiles/run.dir/build

CMakeFiles/run.dir/requires: CMakeFiles/run.dir/src/main.cpp.o.requires

.PHONY : CMakeFiles/run.dir/requires

CMakeFiles/run.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run.dir/clean

CMakeFiles/run.dir/depend:
	cd /home/asy/Documents/sim/simulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/asy/Documents/sim/simulation /home/asy/Documents/sim/simulation /home/asy/Documents/sim/simulation/build /home/asy/Documents/sim/simulation/build /home/asy/Documents/sim/simulation/build/CMakeFiles/run.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/run.dir/depend

