# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.4

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
CMAKE_COMMAND = /opt/cmake-3.4.3/bin/cmake

# The command to remove a file.
RM = /opt/cmake-3.4.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ph123693/Documents/2D/solver2d/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ph123693/Documents/2D/solver2d/src

# Include any dependencies generated for this target.
include libfem/CMakeFiles/fem.dir/depend.make

# Include the progress variables for this target.
include libfem/CMakeFiles/fem.dir/progress.make

# Include the compile flags for this target's objects.
include libfem/CMakeFiles/fem.dir/flags.make

libfem/CMakeFiles/fem.dir/solver.cpp.o: libfem/CMakeFiles/fem.dir/flags.make
libfem/CMakeFiles/fem.dir/solver.cpp.o: libfem/solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ph123693/Documents/2D/solver2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libfem/CMakeFiles/fem.dir/solver.cpp.o"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fem.dir/solver.cpp.o -c /home/ph123693/Documents/2D/solver2d/src/libfem/solver.cpp

libfem/CMakeFiles/fem.dir/solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fem.dir/solver.cpp.i"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ph123693/Documents/2D/solver2d/src/libfem/solver.cpp > CMakeFiles/fem.dir/solver.cpp.i

libfem/CMakeFiles/fem.dir/solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fem.dir/solver.cpp.s"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ph123693/Documents/2D/solver2d/src/libfem/solver.cpp -o CMakeFiles/fem.dir/solver.cpp.s

libfem/CMakeFiles/fem.dir/solver.cpp.o.requires:

.PHONY : libfem/CMakeFiles/fem.dir/solver.cpp.o.requires

libfem/CMakeFiles/fem.dir/solver.cpp.o.provides: libfem/CMakeFiles/fem.dir/solver.cpp.o.requires
	$(MAKE) -f libfem/CMakeFiles/fem.dir/build.make libfem/CMakeFiles/fem.dir/solver.cpp.o.provides.build
.PHONY : libfem/CMakeFiles/fem.dir/solver.cpp.o.provides

libfem/CMakeFiles/fem.dir/solver.cpp.o.provides.build: libfem/CMakeFiles/fem.dir/solver.cpp.o


libfem/CMakeFiles/fem.dir/output_fem.cpp.o: libfem/CMakeFiles/fem.dir/flags.make
libfem/CMakeFiles/fem.dir/output_fem.cpp.o: libfem/output_fem.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ph123693/Documents/2D/solver2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object libfem/CMakeFiles/fem.dir/output_fem.cpp.o"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fem.dir/output_fem.cpp.o -c /home/ph123693/Documents/2D/solver2d/src/libfem/output_fem.cpp

libfem/CMakeFiles/fem.dir/output_fem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fem.dir/output_fem.cpp.i"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ph123693/Documents/2D/solver2d/src/libfem/output_fem.cpp > CMakeFiles/fem.dir/output_fem.cpp.i

libfem/CMakeFiles/fem.dir/output_fem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fem.dir/output_fem.cpp.s"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ph123693/Documents/2D/solver2d/src/libfem/output_fem.cpp -o CMakeFiles/fem.dir/output_fem.cpp.s

libfem/CMakeFiles/fem.dir/output_fem.cpp.o.requires:

.PHONY : libfem/CMakeFiles/fem.dir/output_fem.cpp.o.requires

libfem/CMakeFiles/fem.dir/output_fem.cpp.o.provides: libfem/CMakeFiles/fem.dir/output_fem.cpp.o.requires
	$(MAKE) -f libfem/CMakeFiles/fem.dir/build.make libfem/CMakeFiles/fem.dir/output_fem.cpp.o.provides.build
.PHONY : libfem/CMakeFiles/fem.dir/output_fem.cpp.o.provides

libfem/CMakeFiles/fem.dir/output_fem.cpp.o.provides.build: libfem/CMakeFiles/fem.dir/output_fem.cpp.o


libfem/CMakeFiles/fem.dir/four_face_current.cpp.o: libfem/CMakeFiles/fem.dir/flags.make
libfem/CMakeFiles/fem.dir/four_face_current.cpp.o: libfem/four_face_current.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ph123693/Documents/2D/solver2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object libfem/CMakeFiles/fem.dir/four_face_current.cpp.o"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fem.dir/four_face_current.cpp.o -c /home/ph123693/Documents/2D/solver2d/src/libfem/four_face_current.cpp

libfem/CMakeFiles/fem.dir/four_face_current.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fem.dir/four_face_current.cpp.i"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ph123693/Documents/2D/solver2d/src/libfem/four_face_current.cpp > CMakeFiles/fem.dir/four_face_current.cpp.i

libfem/CMakeFiles/fem.dir/four_face_current.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fem.dir/four_face_current.cpp.s"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ph123693/Documents/2D/solver2d/src/libfem/four_face_current.cpp -o CMakeFiles/fem.dir/four_face_current.cpp.s

libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.requires:

.PHONY : libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.requires

libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.provides: libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.requires
	$(MAKE) -f libfem/CMakeFiles/fem.dir/build.make libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.provides.build
.PHONY : libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.provides

libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.provides.build: libfem/CMakeFiles/fem.dir/four_face_current.cpp.o


libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o: libfem/CMakeFiles/fem.dir/flags.make
libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o: libfem/seven_face_current.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ph123693/Documents/2D/solver2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fem.dir/seven_face_current.cpp.o -c /home/ph123693/Documents/2D/solver2d/src/libfem/seven_face_current.cpp

libfem/CMakeFiles/fem.dir/seven_face_current.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fem.dir/seven_face_current.cpp.i"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ph123693/Documents/2D/solver2d/src/libfem/seven_face_current.cpp > CMakeFiles/fem.dir/seven_face_current.cpp.i

libfem/CMakeFiles/fem.dir/seven_face_current.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fem.dir/seven_face_current.cpp.s"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ph123693/Documents/2D/solver2d/src/libfem/seven_face_current.cpp -o CMakeFiles/fem.dir/seven_face_current.cpp.s

libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.requires:

.PHONY : libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.requires

libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.provides: libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.requires
	$(MAKE) -f libfem/CMakeFiles/fem.dir/build.make libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.provides.build
.PHONY : libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.provides

libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.provides.build: libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o


libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o: libfem/CMakeFiles/fem.dir/flags.make
libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o: libfem/ten_face_current.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ph123693/Documents/2D/solver2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fem.dir/ten_face_current.cpp.o -c /home/ph123693/Documents/2D/solver2d/src/libfem/ten_face_current.cpp

libfem/CMakeFiles/fem.dir/ten_face_current.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fem.dir/ten_face_current.cpp.i"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ph123693/Documents/2D/solver2d/src/libfem/ten_face_current.cpp > CMakeFiles/fem.dir/ten_face_current.cpp.i

libfem/CMakeFiles/fem.dir/ten_face_current.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fem.dir/ten_face_current.cpp.s"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && /opt/gcc/5.2.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ph123693/Documents/2D/solver2d/src/libfem/ten_face_current.cpp -o CMakeFiles/fem.dir/ten_face_current.cpp.s

libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.requires:

.PHONY : libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.requires

libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.provides: libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.requires
	$(MAKE) -f libfem/CMakeFiles/fem.dir/build.make libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.provides.build
.PHONY : libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.provides

libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.provides.build: libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o


# Object files for target fem
fem_OBJECTS = \
"CMakeFiles/fem.dir/solver.cpp.o" \
"CMakeFiles/fem.dir/output_fem.cpp.o" \
"CMakeFiles/fem.dir/four_face_current.cpp.o" \
"CMakeFiles/fem.dir/seven_face_current.cpp.o" \
"CMakeFiles/fem.dir/ten_face_current.cpp.o"

# External object files for target fem
fem_EXTERNAL_OBJECTS =

libfem/libfem.a: libfem/CMakeFiles/fem.dir/solver.cpp.o
libfem/libfem.a: libfem/CMakeFiles/fem.dir/output_fem.cpp.o
libfem/libfem.a: libfem/CMakeFiles/fem.dir/four_face_current.cpp.o
libfem/libfem.a: libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o
libfem/libfem.a: libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o
libfem/libfem.a: libfem/CMakeFiles/fem.dir/build.make
libfem/libfem.a: libfem/CMakeFiles/fem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ph123693/Documents/2D/solver2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library libfem.a"
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && $(CMAKE_COMMAND) -P CMakeFiles/fem.dir/cmake_clean_target.cmake
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fem.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libfem/CMakeFiles/fem.dir/build: libfem/libfem.a

.PHONY : libfem/CMakeFiles/fem.dir/build

libfem/CMakeFiles/fem.dir/requires: libfem/CMakeFiles/fem.dir/solver.cpp.o.requires
libfem/CMakeFiles/fem.dir/requires: libfem/CMakeFiles/fem.dir/output_fem.cpp.o.requires
libfem/CMakeFiles/fem.dir/requires: libfem/CMakeFiles/fem.dir/four_face_current.cpp.o.requires
libfem/CMakeFiles/fem.dir/requires: libfem/CMakeFiles/fem.dir/seven_face_current.cpp.o.requires
libfem/CMakeFiles/fem.dir/requires: libfem/CMakeFiles/fem.dir/ten_face_current.cpp.o.requires

.PHONY : libfem/CMakeFiles/fem.dir/requires

libfem/CMakeFiles/fem.dir/clean:
	cd /home/ph123693/Documents/2D/solver2d/src/libfem && $(CMAKE_COMMAND) -P CMakeFiles/fem.dir/cmake_clean.cmake
.PHONY : libfem/CMakeFiles/fem.dir/clean

libfem/CMakeFiles/fem.dir/depend:
	cd /home/ph123693/Documents/2D/solver2d/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ph123693/Documents/2D/solver2d/src /home/ph123693/Documents/2D/solver2d/src/libfem /home/ph123693/Documents/2D/solver2d/src /home/ph123693/Documents/2D/solver2d/src/libfem /home/ph123693/Documents/2D/solver2d/src/libfem/CMakeFiles/fem.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libfem/CMakeFiles/fem.dir/depend

