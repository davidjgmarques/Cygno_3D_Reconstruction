# Compiler
CXX = c++

# Compiler flags
CXXFLAGS = `root-config --cflags`

# Linker flags
LDFLAGS = `root-config --glibs` -lSpectrum

# Target executable
TARGET = alpha_3d.out

# Source files
SRCS = Alpha_3D_LIME.cpp Track_analyzer/Analyzer.cxx plotting_functions.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link the object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile the source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up the build
clean:
	rm -f $(TARGET) $(OBJS)

# Phony targets
.PHONY: all clean
