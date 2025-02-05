# Project root directory (relative to this Makefile)
ROOT_DIR = ..

# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -I$(ROOT_DIR)

# Default source file (can be overridden)
# If no source file is specified, use main.cpp from parent directory
MAIN_SRC ?= $(ROOT_DIR)/main.cpp

# Source files from root directory
ROOT_SOURCES = $(ROOT_DIR)/AutoDiff_lib.cpp \
              $(ROOT_DIR)/fft.cpp \
              $(ROOT_DIR)/lib_basic.cpp \
              $(ROOT_DIR)/Matrix_lib.cpp \
              $(ROOT_DIR)/Numerical_lib.cpp \
              $(ROOT_DIR)/pch.cpp \
              $(ROOT_DIR)/Vector_lib.cpp

# All source files
ALL_SOURCES = $(ROOT_SOURCES)

# Object files will be created in the current directory
OBJECTS = $(addprefix obj/, $(notdir $(ALL_SOURCES:.cpp=.o)))
MAIN_OBJ = obj/$(notdir $(MAIN_SRC:.cpp=.o))

# Create obj directory
$(shell mkdir -p obj)

# Output executable name
TARGET = matrix

# Default target
all: $(TARGET)

# Main target
$(TARGET): $(OBJECTS) $(MAIN_OBJ)
	$(CXX) $(OBJECTS) $(MAIN_OBJ) -o $(TARGET)

# Compile source files from root directory
obj/%.o: $(ROOT_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile local source files
obj/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Generate dependencies
depend: .depend
.depend: $(ALL_SOURCES) $(MAIN_SRC)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $(ALL_SOURCES) $(MAIN_SRC) | sed 's|^|obj/|' > ./.depend

# Clean target
clean:
	rm -rf obj $(TARGET) .depend

# Help target
help:
	@echo "Usage:"
	@echo "  make                   - Build with default main source ($(MAIN_SRC))"
	@echo "  make MAIN_SRC=path    - Build with specified main source"
	@echo "  make clean            - Remove all built files"
	@echo "  make help             - Show this help message"
	@echo ""
	@echo "Notes:"
	@echo "  - Place this Makefile and your source files in the newsourceq4/ directory"
	@echo "  - The executable will be created in the same directory"
	@echo "  - Object files will be placed in the obj/ subdirectory"

.PHONY: all clean help depend

# Include dependencies if they exist
-include .depend