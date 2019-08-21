################################################################################


TARGET := filter

#-----------------------------------------------------------------------------
# COMPILER AND OPT/DEBUG FLAGS
CC = gcc-9
CXX = g++-9


CXXFLAGS := -g -Wall -O0 --debug -std=c++11 -DDEBUG
# CXXFLAGS := -O3 -g -std=c++11 

# define some options
DEF :=

# linker specific flags
# LDFLAGS := -L/Users/DeeGee/Documents/lide.space/FDR_data/Kalman/eigen/build/include/Eigen
LDFLAGS := -L/usr/local/Cellar/gcc/9.1.0/lib/gcc/9 -lstdc++

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./
OBJ_DIR := ./build

## add the headers to the vpaths
# INC := -I/Users/DeeGee/Documents/lide.space/FDR_data/Kalman/eigen/build/include/Eigen
INC := -I/Users/DeeGee/Documents/lide.space/FDR_data/Kalman/eigen/

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))

## generate object list
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -c $< -o $@
	# $(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################

default: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIB)

test:
	@echo $(SRC)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET)

destroy:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(TARGET)
	rm -f $(OBJ_DIR)/*

info:
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info DEP = $(DEP))

-include $(DEP)
