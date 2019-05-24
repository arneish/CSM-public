CXX = g++
CPPFLAGS += -Iinclude -MMD -MP
CXXFLAGS += -std=c++14 -flto -O3

SRC_DIR = src
OBJ_DIR = build
BIN_DIR = bin
SRC = $(SRC_DIR)/replica.cpp $(SRC_DIR)/algorithm_nosb.cpp $(SRC_DIR)/DFScode.cpp
SRC_WSB = $(SRC_DIR)/replica.cpp $(SRC_DIR)/algorithm_wsb.cpp $(SRC_DIR)/DFScode.cpp
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OBJ_WSB = $(SRC_WSB:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
DEPS = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.d)
DEPS_WSB = $(SRC_WSB:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.d)
EXE = $(BIN_DIR)/csm
EXE_WSB = $(BIN_DIR)/csm_wsb

GCC_CXXFLAGS = -DMESSAGE='"GCC: g++"'
CLANG_CXXFLAGS = -DMESSAGE='"Clang: clang++"'
UNKNOWN_CXXFLAGS = -DMESSAGE='"Unknown compiler"'

ifeq ($(CXX),g++)
  CXXFLAGS += $(GCC_CXXFLAGS)
else ifeq ($(CXX),clang++)
  CXXFLAGS += $(CLANG_CXXFLAGS)
else
  CXXFLAGS += $(UNKNOWN_CXXFLAGS)
endif

.PHONY: all clean 

all nosb: $(EXE)

wsb: $(EXE_WSB)

$(EXE): $(OBJ) 

$(EXE_WSB): $(OBJ_WSB)

$(EXE) $(EXE_WSB): 
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	$(RM) $(BIN_DIR)/* $(OBJ_DIR)/* 

print-%  : ; @echo $* = $($*)

-include $(DEPS)
-include $(DEPS_WSB)