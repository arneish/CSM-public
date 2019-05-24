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

# %.o: %.cpp
# 	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
# $(MAIN): $(OBJS)
	# $(CC) $(CFLAGS) $^ -o $@

# all: $(MAIN)
# 	@echo Compiled csm

# $(MAIN): $(OBJS)
# 	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS)

# .c.o:
# 	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# clean:
# 	$(RM) *.o *~ $(MAIN)

# depend: $(SRCS)
# 	makedepend $(INCLUDES) $^

# all NOSB: csm_NOSB 
# 	$(info Compiled CSM (NO size-bound))
	
# # wsb: csm_wsb 
# # 	$(info Compiled CSM (WITH size-bound))
	
# csm_NOSB: algorithm_NOSB.o DFScode.o replica.o
# 	$(CC) $(CFLAGS) -o csm algorithm_NOSB.o DFScode.o replica.o

# # csm_wsb: algorithm_wsb.o DFScode.o replica.o
# # 	@$(CC) $(CFLAGS) -o csm algorithm_wsb.o DFScode.o replica.o

# algorithm_NOSB.o: algorithm_NOSB.cpp csm.h
# 	$(CC) $(CFLAGS) -c algorithm_NOSB.cpp

# # algorithm_wsb.o: algorithm_wsb.cpp csm.h
# # 	@$(CC) $(CFLAGS) -c algorithm_wsb.cpp

# DFScode.o: DFScode.cpp DFScode.h 
# 	$(CC) $(CFLAGS) -c DFScode.cpp

# replica.o: replica.cpp csm.h
# 	$(CC) $(CFLAGS) -c replica.cpp
# g++ -std=c++11 -o csm -O3 algorithm_NOSB.cpp DFScode.cpp replica.cpp
