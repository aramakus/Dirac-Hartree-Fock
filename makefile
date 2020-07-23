# Create a folder named "build" and adjust compiler according to the platform.
CXX = g++
CXXFLAGS = -O3 -std=c++11
OBJ_DIR = build
SRC_FILES = $(wildcard *.cpp)
OBJ_FILES = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

$(OBJ_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

Atom: $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	-rm $(OBJ_DIR)/*