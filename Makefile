# https://makefiletutorial.com/#makefile-cookbook
# https://spin.atomicobject.com/makefile-c-projects/

CC = gcc
CXX = g++
CXXFLAGS ?= -O2 -Wall -Wno-sign-compare -march=native
# override commandline: specify standard
override CXXFLAGS := -std=c++20 ${CXXFLAGS}

TARGET = main
BUILD_DIR = build
SRCS := main.cpp solvers.cpp local_search.cpp msls_ils.cpp convexity.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.cpp.o)
DEPS := $(OBJS:.cpp.o=.cpp.d)

.phony: all clean run rebuild

all: $(TARGET)

# allows `make -j4`
rebuild: clean
	@$(MAKE) --no-print-directory main
# make all

run: $(TARGET)
	./$(TARGET)
# ./$(TARGET) TSPA.csv
# ./$(TARGET) TSPB.csv

$(BUILD_DIR):
	@mkdir -p $@

$(TARGET): $(OBJS) | $(BUILD_DIR)
	$(CXX) -pipe $(CXXFLAGS) -o $@ $^

# in build/ folder
$(BUILD_DIR)/%.cpp.o: %.cpp | $(BUILD_DIR)
	$(CXX) -pipe $(CXXFLAGS) -MMD -c $< -o $@

-include $(DEPS)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
#	rm -f $(TARGET) $(OBJS) $(DEPS)

# https://stackoverflow.com/questions/16467718/how-to-print-out-a-variable-in-makefile
# print-%  : ; @echo $* = $($*)
