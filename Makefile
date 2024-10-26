# https://makefiletutorial.com/#makefile-cookbook
# https://spin.atomicobject.com/makefile-c-projects/

CC = gcc
CXX = g++
CXXFLAGS ?= -O2 -Wall -Wno-sign-compare
# override commandline: specify standard
override CXXFLAGS := -std=c++20 ${CXXFLAGS}

TARGET = main
BUILD_DIR = build
SRCS := main.cpp solvers.cpp local_search.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.cpp.o)
DEPS := $(OBJS:.cpp.o=.cpp.d)

.phony: all clean main run

all: $(TARGET)

run: $(TARGET)
# ./$(TARGET)
	./$(TARGET) TSPA.csv
	./$(TARGET) TSPB.csv

$(TARGET): $(OBJS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# in build/ folder
$(BUILD_DIR)/%.cpp.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPS)

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
#	rm -f $(TARGET) $(OBJS) $(DEPS)

# https://stackoverflow.com/questions/16467718/how-to-print-out-a-variable-in-makefile
# print-%  : ; @echo $* = $($*)
