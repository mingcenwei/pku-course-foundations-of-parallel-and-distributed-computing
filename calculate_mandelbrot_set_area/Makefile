CXX         := -c++
CXXFLAGS    := \
	-std=gnu++17 \
	-pedantic-errors \
	-Wall \
	-Wextra \
	-Wconversion \
	-Wsign-conversion \
	-fopenmp \
	-fconstexpr-ops-limit=10000000000000 \

LDFLAGS     := \

BUILD       := ./build
OBJ_DIR     := $(BUILD)/objects
APP_DIR     := $(BUILD)/bin
INSTALL_DIR := .
TARGET      := calculate_mandelbrot_set_area
INCLUDE     := -I./include/
SRC         := \
	$(wildcard *.cpp) \
	$(wildcard src/*.cpp) \

OBJECTS     := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: release_with_debug_info install clean
.PHONY: all

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(LDFLAGS)

# .PHONY: all build clean debug release

build: $(APP_DIR)/$(TARGET)
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
.PHONY: build

install:
	@mkdir -p $(INSTALL_DIR)
	@cp -v $(APP_DIR)/$(TARGET) $(INSTALL_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: build
.PHONY: debug

release: CXXFLAGS += -O3
release: build
.PHONY: release

release_with_debug_info: CXXFLAGS += -O3 -g
release_with_debug_info: build
.PHONY: release_with_debug_info

# -@rm -rv $(OBJ_DIR)
# -@rm -rv $(APP_DIR)
clean:
	-@rm -rv $(BUILD)
.PHONY: clean

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... build"
	@echo "... clean"
	@echo "... debug"
	@echo "... help"
	@echo "... install"
	@echo "... release"
	@echo "... release_with_debug_info"
.PHONY: help