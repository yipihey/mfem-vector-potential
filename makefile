# Build the vector-potential study: libvp (src/*.cpp) + experiments/*.cpp -> bin/*
#
#   make                       # library + all experiments
#   make MFEM_INSTALL_DIR=/path/to/install/mfem
#   make bin/exp1_static
#
MFEM_INSTALL_DIR ?= $(CURDIR)/external/install/mfem
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
ifeq (,$(wildcard $(CONFIG_MK)))
$(error $(CONFIG_MK) not found: set MFEM_INSTALL_DIR to an installed MFEM)
endif
include $(CONFIG_MK)

CXX      = $(MFEM_CXX)
# MFEM_FLAGS carries -std=c++17 -O3 -DNDEBUG and the include paths;
# the later -O2 wins.
CXXFLAGS = $(MFEM_FLAGS) -O2 -Wall -Wno-unused-local-typedefs -Isrc

SRCS  = $(wildcard src/*.cpp)
OBJS  = $(patsubst src/%.cpp,build/%.o,$(SRCS))
EXPS  = $(wildcard experiments/*.cpp)
BINS  = $(patsubst experiments/%.cpp,bin/%,$(EXPS))
LIB   = build/libvp.a

.PHONY: all clean
all: $(BINS)

build/%.o: src/%.cpp src/vp_core.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJS)
	ar rcs $@ $^

bin/%: experiments/%.cpp $(LIB) src/vp_core.hpp $(wildcard experiments/*.hpp)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIB) $(MFEM_LIBS)

clean:
	rm -rf build bin
