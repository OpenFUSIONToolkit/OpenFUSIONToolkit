#---------------------------------------
# File make_inc.mk for systems using build_libs.py
#---------------------------------------

# Path to src directory
OFT_BASE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Include build settings
include $(OFT_BASE_DIR)/make_libs.mk

# Setup default options (Optimized+OpenMP)
ifndef BUILD_TYPE
BUILD_TYPE = OPT
endif
ifndef USE_OMP
USE_OMP = 1
endif

# Set compiler options
CC_FLAGS = $(BASE_CFLAGS)
FC_FLAGS = $(BASE_FFLAGS)
ifeq ($(BUILD_TYPE), OPT)
CC_FLAGS += $(OPT_FLAGS)
FC_FLAGS += $(OPT_FLAGS)
else ifeq ($(BUILD_TYPE), DEBUG)
CC_FLAGS += $(DEBUG_FLAGS)
FC_FLAGS += $(DEBUG_FLAGS) $(CHK_FLAGS)
else ifeq ($(BUILD_TYPE), NO_OPT)
CC_FLAGS += $(DEBUG_FLAGS)
FC_FLAGS += $(DEBUG_FLAGS)
else
$(error Invalid "BUILD_TYPE": must be "OPT", "DEBUG" or "NO_OPT")
endif
ifeq ($(USE_OMP), 1)
CC_FLAGS += $(OMP_FLAGS)
FC_FLAGS += $(OMP_FLAGS)
endif
ifeq ($(USE_STACK), 1)
COMP_DEFS += -DOFT_STACK -DOFT_PROFILE
endif
