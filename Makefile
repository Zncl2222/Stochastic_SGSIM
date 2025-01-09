# Detect the operating system
ifeq ($(OS),Windows_NT)
    BUILD_SCRIPT = .\build_scripts\cmake_build.bat
else
    BUILD_SCRIPT = ./build_scripts/cmake_build.sh
    CMD_PREFIX = sh
endif

# Default target
all:
	@echo "Building based on the OS-specific script..."
	$(CMD_PREFIX) $(BUILD_SCRIPT)

# Target for shared library
shared:
	@echo "Building shared library..."
	$(CMD_PREFIX) $(BUILD_SCRIPT) -s

# Clean build directories
clean:
	@echo "Cleaning build directories..."
ifeq ($(OS),Windows_NT)
	@if exist cbuild (rmdir /s /q cbuild)
else
	rm -rf cmake_build
endif
