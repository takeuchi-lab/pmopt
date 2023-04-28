.PHONY: install, all, clean, test

TARGET		:= pmopt
BUILD		:= build
TEST		:= 
TEST_ARGS	:=
CMAKE_NJOBS := 8
MODE		:= Debug

# Define this to suppress some debug info.
# Do NOT change here !!
GNUMAKEFLAGS 	:= --no-print-directory

install:
	cmake -DCMAKE_BUILD_TYPE=$(MODE) -S . -B $(BUILD)
	cmake --build $(BUILD) --target $(TARGET) -j $(CMAKE_NJOBS)
	@ mv $(BUILD)/src/$(TARGET) .


all:  clean test $(TARGET)
	

clean:
	rm -rf build $(TARGET)


test:
	cmake -S . -B $(BUILD)
	cmake --build $(BUILD) $(if $(TEST),--target $(TEST),) -j $(CMAKE_NJOBS)
	ctest --test-dir $(BUILD) --output-on-failure $(if $(TEST),-R $(TEST),) $(TEST_ARGS)