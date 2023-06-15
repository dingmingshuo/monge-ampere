CC = g++
CFLAGS = -W -Wall -Wextra -Wconversion -Werror -std=c++17 -Ofast -march=native -g
INC = -I/usr/local/include -I./include
LIB = -lma -lpthread -fopenmp
LIBMA = -L./$(BUILD_DIR) -lma
TEST_LIB = $(LIBMA) -L/usr/local/lib -lgtest

BUILD_DIR = build
SOURCE_DIR = source
TEST_DIR = test
EXAMPLE_DIR = example

SOURCES = $(wildcard $(SOURCE_DIR)/*.cpp)
OBJECTS = $(addprefix $(BUILD_DIR)/, $(SOURCES:.cpp=.o))

TEST_SOURCES = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS = $(addprefix $(BUILD_DIR)/, $(TEST_SOURCES:.cpp=.o))

EXAMPLE_SOURCES = $(wildcard $(EXAMPLE_DIR)/*.cpp)
EXAMPLE_OBJECTS = $(addprefix $(BUILD_DIR)/, $(EXAMPLE_SOURCES:.cpp=.o))
EXAMPLE_EXE = $(addprefix $(BUILD_DIR)/, $(EXAMPLE_SOURCES:.cpp=))

all: prepare libma alltest $(EXAMPLE_EXE)

prepare:
	@echo "Compiling libma ..."
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/$(SOURCE_DIR)
	mkdir -p $(BUILD_DIR)/$(TEST_DIR)
	mkdir -p $(BUILD_DIR)/$(EXAMPLE_DIR)
	@echo $(SOURCES)
	@echo $(OBJECTS)
	@echo $(TEST_SOURCES)
	@echo $(TEST_OBJECTS)
	@echo $(EXAMPLE_SOURCES)
	@echo $(EXAMPLE_OBJECTS)

libma: prepare $(OBJECTS)
	ar cr $(BUILD_DIR)/libma.a $(OBJECTS)

$(BUILD_DIR)/$(SOURCE_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(CC) $(CFLAGS) $(LIB) $(INC) -c $< -o $@

alltest: prepare libma $(TEST_OBJECTS)
	$(CC) $(CFLAGS) $(INC) -o $(BUILD_DIR)/test_main $(TEST_OBJECTS) $(LIB) $(TEST_LIB)
	./$(BUILD_DIR)/test_main
	
$(BUILD_DIR)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp libma
	$(CC) $(CFLAGS) $(INC) -c $< -o $@ $(LIB) $(TEST_LIB)

$(BUILD_DIR)/$(EXAMPLE_DIR)/%: $(EXAMPLE_DIR)/%.cpp libma
	$(CC) $(CFLAGS) $(INC) $< -o $@ $(LIB) $(LIBMA)

clean:
	rm -rf build

.PHONY: all clean prepare alltest libma