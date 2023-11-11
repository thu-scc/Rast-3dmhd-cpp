CXX = mpic++
CXXFLAGS := -O2 -std=c++20
# CXXFLAGS += -fsanitize=address -fsanitize=undefined
TARGET = 3dmhd
SRCS = 3dmhd.cpp
OBJS = $(SRCS:.cpp=.o)
HDRS = $(wildcard *.h)

NUM_PROC = 24

$(TARGET): $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	mpirun -np $(NUM_PROC) ./$(TARGET)