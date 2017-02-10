CC = icc
CFLAGS  = -std=c++11 -Wall -fast
LDFLAGS = #-opt-report=4 -opt-report-phase=loop,vec -opt-report-file=stderr
HEAD = $(wildcard *.h)
SOURCES = $(wildcard *.cpp) 
EXECUTABLE = FD

OBJECTS = $(SOURCES:.cpp=.o)
$(info Object Files = $(OBJECTS))
$(info header Files = $(HEAD))

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< 
clean:
	rm *.o $(EXECUTABLE)
