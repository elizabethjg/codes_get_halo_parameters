CC = g++
CFLAGS = -O3 -ffast-math -faggressive-loop-optimizations -ftree-vectorize
LDFLAGS = -lgsl -lgslcblas -lm

SOURCES := halo_energy.cpp recentering.cpp get_halo_props_PIC_v2.cpp
OBJECTS := $(SOURCES:.cpp=.o)
TARGET = get_halo_props_PIC_v2

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	$(RM) $(TARGET) $(OBJECTS)
