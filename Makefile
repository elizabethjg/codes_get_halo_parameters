CC = g++
CFLAGS = -O3 -lgsl -lgslcblas -lm

SOURCES := halo_energy.cpp recentering.cpp get_halo_props_PIC_v2.cpp
OBJECTS := $(SOURCES:.cpp=.o)
TARGET = get_halo_props_PIC

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) $(TARGET) $(OBJECTS)
