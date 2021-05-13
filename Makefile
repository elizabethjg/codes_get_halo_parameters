CC = g++
CFLAGS = -O3 -lgsl -lgslcblas -lm
TARGET = get_halo_props_PIC

all:$(TARGET).cpp
	$(CC) $(TARGET).cpp $(CFLAGS) -o $(TARGET)

clean:
	$(RM) $(Target)
