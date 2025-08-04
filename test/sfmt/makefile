CC = gcc
CFLAGS = -Wall -O2
TARGET = main
OBJS = main.o sfmt/SFMT.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

main.o: main.c sfmt/SFMT.h
	$(CC) $(CFLAGS) -DSFMT_MEXP=19937 -c main.c

sfmt/SFMT.o: sfmt/SFMT.c sfmt/SFMT.h sfmt/SFMT-params.h sfmt/SFMT-params19937.h
	$(CC) $(CFLAGS) -DSFMT_MEXP=19937 -c sfmt/SFMT.c -o sfmt/SFMT.o

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean

