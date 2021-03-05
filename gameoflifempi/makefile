
CC = mpicc

CFLAGS  = -g -Wall -std=c99 -lc -D _BSD_SOURCE -lm

TARGET = gameoflife

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c


run: all
	clear
	echo $$(($$(date +%s%N)/1000000))
	mpirun -np 4 ./$(TARGET)
	echo $$(($$(date +%s%N)/1000000))
	
clean:
	$(RM) $(TARGET)
	$(RM) data/*.vti
	$(RM) data/*.pvti
