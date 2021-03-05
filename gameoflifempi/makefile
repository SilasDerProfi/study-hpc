CC = mpicc
TARGET = gameoflifempi

$(TARGET): $(TARGET).c
	$(CC) -o $(TARGET) $(TARGET).c

run: $(TARGET)
	mpirun -np 4 ./$(TARGET)