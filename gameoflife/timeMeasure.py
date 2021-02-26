import os
import re
import csv
import subprocess as sp

thread_sizes = [(1,1), (1,2), (2,1), (2,2), (4,1), (1,4)]
board_sizes = [1024, 2048, 4096]
timestep = 400

for thread_size_x, thread_size_y in thread_sizes:
    for board_size in board_sizes:
        px = thread_size_x
        py = thread_size_y
        nx = board_size // thread_size_x
        ny = board_size // thread_size_y
        
        result = sp.run(["time", "./gameoflife", str(timestep), str(nx), str(ny), str(px), str(py)], stderr=sp.PIPE)
        elapsed = str(result.stderr.split()[2])
        time_search = re.search("([0-9]+):([0-9]+).([0-9]+)", elapsed)
        print(elapsed)
        print(time_search)
        minutes = time_search.group(1)
        seconds = time_search.group(2)
        milis = time_search.group(3)

        print("ITERATION", 0, "Elapsed time:", minutes,"minutes", seconds, "seconds", milis, "miliseconds")
