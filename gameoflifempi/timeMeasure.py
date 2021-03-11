import os
import re
import csv
import subprocess as sp
from datetime import datetime

thread_sizes = [(1,1), (1,2), (2,1), (2,2), (4,1), (1,4)]
board_sizes = [1024, 2048, 4096]
timestep = 654
timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

if not os.path.exists("csv"):
    os.mkdir("csv")

with open("csv/results"+ timestamp +".csv", "w") as out:  
    writer = csv.writer(out)
    writer.writerow(["px", "py", "nx", "ny", "minutes", "seconds", "millis", "total_millis"])

    for board_size in board_sizes:
        for thread_size_x, thread_size_y in thread_sizes:
            
            px = thread_size_x
            py = thread_size_y

            p = px * py

            nx = board_size // thread_size_x
            ny = board_size // thread_size_y
            
            for i in range(5):
                result = sp.run(["time", "mpirun", "-n", str(p), "./gameoflifempi", str(timestep), str(nx), str(ny), str(px), str(py)], stderr=sp.PIPE)
                elapsed = str(result.stderr.split()[2])
                time_search = re.search("([0-9]+):([0-9]+).([0-9]+)", elapsed)
                minutes = time_search.group(1)
                seconds = time_search.group(2)
                millis = time_search.group(3)
                total_millis = int(minutes) * 60000 + int(seconds) * 1000 + int(millis)
                print("ITERATION", i, "Elapsed time:", minutes,"minutes", seconds, "seconds", millis, "miliseconds", total_millis, "Total")
                writer.writerow([px, py, nx, ny, minutes, seconds, millis, total_millis])