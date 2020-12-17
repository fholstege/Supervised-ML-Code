# Saved as calculate_squares.py
from mpi4py import MPI

def square_me(nstart, nstop):
    """
    This function computes all squares from nstart till nstop and returns them in a list.
    """
    squares = []
    for n in range(nstart, nstop):
        squares.append((n, n * n))
    return squares
    
# Define communicator and determine size and rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Divide the work over the all processes
number_of_squares = 1000000
nstart = (number_of_squares // size) * rank
if rank == (size - 1):
    nstop = number_of_squares
else:
    nstop = (number_of_squares // size) * (rank + 1)

# Compute the squares    
print('On rank {} doing {} till {}'.format(rank, nstart, nstop))
result = square_me(nstart, nstop)

# Gather the squares on the root rank
total_squares = comm.gather(result, root=0)
if rank == 0:
    squares = []
    for item in total_squares:
        squares.extend(item)
    squares.sort()

# Write the results to disk
with open('total_squares.txt', 'w') as fp:
    fp.write(str(squares))