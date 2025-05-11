import numpy as np

#mat_data = np.fromfile('mat-d20-b5-p4.bin')
#vec_data = np.fromfile('x-d20.txt.bin')

with open("x-d20.txt.bin", "rb") as f:
    len_vec = np.fromfile(f, dtype=np.int32, count=1)[0]
    vec = np.fromfile(f, dtype=np.float64, count=len_vec)

with open("mat-d20-b5-p4.bin", "rb") as f:
    len_mat = np.fromfile(f, dtype=np.int32, count=1)[0]
    mat = np.fromfile(f, dtype=np.float64, count=len_mat*len_mat)

#vec = np.reshape(vec, (len_mat, 1))
mat = np.reshape(mat, (4, 20, 5))
mat = np.block(list(mat))

print(mat)
print(vec)


print(mat @ vec)
