import numpy as np
def reshape_distances(distances_array):
    distances_array_reshaped = np.zeros(
        (
            len(distances_array),
            12,
            distances_array.shape[2] if len(distances_array.shape) > 2 else 1,
        ),
        dtype=distances_array.dtype,
    )
    for en, x in enumerate(distances_array):
        if len(x.shape) == 1:
            x = x.reshape(-1, 1)
        for v in range(x.shape[1]):
            for e in range(x.shape[0]):
                if e < 4:
                    distances_array_reshaped[en, e, v] = -(np.sum(x[e:5, v]))
                elif e == 4:
                    distances_array_reshaped[en, e, v] = -x[e, v]
                    distances_array_reshaped[en, e + 1, v] = 0
                else:
                    distances_array_reshaped[en, e + 1, v] = (
                        distances_array_reshaped[en, e, v] + x[e, v]
                    )
    return distances_array_reshaped.squeeze()