# -*- coding: utf-8 -*-
"""monte_carlo_wo_list_wo_np_arr.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1krWggGEhp8Irn1dd65oP49mzyPzti6k2
"""

#
#!pip install matplotlib

import numpy as np
from typing import Union, Sequence

# noinspection SpellCheckingInspection
def mod(
    arr: np.ndarray, dim: Union[int, Sequence[int]] = -1, keepdims: bool = True
):
    return np.sum(a=(arr ** 2), axis=dim, keepdims=keepdims) ** 0.5


def norm(
        arr: np.ndarray,
        dim: Union[int, Sequence[int]] = -1,
        fill_value: Union[int, float] = 0
) -> np.ndarray:
    x2 = mod(arr=arr, dim=dim, keepdims=True)
    x1 = arr.astype(x2.dtype)
    result = np.full_like(a=arr, fill_value=fill_value)
    return np.divide(
        x1,
        x2,
        out=result,
        where=(x2 != 0)
    )

def hamming(
        big_j: float,
        big_l: int,
        spin: np.ndarray,
        a_z: float,
        i: Union[int, np.ndarray],
        j: Union[int, np.ndarray],
        k: Union[int, np.ndarray]
) -> Union[float, np.ndarray]:
    h0 = (
        big_j * spin[k, i, j, 0] * (
            spin[k, ((i + 1) % big_l), j, 0] +
            spin[k, ((i - 1) % big_l), j, 0] +
            spin[k, i, ((j + 1) % big_l), 0] +
            spin[k, i, ((j - 1) % big_l), 0]
        )
    ) - (a_z * spin[k, i, j, 2] ** 2)

    h1 = (
        big_j * spin[k, i, j, 1] * (
            spin[k, ((i + 1) % big_l), j, 1] +
            spin[k, ((i - 1) % big_l), j, 1] +
            spin[k, i, ((j + 1) % big_l), 1] +
            spin[k, i, ((j - 1) % big_l), 1]
        )
    ) - (a_z * spin[k, i, j, 2] ** 2)

    h2 = (
        big_j * spin[k, i, j, 2] * (
            spin[k, ((i + 1) % big_l), j, 2] +
            spin[k, ((i - 1) % big_l), j, 2] +
            spin[k, i, ((j + 1) % big_l), 2] +
            spin[k, i, ((j - 1) % big_l), 2]
        )
    ) - (a_z * spin[k, i, j, 2] ** 2)

    big_h = -(h0 + h1 + h2)
    return big_h

def simulate(
        big_j: float, big_l: int, big_n: int, a_z: float, n_eq: int, t_max: int
) -> Sequence[float]:
    # Generate all kt: (t_max,)
    kt = 0.4 + (np.arange(0, t_max) / 10)

    # Randomly initialize spin: (t_max, big_l, big_l, 3)
    spin = np.repeat(
        a=norm(
            arr=np.random.uniform(
                low=-1.0,
                high=1.0,
                size=(1, big_l, big_l, 3)
            ),
            dim=-1,
            fill_value=0
        ),
        repeats=t_max,
        axis=0
    )

    # Generate indexes: ((t_max * big_l * big_l),), ((t_max * big_l * big_l),)
    idx_arr = np.arange(big_l)
    k_arr, i_arr, j_arr = np.meshgrid(np.arange(t_max), idx_arr, idx_arr)
    k_arr, i_arr, j_arr = k_arr.ravel(), i_arr.ravel(), j_arr.ravel()

    # Calculate energy & magnitude: (t_max,), (t_max, 3)
    energy = 0.5 * np.sum(
        hamming(
            big_j=big_j,
            big_l=big_l,
            spin=spin,
            a_z=a_z,
            i=i_arr,
            j=j_arr,
            k=k_arr
        ).reshape((t_max, -1)),
        axis=-1,
        keepdims=False
    )
    mag = np.repeat(
        a=np.expand_dims(
            a=np.sum(a=spin, axis=(1, 2), keepdims=False),
            axis=1
        ),
        repeats=big_n,
        axis=1
    )

    # Randomly generate indexes
    # for `big_n` trials: ((t_max * big_n * big_l * big_l),)
    rnd_i = np.random.randint(
        low=0, high=(big_l - 1), size=(t_max * big_n * (big_l ** 2))
    )
    rnd_j = np.random.randint(
        low=0, high=(big_l - 1), size=(t_max * big_n * (big_l ** 2))
    )
    rnd_k = np.random.randint(
        low=0, high=t_max, size=(t_max * big_n * (big_l ** 2))
    )

    # Preserve spins & energy at indexes:
    # (t_max, big_n, (big_l * big_l), 3), (t_max, big_n, (big_l * big_l))
    spin_initial = spin[rnd_k, rnd_i, rnd_j, :].reshape(
        (t_max, big_n, (big_l ** 2), 3)
    )
    energy_initial = hamming(
        big_j=big_j, big_l=big_l, spin=spin, a_z=a_z, i=rnd_i, j=rnd_j, k=rnd_k
    ).reshape((t_max, big_n, (big_l ** 2)))

    # Randomly modify spin at indexes
    spin[rnd_k, rnd_i, rnd_j, :] = norm(
        np.random.uniform(
            low=-1.0, high=1.0, size=((t_max * big_n * (big_l ** 2)), 3)
        ),
        dim=-1,
        fill_value=0
    )

    # Calculate energy after modifying spins:
    # (t_max, big_n, (big_l * big_l), 3), (t_max, big_n, (big_l * big_l))
    spin_final = spin[rnd_k, rnd_i, rnd_j, :].reshape(
        (t_max, big_n, (big_l ** 2), 3)
    )
    energy_final = hamming(
        big_j=big_j, big_l=big_l, spin=spin, a_z=a_z, i=rnd_i, j=rnd_j, k=rnd_k
    ).reshape((t_max, big_n, (big_l ** 2)))

    # Calculate ??(energy) & ??(magnitude):
    # (t_max, big_n, (big_l ** 2)), (t_max, big_n, (big_l ** 2), 3)
    delta_energy = energy_final - energy_initial
    delta_mag = spin_final - spin_initial

    # Make kt broadcastable and calculate exp(-??(energy) / kt):
    # (t_max, 1, 1), ((t_max, big_n, (big_l * big_l)))
    kt_tensor = np.expand_dims(a=kt, axis=tuple(range(1, delta_energy.ndim)))
    exp_arr = np.exp(-delta_energy / kt_tensor)

    # Find condition: ??(energy) < 0 and stochastically decide otherwise:
    # (t_max, big_n, (big_l ** 2))
    condition = delta_energy <= 0
    complement_condition = np.logical_not(condition)
    condition[complement_condition] = (
        np.random.uniform(
            low=0, high=1, size=complement_condition.sum()
        ) < exp_arr[complement_condition]
    )

    # Set ??(energy) = 0 where condition is not satisfied:
    # (t_max, big_n, (big_l ** 2))
    delta_energy = delta_energy * condition

    # Reduce ??(energy) by summing along trials and experiment dimension
    # and add with energy
    energy = energy + np.sum(
        a=delta_energy,
        axis=tuple(range(1, delta_energy.ndim)),
        keepdims=False
    )

    # Reduce ??(mag) by summing along trials and experiment dimension
    # and add with mag
    mag = mag + np.sum(
        a=delta_mag,
        axis=2,
        keepdims=False
    )

    mag_avg = np.abs(
        np.mean(
            a=mag[:, n_eq:, :],
            axis=1,
            keepdims=False
        )
    )
    mag_avg_total = mod(arr=mag_avg, dim=1, keepdims=False)
    return energy, kt, mag_avg_total

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    big__l = 10
    big__n = 1000
    big__j = 1.0
    a__z = 0.4
    n__eq = 500
    t__max = 11
    a, b, c = simulate(big__j, big__l, big__n, a__z, n__eq, t__max)
    plt.plot(b, a)
    plt.show()
