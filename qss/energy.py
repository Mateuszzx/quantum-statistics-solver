import numpy as np

import chemical_potential as cp
import distributions as dist


def energy(distro: dist.DistType, beta: float | int,
            n: float | np.ndarray, e_k: dist.DistData,
            kx: dist.LaticeData | None = None,
            ky: dist.LaticeData | None = None,
            mu: float | None = None) -> dist.EnergyData:

    if isinstance(n, float):

        if mu is None:
            mu_calc:float = cp.chemical_potential_computation(distro=distro, beta=beta, n=n, e_k=e_k, kx=kx, ky=ky)
        else:
            mu_calc = mu

        if isinstance(e_k, np.ndarray):
            result:np.ndarray = dist.distribution(dist=distro, beta=beta, e_k=e_k, mu=mu_calc)
            result = np.multiply(result, e_k)
            return result.mean()
        if isinstance(e_k,(float, int)):
            result = dist.distribution(dist=distro, beta=beta, e_k=e_k, mu=mu_calc)
            if result is not None:
                return result * e_k
            return None
        raise TypeError("Invalid input energy data")

    return np.array([energy(distro, beta, ni, e_k, kx, ky, mu_i) for ni, mu_i in zip(n,mu)])
