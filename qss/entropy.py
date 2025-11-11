"""Entropy functions"""
import numpy as np

import chemical_potential as cp
import distributions as dist


# Helper function to safe get logarithm
def safe_log(x: np.ndarray) -> np.ndarray:
    """Compute the logarithm safely, avoiding log(0) issues."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x > 0, np.log(x), 0)


def calculate_entropy_formula(nk: np.ndarray, distro: dist.DistType) -> np.ndarray:
        """Helper function to calculate entropy based on distribution."""
        if distro == "FD":
            nk = nk/2
            return 2*(-nk * safe_log(nk)
                    - (1 - nk) * safe_log(1 - nk))

        if distro == "FLB":
            # nk = nk/2
            return (-nk * safe_log(nk)
                    + 2 * nk * safe_log(2)
                    - 2 * (1 - nk/2) * safe_log(1 - nk/2))
        if distro in ["NFL_i", "NFL"]:
            return (-nk  * safe_log(nk)
                    + 2 * nk * safe_log(2)
                    + 2 * (1 - nk / 2) * safe_log(1 - nk / 2)
                    - 2 * (1 - nk ) * safe_log(1 - nk ))
        raise ValueError("Unknown distribution")



def calculate_entropy(n: float | np.ndarray, e_k: np.ndarray, beta: float, distro: dist.DistType,
                    kx: np.ndarray | None = None, ky: np.ndarray | None = None, mu: float | None = None,
                    dim: int = 2) -> float | np.ndarray:
    """
    Calculate the entropy for a given distribution and parameters.

    Args:
        n (Union[float, np.ndarray]): Density value(s).
        e_k (np.ndarray): Wave vector array.
        beta (float): Inverse temperature.
        distro (dist.DistType): Distribution type.
        kx (Optional[np.ndarray]): Optional x-axis grid for 2D.
        ky (Optional[np.ndarray]): Optional y-axis grid for 2D.
        dim (int): Dimension (1 or 2).

    Returns:
        Union[float, np.ndarray]: Calculated entropy.
    """


    if isinstance(n, (int, float)):
        if mu is None:
            mu = cp.chemical_potential_computation(distro=distro, beta=beta, n=n, dim=dim, e_k=e_k, kx=kx, ky=ky)
            print(f'mu: {mu:.3f}, beta: {beta}, distro: {distro}')
            print(f" n_computed: {cp.compute_density(mu=mu, dist_type=distro, beta=beta, e_k=e_k, dim=2, kx = kx, ky = ky):.2f}, n_assigned: {n:.2f}")

        if dim == 2 and kx is not None and ky is not None:
            nk = dist.distribution(dist=distro, e_k=e_k, mu=mu, beta=beta)

        elif dim == 2 and kx is None and ky is None:
            nk = dist.distribution_symmetry_2d(dist=distro, e_k=e_k, mu=mu, beta=beta)

        elif dim == 1:
            nk = dist.distribution(dist=distro, e_k=e_k, mu=mu, beta=beta)
        else:
            raise ValueError("Unknown dimension")

        entropy = calculate_entropy_formula(nk, distro=distro)

        if kx is None or ky is None:
            entropy = cp.left_edge_from_samples_1d(entropy , e_k)
        else:
            entropy = cp.left_edge_from_samples_2d(entropy, kx, ky)

        if mu is None:
            print(f"Entropy: {entropy:.3f} \n")
        return entropy

    if isinstance(n, np.ndarray):
        return np.array([calculate_entropy(n=n_i, e_k=e_k, beta=beta, distro=distro, dim=dim, kx=kx, ky=ky, mu=mu_i) for n_i, mu_i in zip(n,mu)])

    raise ValueError("Invalid density value, should be a float or array of floats.")
