""" Distributions functions"""
from typing import Union, Tuple
import numpy as np

EnergyData = Union[np.ndarray, float, int]
LaticeData = Union[np.ndarray, float, int]
def prepare_2D_lattice_data(kx: LaticeData, ky: LaticeData, t = -1,
                        t_prim = 0.25) -> Tuple[EnergyData,LaticeData,LaticeData]:
    """Prepare 2D lattice data."""
    kx = np.array(kx, dtype=float) if isinstance(kx, (int, float)) else kx
    ky = np.array(ky, dtype=float) if isinstance(ky, (int, float)) else ky
    KX, KY = np.meshgrid(kx, ky)

    # e_k =  2*t*(np.cos(KX)+ np.cos(KY)) + 4*t_prim*np.cos(KY/2)*np.cos(KX/2)
    e_k =  2*t*(np.cos(KX) + np.cos(KY)) + 2*t_prim*(np.cos(KX-KY)+np.cos(KX+KY))
    
    return e_k, KX, KY
