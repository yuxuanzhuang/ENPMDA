"""Top-level package for Ensemble Analysis."""

__author__ = """Yuxuan Zhuang"""
__email__ = 'yuxuan.zhuang@dbb.su.se'
__version__ = '0.1.0'

from .analysis import get_backbonetorsion, get_interatomic_distance, get_atomic_position, get_domain_position, get_domain_interdistance, get_domain_intradistance, get_c_alpha_distance, get_c_alpha_distance_filtered, get_c_alpha_distance_filtered_inverse, get_rmsd_init, get_rmsd_ref, get_pore_hydration