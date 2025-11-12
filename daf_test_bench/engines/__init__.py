"""Engine wrappers for DAF Plant implementations."""

# Import wrappers with graceful fallback if dependencies not installed
_available_wrappers = []

try:
    from .pillar3_wrapper import Pillar3PhysicsWrapper
    _available_wrappers.append("Pillar3PhysicsWrapper")
except ImportError as e:
    import warnings
    warnings.warn(f"Pillar3PhysicsWrapper unavailable: {e}")
    Pillar3PhysicsWrapper = None

try:
    from .floc_pbm_wrapper import FlocPBMWrapper
    _available_wrappers.append("FlocPBMWrapper")
except ImportError as e:
    import warnings
    warnings.warn(f"FlocPBMWrapper unavailable: {e}")
    FlocPBMWrapper = None

__all__ = [name for name in ["Pillar3PhysicsWrapper", "FlocPBMWrapper"]
           if name in _available_wrappers]
