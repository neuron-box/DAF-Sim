"""
DAF Test Bench Framework

A unified testing system for Dissolved Air Flotation (DAF) simulation engines.
Implements the pluggable architecture defined in DAF-TB-HLD-v1.0.

Components:
- IDAFPlant: Abstract interface for all DAF engines
- Engine Wrappers: Adapters for specific DAF Plant implementations
- Test Harness: Main orchestration and benchmarking driver
- Metrics Module: Performance metric collection and reporting
- V&V Module: Verification and validation tools
"""

__version__ = "1.0.0"
__author__ = "DAF-Sim Development Team"

from .interfaces.idaf_plant import IDAFPlant
from .metrics.metrics_collector import MetricsCollector, BenchmarkResult
from .test_harness import TestHarness

__all__ = [
    "IDAFPlant",
    "MetricsCollector",
    "BenchmarkResult",
    "TestHarness",
]
