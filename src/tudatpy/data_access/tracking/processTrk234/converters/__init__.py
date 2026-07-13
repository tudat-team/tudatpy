__all__ = [
    "Converter",
    "RadioBase",
    "RampConverter",
    "OpenRampHandling",
    "DerivedDopplerConverter",
    "DerivedSraRangeConverter",
]


def __getattr__(name):
    if name == "Converter":
        from .converter import Converter

        return Converter
    if name == "RadioBase":
        from .radioBase import RadioBase

        return RadioBase
    if name in {"RampConverter", "OpenRampHandling"}:
        from .ramp import OpenRampHandling, RampConverter

        return {"RampConverter": RampConverter, "OpenRampHandling": OpenRampHandling}[name]
    if name == "DerivedDopplerConverter":
        from .derivedDoppler import DerivedDopplerConverter

        return DerivedDopplerConverter
    if name == "DerivedSraRangeConverter":
        from .derivedSraRange import DerivedSraRangeConverter

        return DerivedSraRangeConverter
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
