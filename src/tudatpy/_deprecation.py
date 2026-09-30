import warnings


def deprecation_warning(old_name: str, new_name: str) -> None:
    warnings.warn(
        f"{old_name} is deprecated. Use {new_name} instead.",
        stacklevel=3,
    )
