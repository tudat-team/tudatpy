from setup import StubGenerator, TUDATPY_ROOT, ChangeDir


with ChangeDir("src"):

    # Ensure that tudatpy is installed
    try:
        from tudatpy import __version__  # type: ignore
    except ImportError:
        raise ImportError("Failed to generate stubs: tudatpy is not installed")

    # Generate stubs
    StubGenerator(clean=True).generate_stubs(TUDATPY_ROOT)
