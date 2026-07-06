# Serialization Design

Status: Updated for second PR

## Purpose

This document defines the serialization approach for Tudat and tudatpy so that C++ serialization, Python pickle support, file formats, versioning, and tests all follow the same rules.

The goal is not to describe every class in detail. The goal is to establish one consistent contract that can be applied across the codebase while the remaining classes are finished.

## Current State

The repository already has a partial cereal-based serialization stack:

- `include/tudat/io/serialization/base.h` provides cereal archive helpers, Eigen support, and binary string helpers.
- `include/tudat/io/serialization.h` is the public umbrella header for the serialization layer.
- `include/tudat/io/observationSerialization.h` is a deprecated compatibility include that forwards to the new serialization base.
- C++ tests already exercise roundtrip serialization for settings, observation containers, and propagation results using cereal binary archives.
- Python bindings already expose pickle support for several classes by reusing cereal through pybind11.
- Python equality is already exposed for a number of serializable classes, which is the right foundation for roundtrip tests.

At the same time, a few pieces are still missing or inconsistent:

- A unified `save` / `load` convention is not yet enforced everywhere: some classes currently implement a single `serialize(archive ar)` method instead.
- Polymorphic registration is currently handled in binding code for some Python-exposed types, but not yet as a clearly centralized per-submodule policy.

That means the project already has working pieces, but not yet a single serialization contract that every developer can follow.

## Design Principles

1. Use cereal as the single serialization library.
2. Every serializable class SHALL implement explicit `save` and `load` methods. A generic `serialize` method must not be used.
3. Use JSON for non-tabulated settings objects.
4. Use regular binary serialization for tabulated settings and other data-heavy objects.
5. No class-level version numbers. Serialization compatibility is determined by the exact commit hash used to produce the file. Incompatible formats should fail loudly on load.
6. Register polymorphic cereal types in one place per submodule.
7. Treat roundtrip survival plus equality as the primary correctness check.
8. Expose the same roundtrip behavior to Python through pybind11.
9. Loading is done through a static method of the class. Saving is done through a method per format type.

### Why `save` and `load` instead of `serialize`

The split between `save` and `load` exists because deserialization sometimes requires a specific member construction order before the rest of the state can be restored — in particular, certain function pointers and callbacks can only be reconstructed during `load` after some members are already in place. A single `serialize` method makes this awkward.

Using `save` / `load` everywhere also ensures consistency across the codebase: polymorphic types require one method signature to be chosen, and having a uniform convention simplifies both implementation and code review. The extra lines of code this requires are minimal and can largely be generated with tooling.

## File Format Policy

The file format choice should follow the nature of the data:

- Non-tabulated settings: JSON.
- Tabulated settings: binary.

The reason for this split is practical:

- JSON is readable, diffable, and suitable for small structured settings.
- Binary is better for large numeric tables, avoids huge text files, and preserves floating-point data without text conversion.

The file format is a persistence concern. It should not be confused with Python pickle. Python pickle is an in-memory transport mechanism and does not need to use the on-disk JSON format.

### Adding A New Archive Type

If we want a new file format such as HDF5, the cereal path is to add a new archive type, register it, and then reuse the same class `save` / `load` methods.

Minimal sketch:

```cpp
// archive definition (very small sketch, not a full HDF5 implementation)
class Hdf5OutputArchive : public cereal::OutputArchive< Hdf5OutputArchive >
{
public:
	explicit Hdf5OutputArchive( H5::H5File& file ) : file_( file ) {}

	template< class T >
	Hdf5OutputArchive& operator()( const T& value )
	{
		// write value to an HDF5 dataset or attribute
		return *this;
	}

private:
	H5::H5File& file_;
};

// make cereal aware that this is an archive type
CEREAL_REGISTER_ARCHIVE( Hdf5OutputArchive )
```

```cpp
// usage stays the same in the serializable class
template < class Archive >
void save( Archive& archive ) const
{
	archive( cereal::make_nvp( "value", value_ ) );
}
```

```cpp
// file writing code chooses the archive type
H5::H5File file( "settings.h5", H5F_ACC_TRUNC );
Hdf5OutputArchive archive( file );
archive( cereal::make_nvp( "MySettings", settings ) );
```

The key idea is that the class serialization does not change when the archive changes. Only the archive wrapper and its registration are new.

## Class Serialization Contract

Each serializable class should follow the same pattern:

- Implement `save` and `load` methods.
- Never use a single generic `serialize` method.
- Every serialized field shall be wrapped in `cereal::make_nvp("field_name", value)` so that all fields are named in the archive. Bare `archive(value)` without a name is not allowed.
- Keep transient or runtime-only members out of the serialized state unless they can be reconstructed safely.
- Loading is exposed via a static method on the class, e.g. `MySettings::load_from_json(path)`.
  For polymorphic hierarchies, calling `BaseType::load_from_json(path)` returns the correct
  derived-class instance automatically thanks to cereal's polymorphic registration (the
  `CEREAL_REGISTER_TYPE` and `CEREAL_REGISTER_POLYMORPHIC_RELATION` declarations tell cereal
  which derived type to instantiate). You may also call `DerivedType::load_from_json(path)`
  directly if you know the exact type and want to be explicit about it.
- Saving is exposed via one method per format, e.g. `save_to_json(path)`, `save_to_bin(path)`, `save_to_yaml(path)`.

This is especially important for classes whose deserialization requires some members to be constructed before the rest of the state can be restored.

### Versioning and Compatibility

Serialized files carry no class version numbers. Compatibility is determined by the exact Tudat commit hash used to create the file.

If a file was created with an incompatible commit — i.e. the on-disk format has changed — loading it should fail with a clear error. The failure should look like:

> "Could not deserialize this object. This file was created with commit `<hash>`."

The Tudat releasing schedule is the mechanism for providing conversion scripts. When a breaking change to the serialized format is introduced between releases, the team may provide a migration script alongside the release notes. There is no built-in version branching inside `load` methods.

### Metadata

Every serialized file should carry metadata that identifies:

- The exact commit hash of Tudat / tudatpy used to generate the file.
- The serialization format family, if needed for future debugging.

This metadata should be stored in the serialized file itself, not in operating-system file properties. Cross-platform file contents are the only reliable place for it.

If the metadata commit hash does not match a known compatible set, the loader reports the mismatch and refuses to proceed.

## Polymorphism

Polymorphic types must be registered in one canonical place per submodule.

That registration layer should contain the cereal type registration and the base-to-derived relationships required for pointer serialization.

Cereal requires polymorphic types to be registered so it knows which derived type to instantiate during deserialization. While cereal does not restrict this registration to a single translation unit, centralizing it per submodule makes the registrations easy to audit and avoids accidentally omitting a type. Spreading registrations across unrelated files (such as binding code) makes it easy to miss a type and hard to review the complete picture.

Recommended layout:

- One dedicated registration header per submodule.
- All derived cereal registrations for that submodule in that header.
- Binding code and any other serialization users include that header.
- Only use a source-file registration TU when there is a strong reason to hide implementation details or control archive-specific linking behavior.

For this codebase, the header-based option is the default because it is easier to audit, easier to include from pybind modules, and less likely to create platform-specific linkage surprises.

Minimal example:

```cpp
// include/tudat/io/serialization/registrations/my_settings_registration.h
#pragma once

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include "my_settings.h"

CEREAL_REGISTER_TYPE( MyDerivedSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( MyBaseSettings, MyDerivedSettings )
```

```cpp
// src/tudatpy/some_module/expose_my_settings.cpp
#include "tudat/io/serialization/registrations/my_settings_registration.h"

// the bindings can now use the registered polymorphic types safely
```

This keeps the registration in one obvious place while allowing the pybind or C++ users of the type to include the same header wherever they need the polymorphic metadata.

## Python Exposure

Python exposure should remain aligned with the C++ contract.

### Persistence Policy

Pickle is a convenience API, not the primary persistence mechanism.

Recommended policy:

- Use pickle for roundtrips, multiprocessing, caching, and tests.
- Use explicit file-oriented save/load methods for durable storage.
- Expose JSON-saving methods for non-tabulated settings.
- Expose binary-saving methods for tabulated or data-heavy objects.

This keeps the Python API honest about what format is being used and avoids presenting pickle as the canonical file format.

### Pickle

Pickle support should be exposed through pybind11.

For value types, pickle can serialize the object state directly.

For polymorphic types, pickle should preserve the dynamic type through a base pointer, not flatten the object to the base class.

The current helper approach is the right shape:

- non-polymorphic types use a value-based pickle helper,
- polymorphic types use a shared-pointer-based pickle helper.

### Equality

The equality operator must be exposed to Python for every class that participates in roundtrip tests.

Recommended policy:

- The public Python API exposes `__eq__`.
- The C++ class exposes `operator==`.
- The class implementation should use a protected `equals` helper where inheritance makes that easier.
- For polymorphic hierarchies, the base class should expose the operators and the derived classes can specialize the equality implementation.

This matches the existing pattern in several parts of the codebase and keeps the test surface compact.

### JSON Settings and Python Pickle

Python pickle should not be coupled to JSON file encoding. Pickle should serialize the in-memory C++ object state directly, using cereal through pybind11. JSON remains the on-disk interchange format for settings files.

This separation keeps pickle simple, avoids forcing JSON parsing into Python roundtrips, and allows future changes to the file format without breaking Python pickle semantics.

If a class needs a Python-visible persistence API beyond pickle, add explicit file export/import functions instead of overloading pickle with file-format concerns.

### Python File API

The preferred Python-facing persistence API should be explicit and format-specific.

Recommended method names:

- `save_to_json(path)` / static `load_from_json(path)` for JSON-backed settings objects.
- `save_to_binary(path)` / static `load_from_binary(path)` for binary-backed objects.
- `save_to_yaml(path)` / static `load_from_yaml(path)` for YAML-backed objects.
- `to_json_string()` / `from_json_string()` only if an in-memory JSON representation is genuinely useful.

The file suffix (`.json`, `.bin`, `.yaml`, etc.) should be appended automatically by the file writer — the caller provides a stem or full path without needing to specify the extension.

Pickle should remain available, but it should not be the method users are told to use when they want to save a file for later reuse.

### Minimal Example

The intended pattern is:

```cpp
// C++ class: only one save/load pair is needed for cereal.
class MySettings
{
public:
	template <class Archive>
	void save( Archive& archive ) const
	{
		archive( cereal::make_nvp( "value", value_ ) );
	}

	template <class Archive>
	void load( Archive& archive )
	{
		archive( cereal::make_nvp( "value", value_ ) );
	}

	// Static factory method for loading
	static MySettings load_from_json( const std::string& path )
	{
		std::ifstream stream( path );
		cereal::JSONInputArchive archive( stream );
		MySettings obj;
		archive( obj );
		return obj;
	}

	static MySettings load_from_binary( const std::string& path )
	{
		std::ifstream stream( path, std::ios::binary );
		cereal::BinaryInputArchive archive( stream );
		MySettings obj;
		archive( obj );
		return obj;
	}

	// Per-format save methods
	void save_to_json( const std::string& path ) const
	{
		std::ofstream stream( path );
		cereal::JSONOutputArchive archive( stream );
		archive( cereal::make_nvp( "MySettings", *this ) );
	}

	void save_to_binary( const std::string& path ) const
	{
		std::ofstream stream( path, std::ios::binary );
		cereal::BinaryOutputArchive archive( stream );
		archive( cereal::make_nvp( "MySettings", *this ) );
	}

	int value_ = 0;
};
```

```cpp
// Python exposure: explicit file methods, implemented as thin wrappers.
py::class_< MySettings, std::shared_ptr< MySettings > >( m, "MySettings" )
	.def( py::init<>() )
	.def( "save_to_json", &MySettings::save_to_json )
	.def( "save_to_binary", &MySettings::save_to_binary )
	.def_static( "load_from_json", &MySettings::load_from_json )
	.def_static( "load_from_binary", &MySettings::load_from_binary )
	.def( py::pickle(
		[]( const MySettings& obj )
		{
			return py::bytes( tudat::serialization::serializeToBinaryString( obj ) );
		},
		[]( py::bytes data )
		{
			return tudat::serialization::deserializeFromBinaryString< MySettings >( data.cast< std::string >() );
		} ) );
```

The rules this example demonstrates are:

- cereal `save` / `load` live on the C++ class.
- Loading is exposed as a static method of the class.
- Saving is exposed as one method per format type.
- pickle stays a separate binary convenience path.

## File Extensions

Recommended policy:

- Use `.json` for non-tabulated settings files.
- Use `.bin` for binary tabulated files.
- Use `.yaml` for YAML settings files.

The suffix should be appended automatically by the file writer. The caller passes a path or stem, and the writer adds the appropriate extension based on the format method called.

If the team wants a more explicit binary extension, a domain-specific suffix such as `.tudatbin` is acceptable, but it is not required.

The important part is consistency across the repository, not the exact naming choice.

## Testing Policy

The primary test criterion should be roundtrip survival.

### C++ Tests

The C++ side should verify:

- The object can be serialized and deserialized.
- The restored object compares equal to the original object.
- The restored object keeps the correct dynamic type when polymorphism is involved.

Where possible, the tests should exercise real constructors and representative sample data instead of mocked internal state.

### Python Tests

The Python side should verify:

- Pickle roundtrip survival.
- Equality after unpickling.
- Concrete type preservation after unpickling for polymorphic classes.

At this stage, Python tests only need to confirm that the object survives the roundtrip. Later functional tests can assert more behavior once the serialization surface is stable.

### Equality Pattern

The recommended test pattern is:

1. Serialize the object.
2. Deserialize it.
3. Compare the result with `==`.
4. For polymorphic types, also verify `type(restored) is type(original)` on the Python side or an equivalent C++ dynamic-cast check.

This keeps the tests simple and makes version mismatches or missing registrations fail immediately.

## Current Implementation Limitations

The current codebase still needs a few follow-up decisions and implementation passes:

- A standardized `save` / `load` pattern is not yet enforced everywhere; some classes use `serialize`.
- Static `load_from_*` and per-format `save_to_*` methods are not yet implemented.
- Type registration is not yet visibly centralized across all submodules.
- Some classes still need a final decision on whether they are serializable at all if they contain callbacks or other non-serializable runtime state.

One concrete example is callback-heavy or function-heavy settings objects. Those should either provide custom serialization logic or be explicitly excluded with a documented reason.

## Rollout Plan

1. Convert all remaining classes that use `serialize` to explicit `save` / `load` methods.
2. Add static `load_from_*` methods and per-format `save_to_*` methods to all serializable classes.
3. Add commit-hash metadata to serialized output.
4. Move polymorphic registration into one header per submodule.
5. Expand C++ roundtrip tests for the remaining classes.
6. Expand Python pickle tests to cover the remaining exposed types.

## Acceptance Criteria

The serialization feature is ready when:

- Every intended serializable class has a `save` / `load` implementation.
- Every intended serializable class exposes static `load_from_*` and per-format `save_to_*` methods.
- Every intended polymorphic type is registered in the correct submodule registration header.
- Every public serializable class has a roundtrip equality test in C++.
- Every Python-exposed serializable class has a pickle roundtrip test.
- File metadata records the exact commit hash used to create the file.
- Loading a file created with an incompatible commit produces a clear error message identifying the source commit.