# Serialization Design

Status: draft for team alignment

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

- A centralized class version manifest does not yet exist in `src/` or `include/`.
- `CEREAL_CLASS_VERSION` declarations are not yet present in the source tree.
- A unified `save` / `load` convention is not yet enforced everywhere: some classes currently implement a single `serialize(archive ar)` method instead.
- Polymorphic registration is currently handled in binding code for some Python-exposed types, but not yet as a clearly centralized per-submodule policy.

That means the project already has working pieces, but not yet a single serialization contract that every developer can follow.

## Design Principles

1. Use cereal as the single serialization library.
2. Every serializable class SHALL implement explicit `save` and `load` methods. A generic `serialize` method must not be used.
3. Use XML for non-tabulated settings objects.
4. Use regular binary serialization for tabulated settings and other data-heavy objects.
5. Every serializable class SHALL declare a version number using `CEREAL_CLASS_VERSION`.
6. The software SHALL maintain a central manifest listing the current serialization version of every serializable class.
7. Register polymorphic cereal types in one place per submodule.
8. Treat roundtrip survival plus equality as the primary correctness check.
9. Expose the same roundtrip behavior to Python through pybind11.

### Why `save` and `load` instead of `serialize`

The split between `save` and `load` exists because deserialization sometimes requires a specific member construction order before the rest of the state can be restored — in particular, certain function pointers and callbacks can only be reconstructed during `load` after some members are already in place. A single `serialize` method makes this awkward.

Using `save` / `load` everywhere also ensures consistency across the codebase: polymorphic types require one method signature to be chosen, and having a uniform convention simplifies both implementation and code review. The extra lines of code this requires are minimal and can largely be generated with tooling.

## File Format Policy

The file format choice should follow the nature of the data:

- Non-tabulated settings: XML.
- Tabulated settings: binary.

The reason for this split is practical:

- XML is readable, diffable, and suitable for small structured settings.
- Binary is better for large numeric tables, avoids huge text files, and preserves floating-point data without text conversion.

The file format is a persistence concern. It should not be confused with Python pickle. Python pickle is an in-memory transport mechanism and does not need to use the on-disk XML format.

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
- Include version-aware branching inside the class methods when the format changes.
- Keep transient or runtime-only members out of the serialized state unless they can be reconstructed safely.

This is especially important for classes whose deserialization requires some members to be constructed before the rest of the state can be restored.

### Versioning

Every serializable class must declare a version identifier using the `CEREAL_CLASS_VERSION` macro.

The version scheme is: `TUDATVERSION::CLASSVERSION.CLASSPATCH`

For example, all new classes at initial release would carry something like `1.0.0::1.0`. The Tudat/tudatpy version prefix ties the class version to the software release it was introduced in, while the class-level `CLASSVERSION.CLASSPATCH` tracks independent evolution of that class's serialized format.

```cpp
// example version declaration
CEREAL_CLASS_VERSION( MySettings, /* version integer registered with cereal */ )
```

The current version for all classes is published in one central manifest (see below). The version number is written into the file metadata and used by class `load` methods to branch when reading older files.

### Central Version Manifest

The software SHALL maintain a single, centralized manifest that lists every serializable class and its current declared version. This makes it immediately visible when a class version changes, which matters because such a change may break compatibility with files users have already saved.

When a version does change, the team must decide whether to provide a migration path or require users to regenerate their files. There is no generic migration mechanism: every version bump that affects the on-disk format will require a custom conversion solution, which may sometimes be straightforward (field reordering) and sometimes require additional user input (new required fields with no default). For this reason the manifest is the canonical place to track and review such changes.

### Metadata

Every serialized file should carry metadata that identifies:

- The object or class version.
- The Tudat / tudatpy version used to generate the file.
- The serialization format family, if needed for future debugging.

This metadata should be stored in the serialized file itself, not in operating-system file properties. Cross-platform file contents are the only reliable place for it.

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
#include <cereal/archives/xml.hpp>
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
- Expose XML-saving methods for non-tabulated settings.
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

### XML Settings and Python Pickle

Python pickle should not be coupled to XML file encoding. Pickle should serialize the in-memory C++ object state directly, using cereal through pybind11. XML remains the on-disk interchange format for settings files.

This separation keeps pickle simple, avoids forcing XML parsing into Python roundtrips, and allows future changes to the file format without breaking Python pickle semantics.

If a class needs a Python-visible persistence API beyond pickle, add explicit file export/import functions instead of overloading pickle with file-format concerns.

### Python File API

The preferred Python-facing persistence API should be explicit and format-specific.

Recommended method names:

- `save_xml(path)` / `load_xml(path)` for XML-backed settings objects.
- `save_binary(path)` / `load_binary(path)` for binary-backed objects.
- `to_xml_string()` / `from_xml_string()` only if an in-memory XML representation is genuinely useful.

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

	int value_ = 0;
};

CEREAL_CLASS_VERSION( MySettings, /* version integer */ )
```

```cpp
// Python exposure: explicit file methods, implemented as thin wrappers.
py::class_< MySettings, std::shared_ptr< MySettings > >( m, "MySettings" )
	.def( py::init<>() )
	.def( "save_xml", []( const MySettings& obj, const std::string& path )
	{
		std::ofstream stream( path );
		cereal::XMLOutputArchive archive( stream );
		archive( cereal::make_nvp( "MySettings", obj ) );
	} )
	.def_static( "load_xml", []( const std::string& path )
	{
		std::ifstream stream( path );
		cereal::XMLInputArchive archive( stream );

		MySettings obj;
		archive( obj );
		return obj;
	} )
	.def( "save_binary", []( const MySettings& obj, const std::string& path )
	{
		std::ofstream stream( path, std::ios::binary );
		cereal::BinaryOutputArchive archive( stream );
		archive( cereal::make_nvp( "MySettings", obj ) );
	} )
	.def_static( "load_binary", []( const std::string& path )
	{
		std::ifstream stream( path, std::ios::binary );
		cereal::BinaryInputArchive archive( stream );

		MySettings obj;
		archive( obj );
		return obj;
	} )
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

The rule this example demonstrates is simple:

- cereal `save` / `load` live on the C++ class,
- Python file persistence is an explicit wrapper around the archive type,
- pickle stays a separate binary convenience path.

## File Extensions

Recommended policy:

- Use `.xml` for non-tabulated settings files.
- Use `.bin` for binary tabulated files.

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

- A central version manifest does not yet exist.
- A standardized `save` / `load` pattern is not yet enforced everywhere; some classes use `serialize`.
- Type registration is not yet visibly centralized across all submodules.
- Some classes still need a final decision on whether they are serializable at all if they contain callbacks or other non-serializable runtime state.

One concrete example is callback-heavy or function-heavy settings objects. Those should either provide custom serialization logic or be explicitly excluded with a documented reason.

## Rollout Plan

1. Add the central version manifest, using the `TUDATVERSION::CLASSVERSION.CLASSPATCH` scheme.
2. Convert all remaining classes that use `serialize` to explicit `save` / `load` methods.
3. Add `CEREAL_CLASS_VERSION` declarations to all serializable classes.
4. Move polymorphic registration into one header per submodule.
5. Add file metadata support for version and Tudat/tudatpy version.
6. Expand C++ roundtrip tests for the remaining classes.
7. Expand Python pickle tests to cover the remaining exposed types.
8. Add migration tests for at least one older file format version per major object family.

## Acceptance Criteria

The serialization feature is ready when:

- Every intended serializable class has a versioned `save` / `load` implementation using `CEREAL_CLASS_VERSION`.
- The central version manifest exists and lists every serializable class.
- Every intended polymorphic type is registered in the correct submodule registration header.
- Every public serializable class has a roundtrip equality test in C++.
- Every Python-exposed serializable class has a pickle roundtrip test.
- Old files can be loaded according to their declared version.
- File metadata records both the file version and the Tudat/tudatpy version used to create the file.