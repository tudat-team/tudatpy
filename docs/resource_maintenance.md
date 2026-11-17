# Maintaining the Zenodo resource records

This guide is for resource maintainers who need to add, replace or remove a
Tudat resource. Downloading resources does not require a Zenodo account; changing
the published records requires edit permission on those records. Request access
from their owners if the editing controls are unavailable.

## Records and installer contract

The authoritative record/file URLs are `RESOURCE_METADATA_URLS` in
[`src/tudatpy/resource_catalog.py`](../src/tudatpy/resource_catalog.py). At the time
of writing they select:

| Resource set | Published record | Manifest | Checksums | Payloads |
|---|---|---|---|---|
| Rolling | [21261530](https://zenodo.org/records/21261530) | `manifest-rolling.txt` | `SHA256SUM-rolling` | `resource-rolling.tar.gz` |
| Static | [21277280](https://zenodo.org/records/21277280) | `manifest-static.txt` | `SHA256SUM-static` | Separate `.tar.gz` archives, listed in the manifest |

Rolling currently contains Earth-orientation and space-weather data. Static
contains atmosphere, gravity, ephemeris, SPICE and other reference data. Follow
the existing manifest grouping. An addition to an existing Earth-orientation
series belongs with the other Earth-orientation files; a directory name alone
is not a rule for assigning all future datasets. Agree on new categories with
the resource maintainers. A resource path must occur only once across both sets.

A Zenodo record is the publication containing a manifest, a checksum file and
one or more tar archives. Individual resource files are **members inside those
archives**. Uploading a loose resource beside an unchanged archive does not
make it installable through the catalog.

## 1. Choose how to publish the change

There are two routes. Choose one before generating the manifest, because its
first line identifies the record containing its archives.

**Update the existing record's files**, if Zenodo permits it. The record number
and DOI stay the same. In the record editor, use **Edit files → Edit published
files** when available, follow the confirmation steps, then save and republish
the coordinated changes. Metadata editing alone does not replace files.
Zenodo's [file-management guide](https://help.zenodo.org/docs/deposit/manage-files/#modify)
limits ordinary minor file corrections to its correction window (currently
30 days); older records normally require a new version. Do not assume that the
word "rolling" grants permanent permission to overwrite published files.

**Publish a new version** for routine changes outside that window or when
preserving the previous file set. Open the existing record, choose **New
version**, then **Import files** to retain unchanged uploaded archives. Do not
create an unrelated upload. Obtain/reserve the new version DOI in the draft
before preparing its manifest. Use the DOI assigned to this version, not the
concept DOI that represents all versions. Follow Zenodo's
[versioning instructions](https://help.zenodo.org/docs/deposit/manage-versions/)
and [DOI reservation instructions](https://help.zenodo.org/docs/deposit/create-new-upload/#metadata).
A new version has a new record number and version DOI; the concept DOI links
the version history.

For either route, keep a local copy of the previous manifest, checksum file
and any archive you change. Use a separate staging directory. Record which
resources changed, their scientific source, coverage, units, format and license.
For derived data, retain the conversion procedure and source version as well.

## 2. Add, replace or remove resource files

Download the selected record's manifest, checksum file and the archive(s) you
need to change. Use the public file links from the record; do not download all
large static archives when changing just one. Retain unchanged archives and
their checksum entries exactly.

Inspect the archive layout with `tar -tzf ARCHIVE.tar.gz`, and extract the
archive into a fresh staging directory. Preserve every unchanged resource in
that archive. The new archive should contain canonical paths relative to the
installed resource root: for example
`earth_orientation/eopc01_polar_motion_1846-1961.txt`, without a leading slash,
`..`, `resource/` or an extra staging-directory prefix. Older archives may have
additional directory prefixes; normalize them against the manifest when
preparing a replacement. Do not rely on the installer's basename fallback to
resolve ambiguous paths.

For an **addition**, put the genuine file at its intended relative path and
add its manifest row under the archive that will contain it. For a
**replacement**, retain its path and replace its bytes and provenance as needed.
For a **removal**, remove the member and its manifest row. Before removing a
resource, check consumers and tests: a newer record must still supply every
file required by the Tudatpy versions that will be directed to it. If an archive
becomes empty, remove its archive-header row and upload entry as well.

For example, when adding a supplied C01 polar-motion file to rolling resources,
stage it under `payload/earth_orientation/`, alongside **all other members** of
that archive. From a working directory containing that `payload/` directory:

```bash
tar -czf resource-rolling.tar.gz -C payload earth_orientation space_weather
tar -tzf resource-rolling.tar.gz
```

This command assumes those are the complete top-level directories of the
selected rolling archive. Check its manifest first; include any other listed
directories too. For static resources, rebuild only the affected archive with
its complete member set, for example `gravity_models.tar.gz` from
`payload/gravity_models/`. Do not change numerical data merely to get a test to
pass. A successful install checks packaging, not scientific validity.

## 3. Update the manifest

Manifests use CSV, with a special first line containing the record's DOI:

```text
https://doi.org/10.5281/zenodo.RECORD_ID
resource-rolling.tar.gz,2026-11-13T00:00:00+00:00,
earth_orientation/FILE.txt,2026-11-13T00:00:00+00:00,https://SOURCE/FILE.txt
```

`RECORD_ID`, `FILE.txt` and `https://SOURCE/FILE.txt` are placeholders. Substitute
the actual publication record, resource path and provenance URL. Use a real
modification timestamp with a timezone; preserve timestamps for unchanged data.

- An archive-header row has three fields: archive filename, timestamp and an
  **empty** third field. Resource rows following it belong to that archive until
  the next archive-header row.
- A resource row has path, timestamp and a **nonempty source URL**. Use CSV
  quoting where necessary. The source URL records provenance; it is not normally
  the URL used to install the archived bytes.
- The installer reads the record number from the first line and constructs
  `https://zenodo.org/api/records/RECORD_ID/files/ARCHIVE/content` for each group.
  In a new-version draft, update this first line even for unchanged imported
  archives. Keep the other record's manifest unchanged if it is not updated.
- For historical `raw.githubusercontent.com/tudat-team/tudat-resources/` source
  URLs containing `/resource/`, the parser derives the install path from that
  URL. Ensure it agrees with the manifest path; copying an unrelated old URL
  can silently select the wrong path.

The current static parser recognizes a legacy malformed star-catalog header
and maps it to `star_catalog_biases.tar.gz`. New publications should spell the
actual uploaded filename correctly instead of reproducing that typo.

## 4. Update and check SHA-256 checksums

`SHA256SUM-rolling` and `SHA256SUM-static` use standard `sha256sum` format:

```text
HEX_DIGEST  resource-rolling.tar.gz
HEX_DIGEST  earth_orientation/FILE.txt
```

Replace `HEX_DIGEST` with the 64-character digest. Include both the complete
compressed archive and each resource member, using canonical resource-relative
paths for members. For each changed archive:

1. Remove its old archive digest and all its old member-digest entries from a
   working copy of the checksum file.
2. Generate entries for the rebuilt archive and its complete current member
   set; merge these with the unchanged archives' entries. Do not append
   conflicting old and new digests. Removed resources must have no stale entry.
3. Do not include the manifest or checksum file in its own checksum list. Finish
   writing the archive before hashing it; any subsequent repack needs a new hash.

For the rolling example, generate the changed-archive entries with this Python
standard-library snippet, run beside the rebuilt archive. Its output is a
**replacement block**, not a complete static checksum file:

```python
import hashlib
import tarfile
from pathlib import Path, PurePosixPath

archive = Path("resource-rolling.tar.gz")
with archive.open("rb") as stream:
    print(hashlib.file_digest(stream, "sha256").hexdigest(), archive.name, sep="  ")
with tarfile.open(archive, "r:gz") as bundle:
    seen = set()
    for member in bundle.getmembers():
        if member.isdir():
            continue
        path = PurePosixPath(member.name)
        if (not member.isfile() or path.is_absolute() or ".." in path.parts
                or str(path) in seen):
            raise ValueError(f"Invalid or duplicate archive member: {member.name}")
        seen.add(str(path))
        with bundle.extractfile(member) as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        print(digest, str(path), sep="  ")
```

Before upload, compare the manifest paths for the changed archive against its
member list in both directions: no missing members and no unlisted files.
Check for duplicate install paths across the rolling and static manifests.
Use `sha256sum --check` for archive digests and, from the staged payload root,
member digests. All hashes must match the exact bytes you will publish.

## 5. Upload and publish a consistent record

In the existing-record editor or new-version draft selected in step 1:

1. Retain/import every unchanged archive needed by the manifest.
2. Remove the previous upload of each changed archive and upload its rebuilt
   replacement under the expected filename. Removing a loose file on the record
   page does not remove a member from a tar archive.
3. Replace the manifest and matching checksum file together. Review the draft's
   complete file list, title/version, description, creators, license and public
   file visibility. Describe additions, replacements and removals.
4. Publish the completed change. Verify the public manifest, checksum file and
   archive download links after publication. A draft DOI is not yet an anonymous
   public download location.

Keep the previous version available for reproducibility. When removing a file
for legal/privacy reasons, deleting it from a new version does not remove old
public copies: contact the record owners and Zenodo support separately.

## 6. Connect the publication to Tudatpy

For a same-record update, the existing record URLs can remain unchanged. For a
new version, update both URLs for the changed record in
`RESOURCE_METADATA_URLS`: its manifest and checksum URL. The manifest's first
line must point to the same new record. Leave the unchanged record pinned.
Do not substitute a concept DOI into this tuple: the code expects downloadable
manifest and checksum files from specific records.

Submit that URL change through the normal Tudatpy review/release process.
Publishing a newer Zenodo version does not redirect existing installations of
Tudatpy that still pin an older record. Coordinate the resource publication with
any code that starts requiring a new file; make the public data available before
releasing that code. This guide documents the current installer, not an
automatic latest-version discovery mechanism.

## 7. Verify through the installer

Use a Tudatpy environment containing the intended `RESOURCE_METADATA_URLS`.
With a source build, select its Python package through the appropriate
`PYTHONPATH`; inspect `tudatpy.resource_catalog.__file__` if unsure. From an empty
working directory, refresh metadata into isolated paths:

```bash
python -m tudatpy.resource_installer --mode manifest --catalog ./check-catalog.csv
python -m tudatpy.resource_installer --mode list --catalog ./check-catalog.csv \
  --files earth_orientation/eopc01_polar_motion_1846-1961.txt --dest ./check-resources
python -m tudatpy.resource_installer --mode update --catalog ./check-catalog.csv \
  --files earth_orientation/eopc01_polar_motion_1846-1961.txt \
  --dest ./check-resources --cache-dir ./check-cache
```

Replace the example selection for other changes. Inspect the refreshed CSV
first and confirm that it contains the exact expected path: `--files` selects
by substring, and an unmatched selection alone is not proof of success.
Manifest mode also writes `resource_hashes.json` beside the custom catalog.
The installer uses those hashes automatically. Confirm the actual installed
file exists, its digest matches, and a small consuming application reads it.
An isolated `--dest` does not by itself redirect the C++ resource lookup;
configure the consuming application for that location or test separately in a
correctly configured resource installation.

For a removal, confirm the path is absent from the new catalog and archive.
Also test updating from an older installation: **the installer does not delete
obsolete local files**. Users who want to remove an old local copy must remove
that specific file themselves after checking their applications. No installer
mode named `remove` or `delete` exists.

A fresh cache demonstrates that the public archive is sufficient, rather than
reusing a local copy. Cached tarballs are named by archive basename. With the
current authoritative archive hashes, the installer checks an existing archive
and discards an outdated cached copy before downloading its replacement.
`--mode update` refreshes the selected archives explicitly. Preserve checksum
verification when troubleshooting a failed download or publication.

Resource maintainers should review the staged bytes, scientific provenance,
manifest grouping and record selection before publication. The command examples
and workflow in this guide were AI-assisted; no production Zenodo record was
modified while preparing this documentation.
