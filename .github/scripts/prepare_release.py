"""Prepare a release of TudatPy, its documentation, package recipe and examples.

Each repository gets two pull requests (PRs): the release into master and an
updated development version into develop. By default, preparation stays local;
--publish uploads the prepared branches and opens the PRs on GitHub.
"""

import argparse
import configparser
import os
from pathlib import Path
import re
import subprocess

OWNER = "tudat-team"
# Prepare examples first so the source and documentation can use the matching examples.
REPOSITORIES = ("tudatpy-examples", "tudatpy", "tudat-space", "tudatpy-feedstock")
SUBMODULES = {"tudatpy": "examples/tudatpy", "tudat-space": "docs/source/examples/tudatpy-examples"}
NOTES = "docs/tudatpy/source/release_notes.rst"


def run(*args, cwd=None):
    """Run a command in the chosen folder and return its text output.

    A failed command stops preparation so later steps cannot use incomplete results.
    """
    return subprocess.check_output([str(arg) for arg in args], cwd=cwd, text=True).strip()


def git(repo, *args):
    """Run a Git command in the given repository."""
    return run("git", *args, cwd=repo)


def config_version(text):
    """Read the current version from the contents of .bumpversion.cfg."""
    config = configparser.RawConfigParser()
    config.read_string(text)
    return config["bumpversion"]["current_version"]


def replace(path, pattern, replacement):
    """Replace one expected setting, preserving the file's line endings.

    The pattern must match exactly once. A changed file layout needs review.
    """
    text, count = re.subn(pattern, replacement, path.read_bytes().decode(), flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"{path}: expected one match for {pattern!r}, found {count}")
    path.write_bytes(text.encode())


def set_version(repo, name, version):
    """Update the version in the configured files and reset the package build number.

    This edits files only. The complete release changes are saved together later.
    """
    if name == "tudatpy-examples":
        (repo / "version").write_text(version + "\n")
    else:
        # Allow existing edits: the develop step has just put back its original settings.
        run(
            "bump-my-version",
            "bump",
            "--config-file",
            ".bumpversion.cfg",
            "--new-version",
            version,
            "--no-commit",
            "--no-tag",
            "--allow-dirty",
            cwd=repo,
        )
        if config_version((repo / ".bumpversion.cfg").read_text()) != version:
            raise ValueError(f"{name}: bump-my-version did not write {version}")
    if name == "tudatpy-feedstock":
        # A new version starts at package build number zero.
        replace(repo / "recipe/meta.yaml", r'^([ \t]*\{% set build = )"[^"\r\n]*"', r'\g<1>"0"')


def make_stable(repo, name):
    """Select released packages for the docs/examples and prepare package publishing."""
    environments = {
        "tudatpy": "docs/tudatpy/environment_readthedocs.yaml",
        "tudatpy-examples": "environment.yaml",
    }
    if name in environments:
        replace(repo / environments[name], r"^[ \t]*-[ \t]+tudat-team/label/dev[ \t]*\r?\n", "")
    if name == "tudatpy-feedstock":
        config = repo / "recipe/conda_build_config.yaml"
        replace(config, r"^([ \t]*-[ \t]+tudat-team )dev([ \t]*\r?)$", r"\g<1>main\2")
        replace(config, r"^([ \t]*-[ \t]+)tudat-team/label/dev, ", r"\1")
        # Only builds of master may upload the released package. PR builds can run checks.
        replace(
            repo / "conda-forge.yml",
            r"^upload_on_branch: develop([ \t]*\r?)$",
            r"upload_on_branch: master\1",
        )
        # Regenerate the package build files from the settings above.
        run("conda", "smithy", "rerender", cwd=repo)


def add_notes(repo, version, notes):
    """Add this version's notes above the older releases, keeping their notes intact."""
    path = repo / NOTES
    heading = "Release notes\n=============\n\n"
    existing = path.read_text()
    title = f"Version {version}"
    # Keep each release listed once, under the existing document title.
    if not existing.startswith(heading) or f"\n{title}\n" in existing:
        raise ValueError(f"Review the heading or existing {version} entry in {NOTES}")
    path.write_text(
        heading
        + title
        + "\n"
        + "-" * len(title)
        + "\n\n"
        + notes
        + "\n\n"
        + existing[len(heading) :]
    )


def check_merge(repo, branch, base):
    """Check which files would result from merging the prepared PR.

    Git calculates this without changing either branch. Stop if the changes
    conflict or the target branch has fixes missing from the prepared release.
    """
    try:
        merged = git(repo, "merge-tree", "--write-tree", f"origin/{base}", branch)
    except subprocess.CalledProcessError as error:
        raise ValueError(
            f"{repo.name}: reconcile {base} into develop before releasing:\n{error.output}"
        ) from error
    if merged != git(repo, "rev-parse", f"{branch}^{{tree}}"):
        raise ValueError(
            f"{repo.name}: {base} has changes missing from {branch}; reconcile them first"
        )


def prepare(directory, version, notes, publish=False):
    """Prepare a release PR and a development-version PR in each repository.

    Prepare and check all eight branches locally first. If publish is True,
    upload them and open the PRs for the developer to review and merge.
    """
    if not re.fullmatch(r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)", version):
        raise ValueError("Version must be X.Y.Z, for example 1.1.0")
    # Accept notes copied from Windows or a text editor, with consistent line endings.
    notes = "\n".join(line.rstrip() for line in notes.splitlines()).strip()
    if not notes:
        raise ValueError("Release notes are required")
    # Work in new copies, leaving the developer's existing checkouts alone.
    directory.mkdir(parents=True, exist_ok=False)
    prepared, examples = [], {}
    for name in REPOSITORIES:
        repo = directory / name
        run(
            "git",
            "clone",
            "--filter=blob:none",
            "--no-checkout",
            f"https://github.com/{OWNER}/{name}.git",
            repo,
        )
        git(repo, "config", "user.name", "GitHub Actions")
        git(repo, "config", "user.email", "actions@github.com")
        # Finish earlier release attempts before preparing another one.
        if (
            run(
                "gh",
                "pr",
                "list",
                "--repo",
                f"{OWNER}/{name}",
                "--state",
                "open",
                "--limit",
                "500",
                "--json",
                "headRefName",
                "--jq",
                '[.[] | select(.headRefName | startswith("release/"))] | length',
            )
            != "0"
        ):
            raise ValueError(f"{name}: finish the open release PRs first")
        # Existing branches or tags may belong to a previous attempt at this version.
        for ref in (
            f"refs/tags/v{version}",
            f"refs/tags/v{version}.dev0",
            f"refs/heads/release/v{version}",
            f"refs/heads/release/develop-v{version}",
        ):
            if git(repo, "ls-remote", "origin", ref):
                raise ValueError(f"{name}: {ref} already exists; inspect the previous attempt")
        if name == "tudatpy-examples":
            try:
                git(repo, "show", "origin/develop:version")
            except subprocess.CalledProcessError as error:
                raise ValueError(
                    "Update examples develop from master and install its release setup first"
                ) from error
        else:
            previous = config_version(git(repo, "show", "origin/master:.bumpversion.cfg"))
            # Compare the three numbers, so 1.10.0 follows 1.9.0.
            if tuple(map(int, version.split("."))) <= tuple(
                map(int, previous.split(".dev")[0].split("."))
            ):
                raise ValueError(f"{name}: release must be newer than master ({previous})")

        # The develop update starts from the release commit, so Git remembers
        # that these release changes have already been included when merging later.
        for base, new_version, branch in (
            ("master", version, f"release/v{version}"),
            ("develop", f"{version}.dev0", f"release/develop-v{version}"),
        ):
            git(
                repo,
                "checkout",
                "-b",
                branch,
                "origin/develop" if base == "master" else f"release/v{version}",
            )
            if base == "develop":
                # Put back develop's package settings and build files before updating its
                # version. Also remove files that were created only for the stable release.
                git(repo, "restore", "--source=origin/develop", "--staged", "--worktree", "--", ".")
            set_version(repo, name, new_version)
            if base == "master":
                make_stable(repo, name)
            if name == "tudatpy":
                add_notes(repo, version, notes)
            git(repo, "add", "--all")
            if name in SUBMODULES:
                # Use the matching examples changes prepared above. The number 160000
                # tells Git that this path points to another repository.
                git(
                    repo,
                    "update-index",
                    "--add",
                    "--cacheinfo",
                    f"160000,{examples[base]},{SUBMODULES[name]}",
                )
            git(repo, "commit", "-m", f"Prepare TudatPy {new_version}")
            check_merge(repo, branch, base)
            if name == "tudatpy-examples":
                examples[base] = git(repo, "rev-parse", "HEAD")
            prepared.append((repo, name, branch, base, new_version))

    # Upload only after all eight branches have passed the preparation checks.
    # The developer reviews and merges the resulting PRs.
    for repo, name, branch, base, new_version in prepared:
        if not publish:
            print(f"Prepared locally: {name} {branch} -> {base}")
            continue
        # Tags are created after merging, when the final release files are known.
        git(
            repo,
            "-c",
            "push.followTags=false",
            "push",
            "--no-follow-tags",
            "origin",
            f"{branch}:refs/heads/{branch}",
        )
        body = "Use a regular merge commit. Finish all eight release PRs before the next release."
        if name == "tudatpy-feedstock":
            # The package recipe needs the new source tag, created after the TudatPy merge.
            body += (
                f" Merge the matching TudatPy PR and wait for tag v{new_version}. "
                "Then mark this PR ready and rerun any early failed Azure check."
            )
        url = run(
            "gh",
            "pr",
            "create",
            "--repo",
            f"{OWNER}/{name}",
            "--base",
            base,
            "--head",
            branch,
            "--title",
            f"Prepare TudatPy {new_version}",
            "--body",
            body,
            *(["--draft"] if name == "tudatpy-feedstock" else []),
        )
        print(url, flush=True)
        if os.environ.get("GITHUB_STEP_SUMMARY"):
            with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as stream:
                stream.write(f"- [{name} -> {base}]({url})\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version")
    parser.add_argument("--notes", type=Path, required=True)
    parser.add_argument("--directory", type=Path, default=Path("release-work"))
    parser.add_argument(
        "--publish",
        action="store_true",
        help="Push branches and create PRs (default: local preparation only)",
    )
    args = parser.parse_args()
    prepare(args.directory.resolve(), args.version, args.notes.read_text(), args.publish)
