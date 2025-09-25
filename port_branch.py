import subprocess
from pathlib import Path
import argparse

if __name__ == "__main__":

    # Get name of tudat branch
    parser = argparse.ArgumentParser()
    parser.add_argument("target_branch", help="Name of your tudat branch")
    parser.add_argument(
        "--reference-branch",
        default="tudat/develop",
        metavar="<branch>",
        help="The original branch in which you wanted to merge",
    )
    args = parser.parse_args()

    # Get files that changed between target branch and develop
    out = subprocess.run(
        [
            "git",
            "diff",
            args.reference_branch,
            args.target_branch,
            "--name-only",
        ],
        capture_output=True,
    )
    out.check_returncode()
    files = [Path(item) for item in out.stdout.decode("utf-8").splitlines()]

    # Filter files
    files_automatic = {}
    files_other = []
    for file in files:

        # try:
        #     name = str(file.parents[-2])
        # except:
        #     files_other.append(file)
        #     print(file)
        #     continue

        match str(file).split("/")[0]:
            case "include":
                files_automatic[file] = file
            case "src":
                updated_path = Path("src/tudat") / file.relative_to("src")
                files_automatic[file] = updated_path
            case "tests":
                updated_path = Path("tests/test_tudat") / file.relative_to(
                    "tests"
                )
                files_automatic[file] = updated_path
            case _:
                files_other.append(file)

    # Show files that will have to be taken care of manually
    if len(files_other) != 0:
        print("The following files will have to be adjusted manually")
        for file in files_other:
            print(file)
        print("")

    # Check-out all the automatic files
    for file in files_automatic:

        out = subprocess.run(
            ["git", "checkout", args.target_branch, "--", str(file)],
        )
        out.check_returncode()
        out = subprocess.run(["mv", file, files_automatic[file]])

    # Unstage
    out = subprocess.run(["git", "restore", "--staged", "."])
    out.check_returncode()
