Making a release
================

A release includes the changes on ``develop`` in TudatPy, Tudat-space,
tudatpy-feedstock (the package recipe), and tudatpy-examples. The same procedure
applies to major, minor and patch releases.

What you do
-----------

1. Choose a version number in the form ``X.Y.Z``, such as ``1.1.0``.
   It must be newer than the version on ``master``.
2. Write the release notes in a local file, for example ``RELEASE_NOTES.rst``.
   Leave out the version heading (for example, ``Version 1.1.0``); the script
   adds it to the published release notes. If you include subheadings,
   underline them with ``~`` characters.
3. Start the release preparation with this terminal command, replacing the
   version and filename as needed. This uses the GitHub command-line tool, ``gh``::

       gh workflow run prepare_release.yml --repo tudat-team/tudatpy --ref develop \
         -f version=1.1.0 -F release_notes=@RELEASE_NOTES.rst

4. Open the **Prepare release** run in TudatPy's GitHub Actions page and follow
   its links to the eight pull requests (PRs). Review their changes. Each of
   the four repositories has one PR into ``master`` for the release and one
   into ``develop`` to set its version to ``X.Y.Z.dev0``.
5. Merge the PRs into ``master`` first, in this order: examples, TudatPy,
   Tudat-space, then feedstock. Click **Merge pull request**, then **Confirm merge**.
   If the button shows another merge method, use its dropdown to select
   **Create a merge commit** first.
   Before merging feedstock, wait for the corresponding tag to appear on the
   `TudatPy tags page <https://github.com/tudat-team/tudatpy/tags>`_
   (``vX.Y.Z`` for ``master``, ``vX.Y.Z.dev0`` for ``develop``).
   Select **Ready for review** on the feedstock PR and wait for the checks to pass.
   Then repeat this order for the PRs into ``develop``.
6. Check that the package is available and both documentation sites show the
   release correctly. Finish all eight PRs before starting another release.

What happens automatically
--------------------------

The script prepares and checks all release changes locally. Once these checks
pass, it sends the changes to new branches on GitHub and opens the eight PRs.
It stops if an earlier release is unfinished or if changes on ``master`` need
to be included in ``develop`` first.

For the release, it sets the version to ``X.Y.Z`` and updates the package
settings. For development, it keeps the development settings and sets the
version to ``X.Y.Z.dev0``. Both websites show the same release announcement on master and develop,
linking to the release notes on the stable API docs site. The two feedstock PRs
start as drafts because they need the TudatPy tags before their packages can be built.

After each release PR is merged, a tag marks that exact version of the files:
``vX.Y.Z`` on ``master`` or ``vX.Y.Z.dev0`` on ``develop``. This applies to all
four repositories.

The weekly process waits while TudatPy or feedstock release PRs are open, or
their development versions differ. Its next run that passes the usual checks
creates ``X.Y.Z.dev1`` in those two repositories, followed by ``.dev2``, and so on.

If preparation stops
--------------------

Read the error in the **Prepare release** run and check which PRs were created.
If none were merged, close any PRs from that attempt and delete only their
release branches before trying again. If any PR was already merged, finish
the remaining release steps manually before starting another release.
