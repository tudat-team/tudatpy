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
   Leave out the version heading; it is added automatically. If you include
   subheadings, underline them with ``~`` characters.
3. Start the release preparation with this terminal command, replacing the
   version and filename as needed. This uses the GitHub command-line tool, ``gh``::

       gh workflow run prepare_release.yml --repo tudat-team/tudatpy --ref develop \
         -f version=1.1.0 -F release_notes=@RELEASE_NOTES.rst

4. Open the **Prepare release** run in TudatPy's GitHub Actions page and follow
   its links to the eight pull requests (PRs). Review their changes. Each of
   the four repositories has one PR into ``master`` for the release and one
   into ``develop`` to set its version to ``X.Y.Z.dev0``.
5. Merge the PRs into ``master`` first, in this order: examples, TudatPy,
   Tudat-space, then feedstock. Choose **Create a merge commit** on GitHub.
   Before merging feedstock, wait for the corresponding TudatPy tag to appear,
   select **Ready for review** on its PR, and wait for the checks to pass.
   Rerun any check that failed because the source tag was not yet available.
   Then repeat this order for the PRs into ``develop``.
6. Check that the package is available and both documentation sites show the
   release correctly. Finish all eight PRs before starting another release.

For a single paragraph of notes, you can also start **Prepare release** from
the GitHub Actions page: select ``develop`` and enter the version and notes.
Use the terminal command above for notes with paragraphs or lists.

What happens automatically
--------------------------

The preparation script checks the version and prepares the changes for all eight
PRs before uploading anything. It stops if an earlier release is unfinished or if
changes on ``master`` need to be included in ``develop`` first.

For the release, it sets the version to ``X.Y.Z`` and updates the package
settings. For development, it keeps the development settings and sets the
version to ``X.Y.Z.dev0``. Both include matching examples and keep a copy of the
notes for published releases. There are no separate notes for weekly development
versions. Both websites show the same release announcement on master and develop,
linking to the release notes on the stable API docs site. The two feedstock PRs
start as drafts because they need the TudatPy tags before their packages can be built.

After each release PR is merged, a tag marks that exact version of the files:
``vX.Y.Z`` on ``master`` or ``vX.Y.Z.dev0`` on ``develop``. This applies to all
four repositories. Ordinary PR merges do not create tags. The existing package
and documentation services then build the release; allow time for them to finish.

The weekly process waits while TudatPy or feedstock release PRs are open, or
their development versions differ. Its next run that passes the usual checks
creates ``X.Y.Z.dev1`` in those two repositories, followed by ``.dev2``, and so on.

If preparation stops
--------------------

Read the error in the **Prepare release** run and check which PRs were created.
If none were merged, close any PRs from that attempt and delete only their
release branches before trying again. If any PR was already merged, finish
the remaining release steps manually before starting another release.
