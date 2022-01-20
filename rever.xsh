$PROJECT = 'Tudat'
$ACTIVITIES = [
              'authors', # authors must have commit before executing Rever, otherwise he/she will be omitted.
              'version_bump',  # Changes the version number in various source files (setup.py, __init__.py, etc)
              'changelog',  # Uses files in the news folder to create a changelog for release
              'tag',  # Creates a tag for the new version number
              'push_tag',  # Pushes the tag up to the $TAG_REMOTE
              'forge', # updates feedstock
               # 'pypi',  # Sends the package to pypi
               # 'conda_forge',  # Creates a PR into your package's feedstock
               ]
$VERSION_BUMP_PATTERNS = [  # These note where/how to find the version numbers
                         ('version', '.*,', "'$VERSION'")
                         ]

# FORGE
$FORGE_SOURCE_URL = 'https://github.com/$GITHUB_ORG/$GITHUB_REPO/archive/$VERSION.tar.gz'
$FORGE_FEEDSTOCK = 'https://github.com/tudat-team/tudat-feedstock.git'
$FORGE_FEEDSTOCK_ORG = 'tudat-team'
$FORGE_PATTERNS = [
                ('meta.yaml','version\s=\s"(.*)"', "version='$VERSION'")
]

$CHANGELOG_FILENAME = 'CHANGELOG.rst'  # Filename for the changelog
$CHANGELOG_TEMPLATE = 'TEMPLATE.rst'  # Filename for the news template
$CHANGELOG_AUTHORS_TITLE = 'Authors'
$CHANGELOG_AUTHORS_FORMAT = '* {name}\n'
$PUSH_TAG_REMOTE = 'https://github.com/tudat-team/tudat.git'  # Repo to push tags to

$GITHUB_ORG = 'tudat-team'  # Github org for Github releases and conda-forge
$GITHUB_REPO = 'tudat'  # Github repo for Github releases  and conda-forge

