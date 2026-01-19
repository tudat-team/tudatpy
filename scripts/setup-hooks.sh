#!/bin/bash
# Setup script to configure git hooks for this repository
#
# Run this once after cloning the repository:
#   ./scripts/setup-hooks.sh

set -e

REPO_ROOT="$(git rev-parse --show-toplevel)"

echo "Configuring git to use project hooks..."
git config core.hooksPath scripts/hooks

echo "Git hooks configured successfully!"
echo "The pre-commit hook will now automatically sync docs/ on each commit."
