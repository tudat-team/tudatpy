#!/bin/sh

git filter-branch -f --env-filter '

HIDDEN_EMAIL="hidden@hidden.com"

if [[ "GIT_AUTHOR_EMAIL" = "$HIDDEN_EMAIL" || "$GIT_COMMITTER_EMAIL" = "$HIDDEN_EMAIL" ]];
then
    if [[ "$GIT_COMMITTER_NAME" = "Guido Holtkamp" || "$GIT_AUTHOR_NAME" = "Guido Holtkamp" ]];
    then
      export GIT_AUTHOR_EMAIL="guido.email@unknown.com"
      export GIT_COMMITTER_EMAIL="guido.email@unknown.com"
    fi
    if [[ "$GIT_COMMITTER_NAME" = "Liban Abdulkadir" || "$GIT_AUTHOR_NAME" = "Liban Abdulkadir" ]];
    then
      export GIT_AUTHOR_EMAIL="liban.email@unknown.com"
      export GIT_COMMITTER_EMAIL="liban.email@unknown.com"
    fi
    if [[ "$GIT_COMMITTER_NAME" = "Dominic Dirkx" || "$GIT_AUTHOR_NAME" = "Dominic Dirkx" ]];
    then
      export GIT_AUTHOR_EMAIL="D.Dirkx@tudelft.nl"
      export GIT_COMMITTER_EMAIL="D.Dirkx@tudelft.nl"
    fi
fi
' --tag-name-filter cat -- --branches --tags
