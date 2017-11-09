#!/bin/sh

PLOTCMD="$(dirname "$0")/plot.py"

PROFILE="${1:-"density"}"
shift
PNGDIR="png/$PROFILE"

ARGFILE="$1"
shift

mkdir -p $PNGDIR

parallel --halt now,fail=1 --eta --arg-file "$ARGFILE" \
    test -f "$PNGDIR/{/.}.png" '||' \
    $PLOTCMD --srcfile '{}' --snkfile "$PNGDIR/{/.}.png" \
    --title "file://title.txt" --profile $PROFILE --nvisu 16
