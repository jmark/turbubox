#!/bin/sh

PLOTCMD="$(dirname "$0")/plot.py"

PROFILE="${1:-"density"}"
shift
PNGDIR="png/$PROFILE"

SPLIT="$1"
shift

mkdir -p $PNGDIR

if test -z "$SPLIT"
then
    parallel --eta \
        test -f "$PNGDIR/{/.}.png" '||' \
        $PLOTCMD '{}' "$PNGDIR/{/.}.png" \
        --title "file://title.txt" --profile $PROFILE \
        ::: $(find data -iname 'chkpt_*.h5' | sort)
else
    cat "$SPLIT" | parallel --halt now,fail=1 --eta \
        test -f "$PNGDIR/{/.}.png" '||' \
        $PLOTCMD '{}' "$PNGDIR/{/.}.png" \
        --title "file://title.txt" --profile $PROFILE
fi
