#! /bin/bash

# Purges Jenkins runs, keeping only the most recent N days of runs.
# A "day" is determined by the date prefix in the directory name (YYYYMMDD).
# All directories from a kept day are preserved (e.g., all 9 optional_features
# variants from that day). Everything older is deleted.
#
# WARNING: will delete shared files in /scratch/groups/mcovert/wc_ecoli/
#
# Usage:
#   ./purge.sh <directory> [keep_days]
#     - directory options include daily_build, with_aa, anaerobic, optional_features
#     - keep_days: number of most recent days to keep (default: 14)

set -eu

# Handle arguments
DIR_TO_PURGE="${1:-}"
KEEP_DAYS="${2:-14}"

if [ -z "$DIR_TO_PURGE" ]; then
	echo "Usage: purge.sh <directory> [keep_days]"
	exit 1
fi

if [ "$DIR_TO_PURGE" != "daily_build" ] && [ "$DIR_TO_PURGE" != "with_aa" ] && [ "$DIR_TO_PURGE" != "anaerobic" ] && [ "$DIR_TO_PURGE" != "optional_features" ]; then
	echo "Incorrect directory passed (got: $DIR_TO_PURGE)"
	exit 1
fi

re='^[0-9]+$'
if ! [[ "$KEEP_DAYS" =~ $re ]] || [ "$KEEP_DAYS" -lt 1 ]; then
	echo "keep_days must be a positive integer (got: $KEEP_DAYS)"
	exit 1
fi

DIR="/scratch/groups/mcovert/wc_ecoli/$DIR_TO_PURGE"

if [ ! -d "$DIR" ]; then
	echo "Directory does not exist: $DIR"
	exit 1
fi

# List all run directories sorted alphabetically (which is chronological
# since directory names start with a timestamp like 20260318.003210).
ALL_DIRS=$(find "$DIR" -maxdepth 1 -mindepth 1 -type d | sort)

if [ -z "$ALL_DIRS" ]; then
	echo "No run directories found in $DIR"
	exit 0
fi

# Extract unique date prefixes (YYYYMMDD) from directory names, sorted
ALL_DATES=$(echo "$ALL_DIRS" | xargs -n1 basename | sed 's/\..*//' | sort -u)
N_DATES=$(echo "$ALL_DATES" | wc -l)

if [ "$N_DATES" -le "$KEEP_DAYS" ]; then
	echo "Nothing to purge in $DIR ($N_DATES days of runs, keeping all)"
	exit 0
fi

# Find the cutoff: dates to remove are everything except the last KEEP_DAYS
DATES_TO_REMOVE=$(echo "$ALL_DATES" | head -n -"$KEEP_DAYS")

# Build list of directories to remove (those whose date prefix is in the remove set)
TO_REMOVE=""
while IFS= read -r dir; do
	date_prefix=$(basename "$dir" | sed 's/\..*//')
	if echo "$DATES_TO_REMOVE" | grep -qx "$date_prefix"; then
		TO_REMOVE="$TO_REMOVE
$dir"
	fi
done <<< "$ALL_DIRS"

TO_REMOVE=$(echo "$TO_REMOVE" | grep -v '^$' || true)

if [ -z "$TO_REMOVE" ]; then
	echo "Nothing to purge in $DIR"
	exit 0
fi

N_REMOVE=$(echo "$TO_REMOVE" | wc -l)
N_TOTAL=$(echo "$ALL_DIRS" | wc -l)
N_DAYS_REMOVE=$(echo "$DATES_TO_REMOVE" | wc -l)

echo "Purging $N_REMOVE directories ($N_DAYS_REMOVE days) of $N_TOTAL total in $DIR (keeping $KEEP_DAYS most recent days)"

while IFS= read -r dir; do
	[ -z "$dir" ] && continue
	echo "  Removing: $(basename "$dir")"
	rm -rf "$dir"
done <<< "$TO_REMOVE"

echo "Done. $N_REMOVE directories removed."
