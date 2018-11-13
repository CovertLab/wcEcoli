#! /bin/bash

# Purges Jenkins runs from a given month while keeping the first run of that month.
# Useful for freeing space on Sherlock while still keeping some history of runs.
# Uses a relative month to now so that it can automatically remove old months when
# scheduled to run in a Jenkins script.
#
# WARNING: will delete shared files in $PI_SCRATCH/wc_ecoli/
#
# Usage:
#   ./purge.sh <directory> <months ago>
#     - directory options include daily_build, with_aa, anaerobic
#     - months ago specifies the month to purge by how many months ago from now it is
#       (eg. 2 when the current month is Nov will remove all but the first day run
#       in Sept)

# Handle arguments
DIR_TO_PURGE="$1"
MONTHS_AGO_TO_PURGE="$2"

if [ $# -ne 2 ]; then
	echo "Need to supply directory and month to purge"
	exit
fi

if [ $DIR_TO_PURGE != "daily_build" ] && [ $DIR_TO_PURGE != "with_aa" ] && [ $DIR_TO_PURGE != "anaerobic" ]; then
	echo "Incorrect directory passed"
	exit
fi

re='^[0-9]+$'
if ! [[ $MONTHS_AGO_TO_PURGE =~ $re ]] || [ $MONTHS_AGO_TO_PURGE -lt 6 ]; then
	echo "Need to supply a month greater than 5"
	exit
fi

DIR="/scratch/PI/mcovert/wc_ecoli/$DIR_TO_PURGE"

# Only deletes files from one month which is selected by MONTHS_AGO_TO_PURGE months prior to today
MONTH=`expr $(date +%m) - $MONTHS_AGO_TO_PURGE`
YEAR=$(date +%Y)
while [ $MONTH -lt 1 ]; do
	MONTH=`expr $MONTH + 12`
	YEAR=`expr $YEAR - 1`
done
DAY=2

# Find the first date of the month that has a file (it won't be removed)
FILES=$(find $DIR -maxdepth 1 -type d -newermt $YEAR-$MONTH-1 ! -newermt $YEAR-$MONTH-$DAY | wc -l)
while [ $FILES -lt 1 ] && [ $DAY -lt 31 ]; do
	DAY=`expr $DAY + 1`
	FILES=$(find $DIR -maxdepth 1 -type d -newermt $YEAR-$MONTH-1 ! -newermt $YEAR-$MONTH-$DAY | wc -l)
done

# Set next month so that only one month gets removed
if [ `expr $MONTH + 1` -gt 12 ]; then
	NEXT_MONTH=1
	NEXT_YEAR=`expr $YEAR + 1`
else
	NEXT_MONTH=`expr $MONTH + 1`
	NEXT_YEAR=$YEAR
fi

# Only remove if less than 30 directories to remove (ie 1 month) as a safety check (might not be necessary)
if [ $(find $DIR -maxdepth 1 -type d -newermt $YEAR-$MONTH-$DAY ! -newermt $NEXT_YEAR-$NEXT_MONTH-1 | wc -l) -lt 30 ]; then
	echo "Purging directories from $MONTH/$YEAR to $NEXT_MONTH/$NEXT_YEAR created on or after day $DAY for $DIR"
	find $DIR -maxdepth 1 -type d -newermt $YEAR-$MONTH-$DAY ! -newermt $NEXT_YEAR-$NEXT_MONTH-1 -printf "%p\n" -exec rm -fr {} \;
else
	echo "More than 30 directories to remove in $DIR between $MONTH/$YEAR and $NEXT_MONTH/$NEXT_YEAR"
fi
