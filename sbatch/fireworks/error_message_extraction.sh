#!/bin/bash

# Define the JSON data directly in the script
JSON_DATA='[
    {
        "fw_id": 290,
        "created_on": "2026-03-03T00:43:10.155482",
        "updated_on": "2026-03-03T06:00:05.170683",
        "state": "FIZZLED",
        "name": "SimulationTask__Var_00__Seed_39__Gen_10__Cell_0"
    },
    {
        "fw_id": 389,
        "created_on": "2026-03-03T00:43:10.151657",
        "updated_on": "2026-03-03T05:11:50.488826",
        "state": "FIZZLED",
        "name": "SimulationTask__Var_00__Seed_31__Gen_7__Cell_0"
    }
]'

# Output CSV file
OUTPUT_CSV="simulation_errors.csv"

# Extract the "name" fields using jq
TASK_NAMES=$(echo "$JSON_DATA" | jq -r '.[].name')

# Write CSV header
echo "variant,seed,generation,error_snippet" > "$OUTPUT_CSV"

# Process each simulation task name
for NAME in $TASK_NAMES; do
    # Find the most recent error file
    ERROR_FILE=$(find .. -name "${NAME}*.error" -type f -exec stat -c '%Y %n' {} \; | sort -nr | awk 'NR==1 {print $2}')

    # Skip if no error file is found
    if [ -z "$ERROR_FILE" ]; then
        continue
    fi

    # Extract variant, seed, and generation from the filename
    if [[ $NAME =~ SimulationTask__Var_([0-9]+)__Seed_([0-9]+)__Gen_([0-9]+)__Cell_[0-9]+ ]]; then
        VARIANT="${BASH_REMATCH[1]}"
        SEED="${BASH_REMATCH[2]}"
        GENERATION="${BASH_REMATCH[3]}"
    else
        continue
    fi

    # Extract the last 250 characters of the error file
    ERROR_SNIPPET=$(tail -c 250 "$ERROR_FILE" | tr '\n' ' ' | tr ',' ' ') # Ensure CSV safety

    # Append to CSV file
    echo "$VARIANT,$SEED,$GENERATION,$ERROR_SNIPPET" >> "$OUTPUT_CSV"
done

echo "Processing complete. Results saved in $OUTPUT_CSV."

