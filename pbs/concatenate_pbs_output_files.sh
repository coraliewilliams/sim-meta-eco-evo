#!/bin/bash

# Request input for array ID, name of the study, and job number range
read -p "Enter the array ID: " ARRAY_ID
read -p "Enter the name of the study: " STUDY_NAME
read -p "Enter the starting job number: " START_JOB_NUM
read -p "Enter the ending job number: " END_JOB_NUM

# Create output files for the concatenated error and observation files
ERROR_OUTPUT="${ARRAY_ID}_error_concatenated.txt"
OBS_OUTPUT="${ARRAY_ID}_obs_concatenated.txt"
MISSING_FILES="${ARRAY_ID}_missing_files.txt"

# Clear output files if they exist
> "$ERROR_OUTPUT"
> "$OBS_OUTPUT"
> "$MISSING_FILES"

# Loop through the specified range of job numbers
for job_num in $(seq $START_JOB_NUM $END_JOB_NUM); do
    ERROR_FILE="${STUDY_NAME}.e${ARRAY_ID}.${job_num}"
    OBS_FILE="${STUDY_NAME}.o${ARRAY_ID}.${job_num}"

    # Flag to check if any file is missing
    missing_flag=0

    # Check if the error file exists
    if [ -f "$ERROR_FILE" ]; then
        # Concatenate the job number and the content of the file into the output file
        echo -e "${job_num}\t$(cat "$ERROR_FILE")" >> "$ERROR_OUTPUT"
    else
        missing_flag=1
    fi

    # Check if the observation file exists
    if [ -f "$OBS_FILE" ]; then
        # Concatenate the job number and the content of the file into the output file
        echo -e "${job_num}\t$(cat "$OBS_FILE")" >> "$OBS_OUTPUT"
    else
        missing_flag=1
    fi

    # If either file is missing, record the job number
    if [ $missing_flag -eq 1 ]; then
        echo "${job_num}" >> "$MISSING_FILES"
    fi
done

echo "Concatenation completed. Missing files logged in $MISSING_FILES."

# Delete all error and observation files in the specified range
rm -f "${STUDY_NAME}.e${ARRAY_ID}."* "${STUDY_NAME}.o${ARRAY_ID}."*

echo "All error and observation files in the specified range have been deleted."
