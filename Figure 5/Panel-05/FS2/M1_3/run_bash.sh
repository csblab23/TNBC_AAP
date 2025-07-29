#!/bin/bash

# Set your notebook filename (without extension)
NOTEBOOK_NAME="your_notebook"

# Convert the Jupyter notebook to a Python script
echo "Converting $NOTEBOOK_NAME.ipynb to Python script..."
jupyter nbconvert --to script "${NOTEBOOK_NAME}.ipynb"

# Run the script 100 times and log output
echo "Running ${NOTEBOOK_NAME}.py 100 times..."
for i in {1..100}; do
    echo "Run $i..." >> your_log.txt
    python "${NOTEBOOK_NAME}.py" "$i" >> your_log.txt 2>&1
done

echo "All runs complete. Output saved to your_log.txt."
