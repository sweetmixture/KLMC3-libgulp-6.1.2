#!/bin/bash

#!/bin/bash

# Define the input file
input_file=$1

# Check if the file exists
if [[ ! -f "$input_file" ]]; then
	echo "File not found!"
	exit 1
fi

# Read the file into an array
readarray -t lines < "$input_file"

# Remove whitespace from each element in the array
for i in "${!lines[@]}"; do
    lines[$i]=$(echo "${lines[$i]}" | tr -d '[:space:]')
done

# Print each element in the array
for word in "${lines[@]}"; do

	echo "  if (associated(${word})) then"
	echo "    deallocate(${word})"
	echo "    nullify(${word})"
	echo "  end if"
done


#count=$1
#
#for (( i=1; i<=${count}; i++ )); do
#
#	echo "	if (associated(p)) then"
#	 echo "    deallocate(p)"
#	 echo "    nullify(p)"
#	 echo "  end if"
#
#done
