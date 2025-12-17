for d in s1-*; do
    # Extract numbers from directory name s1-XX_s2-YY
    s1=$(echo "$d" | sed -E 's/s1-([0-9]+)_s2-([0-9]+)/\1/')
    s2=$(echo "$d" | sed -E 's/s1-([0-9]+)_s2-([0-9]+)/\2/')

    # File path
    file="$d/input.txt"

    # Replace S1= and S2= lines in the file
    sed -i "s/^S1=.*/S1= $s1/" "$file"
    sed -i "s/^S2=.*/S2= $s2/" "$file"
done
