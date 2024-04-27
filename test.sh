#!/bin/sh

rm -rf logFile
touch logFile
timestamp=$(date)
echo "Running tests at $timestamp" >> logFile

green='\033[0;32m'
red='\033[0;31m'
end='\033[0m'

# Run the application
for input in ./tests/input/*.txt
do
    name=$(basename $input)
    echo "--------------------------------" >> logFile
    echo "$name" >> logFile
    mpirun -np 8 ./build/dynamic_scc < $input >> logFile
    #/home/mage/Documents/random-graph-generator/tests/sol < $input >output.txt
    #/home/mage/Documents/random-graph-generator/tests/scc < $input >output.txt
    diff ./tests/output/$name ./output.txt > /dev/null
    if [ $? -eq 0 ]
    then
        # green color
        echo -e "${green} Test $name: passed"
        echo "Test $name passed" >> logFile
    else
        # red color
        echo -e "${red} Test $name: failed"
        echo "Test $name failed" >> logFile
    fi
    echo -e "${end}"
    echo "--------------------------------" >> logFile
done
