Code with implementations of:
"A comparative study of Strength Reduction  and Gravity Increase Methods considering Random Fields in slope stability analysis"

# dcproj

#download/install kdevelop

sudo apt-get install kdevelop

#download/install cmake

sudo apt-get install cmake

#download/install git

sudo apt-get install git

#download/install eigen

sudo apt-get install libeigen3-dev

#clone repository

git clone https://github.com/diogocecilio/dcproj.git

#in kdevelop open/import the cloned project
#Choose cmake configure and create a build directory in a parallel directory




MEMORY LEAK CHECK:

valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
         --log-file=valgrind-out.txt \
         ./dcproj exampleParam1
