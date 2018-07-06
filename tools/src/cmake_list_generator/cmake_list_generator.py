import os
import sys

# Used to generate CMakeLists.txt files
# for IRL subdirectories.
# To use, cd into the directory you
# wish the CMakeList.txt to be generated
# for and run `python cmake_gen.py .`

target = sys.argv[1]

fn = sys.argv[2]
if os.path.exists(fn):
  print(os.path.basename(fn))
else:
  print("ERROR, no file path \n")

cwd_directory = str(os.path.abspath(fn))
check_src = cwd_directory.find("src")
if(check_src == -1):
  print("ERROR, src/ directory does not exist in given path \n")

up_to_source_directory = cwd_directory[:check_src+3]

file_to_write = os.path.join(cwd_directory+"/", "CMakeLists.txt")
with open(file_to_write, "w") as a:
    a.write("#List of files from this directory and its subdirectores."+os.linesep)
    for path, subdirs, files in os.walk(cwd_directory):
       for filename in files:
         if(filename == "CMakeLists.txt"):
           continue
         f = str(os.path.join(path, filename))
         relative_f = f[len(up_to_source_directory):]
         a.write("target_sources("+target+" PRIVATE ${IRL_SOURCE_DIR}"+relative_f+")" + os.linesep)
