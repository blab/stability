import sys
import os
import math

def split_file(file_name, num_split_files):
    #file_name = "test_split_file.txt"
    split_file = open(file_name, 'r')

    # determine number of lines in the file
    num_lines = 0
    for line in split_file:
        num_lines += 1
    #lines_in_each_file = math.ceil(float(num_lines) / float(num_split_files))
    lines_in_each_file = num_lines / num_split_files

    split_file.close()
    split_file = open(file_name, 'r')


    current_line_number = 0
    lines_in_new_file = 0
    root_line = split_file.readline()  # get the root line that needs to be put in each new file
    files_made = 1
    current_file = open("test_split_file_" + str(files_made) + ".txt", 'w')
    current_file.write(root_line)
    for line in split_file:
        if files_made != num_split_files:
            if lines_in_new_file < lines_in_each_file:
                current_file.write(line)
                lines_in_new_file += 1
                if lines_in_new_file == lines_in_each_file:
                    current_file.close()
                    lines_in_new_file = 0
                    if files_made != num_split_files:  # haven't made the desired number of files yet
                        files_made += 1
                        current_file = open("test_split_file_" + str(files_made) + ".txt", 'w')
                        current_file.write(root_line)

        else:  # on the last file, write every line to this file
            current_file.write(line)

def main():
    split_file(file_name, num_split_files)

if __name__ == "__main__":
    file_name = sys.argv[1]
    num_split_files = int(sys.argv[2])
    main()