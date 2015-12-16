import sys
import os
import shutil

'''
Opens the original file given. Splits the file into specified number of files.

'''
def move_foldx_files(destination_path):
    '''
    :param destination_path: the path of where you want to move the files
    :return:moves the files from /src to the destination path
    '''
    essential_files_directory = os.getcwd() + "/src"
    essential_files = os.listdir(essential_files_directory)
    for file in essential_files:
        file_name = os.path.join(essential_files_directory, file)
        if os.path.isfile(file_name):
            shutil.copy(file_name, destination_path)

def make_file(files_made, file_name):
    '''
    :param files_made: number of files already made, used to number new directory and file
    :param file_name: name of file to be created
    :return: returns access to the file that was just created in a new directory
    '''
    folder_name = files_made + "_foldx_split"
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        move_foldx_files(folder_name)
    file = open(folder_name + "/" + file_name, 'w')
    return file

def split_file(file_name, num_split_files):
    '''
    :param file_name: name of file to be split up
    :param num_split_files: total number of files desired
    :return: specified number of directories with split up files and all the files included in /src
    '''
    split_file = open(file_name, 'r')
    # determine number of lines in the file
    num_lines = 0
    for line in split_file:
        num_lines += 1
    lines_in_each_file = num_lines / num_split_files
    split_file.close()
    split_file = open(file_name, 'r')

    lines_in_new_file = 0
    files_made = 1
    current_file = make_file(str(files_made), str(files_made) + "_" + file_name[2:])
    for line in split_file:
        if lines_in_new_file < lines_in_each_file:
            current_file.write(line)
            lines_in_new_file += 1
            if lines_in_new_file == lines_in_each_file:
                lines_in_new_file = 0
                if files_made != num_split_files:  # haven't made the desired number of files yet
                    current_file.close()
                    files_made += 1
                    current_file = make_file(str(files_made), str(files_made) + "_" + file_name[2:])
        else:  # on the last file, write every line to this file
            current_file.write(line)

def main():
    split_file(file_name, num_split_files)

if __name__ == "__main__":
    file_name = sys.argv[1]
    num_split_files = int(sys.argv[2])
    main()