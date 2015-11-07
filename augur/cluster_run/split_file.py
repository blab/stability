import sys
import os
import shutil

'''
Opens the original file given. Splits the file into specified number of files. Keeps the 
first line in the original file in all subsequent split files. 

'''
# moves all files needed to run foldx into directory
def move_foldx_files(destination_path):
    essential_files_directory = os.getcwd() + "/essential_foldx_files"
    essential_files = os.listdir(essential_files_directory)
    for file in essential_files:
        file_name = os.path.join(essential_files_directory, file)
        if os.path.isfile(file_name):
            shutil.copy(file_name, destination_path)

# makes split folder and includes all neccesary files in new folder
def make_file(path_total_folder, file_name):
    folder_name = file_name[:len(file_name)-19] + "_foldx_split"
    current_path = path_total_folder
    check_path = current_path + "/" + folder_name
    if not os.path.exists(check_path):
        os.mkdir(check_path)
        move_foldx_files(check_path)
    file = open(path_total_folder + "/" + folder_name + "/" + file_name, 'w')
    return file

def split_file(file_name, num_split_files):
    folder_name = "/multiple_runs_0"
    total_folder = os.getcwd() + folder_name
    if os.path.exists(total_folder):
        while os.path.exists(total_folder):
            folder_number = int(folder_name[len(folder_name) - 1]) + 1
            folder_name = folder_name[:len(folder_name) - 1] + str(folder_number)
            total_folder = total_folder[:len(total_folder) - 1] + str(folder_number)
    os.mkdir(total_folder)



    split_file = open(file_name, 'r')

    # determine number of lines in the file
    num_lines = 0
    for line in split_file:
        num_lines += 1
    #lines_in_each_file = math.ceil(float(num_lines) / float(num_split_files))
    lines_in_each_file = num_lines / num_split_files

    split_file.close()
    split_file = open(file_name, 'r')


    lines_in_new_file = 0
    root_line = split_file.readline()  # get the root line that needs to be put in each new file
    files_made = 1
    current_file = make_file(total_folder, str(files_made) + "_" + file_name[2:])
    current_file.write(root_line)
    for line in split_file:
        if lines_in_new_file < lines_in_each_file:
            current_file.write(line)
            lines_in_new_file += 1
            if lines_in_new_file == lines_in_each_file:
                current_file.close()
                lines_in_new_file = 0
                if files_made != num_split_files:  # haven't made the desired number of files yet
                    files_made += 1
                    current_file = make_file(total_folder, str(files_made) + "_" + file_name[2:])
                    current_file.write(root_line)
        else:  # on the last file, write every line to this file
            current_file.write(line)

def main():
    split_file(file_name, num_split_files)

if __name__ == "__main__":
    file_name = sys.argv[1]
    num_split_files = int(sys.argv[2])
    main()