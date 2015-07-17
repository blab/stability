import sys
import os

pdb_name = sys.argv[1]
pdb_original_name = pdb_name + "_original.pdb"
pdb_repaired_name = pdb_name + "_repaired.pdb"
pdb_nochain_name = pdb_name + "_repaired_nochain.pdb"
pdb_formatted_name = pdb_name + "_formatted.pdb"

def make_repair_runfile():
    output = pdb_repaired_name[:len(pdb_repaired_name) - 4]  # since it automatically adds on ".pdb"
    runfileName = "repair_runfile.txt"
    runfile = open(runfileName, 'w')
    runfile.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<RepairPDB>%s;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;' % (pdb_original_name, output))

def remove_chain_info():
    repaired_file = open(pdb_repaired_name, 'r')
    nochain_file = open(pdb_nochain_name, 'w')
    for line in repaired_file:
        nochain_line = line
        if len(line) > 21:
            if line[21] == "A" or line[21] == "B":
                nochain_line = line[:21] + " " + line[22:]
        nochain_file.write(nochain_line)
    repaired_file.close()
    nochain_file.close()

def adjust_sites():

    nochain_file = open(pdb_nochain_name, 'r')
    formatted_file = open(pdb_formatted_name, 'w')

    foundTER = False
    offset_number = 0

    # some of the structures are off from the normal sequence numbering we use
    adjust = 0   # normal
    #adjust = 40  # 4WE6

    for line in nochain_file:
        newLine = line
        if len(line) > 26:
            current_line_number = int(line[23:26])
            write_number = current_line_number + adjust
            if foundTER:
                write_number = current_line_number + offset_number + adjust
            newLine = line[:23] + str(write_number) + line[26:]
        if line[0:3] == "TER":
            foundTER = True
            offset_number = current_line_number
        formatted_file.write(newLine)
    formatted_file.close()
    os.remove(pdb_nochain_name)

def main():
    make_repair_runfile()
    os.system("./foldx3b6 -runfile repair_runfile.txt")
    remove_chain_info()
    adjust_sites()

main()