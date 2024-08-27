import copy
import sys


def convert_file(input_file, output_file_phy, output_file_nex):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    species_set, num_part = get_species(lines)
    alignment = get_alignment_dict(species_set, num_part)
    boundry_list = []
    part_count = -1
    for line in lines:
        data = list(line.strip().split(" "))
        if data[0].isnumeric():
            boundry_list.append(int(data[-1]))
            part_count += 1
        if not data[0].isnumeric() and not data[0].isspace() and not len(data[0]) == 0:
            alignment[part_count][data[0]] = data[-1]

    alignment_processed = [str(len(species_set)) + " ", str(sum(boundry_list)) + "\n", "\n"]
    for i in species_set:
        sequence = ""
        for j in range(len(alignment)):
            if i in alignment[j] and len(alignment[j][i]) > 1:
                sequence += alignment[j][i]
            else:
                sequence += "-" * boundry_list[j]
        else:
            sequence += "\n"
            alignment_processed.append(i + " ")
            alignment_processed.append(sequence)

    ############## PHYLIP File Generation ################
    ######################################################
    alignment_merged = "".join(alignment_processed)
    with open(output_file_phy, 'w') as file:
        file.write(alignment_merged)

    ############## NEXUS File Generation ################
    #####################################################

    nexus_str = "#nexus\n"
    nexus_str += "begin sets;" + "\n"
    counter_sum = 1
    part_counter = 1
    for i in boundry_list:
        first_pos = counter_sum
        second_pos = counter_sum + i - 1
        nexus_str += "  charset part{} = {}-{};\n".format(part_counter, first_pos, second_pos)
        counter_sum += i
        part_counter += 1
    # nexus_str += "  charpartition mine ="
    # for i in boundry_list:
    #     nexus_str += " MF:part{}".format(i)
    #     if i < len(boundry_list)-1: nexus_str += ","
    # nexus_str += ";\n"
    nexus_str += "end;" + "\n"

    with open(output_file_nex, 'w') as file:
        file.write(nexus_str)


def get_species(lines):
    species = []
    counter = 0
    for i in lines:
        data = list(i.strip().split(" "))
        if not data[0].isnumeric() and not data[0].isspace() and not len(data[0]) == 0:
            species.append(data[0])
            print(data[0])
        if data[0].isnumeric():
            counter += 1
    return list(set(species)), counter


def get_alignment_dict(species, num_partitions):
    alignment = []
    species_dict = {}
    for i in species:
        species_dict[i] = "_"
    alignment = [copy.deepcopy(species_dict) for i in range(num_partitions)]
    return alignment


if __name__ == '__main__':
    print(sys.argv)
    # input_file = sys.argv[1]
    # output_file_phy = sys.argv[2]
    # output_file_nex = sys.argv[3]

    input_file = "/home/piyumal/PHD/TimeTree/Hessian_Reversible_Models/Empirical data/Alvarez-Carretero_etal_SI/aln/00_step1/alignment_4parts.aln"
    output_file_phy = "/home/piyumal/PHD/TimeTree/Hessian_Reversible_Models/Empirical data/Alvarez-Carretero_etal_SI/aln/00_step1/alignment_4parts_iqtree.phy"
    output_file_nex = "/home/piyumal/PHD/TimeTree/Hessian_Reversible_Models/Empirical data/Alvarez-Carretero_etal_SI/aln/00_step1/alignment_4parts_iqtree.nex"
    convert_file(input_file, output_file_phy, output_file_nex)
    # a = '/home/piyumal/test/IQTREE_file_format_conversion/Xenarthra_5parts.aln'
    # b = '/home/piyumal/test/IQTREE_file_format_conversion/Xenarthra_5parts_test1.phy'
    # c = '/home/piyumal/test/IQTREE_file_format_conversion/Xenarthra_5parts_test1.nex'
    # convert_file(a, b, c)

    print("ALL the files generated")
