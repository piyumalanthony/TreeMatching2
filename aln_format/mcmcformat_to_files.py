
def convert_to_files(input_file, output_file_phy):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    alignment_started = False
    alignment_id = 0
    local_lines = []
    species_count = 0
    counter = 0
    for line in lines:
        data = list(line.strip().split(" "))
        print(data[0])
        if data[0].isnumeric():
            print(data)
            alignment_started = True
            species_count = int(data[0])
            print(f'species count: {species_count}')

        if alignment_started and counter > 2 and data[0].isnumeric():
            alignment_id += 1
            with open(f'{output_file_phy}_part_{alignment_id}.phy', 'w') as file:
                file.write("".join(local_lines))
            local_lines = [line]
            print(f'generated {output_file_phy}_part_{alignment_id}.phy')

        else:
            local_lines.append(line)
            counter += 1
    else:
        with open(f'{output_file_phy}_part_{alignment_id+1}.phy', 'w') as file:
            file.write("".join(local_lines))
        print(f'generated {output_file_phy}_part_{alignment_id+1}.phy')


if __name__ == '__main__':

    input_file = "/home/piyumal/PHD/TimeTree/Hessian_Reversible_Models/Discrepancy_analysis/Mammals_backbone/data/alignment_4parts.aln"
    output_file_phy = "/home/piyumal/PHD/TimeTree/Hessian_Reversible_Models/Discrepancy_analysis/Mammals_backbone/data/mammal_backbone"
    convert_to_files(input_file, output_file_phy)


    print("ALL the files generated")