import os

def blast_identify(blast_input,ouput_path,filter_value):
    """
    filter blast result
    """
    blast_output = os.path.join(ouput_path,'blast_filter.txt')
    with open(blast_input,'r') as f, open(blast_output,'w') as f_out :
        lines = f.readlines()
        my_dict = {}
        for line in lines:
            line = line.strip()
            items = line.split()
            value2 = float(items[2])
            if value2 >= filter_value:
                f_out.write(line + '\n')
    return
