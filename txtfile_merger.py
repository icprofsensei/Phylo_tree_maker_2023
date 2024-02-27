filenames = ['microbiome_txt_files/AIRWAYids.txt', 'microbiome_txt_files/FAECALids.txt', 'microbiome_txt_files/SKINids.txt', 'microbiome_txt_files/URINARYids.txt', 'microbiome_txt_files/VAGINALids.txt']
with open('master', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
                                    