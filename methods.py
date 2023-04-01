import os
import natsort
import yaml
import pandas as pd

#natsort a directory
def nsort(dir, full_path):
    nsort = []
    for filename in os.listdir(dir):
        if full_path == True:
            nsort.append(os.path.join(dir, filename).replace("\\", "/"))
        else:
            nsort.append(filename)
    nsort = natsort.natsorted(nsort)
    return nsort

#creates the config file for the snakefile
def config(outdir, d):
    path = os.path.join(f"{outdir}", "config.yaml").replace("\\", "/")
    with open (path, "w") as outfile:
        yaml.dump(d, outfile)
    return path

#checks to make sure number of reads is even and above 0
def check_reads(size):
    if size == 0 or size % 2 != 0:
        raise ValueError("There must be an even number of reads greater than 0 in the folder.")

#takes all of the annotation output files and creates a FASTA file from it
def fasta(all_fasta, outdir):
    # numerically order files
    num_order = nsort(all_fasta, False)

    repeat_seqs = set()
    uid_tracker = {}
    uid = 1
    protein_tracker = {}

    #loop over the annotations folder and record every line
    for filename in num_order:
        head = ''
        header = True
        with open (os.path.join(all_fasta, filename).replace("\\", "/"), 'r') as file:
            for line in file:
                z = line.split("\t")
                entry = [x.strip() for x in z]
                rpnum = filename.split(".")[0].split("_")[0]

                #ignore header line
                # if "Best_Hit_ARO" not in z:
                if header != True:
                    
                    # arg_tracker.append((z[8].replace(' ', '_'), z[17].replace(' ', '_'), z[16].replace(' ', '_'), z[14], rpnum, z[18], z[19]))
                    #do not add anything with repeat sequences
                    if entry[17] not in repeat_seqs or entry[18] not in repeat_seqs or entry[19] not in repeat_seqs:
                        repeat_seqs.update((entry[17], entry[18], entry[19]))

                        #ensure there exists a protein-DNA sequence combo before adding the entry
                        #number every successful entry add with a universal id and store in dictionary
                        if entry[18] != "" and entry[17] != "":
                            uid_tracker[uid] = (rpnum, entry, entry[18])
                            uid+=1
                        elif entry[19] != "" and entry[17] != "":
                            uid_tracker[uid] = (rpnum, entry, entry[19])
                            uid+=1
                elif header == True:
                    head = z
                header = False

    #create FASTA file
    os.mkdir(outdir)
    with open (os.path.join(outdir, "annotations.FASTA").replace("\\", "/"), 'w+') as f:
        for x in uid_tracker:
            f.write(f">{x}|{uid_tracker[x][0]}\n{uid_tracker[x][1][17]}\n")
            if uid_tracker[x][0] not in protein_tracker:
                protein_tracker[uid_tracker[x][0]] = [[x, uid_tracker[x][1], uid_tracker[x][2]]]
            else:
                protein_tracker[uid_tracker[x][0]].append([x, uid_tracker[x][1], uid_tracker[x][2]])

    #create .faa files          
    os.mkdir(os.path.join(outdir, "protein_files").replace("\\", "/"))
    for x in protein_tracker:
        with open (os.path.join(outdir, f"protein_files/{x}.faa").replace("\\", "/"), 'w+') as f:
            for y in protein_tracker[x]:
                # f.write(f">{x}|{y[0]}\n{y[2]}\n")
                f.write(f">{y[0]}|{x}\n{y[2]}\n")

    return uid_tracker, protein_tracker, os.path.join(outdir, "annotations.FASTA").replace("\\", "/"), os.path.join(outdir, "protein_files").replace("\\", "/"), head

#make the three files for metagenomseq
def counts(uid_tracker, kallisto, shortbred, outdir):
    first_col = []
    outer = []
    rpnums = []
    df = pd.DataFrame()
    df2 = pd.DataFrame()
    checked = False
    gene_count = 1

    num_order = nsort(kallisto, False)

    #iterate through kallisto files and create matrix
    for filename in num_order:
        f = os.path.join(kallisto, filename).replace("\\", "/")
        with open (f, 'r') as file:
            inner = []
            #get number from filename
            inner.append(filename.split("_")[0])
            rpnums.append(filename.split("_")[0])
            for line in file:
                z = line.split("\t")
                count = z[3]
                #ignore header
                if z[3] != "est_counts":
                    inner.append(count)

                #increment gene count if first loop
                if checked == False:
                    #add onto first column of kallisto counts by IDing every gene
                    first_col.append(gene_count)
                    gene_count+=1

            #save each kallisto file into its own list
            outer.append(inner)
            checked = True
    
    #reorder columns numerically
    #if you ever have to take a number out, then q=i+1 doesn't work
    #ex: since rp7 isn't usable, 7 will be indexed over, but 17 will not
    fin = []
    # for i in range(len(outer)):
    #     q = i+1
    #     for y in range(len(outer)):
    #         if int(outer[y][0]) == q:
    #             fin.append(outer[y])
    for i in rpnums:
        for y in range(len(outer)):
            if str(outer[y][0]) == i:
                print(str(outer[y][0]))
                fin.append(outer[y])
    
    del first_col[-1]

    #place rows in df to create columns for counts file
    df["ID"] = first_col

    for x in fin:
        header = f"read_pair_{x[0]}"
        print(header)
        x.pop(0)
        df[header] = x

    #create kallisto counts csv
    df.to_csv(os.path.join(outdir, "kallisto_counts.csv").replace("\\", "/"), sep=',', index=False)
    
    #shortbred extraction
    #get the number of files and number of shortbred entries
    num_order = nsort(shortbred, False)
    max = int(list(uid_tracker)[-1])+1
    # sample_total = len(num_order)

    #created double nested dictionary to house shortbred data per [read pair][gene]
    outer = {}
    for x in rpnums:
        outer[str(x)] = {}
        for y in range(1, max):
            outer[str(x)][str(y)] = 0

    #loop over shortbred files and obtain proper dictionary entries
    for filename in num_order:
        f = os.path.join(shortbred, filename).replace("\\", "/")
        with open (f, 'r') as file:
            for line in file:
                z = line.split("\t")
                if "Family" not in z:
                    vars = z[0].split("|")
                    if len(vars) > 1:
                        outer[str(vars[1])][str(vars[0])] = z[1]

    #create OTU label column
    df2["ID"] = first_col

    #extract assigned, in order dictionary values
    fin = []
    for x in outer:
        inner = []
        inner.append(x)
        for y in outer[x]:
            inner.append(outer[x][y])
        fin.append(inner)
    
    #create columns
    for x in fin:
        header = f"read_pair_{x[0]}"
        x.pop(0)
        df2[header] = x

    #create shortbred counts csv
    df2.to_csv(os.path.join(outdir, "shortbred_counts.csv").replace("\\", "/"), sep=',', index=False)

    return first_col

#create otu table
#uid_tracker[x]: rpnum, entry
#NEED TWO COLUMNS SAYING IF DRUG CLASSES/GENE FAMILIES ARE SINGLE/MULTIPLE
def observations(uid_tracker, head, first_col, treatments, outdir):
    df = pd.DataFrame()
    df["ID"] = first_col
    outer = {}

    #extract user specified sample - treatment info
    treatments_dict = {}
    header = True
    with open (treatments, 'r') as file:
        for line in file:
            z = line.split(",")
            print(z)
            # z = [o.replace("\n", "").replace("\t", "").strip for o in z]
            if header == True:
                for x in z[1:]:
                    head.append(x)
            else:
                treatments_dict[z[0]] = z[1:]
            header = False
    i = 0
    
    head.append("DC_MULTIPLE")
    head.append("AMRGF_MULTIPLE")
    head.append("read_pair_number")

    print(treatments_dict)

    #create a nested dictionary for every otu row
    # for x in range(1, int(list(uid_tracker)[-1])+1):
    #     outer[str(x)] = {}
    for x in head:
        outer[x] = []

    #dictionary the rgi information
    for x in uid_tracker:
        otu_info = uid_tracker[x][1]

        #add on the user supplied treatment information to otu rows
        for q in treatments_dict[uid_tracker[x][0]]:
            otu_info.append(q)
        
        #add on if drug class/amr gene family is multiple
        if ";" in otu_info[14]:
            otu_info.append("MULTIPLE")
        else:
            otu_info.append("SINGLE")
        if ";" in otu_info[16]:
            otu_info.append("MULTIPLE")
        else:
            otu_info.append("SINGLE")
        
        #append read pair number last
        otu_info.append(uid_tracker[x][0])

        # print(otu_info)
        #get index of otu row info to match the corresponding header name in head
        for y in range(len(otu_info)):
            outer[head[y]].append(str(otu_info[y]).replace("\t", "").replace("\n", "").strip())
        
    #extract assigned, in order dictionary values
    for x in outer:
        df[x.strip("\n")] = outer[x]

    df.to_csv(os.path.join(outdir, "observations.csv").replace("\\", "/"), sep=',', index=False)

# def genomad(plasmid, virus, outdir):