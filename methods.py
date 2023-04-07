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
def config(d, name, outdir):
    path = os.path.join(f"{outdir}", f"{name}.yaml").replace("\\", "/")
    with open (path, "w") as outfile:
        yaml.dump(d, outfile)
    return path

#checks to make sure number of reads is even and above 0
def check_reads(size):
    if size == 0 or size % 2 != 0:
        raise ValueError("There must be an even number of reads greater than 0 in the folder.")

#takes all of the annotation output files and creates a dictionary of
#labelled ARGs used to create the main FASTA file and individual .faa filesj
def fasta(all_rgi):
    # numerically order files
    num_order = nsort(all_rgi, False) 
    repeat_seqs = set()
    uid_tracker = {}
    uid = 1
    protein_tracker = {}

    #loop over the RGI output files and record every ARG
    for filename in num_order:
        head = ''
        header = True
        with open (os.path.join(all_rgi, filename).replace("\\", "/"), 'r') as file:
            for line in file:
                #split by tab, remove newlines from every split value, 
                #and extract sample number from filename (number will always be first)
                z = line.split("\t")
                entry = [x.strip() for x in z]
                rpnum = filename.split(".")[0].split("_")[0]

                #ignore header line
                if header != True:
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

                #if in the header, grab all the column names
                elif header == True:
                    head = z
                header = False

    return uid_tracker, protein_tracker, head

#make the kallisto and shortbred abundance matrices
def count_matrices(uid_tracker, kallisto, shortbred, outdir):
    #create outdir
    try:
        os.mkdir(outdir)
    except OSError as error:
        print(error)

    #initliaze the two df's for eventual csv creation
    df = pd.DataFrame()
    df2 = pd.DataFrame()

    first_col = []
    outer = []
    rpnums = []
    checked = False
    gene_count = 1

    #sort the files in numerical order
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
    fin = []
    for i in rpnums:
        for y in range(len(outer)):
            if str(outer[y][0]) == i:
                fin.append(outer[y])
    
    del first_col[-1]

    #place rows in df to create columns for counts file
    df["ARG_ID"] = first_col

    #attach all other header names and values to the dataframe
    for x in fin:
        header = f"read_pair_{x[0]}"
        x.pop(0)
        df[header] = x

    #create kallisto counts csv
    df.to_csv(os.path.join(outdir, "kallisto_counts.csv").replace("\\", "/"), sep=',', index=False)
    
    #shortbred extraction

    #order the files numerically and
    #get the number of files and number of shortbred entries
    num_order = nsort(shortbred, False)
    max = int(list(uid_tracker)[-1])+1

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

    #create id label column
    df2["ARG_ID"] = first_col

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

#create observation matrix
def observations(uid_tracker, head, first_col, treatments, outdir):
    df = pd.DataFrame()
    df["ARG_ID"] = first_col
    outer = {}

    #extract user provided sample metadata
    treatments_dict = {}
    header = True
    with open (treatments, 'r') as file:
        for line in file:
            z = line.split(",")

            #grab header names if on first line
            if header == True:
                for x in z[1:]:
                    head.append(x)
            
            #otherwise take first column of sample ids as keys and metadata info as values
            else:
                treatments_dict[z[0]] = z[1:]
            header = False
    i = 0
    
    #name some pre-established column headers for the matrix
    head.append("DC_MULTIPLE")
    head.append("AMRGF_MULTIPLE")
    head.append("read_pair_number")

    #create an empty list for each column header that requisite info will be appended to
    for x in head:
        outer[x] = []

    #access the ARG tracker and store the info
    for x in uid_tracker:
        arg_info = uid_tracker[x][1]

        #add on the user supplied treatment information to arg rows
        #by accessing the sample that the ARG is from
        for q in treatments_dict[uid_tracker[x][0]]:
            arg_info.append(q)
        
        #add on if drug class/amr gene family is multiple
        if ";" in arg_info[14]:
            arg_info.append("MULTIPLE")
        else:
            arg_info.append("SINGLE")
        if ";" in arg_info[16]:
            arg_info.append("MULTIPLE")
        else:
            arg_info.append("SINGLE")
        
        #append read pair number last
        arg_info.append(uid_tracker[x][0])

        #get index of arg row info to match the corresponding header name in head
        for y in range(len(arg_info)):
            outer[head[y]].append(str(arg_info[y]).replace("\t", "").replace("\n", "").strip())
        
    #extract assigned, in order dictionary values
    #and append the headers and values to the dataframe
    for x in outer:
        df[x.strip("\n")] = outer[x]

    #create the observations matrix
    df.to_csv(os.path.join(outdir, "observations.csv").replace("\\", "/"), sep=',', index=False)

#append genomad data to observations matrix
def genomad(plasmid, virus, outdir):
    plasmid_set = set()
    virus_set = set()

    #iterate over plasmid files and get contig names
    for filename in os.listdir(plasmid):
        f = os.path.join(plasmid, filename)
        rpnum = filename.split("_")[0]
        with open (f, 'r') as file:
            for line in file:
                z = line.split("\t")
                plasmid_set.add((z[0], rpnum))
                

    #iterate over virus files and get contig names
    for filename in os.listdir(virus):
        f = os.path.join(virus, filename)
        rpnum = filename.split("_")[0]
        with open (f, 'r') as file:
            for line in file:
                z = line.split("\t")
                virus_set.add((z[0], rpnum))

    #go through observations file and match 
    #ARG contig name with virus or plasmid contig name
    plasmid_col = []
    virus_col = []
    header = True
    with open (f"{outdir}/observations.csv", 'r') as file:
        for line in file:
            z = line.split(",")
            if header != True:
                q = z[2].split("_")
                del q[-1]
                name = "_".join(q)
                combo = (name, z[-1].replace("\n", ""))
                #check if contig name, rpnum combo exist in either plasmid or virus sets
                if combo in plasmid_set:
                    plasmid_col.append("YES")
                else:
                    plasmid_col.append("NO")
                if combo in virus_set:
                    virus_col.append("YES")
                else:
                    virus_col.append("NO")
                
            header = False

    #add new columns to observations.csv and re-write it
    df = pd.read_csv(f"{outdir}/observations.csv")
    df['PLASMID_MGE'] = plasmid_col
    df['VIRUS_MGE'] = virus_col
    df.to_csv(f"{outdir}/observations.csv", index=False)

#append integron_finder data to observations matrix
def integron(files, outdir):
    return None

def spraynpray():
    return None