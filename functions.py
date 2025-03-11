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
        raise ValueError("There must be an even number of reads greater than 0 in the folder. Ensure the only files in your reads folder are the read pairs.")

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
                    #17 = predicted DNA, 18 = predicted protein, 19 = CARD protein
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
            #remove any stray commas that will interfere with re-reading the observation matrix
            arg_info[y] = arg_info[y].replace(",", "")
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
    all_plasmid = [[], []]
    virus_set = set()
    all_virus = [[], []]

    #iterate over plasmid files and get contig names
    for filename in os.listdir(plasmid):
        f = os.path.join(plasmid, filename)
        rpnum = filename.split("_")[0]
        with open (f, 'r') as file:
            for line in file:
                z = line.split("\t")
                plasmid_set.add((z[0], rpnum))
                all_plasmid[0].append(rpnum)
                all_plasmid[1].append(z[0])
                
    #iterate over virus files and get contig names
    for filename in os.listdir(virus):
        f = os.path.join(virus, filename)
        rpnum = filename.split("_")[0]
        with open (f, 'r') as file:
            for line in file:
                z = line.split("\t")
                virus_set.add((z[0], rpnum))
                all_virus[0].append(rpnum)
                all_virus[1].append(z[0])

    #go through observations file and match ARG contig name with virus or plasmid contig name
    plasmid_col = []
    virus_col = []
    header = True
    rp_index = 0
    with open (f"{outdir}/observations.csv", 'r') as file:
        for line in file:
            z = line.split(",")

            #find the index of the "read_pair_number" column
            if header == True:
                for i in range(len(z)):
                    if z[i].strip('\n') == "read_pair_number":
                        rp_index = i

            #find matching contig names
            elif header != True:
                q = z[2].split("_")
                del q[-1]
                name = "_".join(q)
                combo = (name, z[rp_index].replace("\n", ""))

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

    df2 = pd.DataFrame()
    df2["read_pair_number"] = all_plasmid[0]
    df2["contig"] = all_plasmid[1]
    df2.to_csv(f"{outdir}/all_plasmid_contigs.csv", index=False)

    df3 = pd.DataFrame()
    df3["read_pair_number"] = all_virus[0]
    df3["contig"] = all_virus[1]
    df3.to_csv(f"{outdir}/all_virus_contigs.csv", index=False)

#append integron_finder data to observations matrix
def integron(files, outdir):
    integron_dict = {}
    all = [[], []]

    #get integron contigs and counts
    #if there are no integrons found, then there is a "#" as the first character on the 2nd line
    for filename in os.listdir(files):
        header = True
        f = os.path.join(files, filename)
        rpnum = filename.split("_")[0]
        with open (f, 'r') as file:
            for line in file:
                #if any #'s are in the line it is not tab separated
                z = line.split("\t")
                if len(z) > 1:
                    if header == True:
                        header = False
                    else:
                        if (z[1], rpnum) not in integron_dict:
                            integron_dict[(z[1], rpnum)] = 1
                            all[0].append(rpnum)
                            all[1].append(z[1])
                        else:
                            integron_dict[(z[1], rpnum)] += 1
    
    integron_presence = []
    integron_count = []
    header = True
    rp_index = 0
    with open (f"{outdir}/observations.csv", 'r') as file:
        for line in file:
            z = line.split(",")

            #find the index of the "read_pair_number" column
            if header == True:
                for i in range(len(z)):
                    if z[i].strip('\n') == "read_pair_number":
                        rp_index = i

            #find matching contig names
            elif header != True:
                q = z[2].split("_")
                del q[-1]
                name = "_".join(q)
                combo = (name, z[rp_index].replace("\n", ""))

                #check if contig name, rpnum combo exist for integron presence
                if combo in integron_dict:
                    integron_presence.append("YES")
                    integron_count.append(integron_dict[combo])
                else:
                    integron_presence.append("NO")
                    integron_count.append(0)
        
            header = False

    #add new columns to observations.csv and re-write it
    df = pd.read_csv(f"{outdir}/observations.csv")
    df['INTEGRON'] = integron_presence
    df['INTEGRON_COUNT'] = integron_count
    df.to_csv(f"{outdir}/observations.csv", index=False)

    df2 = pd.DataFrame()
    df2["read_pair_number"] = all[0]
    df2["contig"] = all[1]
    df2.to_csv(f"{outdir}/all_integron_contigs.csv", index=False)

#append BLAST found ICEs, ISs, and transposons to observations matrix
#THE OUTDIR FOLDER WILL ALREADY HAVE BEEN MADE
def blast_mges(folder, outdir):
    #order files in folder
    sorted = nsort(folder, False)

    #all 3 possible MGEs that could be labels in the files
    labels = {"IS": {}, "IC": {}, "Tn": {}} #will house another dictionary for each rpnum
    cols = {"IS": {}, "IC": {}, "Tn": {}} #will be appended to the observations matrix
    all = {"IS": {}, "IC": {}, "Tn": {}} #stores all the rpnum-contig pairs to be used later

    #for things that dont need a rpnum inner dictionary initialize outside of loop
    for label in labels:
        cols[label] = []
        all[label] = [[], []]

    #iterate over files and get contig names
    for filename in sorted:
        f = os.path.join(folder, filename)
        # rpnum = filename.split("_")[-1].split(".")[0]
        rpnum = filename.split("_")[0]

        #create the storage structure for each mge and read pair combo
        #this is the structure that will be checked against the contig
        #names in the observations file
        for label in labels:
            labels[label][rpnum] = set()

        #open each file
        with open (f, 'r') as file:
            #for each MGE found, add the contig to the appropriate set
            for line in file:
                z = line.split("\t")
                contig = z[0]
                mge = z[1][:2]
                
                #only pick mges that match the MGE labels we are looking for
                if mge in list(labels.keys()):
                    labels[mge][rpnum].add(contig)
                    all[mge][0].append(rpnum)
                    all[mge][1].append(contig)
                 
    #go through observations file and check each contig name to see
    #if it is associated with an MGE or not
    header = True
    rp_index = 0
    with open (f"{outdir}/observations.csv", 'r') as file:
        for line in file:
            z = line.split(",")

            #find the index of the "read_pair_number" column
            if header == True:
                for i in range(len(z)):
                    if z[i].strip('\n') == "read_pair_number":
                        rp_index = i

            #find matching contig names
            elif header != True:
                q = z[2].split("_")
                del q[-1]
                contig = "_".join(q)
                rpnum = z[rp_index].replace("\n", "")

                #for every MGE type, check contig against it
                #if it appears, then YES in the column, if not then NO
                for label in labels:
                    if contig in labels[label][rpnum]:
                        cols[label].append("YES")
                    else:
                        cols[label].append("NO")

            header = False

    #add new columns to observations.csv and re-write it
    df = pd.read_csv(f"{outdir}/observations.csv")
    for col in list(cols.keys()):
        #only update that specific MGE column if there were actually any contig matches
        if len(all[col][0]) > 0:
            print(f"MGE: {col}, SIZE OF SECTION: {len(all[col][0])}")
            df[col] = cols[col]
    df.to_csv(f"{outdir}/observations.csv", index=False)

    #stores all the contigs to be accessed later
    for col in list(all.keys()):
        #same as above, only create this file if the MGE actually had any contig matches
        if len(all[col][0]) > 0:
            df2 = pd.DataFrame()
            df2["read_pair_number"] = all[col][0]
            df2["contig"] = all[col][1]
            df2.to_csv(f"{outdir}/all_{col}_contigs.csv", index=False)

#append cat taxa and coverage to observations
def cat(files, outdir):
    #get all the args
    arg_set = set()
    arg_contig_names = []

    add_otu = {"CAT_genus": [], "CAT_species": [], "CAT_taxa_coverage": []}

    header = True
    rp_index = 0
    #open up observation file and obtain contig names
    with open (f"{outdir}/observations.csv", 'r') as file:
        for line in file:
            z = line.split(",")

            #find the index of the "read_pair_number" column
            if header == True:
                for i in range(len(z)):
                    if z[i].strip('\n') == "read_pair_number":
                        rp_index = i
            elif header != True:
                #remove last underscore and number
                q = z[2].split("_")
                del q[-1]
                name = "_".join(q)
                combo = (str(name), str(z[rp_index].replace("\n", "")))
                arg_set.add(combo)
                arg_contig_names.append(combo)
                rpnum = str(z[rp_index].replace("\n", ""))

            header = False

    total_cov_dict = {}
    taxa_total_cov_dict = {}
    contig_num = 1
    snp_set = set()
    search = {}
    #loop over filenames in folder
    for filename in os.listdir(files):
        total_cov = 0
        header = True
        # rpnum = filename.split("_")[1].split(".")[0]
        rpnum = filename.split("_")[0]
        taxa_total_cov_dict[rpnum] = {}
        search[rpnum] = {}
        
        #loop over the CAT files
        with open (os.path.join(files, filename).replace("\\", "/"), 'r') as file:
            for line in file:
                #species will remain NA unless a species is present in the line
                species = "NA"
                genus = "NA"
                z = line.split("\t")
                z[-1] = z[-1].strip()

                #make sure header isn't grabbed
                if header != True:
                    contig_name = z[0]
                    cov = float(z[0].split("_")[-1])
                    total_cov += cov

                    #if there is lineage, look for species name
                    if len(z) > 3:
                        if ";" in z[3]: 
                            entry = z[3].split(";")
                            #grab genus and species entries if they exist
                            for taxa in entry:
                                if taxa[:3] == "g__":
                                    genus = taxa[3:]
                                if taxa[:3] == "s__":
                                    species = taxa[3:]

                    #add taxa to dictionary
                    if species not in taxa_total_cov_dict[rpnum]:
                        taxa_total_cov_dict[rpnum][species] = cov
                    else:
                        taxa_total_cov_dict[rpnum][species] += cov

                    #add stuff for otu column
                    combo = (str(contig_name), str(rpnum))
                    snp_set.add(combo)
                    search[rpnum][contig_name] = (genus, species, cov, contig_name, rpnum)

                header = False
                contig_num+=1
            
            #match rpnum to total coverage
            total_cov_dict[rpnum] = total_cov
    
    argcontigs = set()

    #iterate over collected ARG contig names and find matches to CAT contigs
    for x in arg_contig_names:
        contig_name = x[0]
        rpnum = x[1]
        #if you find correct rpnum+contig name, get info
        if contig_name in search[rpnum]:
            info = search[rpnum][contig_name]
            # print(f"RP: {rpnum} | SPECIES: {info[0]} | CONTIG: {contig_name}")
            add_otu["CAT_genus"].append(info[0])
            add_otu["CAT_species"].append(info[1])
            add_otu["CAT_taxa_coverage"].append(info[2])
            argcontigs.add(contig_name)
        else:
            add_otu["CAT_genus"].append("NA")           
            add_otu["CAT_species"].append("NA")
            add_otu["CAT_taxa_coverage"].append("NA")
    
    df = pd.read_csv(f"{outdir}/observations.csv")
    df['CAT_GENUS'] = add_otu["CAT_genus"]
    df['CAT_SPECIES'] = add_otu["CAT_species"]
    df['CAT_TAXA_COVERAGE'] = add_otu["CAT_taxa_coverage"]
    df.to_csv(f"{outdir}/observations.csv", index=False)

    df2 = pd.DataFrame()
    nonarg_contigs = {"read_pair_number": [], "contig": [], "CAT_genus": [], "CAT_species": [], "CAT_taxa_coverage": []}
    #iterate over all collected contigs and pick out the ones that weren't picked
    # Loop through the dictionary without using .items()
    for rpnum in search:
        contigs = search[rpnum]
        for contig_name in contigs:
            value = contigs[contig_name]
            if (rpnum, contig_name) not in argcontigs:
                nonarg_contigs["contig"].append(value[3])  # Use contig_name, not value[-1]
                nonarg_contigs["read_pair_number"].append(f"read_pair_{value[4]}")
                nonarg_contigs["CAT_genus"].append(value[0])
                nonarg_contigs["CAT_species"].append(value[1])
                nonarg_contigs["CAT_taxa_coverage"].append(value[2])
                argcontigs.add((rpnum, contig_name))
    # Loop through the dictionary to assign column names and values
    for x in nonarg_contigs:
        df2[x.strip("\n")] = nonarg_contigs[x]

    df2.to_csv(f"{outdir}/ww_nonarg_cat.csv", sep=',', index=False)
