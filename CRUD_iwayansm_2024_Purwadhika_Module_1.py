import hashlib, os, pickle, atexit

########### | INISIASI DATABASE | ################################################################# DONE
f_name = 'gen_db.pickle'
gen_db = [{'ID.REF':'','SPECIES':'','DNA SEQUENCE':'','%P/Pi':'','AMINO ACID':''}]

########### | DATABASE | ########################################################################## DONE
if not os.path.exists(f_name):
    with open(f_name, 'wb') as f:
        pickle.dump(gen_db, f)
else:
    with open(f_name, 'rb') as f:
        db = pickle.load(f)
        gen_db = db

########### | INPUT SUBMISSION | ################################################################## DONE
########### DNA SEQUENCE ######### DONE
def AA(dna_seq):
    dna = str(dna_seq)
    mRNA = dna.translate(str.maketrans({'A': 'U', 'C': 'G', 'T': 'A', 'G': 'C'}))
    rRNA = mRNA[::-1]
    codon_library = {
                    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
                    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
                    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
                    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
                    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
                    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
                    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
                    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
                    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
                    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
                    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
                    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
                    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
                    }
    #   Codon is a triplet pair of nucleotides
        #   Clancy, S. & Brown, W. (2008) Translation: DNA to mRNA to Protein. Nature Education 1(1):101
    protein_seq = ''
    if len(rRNA)%3 == 0: 
        for i in range(0, len(rRNA), 3): 
            codon = rRNA[i:i + 3] 
            protein_seq+= codon_library[codon]
    amino_acid = ' \u2666 '.join((f'{protein_seq}')[i:i+3] for i in range(0, len((f'{protein_seq}')), 3))
    return amino_acid
    #   Timin (T) in DNA is replaced with Urasil (U) in RNA, corresponding nitrogenous bases [A:T/U] and [C:G]
        #   Two classes of nitrogenous bases are Purines consist of Adenine (A) and Guanine (G), while Pyrimidines consist of Thymine (T), Uracil (U), and Cytosine (C)
            #   [https://doi.org/10.1016/B978-0-12-802823-0.00001-8]

def prop(dna_seq):
    #   Tm = 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC) -- for sequences longer than 13 nucleotides
        #   wxyz = number of the bases A,T,G,C in the sequences, respectively
            #   Assumptions: Annealing occurs under the standard conditions of 50nM primer, 50mM Na+, and pH 7.0
                #   [http://insilico.ehu.eus/] 2003-2024@ University of the Basque Country
    Tm_comp = str(dna_seq)
    Tm_total = int(len(dna_seq))
    Tm_A = int(Tm_comp.count('A'))
    Tm_G = int(Tm_comp.count('G'))
    Tm_C = int(Tm_comp.count('C'))
    Tm_T = int(Tm_comp.count('T'))

    Tm_percA = Tm_A/Tm_total*100
    Tm_percG = Tm_G/Tm_total*100
    Tm_percC = Tm_C/Tm_total*100
    Tm_percT = Tm_T/Tm_total*100
    Tm_calc = 64.9 + 41*(Tm_G+Tm_C-16.4)/Tm_total

    Tm_point = '{:<84} {:>151}'.format('MELTING TEMPERATURE (Tm): 'f'{Tm_calc:.2f}''\u00B0C\t\t\t\t\tSTANDARD CONDITIONS [\x1B[1m50nM\x1B[0m primer, \x1B[1m50mM\x1B[0m Na+, and pH \x1B[1m7.0\x1B[0m]\n', 'CATEGORY \u21D2 \u2265 60\u00B0C\x1B[1m|HIGH Tm|\x1B[0m, BETWEEN 50-60\u00B0C\x1B[1m|IDEAL Tm|\x1B[0m, \u2264 50\u00B0C\x1B[1m|LOW Tm|\x1B[0m\n')
    Tm_perc = 'COMPOSITION OF NUCLEOTIDES:\n''\t\t\t\t|\t\x1B[1mA\x1B[0m\t|\t\x1B[1mT\x1B[0m\t|\t\x1B[1mC\x1B[0m\t|\t\x1B[1mG\x1B[0m\t|\n'f'\t\t\t\t|\t\u21D3{Tm_percA:.2f}%\t|\t\u21D3{Tm_percG:.2f}%\t|\t\u21D3{Tm_percC:.2f}%\t|\t\u21D3{Tm_percT:.2f}%\t|\n'

    GC = Tm_percG + Tm_percC
    GC_content = '{:<69} {:>69}'.format('GUANINE-CYTOSINE (GC) CONTENT: 'f'{GC:.2f}%', 'CATEGORY \u21D2 \u2265 80%\x1B[1m|HIGH|\x1B[0m, BETWEEN 50-55%\x1B[1m|IDEAL|\x1B[0m, \u2264 40%\x1B[1m|LOW|\x1B[0m\n')

    n_molw = 0.0
    n_weights = {'A':335.2, 'T':326.2, 'C':311.2, 'G':351.2}
    for n in dna_seq:
        n_molw += n_weights[n]
        n_molws = 'MOLECULAR WEIGHT: 'f'{n_molw:.2f} dalton(grams/mole)'
    #   MW = (nA×335.2)+(nC×311.2)+(nG×351.2)+(nT×326.2)+P
        #   where nx is the number of nucleotides of A, C, G, or T in the oligonucleotide and P is equal to −101.0 for dephosphorylated (lacking an end phosphate group) or 40.0 for phosphorylated oligonucleotides.
            #   The following demonstrate how molecular weight, molarity, and nucleic acid length relate to DNA quantity.
                #   [https://doi.org/10.1016/B978-012665751-7/50046-9]
        nprop = f'{Tm_point}\n{Tm_perc}\n{GC_content}\n{n_molws}'
    return nprop

def ppi(dna_seq):
    p_pi_comp = str(dna_seq)
    total = int(len(dna_seq))
    pct_p = p_pi_comp.count ('A') / total * 100 + p_pi_comp.count ('G') / total * 100
    pct_pi = p_pi_comp.count ('T') / total * 100 + p_pi_comp.count ('C') / total * 100
    p_pi = f'{pct_p:.2f}% / {pct_pi:.2f}%'
    return p_pi

def gen_seq():
    print('DNA sequence entry in this program is limited to \x1B[1m33-bp\x1B[0m in length')
    print('DNA sequence consists of four nucleotides \x1B[1mA\x1B[0m, \x1B[1mG\x1B[0m, \x1B[1mT\x1B[0m, \x1B[1mC\x1B[0m')
    print('DNA sequence begins with \x1B[1m5-\x1B[0m and ends with \x1B[1m-3\x1B[0m-ex: \x1B[1m5-ATCGAATGGCCAATAGAATTCCAATAGATAGCG-3\x1B[0m')
    while True:
        dna_seq = input('ENTER THE DNA SEQUENCE (\x1B[1m33-bp\x1B[0m): ')
        n_base = ['A', 'G', 'T', 'C']
        if len(dna_seq) == 33 and dna_seq.isupper():
            if all(x in n_base for x in dna_seq):
                print('\n', f'DNA SEQUENCE IS CONFIRMED\n\t5-\x1B[1m{dna_seq}\x1B[0m-3', '\n')
                return dna_seq
            else:
                print('DNA SEQUENCE CONSISTS OF ONLY FOUR NITROGEN BASES: A, G, T, AND C'.center(127, '-'), '\n')
        else:
            print('DNA SEQUENCE CANNOT BE EMPTY AND MUST BE 33 UPPERCASE LETTERS'.center(127, '-'), '\n')

########### SEQUENCING METHOD #### DONE
def seq_md():
    print('\x1B[1m01\x1B[0m: Sanger Sequencing; \x1B[1m02\x1B[0m: Fragment Analysis; \x1B[1m03\x1B[0m: Next-Generation Sequencing')
    while True:
        sq_md = input('ENTER ONE OF THE FOLLOWING SEQUENCING METHOD |\x1B[1m01\x1B[0m|  |\x1B[1m02\x1B[0m|  |\x1B[1m03\x1B[0m|: ')
        sq_md_opt = ['01', '02', '03']
        if sq_md in sq_md_opt:
            if sq_md == '01':
                print('\n', 'SEQUENCING METHOD\n\t\x1B[1mSANGER SEQUENCING\x1B[0m', '\n')
                break
            elif sq_md == '02':
                print('\n', 'SEQUENCING METHOD\n\t\x1B[1mFRAGMENT ANALYSIS\x1B[0m', '\n')
                break
            elif sq_md == '03':
                print('\n', 'SEQUENCING METHOD\n\t\x1B[1mNEXT-GENERATION SEQUENCING\x1B[0m', '\n')
                break
        else:
            print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')
    return sq_md

########### LOCATIONS ############ DONE
def loc():
    print('\x1B[1mPROVINCE\x1B[0m-only \x1B[1m38\x1B[0m registered province in \x1B[1mIndonesia\x1B[0m | '
          '\x1B[1mREGION\x1B[0m-ex: \x1B[1mUTARA\x1B[0m; \x1B[1mSELATAN\x1B[0m; \x1B[1mTIMUR\x1B[0m; \x1B[1mBARAT\x1B[0m')
    while True:
        try:
            province, region = input('ENTER THE SAMPLING LOCATION \x1B[1m(PROVINCE-REGION)\x1B[0m: ').split('-')
        except ValueError:
            print('ENTER THE PROVINCE AND REGION SEPARATED BY A HYPHEN'.center(127, '-'), '\n')
            continue
        else:
            if province.isupper() and region.isupper():
                prov_list = [
                            'ACEH', 'BALI', 'BANTEN', 'BENGKULU', 'DAERAH ISTIMEWA YOGYAKARTA', 'DIY', 'YOGYAKARTA', 'JAKARTA', 'GORONTALO',
                            'JAMBI', 'JAWA BARAT', 'JAWA TENGAH', 'JAWA TIMUR', 'KALIMANTAN BARAT', 'KALIMANTAN SELATAN', 'MALUKU UTARA',
                            'KALIMANTAN TENGAH', 'KALIMANTAN TIMUR', 'KALIMANTAN UTARA', 'BANGKA BELITUNG', 'KEPULAUAN RIAU', 'LAMPUNG', 'MALUKU',
                            'NUSA TENGGARA BARAT', 'NUSA TENGGARA TIMUR', 'PAPUA', 'PAPUA BARAT', 'PAPUA BARAT DAYA', 'PAPUA PEGUNUNGAN',
                            'PAPUA SELATAN', 'PAPUA TENGAH', 'RIAU', 'SULAWESI BARAT', 'SULAWESI SELATAN', 'SULAWESI TENGAH', 'SULAWESI TENGGARA',
                            'SULAWESI UTARA', 'SUMATERA BARAT', 'SUMATERA SELATAN', 'SUMATERA UTARA', 'SUMATERA', 'JAWA', 'NUSA TENGGARA', 'PAPUA', 'SULAWESI',
                            'KALIMANTAN', 'MALUKU'
                            ]
                reg_list = ['UTARA', 'TIMUR', 'BARAT', 'SELATAN', 'TENGAH', 'TENGGARA', 'BARAT DAYA', 'BARAT LAUT', 'TIMUR LAUT']
                if province in prov_list and region in reg_list:
                    print('\n', f'SAMPLING LOCATION\n\t\x1B[1m{province}\x1B[0m-\x1B[1m{region}\x1B[0m', '\n')
                    return province, region
                else:
                    print('\x1B[1mPROVINCE\x1B[0m OR \x1B[1mREGION\x1B[0m NOT REGISTERED OR FORMATTED INCORRECTLY'.center(143, '-'), '\n')
            else:
                print('PROVINCE AND REGION CANNOT BE EMPTY AND MUST BE IN UPPERCASE LETTERS'.center(127, '-'), '\n')

########### DATES ################ DONE
def date():
    print('Submission date \x1B[1mDD/MM/YYYY\x1B[0m')
    while True:
        try:
            day, month, year = input('ENTER THE SUBMISSION DATE: ').split('/')
        except ValueError:
            print('ENTER THE DAY, MONTH, AND YEAR SEPARATED BY A SLASH'.center(127, '-'), '\n')
            continue
        else:
            if day.isnumeric() and month.isnumeric() and year.isnumeric():
                if 1<=int(day)<=31 and 1<=int(month)<=12 and len(day and month) == 2 and len(year) == 4:
                    print('\n', f'SUBMISSION DATE\n\t\x1B[1m{day[0:2]}/{month[0:2]}/{year[0:4]}\x1B[0m', '\n')
                    return day, month, year
                else:
                    print('\x1B[1mDAY\x1B[0m/\x1B[1mMONTH\x1B[0m/\x1B[1mYEAR\x1B[0m OUT OF RANGE OR FORMATTED INCORRECTLY'.center(151, '-'), '\n')
            else:
                print('SUBMISSION DATE CANNOT BE EMPTY AND MUST BE NUMERIC'.center(127, '-'), '\n')

########### SPECIES ############## DONE
def species_id():
    print('\x1B[1mGenus\x1B[0m is capitalized-ex: \x1B[3mOryza\x1B[0m and \x1B[1mspecies\x1B[0m is lowercase-ex: \x1B[3msativa\x1B[0m')
    while True:
        try:
            genus, species = input('ENTER THE GENUS AND SPECIES SEPARATED BY A SPACE: ').split(' ')
        except ValueError:
            print('ENTER THE GENUS AND SPECIES SEPARATED BY A SPACE'.center(127, '-'), '\n')
            continue
        else:
            if genus.isalpha() and species.isalpha():
                if len(genus) == 1 and genus[0].isupper() and species.islower():
                    print('\n', f'SPECIES NAME\n\t\x1B[1m\x1B[3m{genus}\x1B[0m\x1B[0m \x1B[1m\x1B[3m{species}\x1B[0m\x1B[0m', '\n')
                    return genus, species
                else:
                    if len(genus) <= 2 and genus[0].isupper() and genus[-1].islower() and species.islower():
                        print('\n', f'SPECIES NAME\n\t\x1B[1m\x1B[3m{genus}\x1B[0m\x1B[0m \x1B[1m\x1B[3m{species}\x1B[0m\x1B[0m', '\n')
                        return genus, species
                    else:
                        if genus[0].isupper() and genus[1:-1].islower() and genus[-1].islower() and species.islower():
                            print('\n', f'SPECIES NAME\n\t\x1B[1m\x1B[3m{genus}\x1B[0m\x1B[0m \x1B[1m\x1B[3m{species}\x1B[0m\x1B[0m', '\n')
                            return genus, species
                        else:
                            print('\x1B[1m\x1B[3mGenus\x1B[0m\x1B[0m MUST BE CAPITALIZED AND \x1B[1m\x1B[3mspecies\x1B[0m\x1B[0m MUST BE IN LOWERCASE'.center(159, '-'), '\n')
            else:
                print('GENUS AND SPECIES CANNOT BE EMPTY AND MUST BE LETTERS'.center(127, '-'), '\n')

########### | SAVING DATA | ####################################################################### DONE
########### SAVING ############### DONE
def saving_data(metadata, s_date, s_location, sq_md, species_name, dna_seq, p_pi, nprop, amino_acid):
    gen_db_temp = {'ID.REF':metadata,'DATE':s_date,'LOCATION':s_location,'SEQ.MD':sq_md,'SPECIES':species_name,'DNA SEQUENCE':dna_seq,'%P/Pi':p_pi,'PROPERTIES':nprop,'AMINO ACID':amino_acid}
    gen_db.append(gen_db_temp)

########### CONFIRMATION ######### DONE
def conf(genus, species, day, month, year, province, region, sq_md, dna_seq):
    while True:
        try:
            meta_conf = int(input('|1| SAVE SUBMISSION|>>  |2| EDIT SUBMISSION|>>'.center(62)))
        except ValueError:
            print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'))
            continue
        else:
            if meta_conf == 1:
                s_date = f'{day[0:2]}/{month[0:2]}/{year[0:4]}'
                s_location = f'{province}-{region}'
                species_name = f'{genus} {species}'
                (p_pi) = ppi(dna_seq)
                (nprop) = prop(dna_seq)
                (amino_acid) = AA(dna_seq)
                metadata = f'#{genus[0]}{genus[-1]}{species[-1]}{species[0]}{day[0:2]}{month[0:2]}{year[0:4]}{province[0]}{province[-1]}{region[0]}{region[-1]}{sq_md}'
                while True:
                    d_metadata = [i.setdefault('ID.REF') for i in gen_db]
                    d_dna_seq = [i.setdefault('DNA SEQUENCE') for i in gen_db]
                    if metadata in d_metadata and dna_seq in d_dna_seq:
                        print(f'IDENTICAL RECORD IS FOUND [5-{dna_seq}-3]--{metadata}'.center(127, '-'))
                        print(' | VALIDATE OR CONTACT US THROUGH THE SYSTEM/PERSONALLY, THANK YOU | '.center(127, '|'))
                        break
                    else:
                        saving_data(metadata, s_date, s_location, sq_md, species_name, dna_seq, p_pi, nprop, amino_acid)
                        print(' | SUBMISSION COMPLETED | '.center(127, '|'))
                        print(f' | {metadata} | '.center(127, '|'))
                        print(' | SAVE GENERATED ID.REF FOR FURTHER AND FUTURE ACCESS | '.center(127, '|'))
                        break
                break
            elif meta_conf == 2:
                try:
                    meta_conf_edit = int(input('|1| SPECIES|>>  |2| DATE|>>  |3| LOCATION|>>  |4| SEQUENCING METHOD|>>  |5| DNA SEQUENCE|>>'.center(62)))
                except ValueError:
                    print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'))
                    continue
                if meta_conf_edit == 1:
                    newGenus, newSpecies = species_id()
                    ui_flow()
                    genus, species = newGenus, newSpecies
                elif meta_conf_edit == 2:
                    newDay, newMonth, newYear = date()
                    ui_flow()
                    day, month, year = newDay, newMonth, newYear         
                elif meta_conf_edit == 3:
                    newProvince, newRegion = loc()
                    ui_flow()
                    province, region = newProvince, newRegion
                elif meta_conf_edit == 4:
                    newSq_Md = seq_md()
                    ui_flow()
                    sq_md = newSq_Md
                elif meta_conf_edit == 5:
                    newGen_Seq = gen_seq()
                    ui_flow()
                    dna_seq = newGen_Seq
                else:
                    print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'))
            else:
                print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'))

########### METADATA ############# DONE
def meta_data(genus, species, day, month, year, province, region, sq_md, dna_seq):
    while True:
        frame = '-'
        print(frame.center(127, '-'))
        print(' | SUBMISSION METADATA | '.center(127, '|'))
        print(frame.center(127, '-'))
        print('{:<79} {:>79}'.format(f'\x1B[1m\x1B[3m{genus}\x1B[0m\x1B[0m \x1B[1m\x1B[3m{species}\x1B[0m\x1B[0m', 'SPECIES'))
        print(frame.center(127, '-'))
        print('{:<67} {:>67}'.format(f'\x1B[1m{day[0:2]}/{month[0:2]}/{year[0:4]}\x1B[0m', 'DATE'))
        print(frame.center(127, '-'))
        print('{:<71} {:>71}'.format(f'\x1B[1m{province}\x1B[0m-\x1B[1m{region}\x1B[0m', 'LOCATION'))
        print(frame.center(127, '-'))
        print('{:71} {:>71}'.format(f'\x1B[1m\x1B[1m5-{dna_seq}-3\x1B[0m\x1B[0m', 'DNA SEQUENCE'))
        print(frame.center(127, '-'))
        if sq_md == '01':
            print('{:<67} {:>67}'.format('\x1B[1mSANGER SEQUENCING\x1B[0m', 'SEQUENCING METHOD'))
            break
        elif sq_md == '02':
            print('{:<67} {:>67}'.format('\x1B[1mFRAGMENT ANALYSIS\x1B[0m', 'SEQUENCING METHOD'))
            break
        elif sq_md == '03':
            print('{:<67} {:>67}'.format('\x1B[1mNEXT-GENERATION SEQUENCING\x1B[0m', 'SEQUENCING METHOD'))
            break
        else:
            break
    print(frame.center(127, '-'))
    conf(genus, species, day, month, year, province, region, sq_md, dna_seq)

########### | UTILITY | ########################################################################### DONE
########### UI_FLOW ############## DONE
def ui_flow():
    while True:
        try:
            ui = int(input('|1| CONTINUE|\u2713  |2| ABORT|\u2717'.center(62)))
        except ValueError:
            print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'))
            continue
        else:
            if ui == 1:
                break
            elif ui == 2:
                main_page()
            else:
                print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'))

########### VSQ SELECTED RECORD ## DONE
def vsq_md(id_search):
    for id_ref in gen_db:
        if id_search == id_ref['ID.REF']:
            if id_ref['SEQ.MD'] == '01':
                print('{:<67} {:>59}'.format('SANGER SEQUENCING', 'SEQUENCING METHOD'))
                break
            elif id_ref['SEQ.MD'] == '02':
                print('{:<67} {:>59}'.format('FRAGMENT ANALYSIS', 'SEQUENCING METHOD'))
                break
            elif id_ref['SEQ.MD'] == '03':
                print('{:<67} {:>59}'.format('NEXT-GENERATION SEQUENCING', 'SEQUENCING METHOD'))
                break
            else:
                break

########### VIEW SELECTED RECORD # DONE
def ViewSR(id_search):
    frame = '-'
    for id_ref in gen_db:
        if id_search == id_ref['ID.REF']:
            print(frame.center(127, '-'))
            print(f' | RECORD OF ID.REF{id_search} | '.center(127, '|'))
            print(frame.center(127, '-'))
            print('{:<62} {:>72}'.format(f'\x1B[3m{id_ref['SPECIES']}\x1B[0m', 'SPECIES'))
            print(frame.center(127, '-'))
            print('{:<62} {:>64}'.format(id_ref['DATE'], 'DATE'))
            print(frame.center(127, '-'))
            print('{:<62} {:>64}'.format(id_ref['LOCATION'], 'LOCATION'))
            print(frame.center(127, '-'))
            vsq_md(id_search)
            print(frame.center(127, '-'))
            print('{:62} {:>64}'.format('5-'+id_ref['DNA SEQUENCE']+'-3', 'DNA SEQUENCE'))
            print(frame.center(127, '-'))
            print('{:<62} {:>64}'.format(id_ref['%P/Pi'], '%PURIN / %PIRIMIDIN'))
            print(frame.center(127, '-'))
            print('|||| DNA SEQUENCE PROPERTIES ||||'.center(127, '|'))
            print(frame.center(127, '-'))
            print('{:<62} {:>99}'.format(id_ref['PROPERTIES'], ''))
            print(frame.center(127, '-'))
            print('{:<62} {:>63}'.format(id_ref['AMINO ACID'], 'AMINO ACID'))
            print(frame.center(127, '-'))
            print('|'.center(127, '|'))
            print(frame.center(127, '-'))
            break
    else:
        print(f'RECORD NON-EXISTENT OR NOT REGISTERED FOR ID.REF{id_search}'.center(127, '-'), '\n')
    return id_search

########### ID SEARCH ############ DONE
def d_search():
    while True:
        id_search = input('\nENTER REGISTERED \x1B[1mID.REF\x1B[0m: ')
        ui_flow()
        if len(id_search) == 19 and id_search[0] == '#':
            if id_search[1:5].isalpha() and id_search[1].isupper() and id_search[2:5].islower() and id_search[13:17].isupper():
                if 1<=int(id_search[5:7])<=31 and 1<=int(id_search[7:9])<=12 and len(id_search[5:7] and id_search[7:9]) == 2 and len(id_search[9:13]) == 4 and id_search[17:19].isnumeric():
                    print(f'FETCHING RECORD FOR ID.REF{id_search}'.center(127, '-'), '\n')
                    break
                else:
                    print(f'ID.REF{id_search} IS INVALID OR FORMATTED INCORRECTLY'.center(127, '-'))
            else:
                print(f'ID.REF{id_search} IS INVALID OR FORMATTED INCORRECTLY'.center(127, '-'))
        else:
            print(f'ID.REF{id_search} IS INVALID OR FORMATTED INCORRECTLY'.center(127, '-'))
    return id_search

########### VIEW ENTIRE RECORD ### DONE
def ViewER():
    iterable_gen_db = []
    iterated_gen_db = [i.setdefault('ID.REF') for i in gen_db]
    for i in iterated_gen_db:
        if i != '':
            for keys in ['ID.REF', 'DATE', 'LOCATION', 'SEQ.MD', 'SPECIES', 'DNA SEQUENCE', '%P/Pi', 'PROPERTIES', 'AMINO ACID']:
                keys_values = [i.setdefault(keys) for i in gen_db]
                iterable_gen_db.append(keys_values)
    else:
        gen_db_null = []
        for keys in range(9):
            keys_values = [i for i in gen_db_null]
            iterable_gen_db.append(keys_values)
    iterable_gen_db = zip(iterable_gen_db[0], iterable_gen_db[1], iterable_gen_db[2], iterable_gen_db[3], iterable_gen_db[4], iterable_gen_db[5], iterable_gen_db[6], iterable_gen_db[7], iterable_gen_db[8])
    iterable_gen_db = list(iterable_gen_db)
    iterable_gen_db.sort()
    if len(iterable_gen_db) != 0:
        iterable_gen_db.pop(0)
    k1, k2, k3, k4, k5, k6, k7, k8, k9 = ('ID.REF', 'DATE', 'LOCATION', 'SEQ.MD', 'SPECIES', 'DNA SEQUENCE', '%P/Pi', 'PROPERTIES', 'AMINO ACID')
    frame = '-'
    print(frame.center(175, '-'))
    print('|'+k1.center(21, ' ')+'|'+k5.center(29, ' ')+'|'+k6.center(39, ' ')+'|'+k7.center(17, ' ')+'|'+k9.center(63, ' ')+'|')
    print(frame.center(175, '-'))
    for it in iterable_gen_db:
        v1, v2, v3, v4, v5, v6, v7, v8, v9 = (str(it[0]), it[1], it[2], it[3], f'\x1B[3m{it[4]}\x1B[0m', f'5-{it[5]}-3', str(it[6]), it[7], str(it[8]))
        frame = '-'
        print('|'+v1.center(21, ' ')+'|'+v5.center(37, ' ')+'|'+v6.center(39, ' ')+'|'+v7.center(17, ' ')+'|'+v9.center(63, ' ')+'|')
        print(frame.center(175, '-'))

########### | EXIT_DB | ########################################################################### DONE
########### SAVE-EXIT ############ DONE
def save_exit():
    with open(f_name, 'wb') as f:
        pickle.dump(db, f)
db = gen_db

########### EXIT ################# DONE
def exit_db():
    print('\n',' | EXIT DATABASE | '.center(125, '|'))
    print('\n','THANK YOU FOR YOUR CONTRIBUTION TO THE NATIONAL GENETIC DATABASE SYSTEM OF INDONESIA'.center(125, '-'), '\n')
    while True:
        try:
            exitinp = int(input('\x1B[1mBE ADVISED\x1B[0m ANY CHANGES MADE TO THE \x1B[1mRECORD\x1B[0m WILL BE SAVED BEFORE EXIT\n|1| SAVE AND EXIT|>>  |2| ABORT THE OPERATION|>>'))
        except ValueError:
            print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'), '\n')
            continue
        else:
            if exitinp == 1:
                atexit.register(save_exit)
                print(' | CHANGES TO THE RECORD HAVE BEEN SUCCESSFULLY SAVED TO THE DATABASE | '.center(127, '|'), '\n')
                exit()
            elif exitinp == 2:
                print(' | EXIT OPERATION ABORTED | '.center(127, '|'))
                main_page()
            else:
                print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')

########### | DELETE_REC | ######################################################################## DONE
def withdraw_rec():
    print('\n',' | WITHDRAW RECORD | '.center(125, '|'))
    print('\n','WITHDRAW THE RECORD ON THE GENETIC DATABASE'.center(125, '-'), '\n')
    ViewER()
    while True:
        ui_flow()
        wdw = ViewSR(d_search())
        for id_ref in gen_db:
            if wdw == id_ref['ID.REF']:
                try:
                    wdw_conf = int(input('\x1B[1mBE ADVISED\x1B[0m RECORD WITHDRAWAL IS \x1B[1mPERMANENT\x1B[0m\n|1| WITHDRAW|>>  |2| ABORT THE OPERATION|>>'))   
                except ValueError:
                    print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'), '\n')
                    continue
                else:
                    if wdw_conf == 1:
                        gen_db.remove(id_ref)
                        print(f'THE RECORD FOR ID.REF{wdw} HAS BEEN SUCCESSFULLY WITHDRAWN FROM DATABASE'.center(127, '-'), '\n')
                        ViewER()
                    elif wdw_conf == 2:
                        print(f'WITHDRAWAL PROCESS FOR ID.REF{wdw} HAS BEEN CANCELLED'.center(127, '-'), '\n')
                    else:
                        print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')

########### | REVIEW_REC | ######################################################################## DONE
########### RFLOW ################ DONE
def rflow():
    while True:
        try:
            rflow = int(input('|1| CONTINUE REVIEW|>>  |2| SAVE REVIEW|>>'.center(62)))   
        except ValueError:
            print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'))
            continue
        else:
            if rflow == 1:
                break
            elif rflow == 2:
                main_page()
            else:
                print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'))

########### REVIEW ############### DONE
def review_rec():
    print('\n',' | REVIEW RECORDS | '.center(125, '|'))
    print('\n','REVIEW THE RECORD ON THE GENETIC DATABASE'.center(125, '-'), '\n')
    ViewER()
    while True:
        ui_flow()
        rvw = ViewSR(d_search())
        for id_ref in gen_db:
            if rvw == id_ref['ID.REF']:
                try:
                    rvwd = int(input('\x1B[1mBE ADVISED\x1B[0m A SUBSTITUTE \x1B[1mID.REF\x1B[0m WILL BE GENERATED FOR THE UPDATED RECORD\n|1| REVIEW|>>  |2| ABORT THE OPERATION|>>'))
                except ValueError:
                    print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'), '\n')
                    continue
                else:
                    if rvwd == 1:
                        try:
                            rvw_input = int(input('|1| SPECIES|>>  |2| DATE|>>  |3| LOCATION|>>  |4| SEQUENCING METHOD|>>  |5| DNA SEQUENCE|>>'.center(62)))
                        except ValueError:
                            print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'), '\n')
                            continue
                        else:
                            for v in gen_db:
                                if rvw == v['ID.REF']:
                                    if rvw_input == 1:
                                        revGenus, revSpecies = species_id()
                                        genus, species = revGenus, revSpecies
                                        species_name = f'{genus} {species}'
                                        v['SPECIES'] = species_name
                                        update = list(rvw)
                                        update[1] = genus[0]
                                        update[2] = genus[-1]
                                        update[3] = species[-1]
                                        update[4] = species[0]
                                        meta_rev = ''.join(update)
                                        v['ID.REF'] = f'{meta_rev}'
                                        print(f'A REPLACEMENT ID.REF{meta_rev} HAS BEEN GENERATED FOR ID.REF{rvw}'.center(127, '-'), '\n')
                                        input('PRESS ENTER TO CONTINUE')
                                        print(f'SPECIES NAME FOR ID.REF{rvw} HAS BEEN SUCCESSFULLY UPDATED'.center(127, '-'), '\n')
                                        ViewER()
                                        print()
                                        ViewSR(meta_rev)
                                        rflow()
                                    elif rvw_input == 2:
                                        revDay, revMonth, revYear = date()
                                        day, month, year = revDay, revMonth, revYear
                                        s_date = f'{day[0:2]}/{month[0:2]}/{year[0:4]}'
                                        v['DATE'] = s_date                                        
                                        update = list(rvw)
                                        update[5:7] = day[0:2]
                                        update[7:9] = month[0:2]
                                        update[9:13] = year[0:4]
                                        meta_rev = ''.join(update)
                                        v['ID.REF'] = f'{meta_rev}'
                                        print(f'A REPLACEMENT ID.REF{meta_rev} HAS BEEN GENERATED FOR ID.REF{rvw}'.center(127, '-'), '\n')
                                        input('PRESS ENTER TO CONTINUE')
                                        print(f'DATE FOR ID.REF{rvw} HAS BEEN SUCCESSFULLY UPDATED'.center(127, '-'), '\n')
                                        ViewER()
                                        print()
                                        ViewSR(meta_rev)
                                        rflow()
                                    elif rvw_input == 3:
                                        revProvince, revRegion = loc()
                                        province, region = revProvince, revRegion
                                        s_location = f'{province}-{region}'
                                        v['LOCATION'] = s_location
                                        update = list(rvw)
                                        update[13] = province[0]
                                        update[14] = province[-1]
                                        update[15] = region[0]
                                        update[16] = region[-1]
                                        meta_rev = ''.join(update)
                                        v['ID.REF'] = f'{meta_rev}'
                                        print(f'A REPLACEMENT ID.REF{meta_rev} HAS BEEN GENERATED FOR ID.REF{rvw}'.center(127, '-'), '\n')
                                        input('PRESS ENTER TO CONTINUE')
                                        print(f'LOCATION FOR ID.REF{rvw} HAS BEEN SUCCESSFULLY UPDATED'.center(127, '-'), '\n')
                                        ViewER()
                                        print()
                                        ViewSR(meta_rev)
                                        rflow()
                                    elif rvw_input == 4:
                                        revSq_Md = seq_md()
                                        sq_md = revSq_Md
                                        v['SEQ.MD'] = sq_md
                                        update = list(rvw)
                                        update[17:19] = sq_md[0:2]
                                        meta_rev = ''.join(update)
                                        v['ID.REF'] = f'{meta_rev}'
                                        print(f'A REPLACEMENT ID.REF{meta_rev} HAS BEEN GENERATED FOR ID.REF{rvw}'.center(127, '-'), '\n')
                                        input('PRESS ENTER TO CONTINUE')
                                        print(f'SEQUENCING METHOD FOR ID.REF{rvw} HAS BEEN SUCCESSFULLY UPDATED'.center(127, '-'), '\n')
                                        ViewER()
                                        print()
                                        ViewSR(meta_rev)
                                        rflow()
                                    elif rvw_input == 5:
                                        revDNA_Seq = gen_seq()
                                        dna_seq = revDNA_Seq
                                        v['DNA SEQUENCE'] = dna_seq
                                        p_pi = ppi(dna_seq)
                                        v['%P/Pi'] = p_pi
                                        nprop = prop(dna_seq)
                                        v['PROPERTIES'] = nprop
                                        amino_acid = AA(dna_seq)
                                        v['AMINO ACID'] = amino_acid
                                        print(f'DNA SEQUENCE, %P/Pi, PROPERTIES, AND AMINO ACID FOR ID.REF{rvw} HAS BEEN SUCCESSFULLY UPDATED'.center(127, '-'), '\n')
                                        input('PRESS ENTER TO CONTINUE')
                                        ViewER()
                                        print()
                                        ViewSR(rvw)
                                        rflow()
                                    else:
                                        print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')
                    elif rvw == 2:
                        print(f'REVIEWING PROCESS FOR ID.REF{rvw} HAS BEEN CANCELLED'.center(127, '-'), '\n')
                    else:
                        print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')
            else:
                print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')

########### | FETCH_REC | ######################################################################### DONE
def fetch_rec():
    print('\n',' | FETCH RECORDS | '.center(125, '|'))
    print('\n','VIEW THE RECORDS ON THE GENETIC DATABASE'.center(125, '-'), '\n')
    while True:
        ui_flow()
        try:
            fetch_rec = int(input('|1| VIEW ENTIRE RECORD|>>  |2| VIEW SELECTED RECORD|>>'.center(62)))
        except ValueError:
            print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'), '\n')
            continue
        else:
            if fetch_rec == 1:
                print('\n', 'ENTIRE RECORD OF GENOMIC SEQUENCES IN THE DATABASE'.center(125, '-'))
                ViewER()
            elif fetch_rec == 2:
                print('\n', 'SELECTED GENOMIC SEQUENCE RECORD IN THE DATABASE'.center(125, '-'))
                ViewSR(d_search())
            else:
                print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')

########### | SUBMIT_REC | ######################################################################## DONE
def submit_rec():
    print('\n',' | SUBMIT RECORD | '.center(125, '|'))
    print('\n','START THE SUBMISSION PROCESS BY FILLING OUT THE FOLLOWING INFORMATION ABOUT THE SAMPLE'.center(125, '-'), '\n')
    (genus, species) = species_id()
    (day, month, year) = date()
    (province, region) = loc()
    (sq_md) = seq_md()
    (dna_seq) = gen_seq()
    meta_data(genus, species, day, month, year, province, region, sq_md, dna_seq)

########### | MAIN_PAGE | ######################################################################### DONE
def main_page():
    while True:
        main = ' WELCOME TO THE NATIONAL GENETIC DATABASE SYSTEM OF INDONESIA '
        middle = ' | A STRONG COUNTRY PROTECTS ITS RICHES | '
        under = '='
        print('\n', main.center(125, '#'))
        print(middle.center(127, '\u2605'), '\n')
        print('+++---------------+++----------------+++----------------+++-----------------+++---------------+'.center(127))
        print('|1| SUBMIT RECORD |2| FETCH RECORDS  |3| REVIEW RECORDS |4| WITHDRAW RECORD |5| EXIT DATABASE |'.center(127))
        print('+++---------------+++----------------+++----------------+++-----------------+++---------------+'.center(127), '\n')
        print(under.center(127, '='))
        try:
            main_input = int(input('ENTER THE FOLLOWING OPTION TO PROCEED |1|2|3|4|5|: '))
        except ValueError:
            print(' | INPUT IS NOT RECOGNIZED | '.center(127, '|'))
            continue
        else:
            if main_input == 1:
                ui_flow()
                submit_rec()
            elif main_input == 2:
                ui_flow()
                fetch_rec()
            elif main_input == 3:
                ui_flow()
                review_rec()
            elif main_input == 4:
                ui_flow()
                withdraw_rec()
            elif main_input == 5:
                ui_flow()
                exit_db()
            else:
                print('INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'))

########### | LOGIN | ############################################################################# DONE
res_db = {'record_id':[],
          'record_id_pwd':[]}

def opt():
    while True:
        try:
            opt = int(input('|1| NEXT|>>  |2| BACK|>>'.center(62)))
        except ValueError:
            print('\n', ' | INPUT IS NOT RECOGNIZED | '.center(127, '|'), '\n')
            continue
        else:
            if opt == 1:
                break
            elif opt == 2:
                ui_page()
            else:
                print('\n', 'INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'), '\n')

def id_card_number():
    while True:
        id_pwd = input(' ID CARD NUMBER (16-digit): ')
        if len(id_pwd) == 16 and id_pwd.isnumeric():
            break
        else:
            print('\n', 'ID CARD NUMBER CANNOT BE EMPTY OR FORMATTED INCORRECTLY'.center(127, '-'), '\n')
    return id_pwd

def researcher_id():
    while True:
        id = input(' RESEARCHER ID: ')
        specials = ["'"]
        if len(id) != 0 and id.isupper():
            if any(x in specials for x in id) or all(x.isalpha() or x.isspace() for x in id):
                break
            else:
                print('\n', 'RESEARCHER ID MUST BE IDENTICAL TO ID CARD'.center(127, '-'), '\n')
        else:
            print('\n', 'RESEARCHER ID CANNOT BE EMPTY AND ONLY USE UPPERCASE LETTERS WITH SPACES'.center(127, '-'), '\n')
    return id

def REGISTER():
    print('\n', ' | REGISTER | '.center(127, '|'))
    print('\n', 'WELCOME, PLEASE ENTER LEGITIMATE RESEARCH ID AND ID CARD NUMBER'.center(127, '-'), '\n')
    opt()
    reg_id = researcher_id()
    rec_id = reg_id.encode()
    rec_id_1 = hashlib.md5(rec_id).hexdigest() 
    reg_id_pwd = id_card_number()
    val_id = reg_id_pwd.encode()
    val_id_1 = hashlib.md5(val_id).hexdigest()
    opt()
    while True:
        if any(val_id_1 in val for val in res_db.values()) and not any(rec_id_1 in val for val in res_db.values()):
            print('\n', 'ACCOUNT ALREADY REGISTERED WITH A DIFFERENT RESEARCH ID'.center(127, '-'), '\n')
            opt()
        elif any(rec_id_1 and val_id_1 in val for val in res_db.values()):
            print('\n', f' | {reg_id}, YOUR ACCOUNT IS ALREADY REGISTERED IN OUR DATABASE, PLEASE SIGN IN | '.center(127, '|'))
            return SIGNIN()
        elif any(rec_id_1 in val for val in res_db.values()) and not any(val_id_1 in val for val in res_db.values()):
            print('\n', 'ACCOUNT ALREADY REGISTERED WITH A DIFFERENT ID CARD NUMBER'.center(127, '-'), '\n')
            opt()
        else:
            while True:
                conf_pwd = input(' CONFIRM ID CARD NUMBER (16-digit): ')
                if conf_pwd == reg_id_pwd:
                    rec = conf_pwd.encode()
                    rec_1 = hashlib.md5(rec).hexdigest()
                    res_db.setdefault('record_id',[]).append(rec_id_1)
                    res_db.setdefault('record_id_pwd',[]).append(rec_1)
                    print('\n', ' | REGISTRATION COMPLETE | '.center(127, '|'))
                    #print(res_db)
                    break
                else:
                    print('\n', ' | RE-ENTER ID CARD NUMBER | '.center(127, '|'), '\n')
            return conf_pwd

def SIGNIN():
    print('\n', ' | SIGN IN | '.center(127, '|'))
    print('\n', 'WELCOME, PLEASE ENTER REGISTERED RESEARCH ID AND ID CARD NUMBER'.center(127, '-'), '\n')
    opt()
    sign_id = researcher_id()
    auth_id = sign_id.encode()
    auth_id_1 = hashlib.md5(auth_id).hexdigest()
    sign_id_pwd = id_card_number()
    auth = sign_id_pwd.encode()
    auth_1 = hashlib.md5(auth).hexdigest()
    opt()
    while True:
        if any(auth_id_1 in val for val in res_db.values()) and not any(auth_1 in val for val in res_db.values()):
            print('\n', 'VERIFY RESEARCH ID OR ID CARD NUMBER'.center(127, '-'))
            return SIGNIN()
        elif any(auth_1 in val for val in res_db.values()) and not any(auth_id_1 in val for val in res_db.values()):
            print('\n', 'VERIFY RESEARCH ID OR ID CARD NUMBER'.center(127, '-'))
            return SIGNIN()
        elif any(auth_id_1 and auth_1 in val for val in res_db.values()):
            print('\n', f' | WELCOME BACK {sign_id}, THANK YOU FOR YOUR CONTRIBUTION | '.center(127, '|'))
            main_page()
        else:
            print('\n', ' | ACCOUNT NON-EXISTENT OR NOT REGISTERED | '.center(127, '|'))
            print('\n', ' | REGISTER TO START CONTRIBUTING | '.center(127, '|'))
        return REGISTER()

def ui_page():
    while True:
        ui = ' WELCOME TO THE NATIONAL GENETIC DATABASE SYSTEM OF INDONESIA '
        middle = ' | PLEASE SIGN IN OR REGISTER TO ACCESS | '
        under = '='
        print('\n', ui.center(127, '#'),'\n')
        print('+--------------------------------------+'.center(127, ' '))
        print(middle.center(127, ' '))
        print('+--------------------------------------+'.center(127, ' '))
        print('\t''\t''\t''\t''\t','   |1| SIGN IN  |2| REGISTER  |3| SIGN OUT|')
        print('+--------------------------------------+'.center(127, ' '))
        print('\n', under.center(127, '='))
        try:
            entry = int(input(' ENTER THE FOLLOWING OPTION TO PROCEED |1|2|3|: '))
        except ValueError:
            print('\n', ' | INPUT IS NOT RECOGNIZED | '.center(127, '|'))
            continue
        else:
            if entry == 1:
                SIGNIN()
            elif entry == 2:
                REGISTER()
            elif entry == 3:
                exit()
            else:
                print('\n', 'INVALID ENTRY OR OPTION NOT AVAILABLE'.center(127, '-'))

########### | USER-INTERFACE | #################################################################### DONE
main_page()
ui_page()