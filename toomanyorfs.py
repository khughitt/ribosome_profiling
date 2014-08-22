from Bio import SeqIO
from Bio.Data import CodonTable
# import time
input_file = "lmajor.fasta"
# record = SeqIO.read("saccharomyces_cerevisiae.fasta","fasta")
table = CodonTable.ambiguous_dna_by_id[1]
min_pro_len = 30

def find_orfs_with_trans(seq, trans_table, min_protein_length, chromosome):
    """
    find_orfs_with_trans:  Search a chromosome for a semi-exhaustive set of putative ORFs
    """
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        aa_start = 0
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = trans.find("M", aa_start)
            aa_end = 0
            while (aa_start < trans_len):
                aa_end = trans.find("*", aa_start)
                if (aa_end == -1):
                    aa_end = trans_len
                if (aa_end - aa_start >= min_protein_length):
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    ## $bed_entry = qq"$chr  $nt  $nt_end  ${id}1  0  $plusminus  $nt $nt_end 125,125,125\n";
                    if int(frame) == 0:
                        color = "255,0,0"
                    elif int(frame) == 1:
                        color = "0,255,0"
                    elif int(frame) == 2:
                        color = "0,0,255"
                    else:
                        color = "0,0,0"
                    if strand == 1:
                        pos = "Pos"
                        plusminus = "+"
                    else:
                        pos = "Neg"
                        plusminus = "-"
                    string = "%s\t%s\t%s\t%s%s\t0\t%s\t%s\t%s\t%s\n" % (chromosome, start, end, pos, frame, plusminus, start, end, color)

                    ## print string
                    ## time.sleep(1)
                    if (strand == 1):
                        if (start > 0 and end < seq_len):
                            f_bedfile.write(string)
                    else:
                        if (start > 0 and end < seq_len):
                            r_bedfile.write(string)
                aa_start = trans.find("M", aa_start + 1)
                if (aa_start == -1):
                    break
    return

f_bedfile = open("f_orfs.bed", "w")
r_bedfile = open("r_orfs.bed", "w")
for record in SeqIO.parse(input_file, "fasta"):
    find_orfs_with_trans(record.seq, table, min_pro_len, record.name)
# }py
f_bedfile.close()
r_bedfile.close()

