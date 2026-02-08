from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
print("all successful")

#2. Sequence Quality Analysis

record = SeqIO.read("E5KP27.fasta", "fasta")

print("The Sequence ID:", record.id)
print("The Sequence Description:", record.description)
print("Sequence:", record.seq)
print("Type of the seq record:", type(record))

print("The Length of the Amino acid:", len(record))

# Amino acid composition

sequence = str(record.seq)

aa_count = {}

for aa in sequence:
    aa_count[aa] = aa_count.get(aa, 0) + 1

print("Amino acid composition:")
for aa, count in sorted(aa_count.items()):
    print(f"{aa}: {count}")

    #3. Sequence Filtering & Validation

if len(record.seq) < 50:
    print("Too short sequence")
else:
    print("Sequence is fine")

#4. Homology Search (BLAST)


result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence= record.seq
)

with open("blast_result.xml", "w") as b:
    b.write(result_handle.read())


print("BLAST performed successfully. Results saved to blast_results.xml")

with open("blast_result.xml") as b:
    blast_record = NCBIXML.read(b)

top_hit = blast_record.alignments[0]
top_hsp = top_hit.hsps[0]

#closest homolog

print("Closest homolog:")
print("Protein:", top_hit.hit_def)
print("Length:", top_hit.length)
print("E-value:", top_hsp.expect)
print("Identity:", top_hsp.identities, "/", top_hsp.align_length)
print("Percent identity:", (top_hsp.identities / top_hsp.align_length) * 100)

#Identify Conserved Regions

print("Query sequence:")
print(top_hsp.query)

print("Match:")
print(top_hsp.match)

print("Subject sequence:")
print(top_hsp.sbjct)

#Step 5: Functional Annotation

from Bio import SeqIO
record = SeqIO.read("E5KP27.fasta", "fasta")
print(record.id)
print(record.description)
print(len(record.seq))
print(record.seq)

#Predict: Function, Biological role, Organism relevance

#Functional Prediction
#Sequence annotation, and homology searching indicate that the protein is likely an adenine DNA glycosylase, an essential enzyme in base excision repair (BER) pathway, which plays a key role in excising damaged or mispaired adenines from DNA.

#Biological Role
#The protein functions in maintaining homeostasis of genome, preventing mutagenesis as result of treatment by oxidative DNA damaging reagents. It fixes A:G mismatches and as a result decreases transversions mutations in DNA replication.

#Organismal Relevance
#In human cells, MUTYH is required to repair DNA and guard against mutagenesis. Its mutation is correlated with hereditary colorectal cancer and it has critical implications in cancer prevention and genome integrity.


#Step 6: Biological Interpretation (Research Outcome)

#1)What does this sequence likely do?
#The sequence encodes the enzyme Adenine DNA glycosylase which functions as a crucial component of humans base excision repair system for DNA restoration. The protein detects adenine bases that have formed incorrect pairs with guanine due to oxidative damage to DNA. When the system detects a mismatch it separates the mispaired adenine by breaking the N-glycosidic bond which creates an abasic site in the DNA structure. Other DNA repair enzymes execute processing of the created lesion to achieve restoration of the original DNA sequence. The enzyme maintains genetic integrity by fixing these mismatches. The activity serves as a critical function which safeguards genomic stability. The sequence functions as a DNA repair enzyme which protects human cells from mutagenic damage.

#2) Why do they think so?
#The functional prediction depends on sequence annotation and gene identity information which the FASTA header supplies. The sequence exists as an annotated Adenine DNA glycosylase which the MUTYH gene encodes in the human species. The protein length matches exactly the length of all known MUTYH proteins which exist in existing databases. The sequence exhibits a composition and structural pattern which matches the characteristics of proteins that mediate DNA binding and repair activities. The protein MUTYH serves as a research subject to study its ability to repair oxidative damage to DNA. The biological function of MUTYH remains unchanged throughout different species. The sequence is expected to demonstrate adenine DNA glycosylase activity based on its structure.

#3) What evidence supports their claim?
#BLAST analysis provides the strongest proof of this prediction because it shows high similarity with known human MUTYH proteins. The sequence ID and description make clear that the protein exists as Adenine DNA glycosylase from Homo sapiens. The gene name MUTYH is clearly specified in the annotation. The protein length and conserved reg