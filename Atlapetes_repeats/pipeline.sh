mamba activate repeats
BuildDatabase -name Atlapetes -engine ncbi /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/blancae.p_ctg.fa

RepeatModeler -pa 16 -engine ncbi -database Atlapetes 2>&1 | tee 00_repeatmodeler.log

cat Atlapetes-families.fa | seqkit fx2tab | awk '{ print "Ablancae_"$0 }' | seqkit tab2fx > Atlapetes-families.prefix.fa

cat Atlapetes-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > Atlapetes-families.prefix.fa.known
cat Atlapetes-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > Atlapetes-families.prefix.fa.unknown

# quantify number of classified elements
grep -c ">" Atlapetes-families.prefix.fa.known
# quantify number of unknown elements
grep -c ">" Atlapetes-families.prefix.fa.unknown

mkdir -p  01_simple_out 02_Aves_out 03_known_out 04_unknown_out

#ROUND 1: Annotate/mask the simple repeats
RepeatMasker -pa 16 -a -e ncbi -dir 01_simple_out -noint -xsmall /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/blancae.p_ctg.fa  2>&1 | tee logs/01_simplemask.log

#rename outputs
rename fa simple_mask 01_simple_out/blancae*
rename .masked .masked.fasta 01_simple_out/blancae*

#ROUND 2: Annotate/mask Aves elements sourced from Repbase using output from 1st round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 02_Aves_out -nolow \
        -species Aves 01_simple_out/blancae.p_ctg.simple_mask.masked.fasta 2>&1 | tee logs/02_Avesmask.log

#rename outputs
rename simple_mask.masked.fasta Aves_mask 02_Aves_out/blancae*
rename .masked .masked.fasta 02_Aves_out/blancae*

#ROUND 3: Annotate/mask known elements identified in repeatmodeler
RepeatMasker -pa 16 -a -e ncbi -dir 03_known_out -nolow \
        -lib Atlapetes-families.prefix.fa.known 02_Aves_out/blancae.p_ctg.Aves_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log

#rename outputs
rename Aves_mask.masked.fasta known_mask 03_known_out/blancae*
rename .masked .masked.fasta 03_known_out/blancae*

#ROUND 4: Annotate/mask unknown elements identified in repeatmodeler
RepeatMasker -pa 16 -a -e ncbi -dir 04_unknown_out -nolow \
        -lib Atlapetes-families.prefix.fa.unknown 03_known_out/blancae.p_ctg.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log

#rename outputs
rename known_mask.masked.fasta unknown_mask 04_unknown_out/blancae*
rename .masked .masked.fasta 04_unknown_out/blancae*

#Combine Masking Rounds:
mkdir -p 05_full_out

# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/blancae.p_ctg.simple_mask.cat.gz \
02_Aves_out/blancae.p_ctg.Aves_mask.cat.gz \
03_known_out/blancae.p_ctg.known_mask.cat.gz \
04_unknown_out/blancae.p_ctg.unknown_mask.cat.gz \
> 05_full_out/blancae.p_ctg.full_mask.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/blancae.p_ctg.simple_mask.out \
<(cat 02_Aves_out/blancae.p_ctg.Aves_mask.out | tail -n +4) \
<(cat 03_known_out/blancae.p_ctg.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/blancae.p_ctg.unknown_mask.out | tail -n +4) \
> 05_full_out/blancae.p_ctg.full_mask.out

# copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/blancae.p_ctg.simple_mask.out > 05_full_out/blancae.p_ctg.simple_mask.out

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat 02_Aves_out/blancae.p_ctg.Aves_mask.out \
<(cat 03_known_out/blancae.p_ctg.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/blancae.p_ctg.unknown_mask.out | tail -n +4) \
> 05_full_out/blancae.p_ctg.complex_mask.out

# combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_out/blancae.p_ctg.simple_mask.align \
02_Aves_out/blancae.p_ctg.Aves_mask.align \
03_known_out/blancae.p_ctg.known_mask.align \
04_unknown_out/blancae.p_ctg.unknown_mask.align \
> 05_full_out/blancae.p_ctg.full_mask.align

# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species Aves 05_full_out/blancae.p_ctg.full_mask.cat.gz 2>&1 | tee logs/05_fullmask.log
#outputs .tbl file

# use Daren's custom script to convert .out to .gff3 for all repeats, simple repeats only, and complex repeats only
rmOutToGFF3custom -o 05_full_out/blancae.p_ctg.full_mask.out > 05_full_out/blancae.p_ctg.full_mask.gff3
rmOutToGFF3custom -o 05_full_out/blancae.p_ctg.simple_mask.out > 05_full_out/blancae.p_ctg.simple_mask.gff3
rmOutToGFF3custom -o 05_full_out/blancae.p_ctg.complex_mask.out > 05_full_out/blancae.p_ctg.complex_mask.gff3

# create masked genome FASTA files
# create simple repeat soft-masked genome
bedtools maskfasta -soft -fi /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/blancae.p_ctg.fa -bed 05_full_out/blancae.p_ctg.simple_mask.gff3 \
        -fo 05_full_out/blancae.p_ctg.fasta.simple_mask.soft.fasta
# create complex repeat hard-masked genome
bedtools maskfasta -fi 05_full_out/blancae.p_ctg.fasta.simple_mask.soft.fasta \
        -bed 05_full_out/blancae.p_ctg.complex_mask.gff3 \
        -fo 05_full_out/blancae.p_ctg.simple_mask.soft.complex_mask.hard.fasta
