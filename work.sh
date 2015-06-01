cat mPing_Ping_Pong/mPing_Ping_Pong.fa otherTE/otherTEs.fa otherTE/otherTEs1.fa otherTE/otherTEs2.fa otherRETROs/nonRedundant.retro.fa otherDNATE/OS_NDTE_2010May11.lib.dusted > Rice.TE.fa

echo ""
formatdb -i Rice.TE.fa -p F
blastall -p blastn -i Rice.TE.fa -d Rice.TE.fa -o Rice.TE.self.blast -e 1e-5 &
perl /rhome/cjinfeng/software/bin/blast_parser.pl --tophit 2 --topmatch 1 Rice.TE.self.blast | awk '$1!=$5'| less -S


echo "make repeat short, fast for mapping reads"
python TrimTransposon.py --input Rice.TE.fa

echo "unique short: use unique as very few difference from representative"
formatdb -i Rice.TE.short.fa -p F
blastall -p blastn -i Rice.TE.short.fa -d Rice.TE.short.fa -o Rice.TE.short.self.blast -e 1e-5 &
perl /rhome/cjinfeng/software/bin/blast_parser.pl Rice.TE.short.self.blast > Rice.TE.short.self.blast.table
#3065, use one sequence to represent these 99% of all element
python TrimDuplicate.py --input Rice.TE.short.self.blast.table --fasta Rice.TE.short.fa > Rice.TE.short.unique.list
#3041, use one sequence to represent these 99% similarity in 100 of start and end
python TrimHomoglogous.py --input Rice.TE.short.self.blast.table --fasta Rice.TE.short.fa > Rice.TE.short.representative.list
