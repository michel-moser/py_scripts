'''Usage: python fasta_stats.py yourfile.fasta

provides standart metrics of the fasta file:
Filename
Total Length
Total Ns
Total Length without Ns
GC content

Number of sequences
Max sequence size
Min sequence size 
Mean sequence size

Number of sequences larger than 1 Mb
Number of sequences larger than 100 Kb

N10	L10 
N50 	L50
N90	L90
N95	L95

'''


from Bio import SeqIO
import sys

seq_file = sys.argv[1]
contigsLength = []
sumL = 0
nS = 0
gS = 0
cS = 0
mP = 0
kP = 0
nr_seq = 0

for seq_record in SeqIO.parse(open(seq_file), "fasta"):
    nr_seq += 1
    nS += seq_record.seq.count('N')
    gS += seq_record.seq.count('G')
    cS += seq_record.seq.count('C')
    sumL += len(seq_record.seq)
    contigsLength.append(len(seq_record.seq))
    if len(seq_record.seq) > 1000000:
        mP += 1
    if len(seq_record.seq) > 100000:
        kP += 1
print('\nFilename: '+seq_file+'\n')

print('\nTotal Length: '+'{0:,}'.format(sumL)+ ' bp')
print('Total Ns: '+'{0:,}'.format(nS)+ ' bp')
print('Total None Ns: '+'{0:,}'.format(sumL - nS)+ ' bp\n')

print('GC content: '+ str((float(gS) + float(cS))/(sumL-nS)))
print('Number of scaffolds: ' + '{0:,}'.format(nr_seq)+ '\n')

print('Max. scaffold size: '+ '{0:,}'.format(max(contigsLength))+ ' bp')
print('Min. scaffold size: '+ '{0:,}'.format(min(contigsLength))+ ' bp')

print('Average scaffold size: '+str(np.mean(contigsLength)))

print('Scaffolds > 1 Mb: '+ '{0:,}'.format(mP))
print('Scaffolds > 100 Kb: '+ '{0:,}'.format(kP)+ '\n')

p10_L = sumL * 0.1
p20_L = sumL * 0.2
p30_L = sumL * 0.3
p40_L = sumL * 0.4
p50_L = sumL * 0.5
p60_L = sumL * 0.6
p70_L = sumL * 0.7
p80_L = sumL * 0.8
p90_L = sumL * 0.9
p95_L = sumL * 0.95

contigsLength.sort()
contigsLength.reverse()    

#N10
testSum = 0
L10 = 0
count = 0
for con in contigsLength:
    testSum += con
    count += 1
    if p10_L < testSum:
        L10 = con
        break
N10 = count
del con, testSum, count

#N50
testSum = 0
L50 = 0
count = 0
for con in contigsLength:
    testSum += con
    count += 1
    if p50_L < testSum:
        L50 = con
        break
N50 = count
del con, testSum, count

#N90
testSum = 0
L90 = 0
count = 0
for con in contigsLength:
    testSum += con
    count += 1
    if p90_L < testSum:
        L90 = con
        break
N90 = count
del con, testSum, count

#N95
testSum = 0
L95 = 0
count = 0
for con in contigsLength:
    testSum += con
    count += 1
    if p95_L < testSum:
        L95 = con
        break
N95 = count
del con, testSum, count


print 'N10: ' + str(N10) +'\tL10: ' + '{0:,}'.format(L10)+ ' bp'
print 'N50: ' + str(N50) +'\tL50: ' + '{0:,}'.format(L50)+ ' bp'
print 'N90: ' + str(N90) +'\tL90: ' + '{0:,}'.format(L90)+ ' bp'
print 'N95: ' + str(N95) +'\tL95: ' + '{0:,}'.format(L95)+ ' bp'


