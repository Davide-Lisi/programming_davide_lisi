def score_mat(filename):
    f=open(filename, 'r')
    l=f.readline()
    l=l.rstrip()
    a=int(l.find('A'))
    aa=l[a:a+20]
    D={}
    for i in range(len(aa)):
        l=f.readline()
        l=l.rstrip()
        line=l.split()
        for j in range(len(line)):
            D[ aa[i]+aa[j] ] = int(line[j][:-1])
    f.close()
    return D

def scoring(al1,al2,mat,g):
    s=0
    for i in range(len(al1)):
        if al1[i]=='-' or al2[i]=='-':
            s+=g
        elif al1[i]+al2[i] in mat:
            s+= mat[ al1[i]+al2[i] ]
        else:
            s+= mat[ al2[i]+al1[i] ]
    return s

def read_alignments(filename):
    f=open(filename, 'r')
    l=f.readline()
    l=l.rstrip()
    al1=''
    al2=''
    Alignments=[]
    while len(Alignments) < 3:
        if l[0]=='>':
            l=f.readline()
            l=l.rstrip()
        else:
            if al1=='':
                al1=l
                l=f.readline()
                l=l.rstrip()
            else:
                al2=l
                Alignments.append([al1,al2])
                al1=''
                al2=''
                l=f.readline()
                l=l.rstrip()
    return Alignments


blosum=score_mat('BLOSUM62.txt')
pam=score_mat('PAM250.txt')
g=-2
A=read_alignments('alignments.fasta.txt')

for i in A:
    sp=scoring(i[0],i[1],pam,g)
    sb=scoring(i[0],i[1],blosum,g)
    print(i[0])
    print(i[1])
    print('PAM250=   ',sp)
    print('BLOSUM62= ',sb)
    print('\n')

#this section prints the scores with the affine gap penalty

def scoring_affine(al1,al2,mat):
    s=0
    for i in range(len(al1)):
        if al1[i]=='-' or al2[i]=='-':
            if i==0:
                s-=2
            else:
                if al1[i]=='-':
                    if al1[i-1]=='-':
                        s-=0.5
                    else:
                        s-=2
                else:
                    if al2[i-1]=='-':
                        s-=0.5
                    else:
                        s-=2
        elif al1[i]+al2[i] in mat:
            s+= mat[ al1[i]+al2[i] ]
        else:
            s+= mat[ al2[i]+al1[i] ]
    return s

for i in A:
    sp=scoring_affine(i[0],i[1],pam)
    sb=scoring_affine(i[0],i[1],blosum)
    print(i[0])
    print(i[1])
    print('PAM250=   ',sp)
    print('BLOSUM62= ',sb)
    print('\n')

