import input_data
s1= input_data.seq1
s2= input_data.seq2
ScoreMat= input_data.BLOSUM52
g=-2

def needleman_wunsch(s1,s2,sm,g):
    M=len(s1)+1
    N=len(s2)+1
    F=[[0]*N for x in range(M)]
    P=[['o']*N for x in range(M)]
    for i in range(1,N):
        F[0][i]=i*g
        P[0][i]='l'
    for i in range(1,M):
        F[i][0]=i*g
        P[i][0]='u'
    for i in range(1,M):
        for j in range(1,N):
            sL=F[i][j-1]+g
            sU=F[i-1][j]+g
            sD=F[i-1][j-1]+ sm[s1[i-1]+s2[j-1]]
            F[i][j]=max(sL,sU,sD)
            if F[i][j]==sD:
                P[i][j]='d'
            elif F[i][j]==sU:
                P[i][j]='u'
            else:
                P[i][j]='l'
    return [F,P]

def alignment(FP,s1,s2):
    P=FP[1]
    F=FP[0]
    score=F[-1][-1]
    trace=''
    al1=''
    al2=''
    i=len(F)-1
    j=len(F[i])-1
    while i!=0 or j!=0:
        trace+=P[i][j]
        if trace[-1]=='d':
            i-=1
            j-=1
        elif trace[-1]=='u':
            i-=1
        else:
            j-=1
    trace=trace[::-1]

    for x in trace:
        if x=='d':
            al1+=s1[i]
            al2+=s2[j]
            i+=1
            j+=1
        elif x=='u':
            al1+=s1[i]
            al2+='-'
            i+=1
        else:
            al1+='-'
            al2+=s2[j]
            j+=1
    return [al1,al2,score]

NW_matrices= needleman_wunsch(s1,s2,ScoreMat,g)
Align=alignment(NW_matrices,s1,s2)
print('Optimal global alignment according to N&W algorithm:')
print(Align[0])
print(Align[1])
print('Alignment Score=',Align[2])
