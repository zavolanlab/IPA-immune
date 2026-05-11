BEGIN {
   OFS="\t";
   mapqU=255;
}
{
if (substr($1,1,1)!="@") {

    # not used
    #m=and($2,0x80)/0x80+1;
    #m = int($2/0x80)%2 + 1; #in case and() is not available
    
    if ($1!=readNameOld) delete readSJs;
    readNameOld=$1;

    n=split($6,L,/[A-Z]/)-1;
    split($6,C,/[0-9]*/);
    t=1;g=$4;
    for (k=1;k<=n;k++) {#scan through CIGAR operations
        if (C[k+1]=="S" || C[k+1]=="I") {
           t+=L[k];
        } else if (C[k+1]=="D") {
           g+=L[k];
        } else if (C[k+1]=="N") {
           readStrand=and($2,0x10)/0x10+0; # read reverse strand
           mateStrand=and($2,0x20)/0x20+0; # mate reverse strand
           first=and($2,0x40)/0x40+0; # first in pair
           second=and($2,0x80)/0x80+0; # second in pair
           sj1=$3 "\t" g "\t" g+L[k]-1 "\t" readStrand "\t" mateStrand "\t" first "\t" second; #include strand in the key
           readSJs[sj1]++;
           if (readSJs[sj1]==1) {#only count this junction if it has not been counted for the same read
               SJ[sj1]=1;
               #if ($5>=mapqU) {
                   #SJu[sj1]++;
               #} else {
                   #SJm[sj1]++;
               #};
           };
           g+=L[k];
        } else { # M operation
           g+=L[k];
           t+=L[k];
        };
    };
};
};
END {

for (ii in SJ) {
    print ii;
    #print ii, SJu[ii]+0, SJm[ii]+0;
};

};
