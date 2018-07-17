#!/usr/bin/awk -f
{
    x1  = $1;
    G1 += log(x1);  
    N1++;	
    x2  = $2;
    G2 += log(x2);  
    N2++;	
    x3  = $3;
    G3 += log(x3);  
    N3++;	
    x4  = $4;
    G4 += log(x4);  
    N4++;	
    x5  = $5;
    G5 += log(x5);  
    N5++;	
    x6  = $6;
    G6 += log(x6);  
    N6++;	
    x7  = $7;
    G7 += log(x7);  
    N7++;	
    x8  = $8;
    G8 += log(x8);  
    N8++;	
    x9  = $9;
    G9 += log(x9);  
    N9++;	
}
 
END {
    print exp(G1/N1),exp(G2/N2),exp(G3/N3),exp(G4/N4),exp(G5/N5),exp(G6/N6),exp(G7/N7),exp(G8/N8),exp(G9/N9);
}

