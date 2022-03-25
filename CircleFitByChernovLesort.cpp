int CircleFitByChernovLesort (Data& data, Circle& circleIni, reals LambdaIni, Circle& circle) 
/*                            <------------------ Input ------------------->  <-- Output -->

       Geometric circle fit to a given set of data points (in 2D)
		
       Input:  data     - the class of data (contains the given points):
		
	       data.n   - the number of data points
	       data.X[] - the array of X-coordinates
	       data.Y[] - the array of Y-coordinates
		          
               circleIni - parameters of the initial circle ("initial guess")
		        
	       circleIni.a - the X-coordinate of the center of the initial circle
	       circleIni.b - the Y-coordinate of the center of the initial circle
	       circleIni.r - the radius of the initial circle
		        
	       LambdaIni - the initial value of the control parameter "lambda"
	                   for the Levenberg-Marquardt procedure
	                   (common choice is a small positive number, e.g. 0.001)
		        
       Output:
	       integer function value is a code:
	                  0:  normal termination, the best fitting circle is 
	                      successfully found
	                  1:  the number of outer iterations exceeds the limit (99)
	                      (indicator of a possible divergence)
	                  2:  the number of inner iterations exceeds the limit (99)
	                      (another indicator of a possible divergence)
	                  3:  the coordinates of the center are too large
	                      (a strong indicator of divergence)
		                   
	       circle - parameters of the fitting circle ("best fit")
		        
	       circle.a - the X-coordinate of the center of the fitting circle
	       circle.b - the Y-coordinate of the center of the fitting circle
 	       circle.r - the radius of the fitting circle
 	       circle.s - the root mean square error (the estimate of sigma)
 	       circle.i - the total number of outer iterations (updating the parameters)
 	       circle.j - the total number of inner iterations (adjusting lambda)
 		        
       Algorithm by Nikolai Chernov and Claire Lesort
                         
       See a detailed description in the journal paper:
       
       N. Chernov and C. Lesort, "Least squares fitting of circles"
          in J. Math. Imag. Vision, volume 23, (2005), pages 239-251.
          
       the algorithm is designed to converge from any initial guess,
       but it is complicated and generally very slow
         
		Nikolai Chernov,  February 2014
*/
{
    int code,i,iter,inner,IterMAX=999;
    
    reals factorUp=10.,factorDown=0.04,lambda;
    reals Aold,Fold,Told,Anew,Fnew,Tnew,DD,H,aabb;
    reals Xi,Yi,Zi,Ui,Vi,Gi,CT,ST,D,ADF,SQ,DEN,FACT,DGDAi,DGDFi,DGDTi;
    reals H11,H12,H13,H22,H23,H33,F1,F2,F3,dA,dF,dT;
    reals epsilon=3.e-8;  
    reals G11,G22,G33,G12,G13,G23,D1,D2,D3;
    reals Xshift=0.,Yshift=0.,dX=One,dY=0.,aTemp,bTemp,rTemp;
    
    Circle Old,New;
    
//       starting with the given initial circle (initial guess)
    
    New = circleIni;
    
//       compute the root-mean-square error via function Sigma; see Utilities.cpp

    New.s = Sigma(data,New);
    
    Anew = One/Two/New.r;
    aabb = New.a*New.a + New.b*New.b;
    Fnew = (aabb - New.r*New.r)*Anew;
    Tnew = acos(-New.a/sqrt(aabb));
    if (New.b > 0.) Tnew = Two*Pi - Tnew;
    
    if (One+Four*Anew*Fnew < epsilon) 
    {
        Xshift += dX;
        Yshift += dY;
        
        New.a += dX;
        New.b += dY;
        aabb = New.a*New.a + New.b*New.b;
        Fnew = (aabb - New.r*New.r)*Anew;
        Tnew = acos(-New.a/sqrt(aabb));
        if (New.b > 0.) Tnew = Two*Pi - Tnew;
    }
    
//       initializing lambda, iteration counters, and the exit code
    
    lambda = LambdaIni;
    iter = inner = code = 0;
    
NextIteration:
    
    Aold = Anew;
    Fold = Fnew;
    Told = Tnew;
    Old = New;
    
    if (++iter > IterMAX) {code = 1;  goto enough;}
    
//       computing moments
    
shiftXY:
	
    DD = One + Four*Aold*Fold;
    D = sqrt(DD);
    CT = cos(Told);
    ST = sin(Told);
    
    H11=H12=H13=H22=H23=H33=F1=F2=F3=0.;
    
    for (i=0; i<data.n; i++)
    {
        Xi = data.X[i] + Xshift;
        Yi = data.Y[i] + Yshift;
        Zi = Xi*Xi + Yi*Yi;
        Ui = Xi*CT + Yi*ST;
        Vi =-Xi*ST + Yi*CT;
        
        ADF = Aold*Zi + D*Ui + Fold;
        SQ = sqrt(Four*Aold*ADF + One);
        DEN = SQ + One;
        Gi = Two*ADF/DEN;
        FACT = Two/DEN*(One - Aold*Gi/SQ);
        DGDAi = FACT*(Zi + Two*Fold*Ui/D) - Gi*Gi/SQ;
        DGDFi = FACT*(Two*Aold*Ui/D + One);
        DGDTi = FACT*D*Vi;
        
        H11 += DGDAi*DGDAi;
        H12 += DGDAi*DGDFi;
        H13 += DGDAi*DGDTi;
        H22 += DGDFi*DGDFi;
        H23 += DGDFi*DGDTi;
        H33 += DGDTi*DGDTi;
        
        F1 += Gi*DGDAi;
        F2 += Gi*DGDFi;
        F3 += Gi*DGDTi;
    }
    Old.g = New.g = sqrt(F1*F1 + F2*F2 + F3*F3);
    
try_again:
	
//        Cholesky decomposition
	
    G11 = sqrt(H11 + lambda);
    G12 = H12/G11;
    G13 = H13/G11;
    G22 = sqrt(H22 + lambda - G12*G12);
    G23 = (H23 - G12*G13)/G22;
    G33 = sqrt(H33 + lambda - G13*G13 - G23*G23);
    
    D1 = F1/G11;
    D2 = (F2 - G12*D1)/G22;
    D3 = (F3 - G13*D1 - G23*D2)/G33;
    
    dT = D3/G33;
    dF = (D2 - G23*dT)/G22;
    dA = (D1 - G12*dF - G13*dT)/G11;
    
//       updating the parameters
    
    Anew = Aold - dA;
    Fnew = Fold - dF;
    Tnew = Told - dT;
    
    if (One+Four*Anew*Fnew < epsilon) 
    {
        Xshift += dX;
        Yshift += dY;
        
        H = sqrt(One+Four*Aold*Fold);
        aTemp = -H*cos(Told)/(Aold+Aold) + dX;
        bTemp = -H*sin(Told)/(Aold+Aold) + dY;
        rTemp = One/abs(Aold+Aold);
        
        Aold = One/(rTemp + rTemp);
        aabb = aTemp*aTemp + bTemp*bTemp;
        Fold = (aabb - rTemp*rTemp)*Aold;
        Told = acos(-aTemp/sqrt(aabb));
        if (bTemp > 0.) Told = Two*Pi - Told;
        
        lambda *= factorUp;
        inner++;
        goto shiftXY;
    }
    
    H = sqrt(One+Four*Anew*Fnew);
    New.a = -H*cos(Tnew)/(Anew+Anew) - Xshift;
    New.b = -H*sin(Tnew)/(Anew+Anew) - Yshift;
    New.r = One/abs(Anew+Anew);
    New.s = Sigma(data,New);
    
    if ((abs(New.a-Old.a) + abs(New.b-Old.b) + abs(New.r-Old.r))/(One + Old.r) < epsilon) goto enough;
    
//       check if improvement is gained
    
    if (New.s < Old.s)    //   yes, improvement
    {
    	lambda *= factorDown;
        goto NextIteration;
    }
    else                       //   no improvement
    {
        if (++inner > IterMAX) {code = 2;  goto enough;}
        lambda *= factorUp;
        goto try_again;
    }
    
    //       exit
    
enough:
	
    Old.i = iter;
    Old.j = inner;
    
    circle = Old;
    
    return code;
}

//****************** Perturb *********************************

Circle Perturb (Circle& New, Circle& Old, reals range)
{
    Circle Perturbed;

    if (range==0.) return New;

    Perturbed.a = New.a + (New.a - Old.a)*(range*rand()/RAND_MAX - range/Two);
    Perturbed.b = New.b + (New.b - Old.b)*(range*rand()/RAND_MAX - range/Two);
    Perturbed.r = New.r + (New.r - Old.r)*(range*rand()/RAND_MAX - range/Two);

    return Perturbed;
}

