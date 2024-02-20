/*You can include any C libraries that you normally use*/
#include "math.h"
#include "mex.h"   /*This C library is required*/

/*
This code implements the maximum asymmetry WASSR B0map algorithm
function [freqWMap] = wassr_b0map(ppmlist, posimage, negimage, mask4all
 */

#define PI 3.141592653589793

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;

    u = mxCalloc( n, sizeof(double) );
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else 
    {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) 
    {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
        qn=un=0.0;
	else
    {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	mxFree(u);
}
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
	int klo,khi,k;
	float h,b,a;

	klo=0;
	khi=n-1;
	while (khi-klo > 1) 
    {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double maxSymWASSR(double xval, double *xsp1, int nsp1, double *zsp1, double *xsp2, int nsp2, double *zsp2, double xstep2)
{
    double SSE,*SSE_1, freqWJ,cf;
    int count,j,index;
    double A1, A2, A3;

    cf = 2*xval;

    SSE_1 = mxCalloc( nsp1, sizeof(double) );
    count = 0;
    
#ifdef SKIP
/* Matlab code from Anup */    
for i = 1:W
    for j = 1:H
        if(mask4all(i,j)>0)
            zSpectW = squeeze(double(ImgfreqW1Smooth(i,j, :)))';
            Ref_zSpectW = zSpectW(nufreqW:-1:1);  %%%% reflect zSpect about centeral freq or ref freq
            IntSp_RefzSpectW = spline(tW, Ref_zSpectW(:),tsW);
            [minVal, minIndex] = min(zSpectW);
            start_c = -freqW(minIndex);  %%careful (-ve sign)--initial guess
            c_bnd = (start_c-0.2): IntStep : (start_c +0.2); % in ppm
            SSE_W = zeros(1, length(c_bnd));
            
            for k = 1: length(c_bnd)                
                SSE_W(k) = maxSymWASSR(c_bnd(k), freqW, nufreqW, zSpectW, IntSp_RefzSpectW, tsNew, f_N); 
            end
            [minSSE, minSSEIndex] =  min(SSE_W);
            EstiB0shift(i,j) = -c_bnd(minSSEIndex);
        else 
            EstiB0shift(i,j) = 0;   
        end
    end
end

function SSE = maxSymWASSR(c_bnd, freqW, nufreqW, vect1, ImgIntSp, tsNew, f_N)  
    SSE = 1000000.0;
    cW =2*c_bnd; 
    SSE_1 = zeros(1, nufreqW);
    count = 0;
    if(cW<0)
         for j = 1: nufreqW
             freqWJ = (freqW(j) + cW);
            if(freqWJ >freqW(1) && freqWJ <(f_N + cW))
               count = j;
               index = tsNew == int32(10000*freqWJ);
               SSE_1(j) = (vect1(j) - ImgIntSp(index)).^2; 
            end 
         end
    else 
        for j = 1: nufreqW
            freqWJ = (freqW(j) + cW);
            if(freqWJ >(freqW(1)+ cW) && freqWJ <f_N)
                count = j; 
                index = tsNew == int32(10000*freqWJ);
                SSE_1(j) = (vect1(j) - ImgIntSp(index)).^2; 
            end     
        end  
    end
    if (count>0) 
        SSE = sum(SSE_1);  
    end
end
            
#endif
    
/*  Begin Kim Paper - Magnetic Resonance in Medicine 61:1441-1450 (2009)*/
    
    for (j = 0; j < nsp1; j++)
    {
        A1 = xsp1[j] + cf;
        if (cf < 0)
        {
            A2 = xsp1[0];
            A3 = cf + xsp1[nsp1-1];
        }
        else
        {
            A2 = cf + xsp1[0];
            A3 = xsp1[nsp1-1];
        }
        if ( A1 > A2 && A1 < A3 )
        {
            index = (int)( (A1 - xsp2[0])/xstep2 );
            SSE_1[count] = pow( (zsp1[j] - zsp2[index]),2 );
            count++;
        }
    }
    
/*  End  Kim WASSR Paper- Magnetic Resonance in Medicine 61:1441?1450 (2009) */
    
    
    if (count > 0)
    {
        SSE = 0.0;
        for (j = 0; j < count; j++)
            SSE += SSE_1[j];
    }
    else
        SSE = 1.0e38;
    mxFree (SSE_1);
    return(SSE);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Inside mexFunction---*/
    /*Declarations*/
    const mxArray *ppmptr; /*The pointer to the ppm list array */
    const mxArray *pimgptr; /*The pointer to the pos images corresponding to ppm list */
    const mxArray *nimgptr; /*The pointer to the neg images corresponding to ppm list */
    const mxArray *maskptr; /*The pointer to the mask array */
    const mxArray *arg4; /*The pointer to the b0ppmstep array */
    
    double *ppmlist, *pimg, *nimg, *mask, *ppmlistasc, *b0ppmstepptr, b0ppmstep; /* Pointers to ppmlist,W,H,pos and negative images, and mask arrays */
    mwSize intfactor, nppm, ndim, nelem, ielem, ippm1, ippm2, ndimimg, nelemimg;
    const mwSize *dim_array_mask, *dim_array_image;
    mwSize i1, i2,nzsp1, nzsp2, rangezsp2, minzsp2index, minSSEindex, W, H;
    mwSize *psort;
    double maxppm, minppm, stepppm, stepppm2, stepppm2inv, ftemp1, ftemp2, maxpimg, maxnimg;
    double *zspp1, *zspn1, *zsp1, *xzsp, *xzsp2, *zsp2, *xzsp2r, *zsp2r, *y2;
    double *B0map, zspymin, zspxmin;
    double temp1,temp2;
    
        
    ppmptr = prhs[0];
    ppmlist = mxGetPr(ppmptr);
    nppm=mxGetNumberOfElements(ppmptr);
    pimgptr = prhs[1];
    pimg = mxGetPr(pimgptr);
    nimgptr = prhs[2];
    nimg = mxGetPr(nimgptr);
    nelemimg = mxGetNumberOfElements(nimgptr);
    ndimimg=mxGetNumberOfDimensions(nimgptr);
    dim_array_image=mxGetDimensions(nimgptr);
    mexPrintf("Himg %i Wimg %i Zimg %i\n",dim_array_image[0],dim_array_image[1],dim_array_image[2]);
    
    maskptr = prhs[3];
    mask = mxGetPr(maskptr);
    /* Get the dimensions in the mask array */
    nelem = mxGetNumberOfElements(maskptr);
    ndim=mxGetNumberOfDimensions(maskptr);
    dim_array_mask=mxGetDimensions(maskptr);
    H = dim_array_mask[0];
    W = dim_array_mask[1];
    mexPrintf("H %i W %i \n",H,W);
    
    if (nrhs > 4)
    {
        arg4 = prhs[4];
        b0ppmstepptr = mxGetPr(arg4);
        b0ppmstep = b0ppmstepptr[0];
    }
    else
        b0ppmstep = 0.005;
    
    maxnimg = -1.0E38;
    maxpimg = -1.0E38;
    for ( i1 = 0; i1 < nelemimg; i1++)
    {
        if (pimg[i1] > maxpimg)
            maxpimg = pimg[i1];
        if (nimg[i1] > maxnimg)
            maxnimg = nimg[i1];
    }
    mexPrintf("maxpimg %f maxnimg %f \n", maxpimg,maxnimg);
    
    plhs[0]=mxCreateDoubleMatrix(H,W,mxREAL);
    B0map=mxGetPr(plhs[0]);
    
    
    ppmlistasc = mxCalloc(nppm, sizeof(double));
    psort      = mxCalloc(nppm, sizeof(int));
    maxppm = -10000.0;
    minppm = 10000.0;
    for (ippm1 =0; ippm1 < nppm; ippm1++)
    {
        if (ppmlist[ippm1] > maxppm)
            maxppm = ppmlist[ippm1];
        if (ppmlist[ippm1] < minppm)
            minppm = ppmlist[ippm1];       
    }
    stepppm = (maxppm - minppm)/(nppm-1);
    for (ippm1 =0; ippm1 < nppm; ippm1++)
    {
        ppmlistasc[ippm1] = minppm + ippm1*stepppm;
    }
    for (ippm1 = 0; ippm1 < nppm; ippm1++)
    {
        psort[ippm1] = ippm1;
        for (ippm2 = 0; ippm2 < nppm; ippm2++)
        {
            if ( fabs(ppmlistasc[ippm2] - ppmlist[ippm1]) < 0.001 )
                psort[ippm1] = ippm2;
        }
    }
        
    
    intfactor = floor(stepppm/b0ppmstep);
    nzsp1 = 2*nppm-1;    
    nzsp2 = intfactor*(nzsp1-1)+1;
    stepppm2 = (stepppm/intfactor);
    stepppm2inv = (intfactor/stepppm);
    rangezsp2 = floor( 0.2*stepppm2inv );
	
	/* NW interpolated z-spectrum as output*/
	mwSize outndims = 3;
	const mwSize outdims[]={H,W,nzsp2};
	/*
	outdims[0] = H;
	outdims[1] = W;
	outdims[2] = nzsp2;
	*/
	plhs[1] = mxCreateNumericArray(outndims,outdims,mxDOUBLE_CLASS,mxREAL);
	double *Zspec_int;
	Zspec_int = mxGetPr(plhs[1]);
	
	/* NW interpolated ppm as output */
	plhs[2]=mxCreateDoubleMatrix(nzsp2,1,mxREAL);
	double *ppm_int;
    ppm_int=mxGetPr(plhs[2]);
	
	
        
    zspp1  = mxCalloc(nppm, sizeof(double));
    zspn1  = mxCalloc(nppm, sizeof(double));
    zsp1   = mxCalloc( nzsp1, sizeof(double) );  /* raw z spectrum sorted */
    xzsp   = mxCalloc( nzsp1, sizeof(double) );  /* raw z spectrum  cordinate */
    y2     = mxCalloc( nzsp1, sizeof(double) );  /* storage for second derivative */
    
    xzsp2  = mxCalloc( nzsp2, sizeof(double) );  /* Interpolated Zspetrum ppm coordinates */
    xzsp2r = mxCalloc( nzsp2, sizeof(double) );  /* Interpolated reflected Zspetrum ppm coordinates */
    zsp2   = mxCalloc( nzsp2, sizeof(double) );  /* Interpolated Zspetrum */
    zsp2r  = mxCalloc( nzsp2, sizeof(double) );  /* Interpolated reflected Zspetrum */
    
    for (ippm1 = 0; ippm1 < nppm; ippm1++) 
    {
        ippm2 = ippm1 + nppm -1;
        xzsp[ippm1] = -ppmlist[nppm-1-psort[ippm1]];                
        xzsp[ippm2] = ppmlist[psort[ippm1]];
        mexPrintf("index = %i psort = %i xzsp1 = %.4f    index = %i xzsp1 = %.4f\n",ippm1,psort[ippm1],xzsp[ippm1],ippm2,xzsp[ippm2]);
    }
    
    for (i2 =0; i2< nzsp2; i2++) 
    {
        xzsp2[i2] = -maxppm + i2 * stepppm2;
        xzsp2r[nzsp2-1-i2] = xzsp2[i2];
		/*NW*/
		ppm_int[i2] = xzsp2[i2];
    }
    
    for (ielem = 0; ielem < nelem; ielem++) 
    {
        B0map[ielem] = 0.0;
        if (mask[ielem] > 0) 
        {
            /* Create z spectrum from -maxppm to +maxppm; */
            for  (ippm1 = 0; ippm1 < nppm; ippm1++) 
            {
                i1 = ippm1*nelem + ielem;
                zspp1[ippm1] = pimg[i1];
                zspn1[ippm1] = nimg[i1];
            }
            for (ippm1 = 0; ippm1 < nppm; ippm1++) 
            {
                ippm2 = ippm1 + nppm -1;
                zsp1[ippm1] = zspn1[nppm-1-psort[ippm1]];                
                zsp1[ippm2] = zspp1[psort[ippm1]];
            }
            
            /* Find zspectrum minimum position and value */
            zspymin = 1.0e38;
            for (ippm1 = 0; ippm1 < nzsp1; ippm1++) 
            {
                if (zsp1[ippm1] < zspymin) 
                {
                    zspymin     = zsp1[ippm1];
                    zspxmin     = xzsp[ippm1];
                }
            }
            
            /* Call spline to get second derivatives */
            spline(xzsp, zsp1, nzsp1, zsp1[0], zsp1[nzsp1-1], y2);
            /* Call splint  for interpolations and find index for zsp2 minimum */
            zspymin = 1.0e38;
            for (i2 =0; i2< nzsp2; i2++) 
            {
                splint(xzsp, zsp1, y2, nzsp1, xzsp2[i2], &zsp2[i2]);
                if (zsp2[i2] < zspymin) 
                {
                    zspymin     = zsp2[i2];
                    zspxmin     = xzsp2[i2];
                }
                
                zsp2r[i2] = zsp2[nzsp2-1-i2];
				
				/*NW output variable*/
				/*Zspec_int[ielem*nzsp2 + i2] = zsp2[i2];*/
				Zspec_int[ielem + i2*nelem] = zsp2[i2];

            }
            /* Apply maximum symmetry algorithm to find B0 map position */
            temp1 = -zspxmin;
            ftemp2 = 1.0E30;
            minSSEindex = -rangezsp2;
            for (i1 = -rangezsp2; i1 <= rangezsp2; i1++) 
            {
                temp2 = temp1 + (i1 * stepppm2);
                ftemp1 = maxSymWASSR(temp2, xzsp, nzsp1, zsp1, xzsp2, nzsp2, zsp2r,stepppm2);
                if (ftemp1 < ftemp2) 
                {
                    minSSEindex = i1;
                    ftemp2 = ftemp1;
                }
            }
            B0map[ielem] = -( temp1 + (minSSEindex * stepppm2) );  
        }
    }
    
    mxFree(psort);
    mxFree(ppmlistasc);
    mxFree(zspp1);
    mxFree(zspn1);
    mxFree(y2);
    mxFree(zsp1);
    mxFree(zsp2);
    mxFree(zsp2r);
    mxFree(xzsp);
    mxFree(xzsp2);
    mxFree(xzsp2r);
    
    return;
}
