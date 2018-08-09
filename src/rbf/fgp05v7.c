# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>


// compile this as gcc -mcmodel=medium -o fgp05v6 fgp05v6.c -lm

int ell, printout;
double err, alphaf;
FILE *g33, *g;

double min ( double fd[], int nd)
{
	int i;
	double ymin = 1e30;
	
	for ( i = 0; i < nd; i++ )
	{
		if (fd[i] < ymin)
		{
			ymin = fd[i];
		}
	}
	
	return ymin;
  
}

double max ( double fd[], int nd)
{
	int i;
	double ymax = -1e30;
	
	for ( i = 0; i < nd; i++ )
	{
		if (fd[i] > ymax)
		{
			ymax = fd[i];
		}
	}
	
	return ymax;
  
}

double *array;

int cmp(const void *a, const void *b)
{
    int ia = *(int *)a;
    int ib = *(int *)b;
    return array[ia] < array[ib] ? -1 : array[ia] > array[ib];
}
/*******************************************************************************/

double *lu_doolittle(int d, double** D, double rhs[]) 
{
	// http://www.sci.utah.edu/~wallstedt/LU.cpp
	// for general matricies
	
	int i, j, k, p;
	double sum, *x, temp;
	//double y[d];
	
	x = ( double * ) malloc ( d * sizeof ( double ) );

	
	for(k = 0; k < d; k++)
	{
		
		for(j = k ; j < d; j++)
		{
			sum = 0.0;
			
			for( p = 0; p < k; p++)
			{
				sum += D[k][p] * D[p][j];
			}
			
			D[k][j] = D[k][j] - sum; 
		}
		
		temp = 1.0 / D[k][k];
		
		for(i = k + 1; i < d; i++)
		{
			sum = 0.0;
			
			for(p = 0; p < k; p++)
			{
				sum += D[i][p] * D[p][k];
			}
			
			D[i][k] = (D[i][k] - sum) * temp;
		}
	}
	
	
	for( i = 0; i < d; i++)
	{
		sum = 0.0;
		
		for(k = 0; k < i; k++)
		{
			sum += D[i][k] * x[k];
		}
	
		x[i] = rhs[i] - sum; 
	}
  

	for( i = d - 1; i >= 0; --i)
	{
		sum = 0.0;
		
		for(k = i + 1; k < d; k++)
		{
			sum += D[i][k] * x[k];
		}
		
		x[i] = (x[i] - sum) / D[i][i];
	}
	
	return x;
	
}


/******************************************************************************/

void Step3_FGP05 ( int omega[], int nd, double** phi, int q, int m, int lset[])
{
	int i2, jj, j2;
	double dist2ell[nd - m];
	
	
	// for some reason this next section from step3_fp05 function seems to solve all our problems.
	// but why do I have to delete the rest of the code? Original code is commented out in fgp05_v2.
		
	int index[nd - m];   
	
	for(i2 = 0; i2 < nd - m; i2++)
	{
		index[i2] = i2;
	}
	
	ell = omega[m];				//
	
	//for ( i2 = 0; i2 < nd; i2++ )
	//{
	//	
	//	for ( j2 = 0; j2 < nd; j2++ )
	//	{
	//		printf("phi: %d %d %4.16f\n",i2,j2,phi[i2][j2]);  //why is phi different here?????????
	//	}
	//}
	
	if ( nd - m > q )					//if N-m+1 > q, we use +0 because m starts at 0 in the C code, not 1 as in the matlab code
	{
		//printf("flag =%d %d %d %d\n",flag, nd, i, nd - i + 1);
		
			for (j2 = m; j2 < nd; j2++ )
			{
				jj = omega[j2];
				dist2ell[j2 - m] = phi[ell - 1][jj - 1];  
				
				//printf("dist2ell: %d %d %4.16f\n",ell-1,jj - 1, phi[ell-1][jj - 1]);
			}
			
			//sorting starts
			
			array = dist2ell;
			qsort(index, nd - m, sizeof(*index), cmp);				//dist2ell is in the opposite order in the matlab code after sorting
			
			//sorting ends; sorted array is dist2ell[index[i]]
			
			for(i2 = 0; i2 < q; i2++)
			{
				lset[i2] = omega[index[i2] + m];   // Lset(1:q)=Omega1(PermIndx(1:q)+m-1);
				//printf("lset[i]=: %d %d %d %d %d %d\n",i, m, i+m, index[i+m], lset[i], omega1[index[i+m]]);
				//printf("lset[i]=: %d %d %lf\n",lset[i2], omega[i2], dist2ell[index[i2]]);
			}
			
	}
	else
	{
		for(i2 = 0;i2 < nd-m; i2++) 
		{ 
			lset[i2] = omega[m + i2];				//Lset(1:N-m+1)=Omega1(m:N);
			//printf("lset: %d %d\n",i,lset[i]);
		}
		
	}
	
	// the old code, which fails for some random numbers for unknown reasons is in fgp05_v2.

}


/******************************************************************************/

double *Step4_FGP05 (int pow2, double r0, int lset[], int nd, double rset[][nd], int dim)
{
	int pow3 = pow2 + 1;

	double *zeta;
	//double sysmx2[pow3][pow3];	//SysMx=zeros(nq+1,nq+1);
	double rhs[pow3], r, r00;
	int ii, jj, i, j, k, l, aa;
	
	
	double** sysmx2 = malloc(pow3 * sizeof(double*));    // allocate the rows

	for (aa = 0; aa < pow3; aa++)
	{
		sysmx2[aa] = malloc(pow3 * sizeof(double));    // allocate the columns
	}
	
	
	r00 = r0 * r0;
	
	for ( i = 0; i < pow3; i++ )		
	{
		for ( j = 0; j < pow3; j++ )	
		{
				sysmx2[i][j] = 0.0;			// if I don't do this the calculation fails
		}
	}
	
	for ( i = 0; i < pow2; i++ )		//for i=1:nq
	{
		ii = lset[i] - 1;						//i_=Lset(i); 
		//printf("ii = %d \n", ii);
		for ( j = i + 1; j < pow2; j++ )	//for j=i+1:nq
		{
			jj = lset[j] - 1;					//j_=Lset(j);
			//printf("jj = %d \n", jj);
			
			r = 0.0;
			
			for ( l = 0; l < dim; l++ )				//SysMx(i,j)=sqrt(dot(ri-rj,ri-rj)+c*c);
			{
				r += (rset[l][ii] - rset[l][jj]) * (rset[l][ii] - rset[l][jj]);
			}
			
			sysmx2[i][j] = sqrt(r + r00);
			
			
			//printf("sysmx2 = %lf\n", sysmx2[i][j]);
		}
	}
	
	
	for ( i = 0; i < pow2; i++ ) 		//SysMx(nq+1,1:nq)=1;
	{
		sysmx2[i][pow2] = 1.0;			
	}
	
	
	for ( i = 0; i < pow3; i++ )		//SysMx=SysMx+SysMx';
	{
		rhs[i] = 0.0;
		for ( j = 0; j < pow3; j++ )
		{
			sysmx2[i][j] = sysmx2[i][j] + sysmx2[j][i];
		}
	}
	
	for ( i = 0; i < pow2; i++ )			//for i=1:nq
	{
		sysmx2[i][i] = r0;		//SysMx(i,i)=abs(c);
	}
	
	rhs[0] = 1.0;
	
	
	//zeta = lu_cholesky(pow3, sysmx2, rhs);
	zeta = lu_doolittle(pow3, sysmx2, rhs);
	
	//for ( i = 0; i <= pow2; i++ )		
	//{
	//	printf("%d zeta= %lf \n", i,zeta[i]);
	//}
	//printf("\n\n");
	
	for (aa = 0; aa < pow3; ++aa)
	{
		free(sysmx2[aa]);    // this frees the columns
	}
	
	free(sysmx2);    // this frees the rows
	
	return zeta;
}


/******************************************************************************/

void Step5_FGP05 (double zeta[], int lset[], double res[], int pow3, int nd, double tau1[])
{
	//static double tau1[nd];
	double sum, myu;
	int ii, jj, i, j;
	
	for (i = 0; i < nd; i++) //for j=1:pow
	{
		tau1[i] = 0.0;
	}
	
	
	sum = 0.0;
	
	for (i = 0; i < pow3; i++) //for i=1:pow
	{
		ii = lset[i];		//i_=Lset(i);
		sum = sum + zeta[i] * res[ii - 1];					//sum=sum+zeta(i,1)*res(i_);
	}
	
	myu = sum / zeta[0];		//myu=sum/zeta(1,1);
	
	//printf("myu= %lf \n", myu);
	
	for (j = 0; j < pow3; j++) //for j=1:pow
	{
		jj = lset[j];			//j_=Lset(j);
		tau1[jj - 1] = myu * zeta[j];	//    tau1(j_)=myu*zeta(j,1);
		//printf("jj=%d tau1=%lf\n", jj, tau1[jj - 1]);
	}
	
}


/******************************************************************************/

void Step6_FGP05 (double tau[] , double delta[], int nd, double d[])
{
	double sum1, sum2, beta;
	int i;
	
	sum1 = 0.0;
	sum2 = 0.0;
	
	for (i = 0; i < nd; i++)				//beta=dot(tau,d)/dot(delta,d);
	{
		sum1 += tau[i] * d[i];
		sum2 += delta[i] * d[i];	
	}
	
	beta = sum1 / sum2;
	
	for (i = 0; i < nd; i++)		//delta1=tau-beta*delta;
	{
		delta[i] = tau[i] - beta * delta[i];
	}
	
}


/******************************************************************************/
/******************************************************************************/


double rbf_interp (double r0, int dim, int nd, double rset[][nd], double w[], double alpha, 
	int ni, double xi[][ni], double yi[], double yio[])
/*
  Purpose:

*/
{
	int i;
	int j;
	int k, l;
	double r, r00;
	double vv[nd];
	
	r00 = r0 * r0;
	double ave = 0.0;
	
	//FILE *g33;
	//g33=fopen("interpolated.dat", "a");
	
	
	clock_t begin4 = clock();
	
	for ( k = 0; k < ni; k++ )
	{
		//printf("aaa \n" );
		for ( i = 0; i < nd; i++ )
		{				
			r = 0.0;
			for ( l = 0; l < dim; l++ )
			{
				r += (xi[l][k] - rset[l][i]) * (xi[l][k] - rset[l][i]);
			}
			
			//r = sqrt(r);
			vv[i] = sqrt(r + r00);			//;v[i,j]=sqrt(r[j]*r[j]+r0*r0)
			
		}
	 
		yi[k] = 0.0;
		
		for ( i = 0; i < nd; i++ )
		{
			yi[k] = yi[k] + vv[i] * w[i];
		}
		
		yi[k] = yi[k] + alpha;
		//ave += fabs(yio[k] - yi[k]);
		

	}	
		
	clock_t end4 = clock();
	double time_spent4 = (double)(end4 - begin4) / CLOCKS_PER_SEC;
	
	if (printout == 2 || printout == 5)
	{
		printf("time4: %lf \n", time_spent4);
	
		for ( k = 0; k < ni; k++ )
		{
			printf("xi = %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f , %d , r0 = %.3f , yi = %.6f \n", 
			xi[0][k], xi[1][k], xi[2][k], xi[3][k], xi[4][k], xi[5][k], k + 1, r0, yi[k]);
		}
	}
	
	if (printout == 5)
	{
		fprintf(g33, "time4: %lf \n", time_spent4);
	
		for ( k = 0; k < ni; k++ )
		{
			fprintf(g33, "xi = %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f , %d , r0 = %.3f , yi = %.6f \n", 
			xi[0][k], xi[1][k], xi[2][k], xi[3][k], xi[4][k], xi[5][k], k + 1, r0, yi[k]);
		}
	}
	
	fflush(g33);
	
	//printf("xi = %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f , 5 , r0 = %.3f , yi = %.6f , yio = %.6f , diff = %.6f \n", 
		//xi[0][0], xi[1][0], xi[2][0], xi[3][0], xi[4][0], xi[5][0], r0, yi[0], yio[0], ave / ni);
	
}

/******************************************************************************/


double rbf_weight (double r0, int dim, int nd, double rset[][nd], double** phi, 
	double** sysmx, int freqnum, double** input, int q, double tol, int maxiter, double w[], 
	double** input2, int ni, double xi[][ni], double yi[], double yio[], double** yif)
/*
  Purpose:

    RBF_WEIGHT computes weights for radial basis function interpolation.

  Discussion:

	FGP05 algorithm

  Parameters:

  r0 - the parameter of the multiquadric kernel phi(r)=sqrt(r^2+r0^2);
% xd - array of sources of size [dim,nd], where dim is the space dimensionality, nd is the number of points)
% fd - array of right hand side of size [nd,1] 
% q - the power of the L-sets
% tol - tolerance
% maxiter - max number of iterations
% dim - no of dimensions
% printscreen=0 - no screen output, otherwise outputs error for each iteration step
% values:
% w - array of size [N,1] of coefficients lambda (solution)
% alpha - constant of solution 
% err - the error of the solution (in terms of the Linf of the residual)
% k - number of iterations
% errs - array of size [k,1] the Linf-norm residual error before each iteration
% example of call:
% [lambdas,alpha,err,k,errs]=FGP05_iterator(0.1,rand(2,100),rand(100,1),30,1e-10,100,1);

    Output, double w[nd], the weights.
*/
{
	int i, ells[nd], lset[nd], aa;
	//int lsets[q][nd]; 
	//double  zetas[q][nd];
	int j, l, pow2, pow3, k, ell1, omega[nd], lsetspow[nd];
	double r, r00, temp;
	double *zeta;
	//double zeta[nd];
	double res[nd], res1[nd], yd[nd];
	double tau[nd], tau1[nd], delta[nd], d[nd];
	double ymin, ymax, alpha, err, errs[maxiter];
	//double lambdas[nd]; 
	
	int** lsets = malloc(q * sizeof(int*));    // allocate the rows
	double** zetas = malloc(q * sizeof(double*));
	
	for (aa = 0; aa < q; ++aa)
	{
		lsets[aa] = malloc(nd * sizeof(int));    // allocate the columns
		zetas[aa] = malloc(nd * sizeof(double)); 
	}
	
	
	// Solve for the weights.

	// Step 1 started
	
	r00 = r0 * r0;
	/*ymin = min(fd, nd);
	ymax = max(fd, nd);
	alpha = 0.5 * (ymin + ymax);
	//printf("alpha1 = %lf\n", alpha);
	for ( i = 0; i < nd; i++ )
	{
		res[i] = fd[i] - alpha;
		w[i] = 0.0;
	}*/

	// Step 1 ended


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Intitial conditioning started. This needs to be done only once! The init_fgp05.c can do this, but then 
// the output of that code needs to be read.

	
	
	
	int pr = nd / 10;
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("Initial conditioning has started (Preconditioning #1.).\n");
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "Initial conditioning has started (Preconditioning #1.).\n");
		fflush(g33);
	}
	
	clock_t begin = clock();
	
	for ( i = 0; i < nd; i++ )
	{
		//if(i%pr==0)
		//{
		//	printf("i1 = %d \n", i);
		//}
		//ri[i] = rset[0][i];										//ri(:,1)=rset(:,i); I think this is not needed
		for ( j = 0; j < i; j++ )									//for ( j = 0; j < i; j++ )
		{
			//rj[j] = rset[0][j];									//rj(:,1)=rset(:,j); I think this is not needed
			r = 0.0;
			for ( l = 0; l < dim; l++ )
			{
				r += (rset[l][i] - rset[l][j]) * (rset[l][i] - rset[l][j]);
			}
			//printf("r: %lf\n", r);
			r = sqrt(r);
			phi[j][i] = r;
			sysmx[j][i] = sqrt(r * r + r00);
		}
	}
	
	
	for ( i = 0; i < nd; i++ )
	{
		
		//if(i%pr==0)
		//{
		//	printf("i2 = %d \n", i);
		//}
		
		for ( j = 0; j < nd; j++ )
		{
			phi[i][j] = phi[i][j] + phi[j][i];
			
			//printf("phi: %d %d %4.16f\n",i,j,phi[i][j]);
		}
	}

		
	for ( i = 0; i < nd; i++ )
	{
		sysmx[i][nd] = 1.0;
	}
	
	
	for ( i = 0; i <= nd; i++ )
	{
		//if(i%pr==0)
		//{
		//	printf("i3 = %d \n", i);
		//}
		
		for ( j = 0; j <= nd; j++ )
		{
			sysmx[i][j] = sysmx[i][j] + sysmx[j][i];
		}
	}
	
	
	for ( i = 0; i < nd; i++ )
	{
		//sysmx[i][i] = fabs(r0);		// r0 is always positive for APOGEE
		sysmx[i][i] = r0;
	}
	
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("time: %lf\n", time_spent);
		printf("Intitial conditioning has ended (Preconditioning #1.).\n");
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "time: %lf\n", time_spent);
		fprintf(g33, "Intitial conditioning has ended (Preconditioning #1.).\n");
		fflush(g33);
	}
	
	// Intitial conditioning ended.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("Computing lsets has started (Preconditioning #2.).\n");
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "Computing lsets has started (Preconditioning #2.).\n");
		fflush(g33);
	}
	
	clock_t begin2 = clock();
	
	//lset computation started here, this is the longest part
	
	//omega[0] = rand() % nd + 1;		// to get random numbers

	for (i = 0; i < nd; i++)			// need to start with 1 if we use rand(), 
	{									
		//int rnd = rand() % nd + 1;
		//
		//for(j = 0; j < i; j++)
		//{
		//	if(rnd==omega[j])
		//	{
		//		i--;
		//		break;
		//	}
		//}
		//
		//if(j >= i)
		//{
		//	omega[i]=rnd;
		//}
		
		omega[i] = i + 1;				// this is if we want increasing numbers
		
		//printf("omega = %d\n", omega[i]);
		
	}
	
	
	for ( i = 0; i < nd - 1; i++ )
	{
		
		if(i%pr==0 && (printout == 1 || printout == 2 || printout == 4 || printout == 5))
		{
			printf("i_ls = %d \n", i);
		}
		
		if (i%pr==0 && (printout == 3 || printout == 4 || printout == 5))
		{
			fprintf(g33, "i_ls = %d \n", i);
			fflush(g33);
		}
		
		Step3_FGP05(omega, nd, phi, q, i, lset);
		
		
		//printf("ell3 = %d \n", ell);
		ells[i] = ell;									//ells(m)=ell;
		//printf("ells = %d\n",ells[i]);
		
		if ( nd - i > q )					// pow=length(Lset); this sets up pow2
		{
			pow2 = q;
		}
		else
		{
			pow2 = nd - i;
		}
		//printf("pow2 = %d\n",pow2);
		
		ell1 = ell - 1;
		
		lsetspow[ell1] = pow2;			//Lsetspow(ell)=pow;
		//printf("lsetspow = %d\n",lsetspow[ell - 1]);
		
		for ( j = 0; j < pow2; j++ )	//Lsets(ell,1:pow)=Lset(1,1:pow); 
		{			
			lsets[j][ell1] = lset[j];			
			//printf("lset= %d\n", lset[j]);
		}
		
		zeta=Step4_FGP05(pow2, r0, lset, nd, rset, dim);   		// zeta=Step4_FGP05(c,Lset,rset);
		
		for ( j = 0; j < pow2; j++ )	//zetas(ell,1:pow)=zeta(1,1:pow); 
		{							
			zetas[j][ell1] = zeta[j];	
			//printf("zetas= %lf %lf\n", zetas[ell - 1][j],zeta[j]);
		}
		
	}
	
	clock_t end2 = clock();
	double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("time2: %lf \n", time_spent2);
		printf("Computing lsets has ended (Preconditioning #2.).\n");
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "time2: %lf \n", time_spent2);
		fprintf(g33, "Computing lsets has ended (Preconditioning #2.).\n");
		fflush(g33);
	}
	
	//lset computation ended
	
	/*for ( i = 0; i < q; i++ )
	{
		for ( j = 0; j < nd; j++ )
		{
		printf("lset= %d\n",lsets[i][j]);
		//printf("zetas= %.16f\n", zetas[i][j]);
		//	printf("zetas= %.16f\n", zetas[0][j]);
		}	
	
	}*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//FILE* g3;
	//g3=fopen("interpolated.dat", "a");
	//system("rm interpolated.dat");
	
	// iterations start
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("Computing weights has started... \n");	
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "Computing weights has started... \n");	
		fflush(g33);
	}
				
	int kk;
	
	int pr2 = freqnum / 100;
	
	for(kk = dim; kk < freqnum; kk++)			//for(kk = dim; kk < freqnum; kk++)
	{
		
		if (printout == 2 || printout == 5)
		{
			printf("******************************************************************************************************* \n");
			printf("frequency: %d \n", kk - dim);
		}
		else if (kk%pr2==0 && (printout == 1 || printout == 4))
		{
			printf("frequency: %d \n", kk - dim);
		}
		
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33, "******************************************************************************************************* \n");
			fprintf(g33, "frequency: %d \n", kk - dim);
			fflush(g33);
		}
		
		for(j = 0;j < nd; j++)
		{
			yd[j] = input[j][kk];
		}
		
		ymin = min(yd, nd);
		ymax = max(yd, nd);
		alpha = 0.5 * (ymin + ymax);
		//printf("alpha1 = %lf\n", alpha);
		for ( i = 0; i < nd; i++ )
		{
			res[i] = yd[i] - alpha;
			w[i] = 0.0;
		}
		
		
		
		if (fabs(max(res, nd)) > fabs(min(res, nd)))	//err=max(abs(max(res)),abs(min(res)));
		{
			err = fabs(max(res, nd));
		}
		else
		{
			err = fabs(min(res, nd));
		}
		//printf("%lf\n", err);
		
		k = 0; 
		double gamma, sum1, sum2, newconst1, newconst2, newconst;
		
		if (printout == 2 || printout == 5)
		{
			printf("Iterating... \n");
		}
		
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33, "Iterating... \n");
			fflush(g33);
		}
		
		clock_t begin3 = clock();
		
		//printf("%g %lf %d\n", err, tol, maxiter);
		while (err > tol & k < maxiter) //while err > tol & k < Nitermax
		{
			k++;						//k=k+1;
			
			if (printout == 2 || printout == 5)
			{
				printf("k = %d , err = %.8f , tol = %.8f , alpha = %.8f , w[0] = %.8f\n", k, err, tol, alpha, w[0]);
			}
			
			if (printout == 3 || printout == 4 || printout == 5)
			{
				fprintf(g33, "k = %d , err = %.8f , tol = %.8f , alpha = %.8f , w[0] = %.8f\n", k, err, tol, alpha, w[0]);\
				fflush(g33);
			}
			
			errs[k] = err;	//errs(k)=err;
		
			for (i = 0; i < nd; i++)
			{
				tau[i] = 0.0;			//tau=zeros(N,1);//tau1=zeros(N,1);
				//d[i] = 0.0;			
			}
			
			for (l = 0; l < nd - 1; l++)	// for m=1:N-1
			{
				ell = ells[l] - 1;
				//printf("ell %d\n", ell);
				pow3 = lsetspow[ell];							//pow=Lsetspow(ell);
				//printf("pow= %d\n", pow3);
				
						
				for (j = 0; j < pow3; j++)	//zeta(1,1:pow)=zetas(ell,1:pow);
				{
					zeta[j] = zetas[j][ell];		
					lset[j] = lsets[j][ell];					
					//printf("zeta= %lf \n", zetas[j][ell]);
				}
				
				Step5_FGP05(zeta, lset, res, pow3, nd, tau1);			//tau1=Step5_FGP05(zeta1,Lset,res,pow);
					
				for (i = 0; i < nd; i++)			// tau=tau+tau1;
				{
					tau[i] = tau[i] + tau1[i];
					//printf("tau= %lf \n", tau[i]);
					//printf("res= %lf \n", res[i]);
				}
				
			}
			
			//for (i = 0; i < nd; i++)			
			//{
				//printf("tau= %lf \n", tau[i]);
				//printf("res= %lf \n", res[i]);
				
			//}
			
			
			if (k == 1)
			{
				for ( i = 0; i < nd; i++ )			//delta=tau;
				{
					delta[i] = tau[i];
					//printf("delta= %lf \n", delta[i]);
				}
			}
			else
			{
				Step6_FGP05(tau, delta, nd, d);				//delta=Step6_FGP05(tau,delta,d);
			}
			
			
			for ( i = 0; i < nd; i++ )					//d=Step7_FGP05(SysMx,delta);
			{	
				//d[i] = 0.0;	
				temp = 0.0;
				for ( j = 0; j < nd; j++ )	
				{
					temp += sysmx[i][j] * delta[j];
					//printf("sysmx= %lf %lf %lf\n", sysmx[i][j], delta[j], d[i]);
				}
				
				d[i] = temp;
				//printf("d= %lf \n", d[i]);
			}		
			
		
			//gamma = 0.0;
			sum1 = 0.0;
			sum2 = 0.0;
		
			for (i = 0; i < nd; i++)		//gamma=dot(delta,res)/dot(delta,d);
			{
				sum1 += delta[i] * res[i];
				sum2 += delta[i] * d[i];
				//printf("d = %lf\n", d[i]);
			}
			
			gamma = sum1 / sum2;
			
			//printf("gamma = %lf\n", gamma);
		
		
			for (i = 0; i < nd; i++)		//res1=res-gamma*d;
			{
				res[i] = res[i] - gamma * d[i];
				//printf("res = %lf\n", res[i]);
			}		
			
			newconst1=max(res, nd);			//newconst1=max(res1);
			newconst2=min(res, nd);			//newconst2=min(res1);
		
			if (fabs(newconst1) > fabs(newconst2)) //err=max(abs(newconst1),abs(newconst2));
			{
				err = fabs(newconst1);
			}
			else
			{
				err = fabs(newconst2);
			}
		
			//printf("err1 = %lf\n", err);
		
			newconst = 0.5 * (newconst1 + newconst2); //newconst=.5*(newconst1+newconst2);
			//printf("newconst = %lf\n", newconst);
			alpha = alpha + newconst; //alpha1=alpha+newconst;
			
			//printf("alpha = %lf\n", alpha);
			
			for (i = 0; i < nd; i++)	//res1=res1-newconst;
			{
				res[i] = res[i] - newconst;
				//printf("res = %lf\n", res[i]);
			}
		
			//lambdas1=lambdas+gamma*delta;
	
			for (i = 0; i < nd; i++)
			{
				w[i] = w[i] + gamma * delta[i];
				//printf("w = %lf\n", w[i]);
			}
			//printf("%le %le %d %d\n", err, tol, k, maxiter);
		}
		
		alphaf = alpha;
		
		clock_t end3 = clock();
		double time_spent3 = (double)(end3 - begin3) / CLOCKS_PER_SEC;
		
		if (printout == 2 || printout == 5)
		{
			printf("time3: %lf , q = %d , r0 = %lf , iter = %d\n", time_spent3, q, r0, k);
		}
		
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33, "time3: %lf , q = %d , r0 = %lf , iter = %d\n", time_spent3, q, r0, k);
		}
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//interpolating here
	
		if (printout == 2 || printout == 5)
		{
			printf("Interpolating...\n");
		}
		
		if (printout == 5)
		{
			fprintf(g33, "Interpolating...\n");
			fflush(g33);
		}
		
		for(j = 0;j < ni; j++)
		{
			yio[j] = input2[j][kk];
		}
	
		
		rbf_interp (r0, dim, nd, rset, w, alphaf, ni, xi, yi, yio);
	
		
		for (j = 0; j < ni; j++ )
		{
			yif[j][kk - dim] = yi[j];
			//printf("%lf ", yif[j][kk - dim] );
		}
		

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	}
	
	
	// saving the final spectra
	
	//FILE *g;
	//g = fopen("outputspectra.dat", "w+");	
	
	int aa1 = freqnum - dim;
	
	for (i = 0; i < ni; i++ )
	{
		for (j = 0; j < aa1; j++ )
		{
			fprintf(g, "%le ", yif[i][j]);
			//printf("%lf\n ", yif[i][j] );
		}
		
		fprintf(g, "\n");
	}
	
	//fclose(g);
	
	//fclose(g3);
	
	/*for (aa = 0; aa < q; ++aa)
	{
		free(lsets[aa]);    // this frees the columns
		free(zetas[aa]);  
	}
	free(lsets);    // this frees the rows
	free(zetas);
	free(zeta);*/
	
	// iteration ends
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
	/*for ( i = 0; i < nd; i=i+1 )
	{
	//printf("x1 = %lf , x2 = %lf , x3 = %lf , x4 = %lf , x5 = %lf , x6 = %lf , yd = %lf , res = %lf, lambdas = %.6f\n", rset[0][i], rset[1][i], rset[2][i], rset[3][i], rset[4][i], rset[5][i], fd[i], phi[1][i], lambdas[i]);
	//printf("x1 = %3.2f , yd = %lf , phi = %lf, w = %.6f\n", 
	//rset[0][i], fd[i] , phi[0][i], w[i]);
	}*/
	
}


void statistics (double r0, int dim, int nd, int ni, double rset[][nd], double xi[][ni])
{
	int i, j, k, l;
	
	// computing the average distance 
	
	double r = 0.0, rden = 0.0, r2 = 0.0; 
	
	for ( i = 0; i < nd; i++ )
	{
		for ( j = 0; j < nd; j++ )								
		{								
			r = 0.0;
			for ( l = 0; l < dim; l++ )
			{
				r += (rset[l][i] - rset[l][j]) * (rset[l][i] - rset[l][j]);
			}
			
			r2 += sqrt(r);
			//printf("r: %lf\n", r2);
		}
	}
	
	//r = sqrt(r) / nd;
	r2 = r2 / nd / nd;
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("\n");
		printf("The average distances is: R_ave = %lf\n\n", r2);
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "\n");
		fprintf(g33, "The average distances is: R_ave = %lf\n\n", r2);
	}
	
	rden = 1.0 / r2;
	
	double maxdimvalues[dim], mindimvalues[dim];
	
	for ( l = 0; l < dim; l++ )
	{
		maxdimvalues[l] = max(rset[l], nd);
		mindimvalues[l] = min(rset[l], nd);
	
		if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
		{
			printf("Min and max values in dimension #%d = %lf %lf \n", l+1, mindimvalues[l], maxdimvalues[l]);
		}
		
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33, "Min and max values in dimension #%d = %lf %lf \n", l+1, mindimvalues[l], maxdimvalues[l]);
		}
	}
	
	fflush(g33);
	
	int aaa = 0, bad = 0, warn = 0;
	
	for ( j = 0; j < ni; j++ )
	{
		aaa = 0;
		
		for ( l = 0; l < dim; l++ )
		{
			if (xi[l][j] > maxdimvalues[l])
			{
				if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
				{
					printf("Extrapolation of point #%d in dimension #%d (%lf > %lf)\n", j+1, l+1, xi[l][j], maxdimvalues[l]);
				}
				if (printout == 3 || printout == 4 || printout == 5)
				{
					fprintf(g33, "Extrapolation of point #%d in dimension #%d (%lf > %lf)\n", j+1, l+1, xi[l][j], maxdimvalues[l]);
				}
				
				aaa++;
			}
			if (xi[l][j] < mindimvalues[l])
			{
				if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
				{
					printf("Extrapolation of point #%d in dimension #%d (%lf < %lf)\n", j+1, l+1, xi[l][j], mindimvalues[l]);
				}
				if (printout == 3 || printout == 4 || printout == 5)
				{
					fprintf(g33, "Extrapolation of point #%d in dimension #%d (%lf < %lf)\n", j+1, l+1, xi[l][j], mindimvalues[l]);
				}
				
				aaa++;
			}
		}
		
		if (aaa > 0)
		{
			warn++;
		}
		
		if (aaa == dim)
		{
			bad = 1;
	
			if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
			{
				printf("ERROR: Extrapolation of point #%d in all dimensions!!!\n", j+1);
			}
			if (printout == 3 || printout == 4 || printout == 5)
			{
				fprintf(g33, "ERROR: Extrapolation of point #%d in all dimensions!!!\n", j+1);
			}
		}
		
	}
	
	fflush(g33);
	
	// Density estimates
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("\n");
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "\n");
	}
	
	fflush(g33);
	
	r2 = 0.0;
	
	for ( i = 0; i < ni; i++ )
	{
		
		for ( j = 0; j < nd; j++ )								
		{								
			r = 0.0;
			for ( l = 0; l < dim; l++ )
			{
				r += (xi[l][i] - rset[l][j]) * (xi[l][i] - rset[l][j]);
			}
			//printf("r: %lf\n", r);
			r2 += sqrt(r);
		}
		
		r2 = 1.0 / (r2 / nd);
		
		if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
		{
			printf("Density around unknown point #%d is %lf\n", i+1, r2);
		}
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33, "Density around unknown point #%d is %lf\n", i+1, r2);
			fflush(g33);
		}
	}
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("Average density using the known points is %lf\n", rden);
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "Average density using the known points is %lf\n", rden);
	}
	
	// final say
	
	if (warn == 0 && (printout == 1 || printout == 2 || printout == 4 || printout == 5))
	{
		printf("\nThere were no points to be extrapolated in any dimensions, everything is good to go.\n\n");
	}
	
	if (warn == 0 && (printout == 3 || printout == 4 || printout == 5))
	{
		fprintf(g33, "\nThere were no points to be extrapolated in any dimensions, everything is good to go.\n\n");
	}
	
	
	
	if (warn > 0 && (printout == 1 || printout == 2 || printout == 4 || printout == 5))
	{
		printf("\nThere were %d points to be extrapolated in at least one dimension.\n", warn);
	}
	
	if (warn > 0 && (printout == 3 || printout == 4 || printout == 5))
	{
		fprintf(g33, "\nThere were %d points to be extrapolated in at least one dimension.\n", warn);
	}
	
	
	
	if (bad == 1 && (printout == 1 || printout == 2 || printout == 4 || printout == 5))
	{
		printf("\nThere were points to be extrapolated in all dimensions, program quits!\n\n");
		exit(0);
	}
	
	if (warn > 0 && (printout == 3 || printout == 4 || printout == 5))
	{
		fprintf(g33, "\nThere were points to be extrapolated in all dimensions, program quits!\n\n");
		fflush(g33);
		exit(0);
	}
	
	
	
	if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
	{
		printf("\n");
	}
	
	if (printout == 3 || printout == 4 || printout == 5)
	{
		fprintf(g33, "\n");
	}
}



main(int argc, char *argv[])
{

	if(argv[1]==NULL || argv[2]==NULL || argv[3]==NULL || argv[4]==NULL || argv[5]==NULL || argv[6]==NULL ||
		argv[7]==NULL || argv[8]==NULL || argv[9]==NULL || argv[10]==NULL || argv[11]==NULL)
	{
		printf("\n ./rbf inputname deletedname r0 nd q dim maxiter tol freqnum ni\n\n");
		printf("\n ./rbf final_657.dat deleted_657_carlos.dat 1.6 657 657 6 1000 1e-5 3034 72 1\n\n");
		
		
		printf("inputname: the name of the file with the known positions with frequencies.\n");
		printf("deletedname: the name of the file with the unknown positions.\n");
		printf("r0: scale factor (usually 1-2 * r_average).\n");
		printf("nd: number of known data points (the number of rows in inputname).\n");
		printf("q: the power of the L-sets (usually the same as nd if you have thousands of frequencies).\n");
		printf("dim: the number of dimensions.\n");
		printf("maxiter: maximum number of iterations.\n");
		printf("tol: tolerance of determining weights.\n");
		printf("freqnum: number of frequencies + the dimensions (the number of columns in inputname).\n");
		printf("ni: number of unknown data points (the number of rows in deletedname).\n");
		printf("printout: 0 no printout, no log file\n");
		printf("          1 moderate printout, no log file.\n");
		printf("          2 detailed printout, no log file.\n");
		printf("          3 no printout, detailed log file\n");
		printf("          4 moderate printout, detailed log file\n");
		printf("          5 detailed printout, very detailed log file\n\n");
		//printf("          6 very detailed printout, create log file, this is for debugging only\n\n");
		
		printf("The output interpolated spectra are in the outputspectra.dat file.\n\n");
		
		
		exit -1;
	}

	else
	{

		//int debug = 0;
		int i, j, l, dim, q, maxiter;
		int ni, nd, aa, sw1, freqnum;
		//double *xi, *yi; 
		double r0, tol;
		
		
		char* inputname = argv[1];  //final_3687.dat
		char* delname = argv[2];	//deleted_3687.dat
                char file[64];
		
		r0 = atof(argv[3]); 		/*scale factor*/
		nd = atoi(argv[4]);		/*no. of data points*/
		q = atoi(argv[5]);		/*the power of the L-sets*/
		dim = atoi(argv[6]);		/*no. of dimensions*/
		maxiter = atoi(argv[7]);	/*maximum no of iterations*/
		tol = atof(argv[8]);	 /*tolerance*/
		freqnum = atoi(argv[9]);	 /*no. of frequencies including the dimensions (3034 for me)*/
		ni = atoi(argv[10]);		/*no. of interpolation points*/	
		printout = atoi(argv[11]); //printout: 0 no printout, 1 moderate printout, 2 detailed printout
		
		
		FILE *f, *f2;
		
		//g33=fopen("logfile.dat", "wb");
                strcpy(file,inputname);
                strcat(file,".log");
		g33=fopen(file,"wb");
                strcpy(file,delname);
                strcat(file,".filled");
	        g = fopen(file, "w+");	
		f = fopen(inputname, "rt");			// file with known points
		f2 = fopen(delname, "rt");		// file with unknown points
		//f = fopen("aaa2", "rt");			// file with known points
		//f2 = fopen("aaa3", "rt");		// file with unknown points
		
			
		clock_t begin5 = clock();
		
		
		double** phi = malloc(nd * sizeof(double*));    // allocate the rows
		double** sysmx = malloc((nd + 1) * sizeof(double*));    // allocate the rows

		for (aa = 0; aa < nd; aa++)
		{
			phi[aa] = malloc(nd * sizeof(double));    // allocate the columns
		}
		for (aa = 0; aa <= nd; aa++)
		{
			sysmx[aa] = malloc((nd + 1) * sizeof(double));    // allocate the columns
		}
		
		double** yif = malloc(ni * sizeof(double*));    // allocate the rows
		
		for (aa = 0; aa < ni; aa++)
		{
			yif[aa] = malloc((freqnum - dim) * sizeof(double));    // allocate the columns
		}
		
		
		
		//static double phi[15625][15625], sysmx[15626][15626];
		
		double xi[dim][ni], yi[ni], yio[ni];
		
		double rset[dim][nd]; 		/*known points*/
		double yd[nd]; 				/*y values at known points*/
		double w[nd]; 				// weights
		
		//g33=fopen("interpolated.dat", "a");
		//srand(time(NULL));
		
		
		double v;
		
		//freqnum = 3034;
		//freqnum = 2;
		
		double** input = malloc(nd * sizeof(double*));    // allocate the rows
		double** input2 = malloc(ni * sizeof(double*));    // allocate the rows

		for (aa = 0; aa < nd; aa++)
		{
			input[aa] = malloc(freqnum * sizeof(double));    // allocate the columns
		}
		
		for (aa = 0; aa < ni; aa++)
		{
			input2[aa] = malloc(freqnum * sizeof(double));    // allocate the columns
		}
		
		if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
		{
			printf("Reading input with known positions...\n");
			
		}
		
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33,"Reading input with known positions...\n");
			fflush(g33);
		}
		
		
		for(i = 0; i < nd; i++)
		{
			
			for(j = 0;j < freqnum; j++)
			{
				//while (fscanf(infile, "%lf", &v) == 1)           //scanning for a readable  value in the file
				//{
				fscanf( f, "%lf", &v ); 
				input[i][j] = v;                                 //saving the integer value in the 2d array defined
					
				//}
			
			}	
		}
		fclose(f);
		
		for(i = 0; i < dim; i++)
		{
			for(j = 0;j < nd; j++)
			{
				rset[i][j] = input[j][i];
			}
		}
		
		
		
		if (printout == 1 || printout == 2 || printout == 4 || printout == 5)
		{
			printf("Reading input with unknown positions...\n");
		}
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33,"Reading input with unknown positions...\n");
			fflush(g33);
		}
		
		for(i = 0; i < ni; i++)
		{
			
			for(j = 0; j < dim; j++)				//< freqnum if we want to read the original fluxes too, < dim if we dont.
			{
				//while (fscanf(infile, "%lf", &v) == 1)           //scanning for a readable  value in the file
				//{
				fscanf( f2, "%lf", &v ); 
				input2[i][j] = v;                                 //saving the integer value in the 2d array defined
					
				//}
			
			}	
		}
		fclose(f2);
		
		for(i = 0; i < dim; i++)
		{
			for(j = 0; j < ni; j++)
			{
				xi[i][j] = input2[j][i];
			}
		}
		
		
		
		// the main loop going through all the frequencies is in rbf_weight
		
		statistics(r0, dim, nd, ni, rset, xi);
		
		rbf_weight (r0, dim, nd, rset, phi, sysmx, freqnum, input, q, tol, maxiter, w, input2, ni, xi, yi, yio, yif);
		
		
		
		clock_t end5 = clock();
		double time_spent5 = (double)(end5 - begin5) / CLOCKS_PER_SEC;
			
		printf("Final time5: %lf \n", time_spent5);
		
		if (printout == 3 || printout == 4 || printout == 5)
		{
			fprintf(g33, "Final time5: %lf \n", time_spent5);
		}
		fflush(g33);
		
	}
}
