#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "mex.h" 

#define dotProduct(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

struct simplex {
	int    nvrtx;       /**< Number of simplex's vertices. 			*/
	double vrtx[4][3]; /**< Coordinates of simplex's vertices. 		*/
	int    wids[4];    /**< Label of the simplex's vertices. 			*/
	double lambda[4];    /**< Barycentric coordiantes for each vertex.  */
};

struct body {
    int numpoints;    /**< Number of points defining the body.            */
    double** coord;  /**< Pointer to pointer to the points' coordinates. */
    //double  s[3];    /**< Support mapping computed last.                 */
};

inline int CompareSign(double a, double b) {
	int flag;
	if ((a > 0) && (b > 0)) {
		flag = 1;
	}
	else if((a < 0) && (b < 0)) {
		flag = 1;
	}
	else {
		flag = 0;
	}
	return flag;
}

inline double FindDistance(struct simplex* s) {
    double d;
    double v[3] = { 0,0,0 };
    for (int ii = 0; ii < s->nvrtx; ++ii) {
        v[0] = v[0] + s->lambda[ii] * s->vrtx[ii][0];
        v[1] = v[1] + s->lambda[ii] * s->vrtx[ii][1];
        v[2] = v[2] + s->lambda[ii] * s->vrtx[ii][2];
    }
    d = sqrt(dotProduct(v, v));
    return d;
}

//int readinput(const char* inputfile, double*** pts, int* out) {
//    int npoints = 0;
//    int idx = 0;
//    FILE* fp;
//    /* Open file. */
//
//    if ((fopen_s(&fp, inputfile, "r")) !=0 ) {
//        fprintf(stdout, "ERROR: input file %s not found!\n", inputfile);
//        return 1;
//    }
//
//    /* Read number of input vertices. */
//    if (fscanf_s(fp, "%d", &npoints) != 1)
//        return 1;
//
//    /* Allocate memory. */
//    double** arr = (double**)malloc(npoints * sizeof(double*));
//    for (int i = 0; i < npoints; i++)
//        arr[i] = (double*)malloc(3 * sizeof(double));
//
//    /* Read and store vertices' coordinates. */
//    for (idx = 0; idx < npoints; idx++)
//    {
//        if (fscanf_s(fp, "%lf %lf %lf\n", &arr[idx][0], &arr[idx][1], &arr[idx][2]) != 3)
//            return 1;
//    }
//
//    /* Close file. */
//    fclose(fp);
//    /* Pass pointers. */
//    *pts = arr;
//    *out = idx;
//    return (0);
//}

inline static void S1D(struct simplex* s) {
    int    i = 0;
    int    indexI = 1;
    int    FacetsTest[2];
    double det_ap = 0;
    double det_pb = 0;
    double pt = 0;
    double leng = 0;
    double a[3];
    double b[3];
    double t[3];
    double nu_fabs[3];

    for (i = 0; i < 3; ++i) {
        b[i] = s->vrtx[0][i];
        a[i] = s->vrtx[1][i];
        t[i] = b[i] - a[i];
        leng += t[i];
        nu_fabs[i] = fabs(t[i]);
    }

    if (nu_fabs[0] > nu_fabs[1]) {
        if (nu_fabs[0] > nu_fabs[2])
            indexI = 0;
        else
            indexI = 2;
    }
    else if (nu_fabs[0] < nu_fabs[1]) {
        if (nu_fabs[1] > nu_fabs[2])
            indexI = 1;
        else
            indexI = 2;
    }
    else if (nu_fabs[0] < nu_fabs[2]) {
        indexI = 2;
    }
    else if (nu_fabs[1] < nu_fabs[2]) {
        indexI = 2;
    }

    /* Project origin onto the 1D simplex - line */
    pt = dotProduct(b, t) / (dotProduct(t, t)) * (a[indexI] - b[indexI]) + b[indexI];

    /* Compute signed determinants */
    det_ap = a[indexI] - pt;
    det_pb = pt - b[indexI];

    /* Compare signs of AB and auxiliary simplices */
    FacetsTest[0] = CompareSign(t[indexI], -1 * det_ap);
    FacetsTest[1] = CompareSign(t[indexI], -1 * det_pb);

    if (FacetsTest[0] + FacetsTest[1] == 2) {
        /* The origin is between A and B */
        s->lambda[0] = det_ap * -1.0 / t[indexI];
        s->lambda[1] = 1 - s->lambda[0];
        s->wids[0] = 0;
        s->wids[1] = 1;
        s->nvrtx = 2;
    }
    else if (FacetsTest[0] == 0) {
        /* The origin is beyond A */
        s->lambda[0] = 1;
        s->wids[0] = 0;
        s->nvrtx = 1;
        for (i = 0; i < 3; ++i) {
            s->vrtx[0][i] = s->vrtx[1][i];
        }
    }
    else {
        /* The origin is behind B */
        s->lambda[0] = 1;
        s->wids[0] = 1;
        s->nvrtx = 1;
    }
}

inline static void S2D(struct simplex* s) {

    double    u[3];
    double    v[3];
    double    n[3];
    double    p[3];
    double    s1[3];
    double    p_coef;
    double    mu[3];
    int       indexI[2] = { 1,2 };
    int       indexJ = 0;
    double    mu_max = 0;
    double    C[3];
    int       FacesTest[3];

    for (int ii = 0; ii < 3; ++ii) {
        s1[ii] = s->vrtx[0][ii];
        u[ii] = s->vrtx[1][ii] - s->vrtx[0][ii];
        v[ii] = s->vrtx[2][ii] - s->vrtx[1][ii];
    }
    n[0] = u[1] * v[2] - u[2] * v[1];
    n[1] = u[2] * v[0] - u[0] * v[2];
    n[2] = u[0] * v[1] - u[1] * v[0];
    p_coef = dotProduct(s1, n) / dotProduct(n, n);
    p[0] = p_coef * n[0];
    p[1] = p_coef * n[1];
    p[2] = p_coef * n[2];

    mu[0] = s->vrtx[0][1] * s->vrtx[1][2] + s->vrtx[1][1] * s->vrtx[2][2] + s->vrtx[2][1] * s->vrtx[0][2]
          - s->vrtx[0][2] * s->vrtx[1][1] - s->vrtx[1][2] * s->vrtx[2][1] - s->vrtx[2][2] * s->vrtx[0][1];
    mu[1] = s->vrtx[0][0] * s->vrtx[1][2] + s->vrtx[1][0] * s->vrtx[2][2] + s->vrtx[2][0] * s->vrtx[0][2]
          - s->vrtx[0][2] * s->vrtx[1][0] - s->vrtx[1][2] * s->vrtx[2][0] - s->vrtx[2][2] * s->vrtx[0][0];
    mu[2] = s->vrtx[0][0] * s->vrtx[1][1] + s->vrtx[1][0] * s->vrtx[2][1] + s->vrtx[2][0] * s->vrtx[0][1]
          - s->vrtx[0][1] * s->vrtx[1][0] - s->vrtx[1][1] * s->vrtx[2][0] - s->vrtx[2][1] * s->vrtx[0][0];

    for (int ii = 0; ii < 3; ++ii) {
        if (fabs(mu[ii]) > fabs(mu_max)) {
            mu_max = mu[ii];
            indexJ = ii;
        }
    }

    int jj = 0;
    for (int ii = 0; ii < 3; ++ii) {
        if (indexJ != ii) {
            indexI[jj] = ii;
            jj++;
        }
    }

    C[0] = p[indexI[0]] * s->vrtx[1][indexI[1]] + s->vrtx[1][indexI[0]] * s->vrtx[2][indexI[1]] + s->vrtx[2][indexI[0]] * p[indexI[1]]
         - p[indexI[1]] * s->vrtx[1][indexI[0]] - s->vrtx[1][indexI[1]] * s->vrtx[2][indexI[0]] - s->vrtx[2][indexI[1]] * p[indexI[0]];
    C[1] = s->vrtx[0][indexI[0]] * p[indexI[1]] + p[indexI[0]] * s->vrtx[2][indexI[1]] + s->vrtx[2][indexI[0]] * s->vrtx[0][indexI[1]]
         - s->vrtx[0][indexI[1]] * p[indexI[0]] - p[indexI[1]] * s->vrtx[2][indexI[0]] - s->vrtx[2][indexI[1]] * s->vrtx[0][indexI[0]];
    C[2] = s->vrtx[0][indexI[0]] * s->vrtx[1][indexI[1]] + s->vrtx[1][indexI[0]] * p[indexI[1]] + p[indexI[0]] * s->vrtx[0][indexI[1]]
         - s->vrtx[0][indexI[1]] * s->vrtx[1][indexI[0]] - s->vrtx[1][indexI[1]] * p[indexI[0]] - p[indexI[1]] * s->vrtx[0][indexI[0]];

    FacesTest[0] = CompareSign(mu_max, C[0]);
    FacesTest[1] = CompareSign(mu_max, C[1]);
    FacesTest[2] = CompareSign(mu_max, C[2]);

    if (FacesTest[0] + FacesTest[1] + FacesTest[2] == 3) {
        for (int ii = 0; ii < 3; ++ii) {
            s->lambda[ii] = C[ii] / mu_max;
        }
        s->nvrtx = 3;
    }
    else {
        double d_ = INFINITY;
        double d;
        struct simplex* smin = NULL;
        struct simplex stmp1, stmp2, stmp3;
        struct simplex* stmp[3] = { &stmp1, &stmp2, &stmp3 };
        for (int ii = 0; ii < 3; ++ii) {
            if (FacesTest[ii] == 0) {
                stmp[ii]->nvrtx = 2;
                int kk = 0;
                for (int jj = 0; jj < 3; ++jj) {
                    if (jj != ii) {
                        stmp[ii]->vrtx[kk][0] = s->vrtx[jj][0];
                        stmp[ii]->vrtx[kk][1] = s->vrtx[jj][1];
                        stmp[ii]->vrtx[kk][2] = s->vrtx[jj][2];
                        kk++;
                    }
                }
                S1D(stmp[ii]);
                d = FindDistance(stmp[ii]);
                if (d < d_) {
                    smin = stmp[ii];
                    d_ = d;
                }
            }
        }
        s->nvrtx = smin->nvrtx;
        for (int ii = 0; ii < s->nvrtx; ++ii) {
            s->wids[ii] = smin->wids[ii];
            s->lambda[ii] = smin->lambda[ii];
            s->vrtx[ii][0] = smin->vrtx[ii][0];
            s->vrtx[ii][1] = smin->vrtx[ii][1];
            s->vrtx[ii][2] = smin->vrtx[ii][2];
        }
    }
}

inline static void S3D(struct simplex* s) {
    double    C[4];
    double    detM;
    int       FacesTest[4];

    C[0] = - s->vrtx[1][0] * s->vrtx[2][1] * s->vrtx[3][2] - s->vrtx[2][0] * s->vrtx[3][1] * s->vrtx[1][2] - s->vrtx[3][0] * s->vrtx[1][1] * s->vrtx[2][2]
           + s->vrtx[1][2] * s->vrtx[2][1] * s->vrtx[3][0] + s->vrtx[2][2] * s->vrtx[3][1] * s->vrtx[1][0] + s->vrtx[3][2] * s->vrtx[1][1] * s->vrtx[2][0];

    C[1] =   s->vrtx[0][0] * s->vrtx[2][1] * s->vrtx[3][2] + s->vrtx[2][0] * s->vrtx[3][1] * s->vrtx[0][2] + s->vrtx[3][0] * s->vrtx[0][1] * s->vrtx[2][2]
           - s->vrtx[0][2] * s->vrtx[2][1] * s->vrtx[3][0] - s->vrtx[2][2] * s->vrtx[3][1] * s->vrtx[0][0] - s->vrtx[3][2] * s->vrtx[0][1] * s->vrtx[2][0];

    C[2] = - s->vrtx[0][0] * s->vrtx[1][1] * s->vrtx[3][2] - s->vrtx[1][0] * s->vrtx[3][1] * s->vrtx[0][2] - s->vrtx[3][0] * s->vrtx[0][1] * s->vrtx[1][2]
           + s->vrtx[0][2] * s->vrtx[1][1] * s->vrtx[3][0] + s->vrtx[1][2] * s->vrtx[3][1] * s->vrtx[0][0] + s->vrtx[3][2] * s->vrtx[0][1] * s->vrtx[1][0];

    C[3] =   s->vrtx[0][0] * s->vrtx[1][1] * s->vrtx[2][2] + s->vrtx[1][0] * s->vrtx[2][1] * s->vrtx[0][2] + s->vrtx[2][0] * s->vrtx[0][1] * s->vrtx[1][2]
           - s->vrtx[0][2] * s->vrtx[1][1] * s->vrtx[2][0] - s->vrtx[1][2] * s->vrtx[2][1] * s->vrtx[0][0] - s->vrtx[2][2] * s->vrtx[0][1] * s->vrtx[1][0];
    
    detM = C[0] + C[1] + C[2] + C[3];

    FacesTest[0] = CompareSign(detM, C[0]);
    FacesTest[1] = CompareSign(detM, C[1]);
    FacesTest[2] = CompareSign(detM, C[2]);
    FacesTest[3] = CompareSign(detM, C[3]);

    if (FacesTest[0] + FacesTest[1] + FacesTest[2] + FacesTest[3] == 4) {
        for (int ii = 0; ii < 4; ++ii) {
            s->lambda[ii] = C[ii] / detM;
        }
    }
    else {
        double d_ = INFINITY;
        double d;
        struct simplex* smin = NULL;
        struct simplex stmp1, stmp2, stmp3, stmp4;
        struct simplex* stmp[4] = { &stmp1, &stmp2, &stmp3, &stmp4 };
        for (int ii = 0; ii < 4; ++ii) {
            if (FacesTest[ii] == 0) {
                int kk = 0;
                for (int jj = 0; jj < 4; ++jj) {
                    if (jj != ii) {
                        stmp[ii]->vrtx[kk][0] = s->vrtx[jj][0];
                        stmp[ii]->vrtx[kk][1] = s->vrtx[jj][1];
                        stmp[ii]->vrtx[kk][2] = s->vrtx[jj][2];
                        kk++;
                    }
                }
                S2D(stmp[ii]);
                d = FindDistance(stmp[ii]);
                if (d < d_) {
                    smin = stmp[ii];
                    d_ = d;
                }
            }
        }
        s->nvrtx = smin->nvrtx;
        for (int ii = 0; ii < s->nvrtx; ++ii) {
            s->wids[ii] = smin->wids[ii];
            s->lambda[ii] = smin->lambda[ii];
            s->vrtx[ii][0] = smin->vrtx[ii][0];
            s->vrtx[ii][1] = smin->vrtx[ii][1];
            s->vrtx[ii][2] = smin->vrtx[ii][2];
        }
    }
}

inline static void subalgorithm(struct simplex* s) {
    if (s->nvrtx == 4) {
        S3D(s);
    }
    else if (s->nvrtx == 3) {
        S2D(s);
    }
    else if (s->nvrtx == 2) {
        S1D(s);
    }
    else if (s->nvrtx == 1) {
        s->lambda[0] = 1;
    }
}

double gjk(struct body* bd) {
    struct simplex  s;
    double Wk[3];
    double v[3] = { bd->coord[0][0] ,bd->coord[0][1] ,bd->coord[0][2] };
    double distance = sqrt(dotProduct(v, v));
    s.nvrtx = 1;
    s.vrtx[0][0] = bd->coord[0][0];
    s.vrtx[0][1] = bd->coord[0][1];
    s.vrtx[0][2] = bd->coord[0][2];

    while (dotProduct(v, v) >= 1e-4 && s.nvrtx < 4) {

        // find wk
        double support = -INFINITY;
        double support_temp;
        for (int ii = 0; ii < bd->numpoints; ++ii) {
            support_temp = -bd->coord[ii][0] * v[0] - bd->coord[ii][1] * v[1] - bd->coord[ii][2] * v[2];
            if (support_temp > support) {
                support = support_temp;
                Wk[0] = bd->coord[ii][0];
                Wk[1] = bd->coord[ii][1];
                Wk[2] = bd->coord[ii][2];
            }
        }

        // exit condition
        int exit_flag = 0;
        for (int ii = 0; ii < s.nvrtx; ++ii) {
            if ((s.vrtx[ii][0] == Wk[0]) && (s.vrtx[ii][1] == Wk[1]) && (s.vrtx[ii][2] == Wk[2])) {
                exit_flag == 1;
            }
        }
        if ((dotProduct(v, v) - dotProduct(v, Wk)) <= (1e-8 * dotProduct(v, v))) {
            exit_flag = 1;
        }
        if (exit_flag == 1) {
            return distance;
        }
        
        s.vrtx[s.nvrtx][0] = Wk[0];
        s.vrtx[s.nvrtx][1] = Wk[1];
        s.vrtx[s.nvrtx][2] = Wk[2];
        s.nvrtx = s.nvrtx + 1;

        subalgorithm(&s);

        v[0] = v[1] = v[2] = 0;
        for (int ii = 0; ii < s.nvrtx; ++ii) {
            v[0] = v[0] + s.lambda[ii] * s.vrtx[ii][0];
            v[1] = v[1] + s.lambda[ii] * s.vrtx[ii][1];
            v[2] = v[2] + s.lambda[ii] * s.vrtx[ii][2];
        }
        distance = sqrt(dotProduct(v, v));   
    }
    return distance;
}


//int main() {
//    /* Squared distance computed by openGJK.                                 */
//    double dd;
//    /* Structure of simplex used by openGJK.                                 */
//    struct simplex  s;
//    struct body bd;
//    int nvrtx;
//    double distance;
//    double(**vrtx) = NULL;
//    /* Specify name of input files for body 1 and body 2, respectively.      */
//    char   inputfileA[40] = "userP.dat";
//    /* For importing openGJK this is Step 2: adapt the data structure for the
//     * two bodies that will be passed to the GJK procedure. */
//
//     /* Import coordinates of object 1. */
//    if (readinput(inputfileA, &vrtx,&nvrtx))
//        return (1);
//    bd.coord = vrtx;
//    bd.numpoints = nvrtx;
//    distance = gjk(&bd);
//    //double d = FindDistance(&s);
//    printf("%lf", distance);
//
//    double a = 1;
//    double b = 0;
//    double c = INFINITY;
//
//}

void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{

    double* inCoord;
    size_t  nCoord;
    int     i;
    double* distance;
    int     c = 3;
    int     count = 0;
    double** arr;

    /**************** PARSE INPUTS AND OUTPUTS **********************/
    /*----------------------------------------------------------------*/
    /* Examine input (right-hand-side) arguments. */
    //if (nrhs != 1) {
    //    mexErrMsgIdAndTxt("MyToolbox:gjk:nrhs", "One inputs required.");
    //}
    ///* Examine output (left-hand-side) arguments. */
    //if (nlhs != 1) {
    //    mexErrMsgIdAndTxt("MyToolbox:gjk:nlhs", "One output required.");
    //}

    ///* make sure the input argument are any numerical type */
    //if (!mxIsNumeric(prhs[0])) {
    //    mexErrMsgIdAndTxt("MyToolbox:gjk:notNumeric", "Input matrix must be type numeric.");
    //}


    ///* make sure the input argument have 3 columns */
    //if (mxGetM(prhs[0]) != 3) {
    //    mexErrMsgIdAndTxt("MyToolbox:gjk:notColumnVector", "Input must have 3 columns.");
    //}


    /*----------------------------------------------------------------*/
    /* CREATE DATA COMPATIBLE WITH MATALB  */

    /* create a pointer to the real data in the input matrix  */
    inCoord = mxGetPr(prhs[0]);

    /* get the length of each input vector */
    nCoord = mxGetN(prhs[0]);

    /* Create output */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

    /* get a pointer to the real data in the output matrix */
    distance = mxGetPr(plhs[0]);

    /* Copy data from Matlab's vectors into two new arrays */
    arr = (double**)mxMalloc(sizeof(double*) * (int)nCoord);

    for (i = 0; i < nCoord; i++)
        arr[i] = &inCoord[i * 3];

    /*----------------------------------------------------------------*/
    /* POPULATE BODIES' STRUCTURES  */

    struct body     bd; /* Structure of body A */

    /* Assign number of vertices to each body */
    bd.numpoints = (int)nCoord;

    bd.coord = arr;

    /*----------------------------------------------------------------*/
    /*CALL COMPUTATIONAL ROUTINE  */

    /* Compute squared distance using GJK algorithm */
    distance[0] = gjk(&bd);

    mxFree(arr);

}