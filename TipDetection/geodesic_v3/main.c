#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include <omp.h>
#include "mex.h"
//#include <time.h>

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

static inline int CompareSign(double a, double b) {
    int flag;
    if ((a > 0) && (b > 0)) {
        flag = 1;
    }
    else if ((a < 0) && (b < 0)) {
        flag = 1;
    }
    else {
        flag = 0;
    }
    return flag;
}


static inline double FindDistance(struct simplex* s) {
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

    C[0] = -s->vrtx[1][0] * s->vrtx[2][1] * s->vrtx[3][2] - s->vrtx[2][0] * s->vrtx[3][1] * s->vrtx[1][2] - s->vrtx[3][0] * s->vrtx[1][1] * s->vrtx[2][2]
           + s->vrtx[1][2] * s->vrtx[2][1] * s->vrtx[3][0] + s->vrtx[2][2] * s->vrtx[3][1] * s->vrtx[1][0] + s->vrtx[3][2] * s->vrtx[1][1] * s->vrtx[2][0];

    C[1] = s->vrtx[0][0] * s->vrtx[2][1] * s->vrtx[3][2] + s->vrtx[2][0] * s->vrtx[3][1] * s->vrtx[0][2] + s->vrtx[3][0] * s->vrtx[0][1] * s->vrtx[2][2]
           - s->vrtx[0][2] * s->vrtx[2][1] * s->vrtx[3][0] - s->vrtx[2][2] * s->vrtx[3][1] * s->vrtx[0][0] - s->vrtx[3][2] * s->vrtx[0][1] * s->vrtx[2][0];

    C[2] = -s->vrtx[0][0] * s->vrtx[1][1] * s->vrtx[3][2] - s->vrtx[1][0] * s->vrtx[3][1] * s->vrtx[0][2] - s->vrtx[3][0] * s->vrtx[0][1] * s->vrtx[1][2]
           + s->vrtx[0][2] * s->vrtx[1][1] * s->vrtx[3][0] + s->vrtx[1][2] * s->vrtx[3][1] * s->vrtx[0][0] + s->vrtx[3][2] * s->vrtx[0][1] * s->vrtx[1][0];

    C[3] = s->vrtx[0][0] * s->vrtx[1][1] * s->vrtx[2][2] + s->vrtx[1][0] * s->vrtx[2][1] * s->vrtx[0][2] + s->vrtx[2][0] * s->vrtx[0][1] * s->vrtx[1][2]
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
                exit_flag = 1;
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

//int readinput(char* filename, int* num, int**** data) {
//    FILE* fp;
//    int num_;
//    int x_size = 512;
//    int y_size = 511;
//    int z_size = 61;
//
////    int a = fopen_s(&fp, filename, "r");
//    fp = fopen(filename, "r");
//    fscanf(fp, "%d", &num_);
//
//    int*** arr = (int***)malloc(x_size * sizeof(int**));
//    for (int ii = 0; ii < x_size; ii++) {
//        arr[ii] = (int**)malloc(y_size * sizeof(int*));
//    }
//    for (int ii = 0; ii < x_size; ii++) {
//        for (int jj = 0; jj < y_size; jj++) {
//            arr[ii][jj] = (int*)malloc(z_size * sizeof(int));
//        }
//    }
//
//    for (int kk = 0; kk < z_size; kk++) {
//        for (int ii = 0; ii < x_size; ii++) {
//            for (int jj = 0; jj < y_size; jj++) {
//                fscanf(fp, "%d", &arr[ii][jj][kk]);
//            }
//        }
//    }
//
//    fclose(fp);
//    *num = num_;
//    *data = arr;
//    return (0);
//}
//
//
//int main() {
//
//    int*** fore_all = NULL;
//    int*** bound_all = NULL;
//    int node_num;
//    int bound_num;
//
//    char filename[50] = "/media/wang/D/geodesic/fore_all.dat";
//    readinput(filename, &node_num, &fore_all);
//    char filename2[50] = "/media/wang/D/geodesic/bound_all.dat";
//    readinput(filename2, &bound_num, &bound_all);
//
//
//    int x_size = 512, y_size = 511, z_size = 61;
//    int num_threads = 1;
//
//    float scale[5] = {3,4,5,6,7};
//    int scale_num = 5;
//    float resolution[3] = {1,1,1};
//    double* dist_map = (double*)malloc(scale_num * x_size * y_size * z_size * sizeof(double));
//
//    // initialization
//    int win_size[3];
//    int step[3];
//    double scale_max = 0;
//    for (int ii = 0; ii < scale_num; ii++) {
//        if (scale_max < scale[ii]) {
//            scale_max = scale[ii];
//        }
//    }
//    for (int ii = 0; ii < 3; ii++) {
//        step[ii] = (int)ceil(scale_max/resolution[ii]);
//        win_size[ii] = 2 * step[ii] + 1;
//    }
//
//
//    int* node_edge = (int*)malloc((node_num + 1) * 26 * sizeof(int));
//    int* node_capacity = (int*)malloc((node_num + 1) * 26 * sizeof(int));
//    int* node_edge_num = (int*)malloc((node_num + 1) * sizeof(int));
//    memset(node_edge_num, 0, (node_num + 1) * sizeof(int));
//
//    int bound_list_num = 0;
//    int** bound_list = (int**)malloc(3 * sizeof(int*));
//    for (int ii = 0; ii < 3; ii++) {
//        bound_list[ii] = (int*)malloc(bound_num * sizeof(int));
//    }
//
//    for (int ii = 0; ii < x_size; ii++) {
//        for (int jj = 0; jj < y_size; jj++) {
//            for (int kk = 0; kk < z_size; kk++) {
//                if (bound_all[ii][jj][kk] == 1) {
//                    bound_list[0][bound_list_num] = ii;
//                    bound_list[1][bound_list_num] = jj;
//                    bound_list[2][bound_list_num] = kk;
//                    bound_list_num = bound_list_num + 1;
//                }
//            }
//        }
//    }
//
//    float c[26];
//    float inf = 1e30;
//    int direction[26][4];
//    int direction_ind = 0;
//    for (int xx = -1; xx < 2; xx++) {
//        for (int yy = -1; yy < 2; yy++) {
//            for (int zz = -1; zz < 2; zz++) {
//                if (abs(xx) + abs(yy) + abs(zz) > 0) {
//                    direction[direction_ind][0] = xx;
//                    direction[direction_ind][1] = yy;
//                    direction[direction_ind][2] = zz;
//                    c[direction_ind] = sqrt(abs(xx)*resolution[0]*resolution[0] + abs(yy)*resolution[1]*resolution[1] + abs(zz)*resolution[2]*resolution[2]);
//                    direction_ind++;
//                }
//            }
//        }
//    }
//
//    omp_set_num_threads(num_threads);
//#pragma omp parallel
//    {
//        int node_s, node_t;
//        int xx, yy, zz, ii;
//#pragma omp for
//        for (ii = 0; ii < x_size; ii++) {
//            for (int jj = 0; jj < y_size; jj++) {
//                for (int kk = 0; kk < z_size; kk++) {
//                    if (fore_all[ii][jj][kk] > 0) {
//                        node_s = fore_all[ii][jj][kk];
//                        for (int dd = 0; dd < 26; dd++) {
//                            xx = ii + direction[dd][0];
//                            yy = jj + direction[dd][1];
//                            zz = kk + direction[dd][2];
//                            if (xx >= 0 && xx < x_size && yy >= 0 && yy < y_size && zz >= 0 && zz < z_size) {
//                                if (fore_all[xx][yy][zz] > 0) {
//                                    node_t = fore_all[xx][yy][zz];
//                                    node_edge[node_s * 26 + node_edge_num[node_s]] = node_t;
//                                    node_capacity[node_s * 26 + node_edge_num[node_s]] = dd;
//                                    node_edge_num[node_s] = node_edge_num[node_s] + 1;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//        struct body bd;
//        bd.numpoints = 0;
//        bd.coord = (double**)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(double*));
//        for (int ii = 0; ii < win_size[0] * win_size[1] * win_size[2]; ii++) {
//            bd.coord[ii] = (double*)malloc(3 * sizeof(double));
//        }
//
//        int** node_ind = (int**)malloc(3 * sizeof(int*));
//        for (int ii = 0; ii < 3; ii++) {
//            node_ind[ii] = (int*)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(int));
//        }
//        int* node_list = (int*)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(int));
//
//        float* node_dist = (float*)malloc((node_num + 1) * sizeof(float));
//        memset(node_dist, 127, (node_num + 1) * sizeof(float));
//        int* node_flag = (int*)malloc((node_num + 1) * sizeof(int));
//        memset(node_flag, 0, (node_num + 1) * sizeof(int));
//
//        int** edge_queue = (int**)malloc(2 * sizeof(int*));
//        float** edge_queue_time = (float**)malloc(2 * sizeof(float*));
//        for (int ii = 0; ii < 2; ii++) {
//            edge_queue[ii] = (int*)malloc(26 * win_size[0] * win_size[1] * win_size[2] * sizeof(int));
//            edge_queue_time[ii] = (float*)malloc(26 * win_size[0] * win_size[1] * win_size[2] * sizeof(float));
//        }
//        int* node_queue = (int*)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(int));
//
//        int bound_ind;
//#pragma omp for
//        for (bound_ind = 0; bound_ind < bound_num; bound_ind++) {
//            //printf("%d\n", bound_ind);
//            int node_list_num = 0;
//            for (int ii = 0; ii < win_size[0]; ii++) {
//                for (int jj = 0; jj < win_size[1]; jj++) {
//                    for (int kk = 0; kk < win_size[2]; kk++) {
//                        node_s = fore_all[bound_list[0][bound_ind] - step[0] + ii]
//                        [bound_list[1][bound_ind] - step[1] + jj][bound_list[2][bound_ind] - step[2] + kk];
//                        if (node_s > 0) {
//                            node_ind[0][node_list_num] = ii;
//                            node_ind[1][node_list_num] = jj;
//                            node_ind[2][node_list_num] = kk;
//                            node_list[node_list_num] = node_s;
//                            node_flag[node_s] = 1;
//                            node_list_num = node_list_num + 1;
//                        }
//                    }
//                }
//            }
//
//
//            int edgeind, edgeind_new, nodeind;
//
//            float current_time = 0;
//            float fire_time;
//            edgeind = 0;
//            node_s = fore_all[bound_list[0][bound_ind]][bound_list[1][bound_ind]][bound_list[2][bound_ind]];
//            node_dist[node_s] = 0;
//            for (int ee = 0; ee < node_edge_num[node_s]; ee++) {
//                edge_queue[0][edgeind] = node_s;
//                edge_queue[1][edgeind] = node_edge[node_s * 26 + ee];
//                edge_queue_time[0][edgeind] = 0;
//                edge_queue_time[1][edgeind] = c[node_capacity[node_s * 26 + ee]];
//                edgeind++;
//            }
//
//            while (edgeind > 0) {
//                current_time++;
//                edgeind_new = 0;
//                nodeind = 0;
//                for (int ii = 0; ii < edgeind; ii++) {
//                    edge_queue_time[0][ii] = edge_queue_time[0][ii] + 1;
//                    if (edge_queue_time[0][ii] < edge_queue_time[1][ii]) {
//                        edge_queue[0][edgeind_new] = edge_queue[0][ii];
//                        edge_queue[1][edgeind_new] = edge_queue[1][ii];
//                        edge_queue_time[0][edgeind_new] = edge_queue_time[0][ii];
//                        edge_queue_time[1][edgeind_new] = edge_queue_time[1][ii];
//                        edgeind_new++;
//                    }
//                    else
//                    {
//                        fire_time = current_time - (edge_queue_time[0][ii] - edge_queue_time[1][ii]);
//                        node_t = edge_queue[1][ii];
//                        if (node_dist[node_t] > inf) {
//                            node_queue[nodeind] = node_t;
//                            node_dist[node_t] = fire_time;
//                            nodeind++;
//                        }
//                        else if (node_dist[node_t] > fire_time) {
//                            node_dist[node_t] = fire_time;
//                        }
//                    }
//                }
//                for (int ii = 0; ii < nodeind; ii++) {
//                    node_s = node_queue[ii];
//                    for (int ee = 0; ee < node_edge_num[node_s]; ee++) {
//                        node_t = node_edge[node_s * 26 + ee];
//                        if (node_flag[node_t] == 1 && node_dist[node_t] > inf) {
//                            edge_queue[0][edgeind_new] = node_s;
//                            edge_queue[1][edgeind_new] = node_t;
//                            edge_queue_time[0][edgeind_new] = current_time - node_dist[node_s];
//                            edge_queue_time[1][edgeind_new] = c[node_capacity[node_s * 26 + ee]];
//                            edgeind_new++;
//                        }
//                    }
//                }
//                edgeind = edgeind_new;
//            }
//
//            for (int ss = 0; ss < scale_num; ss++) {
//                bd.numpoints = 0;
//                for (int ii = 0; ii < node_list_num; ii++) {
//                    if (node_dist[node_list[ii]] > scale[ss] - 0.5 && node_dist[node_list[ii]] < scale[ss] + 0.5) {
//                        bd.coord[bd.numpoints][0] = (double) node_ind[0][ii] - step[0];
//                        bd.coord[bd.numpoints][1] = (double) node_ind[1][ii] - step[1];
//                        bd.coord[bd.numpoints][2] = (double) node_ind[2][ii] - step[2];
//                        bd.numpoints++;
//                    }
//                }
//                if (bd.numpoints > 0) {
//                    double distance = gjk(&bd);
//                    dist_map[ss * x_size * y_size * z_size + bound_list[2][bound_ind] * x_size * y_size +
//                             bound_list[1][bound_ind] * x_size + bound_list[0][bound_ind]] = distance;
//                }
//            }
//            for (int ii = 0; ii < node_list_num; ii++) {
//                node_flag[node_list[ii]] = 0;
//                node_dist[node_list[ii]] = 1e38;
//            }
//
//        }
//
//        for (int ii = 0; ii < win_size[0] * win_size[1] * win_size[2]; ii++) {
//            free(bd.coord[ii]);
//        }
//        free(bd.coord);
//        for (int ii = 0; ii < 3; ii++) {
//            free(node_ind[ii]);
//        }
//        free(node_ind);
//        free(node_list);
//
//        free(node_dist);
//        free(node_flag);
//
//        for (int ii = 0; ii < 2; ii++) {
//            free(edge_queue[ii]);
//            free(edge_queue_time[ii]);
//        }
//        free(edge_queue);
//        free(edge_queue_time);
//        free(node_queue);
//    }
//    free(node_edge);
//    free(node_capacity);
//    free(node_edge_num);
//
//    for (int ii = 0; ii < 3; ii++) {
//        free(bound_list[ii]);
//    }
//    free(bound_list);
//
//    for (int ii = 0; ii < x_size; ii++) {
//        for (int jj = 0; jj < y_size; jj++) {
//            free(fore_all[ii][jj]);
//            free(bound_all[ii][jj]);
//        }
//    }
//    for (int ii = 0; ii < x_size; ii++) {
//        free(fore_all[ii]);
//        free(bound_all[ii]);
//    }
//    free(fore_all);
//    free(bound_all);
//
//    int a = 1;
//}



void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {
    //dist_map = (node_num, fore_all, bound_num, bound_all, scale, resolution, num_threads)

    int node_num;
    int bound_num;
    int x_size, y_size, z_size;
    int num_threads;
    int scale_num;

    double* pr_0 = mxGetPr(prhs[0]);
    node_num = (int)pr_0[0];
    double* pr_2 = mxGetPr(prhs[2]);
    bound_num = (int)pr_2[0];
    double* pr_4 = mxGetPr(prhs[4]);
    const mwSize* pr_4_dims = mxGetDimensions(prhs[4]);
    scale_num = (int)pr_4_dims[1];
    float* scale = (float*)malloc(scale_num * sizeof(float*));
    for (int ii = 0; ii < scale_num; ii++) {
        scale[ii] = (float) pr_4[ii];
    }
    double * pr_5 = mxGetPr(prhs[5]);
    float resolution[3] = {(float)pr_5[0],(float)pr_5[1],(float)pr_5[2]};
    double* pr_6 = mxGetPr(prhs[6]);
    num_threads = (int)pr_6[0];
    int max_num_threads = omp_get_num_procs();
    if (num_threads > max_num_threads) {
        printf("Num of threads is too large!");
        return 0;
    }

    const mwSize* dims = mxGetDimensions(prhs[1]);
    x_size = (int)dims[0];
    y_size = (int)dims[1];
    z_size = (int)dims[2];
    double* fore_all_pr = mxGetPr(prhs[1]);
    double* bound_all_pr = mxGetPr(prhs[3]);
    const mwSize dims_output[4] = { x_size, y_size, z_size, scale_num};
    plhs[0] = mxCreateNumericArray(4, dims_output, mxDOUBLE_CLASS, mxREAL);
    double* dist_map = mxGetPr(plhs[0]);

    int*** fore_all = (int***)malloc(x_size * sizeof(int**));
    int*** bound_all = (int***)malloc(x_size * sizeof(int**));
    for (int ii = 0; ii < x_size; ii++) {
        fore_all[ii] = (int**)malloc(y_size * sizeof(int*));
        bound_all[ii] = (int**)malloc(y_size * sizeof(int*));
    }
    for (int ii = 0; ii < x_size; ii++) {
        for (int jj = 0; jj < y_size; jj++) {
            fore_all[ii][jj] = (int*)malloc(z_size * sizeof(int));
            bound_all[ii][jj] = (int*)malloc(z_size * sizeof(int));
        }
    }
    for (int ii = 0; ii < x_size; ii++) {
        for (int jj = 0; jj < y_size; jj++) {
            for (int kk = 0; kk < z_size; kk++) {
                fore_all[ii][jj][kk] = (int)fore_all_pr[kk * x_size * y_size + jj * x_size + ii];
                bound_all[ii][jj][kk] = (int)bound_all_pr[kk * x_size * y_size + jj * x_size + ii];
                for (int ss = 0; ss < scale_num; ss++) {
                    dist_map[ss * x_size * y_size * z_size + kk * x_size * y_size + jj * x_size + ii] = NAN;
                }
            }
        }
    }



    // initialization
    int win_size[3];
    int step[3];
    double scale_max = 0;
    for (int ii = 0; ii < scale_num; ii++) {
        if (scale_max < scale[ii]) {
            scale_max = scale[ii];
        }
    }
    for (int ii = 0; ii < 3; ii++) {
        step[ii] = (int)ceil(scale_max/resolution[ii]);
        win_size[ii] = 2 * step[ii] + 1;
    }


    int* node_edge = (int*)malloc((node_num + 1) * 26 * sizeof(int));
    int* node_capacity = (int*)malloc((node_num + 1) * 26 * sizeof(int));
    int* node_edge_num = (int*)malloc((node_num + 1) * sizeof(int));
    memset(node_edge_num, 0, (node_num + 1) * sizeof(int));

    int bound_list_num = 0;
    int** bound_list = (int**)malloc(3 * sizeof(int*));
    for (int ii = 0; ii < 3; ii++) {
        bound_list[ii] = (int*)malloc(bound_num * sizeof(int));
        }

    for (int ii = 0; ii < x_size; ii++) {
        for (int jj = 0; jj < y_size; jj++) {
            for (int kk = 0; kk < z_size; kk++) {
                if (bound_all[ii][jj][kk] == 1) {
                    bound_list[0][bound_list_num] = ii;
                    bound_list[1][bound_list_num] = jj;
                    bound_list[2][bound_list_num] = kk;
                    bound_list_num = bound_list_num + 1;
                }
            }
        }
    }

    float c[26];
    float inf = 1e30;
    int direction[26][4];
    int direction_ind = 0;
    for (int xx = -1; xx < 2; xx++) {
        for (int yy = -1; yy < 2; yy++) {
            for (int zz = -1; zz < 2; zz++) {
                if (abs(xx) + abs(yy) + abs(zz) > 0) {
                    direction[direction_ind][0] = xx;
                    direction[direction_ind][1] = yy;
                    direction[direction_ind][2] = zz;
                    c[direction_ind] = sqrt(abs(xx)*resolution[0]*resolution[0] + abs(yy)*resolution[1]*resolution[1] + abs(zz)*resolution[2]*resolution[2]);
                    direction_ind++;
                }
            }
        }
    }

    omp_set_num_threads(num_threads);
#pragma omp parallel
    {
        int node_s, node_t;
        int xx, yy, zz, ii;
#pragma omp for
        for (ii = 0; ii < x_size; ii++) {
            for (int jj = 0; jj < y_size; jj++) {
                for (int kk = 0; kk < z_size; kk++) {
                    if (fore_all[ii][jj][kk] > 0) {
                        node_s = fore_all[ii][jj][kk];
                        for (int dd = 0; dd < 26; dd++) {
                            xx = ii + direction[dd][0];
                            yy = jj + direction[dd][1];
                            zz = kk + direction[dd][2];
                            if (xx >= 0 && xx < x_size && yy >= 0 && yy < y_size && zz >= 0 && zz < z_size) {
                                if (fore_all[xx][yy][zz] > 0) {
                                    node_t = fore_all[xx][yy][zz];
                                    node_edge[node_s * 26 + node_edge_num[node_s]] = node_t;
                                    node_capacity[node_s * 26 + node_edge_num[node_s]] = dd;
                                    node_edge_num[node_s] = node_edge_num[node_s] + 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        struct body bd;
        bd.numpoints = 0;
        bd.coord = (double**)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(double*));
        for (int ii = 0; ii < win_size[0] * win_size[1] * win_size[2]; ii++) {
            bd.coord[ii] = (double*)malloc(3 * sizeof(double));
        }

        int** node_ind = (int**)malloc(3 * sizeof(int*));
        for (int ii = 0; ii < 3; ii++) {
            node_ind[ii] = (int*)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(int));
        }
        int* node_list = (int*)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(int));

        float* node_dist = (float*)malloc((node_num + 1) * sizeof(float));
        memset(node_dist, 127, (node_num + 1) * sizeof(float));
        int* node_flag = (int*)malloc((node_num + 1) * sizeof(int));
        memset(node_flag, 0, (node_num + 1) * sizeof(int));

        int** edge_queue = (int**)malloc(2 * sizeof(int*));
        float** edge_queue_time = (float**)malloc(2 * sizeof(float*));
        for (int ii = 0; ii < 2; ii++) {
            edge_queue[ii] = (int*)malloc(26 * win_size[0] * win_size[1] * win_size[2] * sizeof(int));
            edge_queue_time[ii] = (float*)malloc(26 * win_size[0] * win_size[1] * win_size[2] * sizeof(float));
        }
        int* node_queue = (int*)malloc(win_size[0] * win_size[1] * win_size[2] * sizeof(int));

        int bound_ind;
#pragma omp for
        for (bound_ind = 0; bound_ind < bound_num; bound_ind++) {
            //printf("%d\n", bound_ind);
            int node_list_num = 0;
            for (int ii = 0; ii < win_size[0]; ii++) {
                for (int jj = 0; jj < win_size[1]; jj++) {
                    for (int kk = 0; kk < win_size[2]; kk++) {
                        node_s = fore_all[bound_list[0][bound_ind] - step[0] + ii]
                        [bound_list[1][bound_ind] - step[1] + jj][bound_list[2][bound_ind] - step[2] + kk];
                        if (node_s > 0) {
                            node_ind[0][node_list_num] = ii;
                            node_ind[1][node_list_num] = jj;
                            node_ind[2][node_list_num] = kk;
                            node_list[node_list_num] = node_s;
                            node_flag[node_s] = 1;
                            node_list_num = node_list_num + 1;
                        }
                    }
                }
            }


            int edgeind, edgeind_new, nodeind;

            float current_time = 0;
            float fire_time;
            edgeind = 0;
            node_s = fore_all[bound_list[0][bound_ind]][bound_list[1][bound_ind]][bound_list[2][bound_ind]];
            node_dist[node_s] = 0;
            for (int ee = 0; ee < node_edge_num[node_s]; ee++) {
                edge_queue[0][edgeind] = node_s;
                edge_queue[1][edgeind] = node_edge[node_s * 26 + ee];
                edge_queue_time[0][edgeind] = 0;
                edge_queue_time[1][edgeind] = c[node_capacity[node_s * 26 + ee]];
                edgeind++;
            }

            while (edgeind > 0) {
                current_time++;
                edgeind_new = 0;
                nodeind = 0;
                for (int ii = 0; ii < edgeind; ii++) {
                    edge_queue_time[0][ii] = edge_queue_time[0][ii] + 1;
                    if (edge_queue_time[0][ii] < edge_queue_time[1][ii]) {
                        edge_queue[0][edgeind_new] = edge_queue[0][ii];
                        edge_queue[1][edgeind_new] = edge_queue[1][ii];
                        edge_queue_time[0][edgeind_new] = edge_queue_time[0][ii];
                        edge_queue_time[1][edgeind_new] = edge_queue_time[1][ii];
                        edgeind_new++;
                    }
                    else
                    {
                        fire_time = current_time - (edge_queue_time[0][ii] - edge_queue_time[1][ii]);
                        node_t = edge_queue[1][ii];
                        if (node_dist[node_t] > inf) {
                            node_queue[nodeind] = node_t;
                            node_dist[node_t] = fire_time;
                            nodeind++;
                        }
                        else if (node_dist[node_t] > fire_time) {
                            node_dist[node_t] = fire_time;
                        }
                    }
                }
                for (int ii = 0; ii < nodeind; ii++) {
                    node_s = node_queue[ii];
                    for (int ee = 0; ee < node_edge_num[node_s]; ee++) {
                        node_t = node_edge[node_s * 26 + ee];
                        if (node_flag[node_t] == 1 && node_dist[node_t] > inf) {
                            edge_queue[0][edgeind_new] = node_s;
                            edge_queue[1][edgeind_new] = node_t;
                            edge_queue_time[0][edgeind_new] = current_time - node_dist[node_s];
                            edge_queue_time[1][edgeind_new] = c[node_capacity[node_s * 26 + ee]];
                            edgeind_new++;
                        }
                    }
                }
                edgeind = edgeind_new;
            }

            for (int ss = 0; ss < scale_num; ss++) {
                bd.numpoints = 0;
                for (int ii = 0; ii < node_list_num; ii++) {
                    if (node_dist[node_list[ii]] > scale[ss] - 0.5 && node_dist[node_list[ii]] < scale[ss] + 0.5) {
                        bd.coord[bd.numpoints][0] = (double) node_ind[0][ii] - step[0];
                        bd.coord[bd.numpoints][1] = (double) node_ind[1][ii] - step[1];
                        bd.coord[bd.numpoints][2] = (double) node_ind[2][ii] - step[2];
                        bd.numpoints++;
                    }
                }
                if (bd.numpoints > 0) {
                    double distance = gjk(&bd);
                    dist_map[ss * x_size * y_size * z_size + bound_list[2][bound_ind] * x_size * y_size +
                             bound_list[1][bound_ind] * x_size + bound_list[0][bound_ind]] = distance;
                }
            }
            for (int ii = 0; ii < node_list_num; ii++) {
                node_flag[node_list[ii]] = 0;
                node_dist[node_list[ii]] = 1e38;
            }

        }

        for (int ii = 0; ii < win_size[0] * win_size[1] * win_size[2]; ii++) {
            free(bd.coord[ii]);
        }
        free(bd.coord);
        for (int ii = 0; ii < 3; ii++) {
            free(node_ind[ii]);
        }
        free(node_ind);
        free(node_list);

        free(node_dist);
        free(node_flag);

        for (int ii = 0; ii < 2; ii++) {
            free(edge_queue[ii]);
            free(edge_queue_time[ii]);
        }
        free(edge_queue);
        free(edge_queue_time);
        free(node_queue);
    }
    free(node_edge);
    free(node_capacity);
    free(node_edge_num);

    for (int ii = 0; ii < 3; ii++) {
        free(bound_list[ii]);
    }
    free(bound_list);

    for (int ii = 0; ii < x_size; ii++) {
        for (int jj = 0; jj < y_size; jj++) {
            free(fore_all[ii][jj]);
            free(bound_all[ii][jj]);
        }
    }
    for (int ii = 0; ii < x_size; ii++) {
        free(fore_all[ii]);
        free(bound_all[ii]);
    }
    free(fore_all);
    free(bound_all);
    free(scale);

}