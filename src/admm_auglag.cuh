__device__
double eval_f_kernel(int n, double scale, double *x, double *param,
                     double YffR, double YffI,
                     double YftR, double YftI,
                     double YttR, double YttI,
                     double YtfR, double YtfI)
{
    int I = blockIdx.x;
    double f = 0, c1, c2, c3, c4, c5, c6, raug;

    int start = 31*I;
    for (int j = 0; j < 6; j++) {
        f += param[start + j]*x[j];
    }
    f += param[start + 24]*x[8];
    f += param[start + 25]*x[9];
    for (int j = 0; j < 6; j++) {
        f += 0.5*(param[start + 6+j]*(x[j] - param[start + 12+j])*(x[j] - param[start + 12+j]));
    }
    f += 0.5*(param[start + 26]*((x[8] - param[start + 28])*(x[8] - param[start + 28])));
    f += 0.5*(param[start + 27]*((x[9] - param[start + 29])*(x[9] - param[start + 29])));

    c1 = (x[0] - (YffR*x[4] + YftR*x[6] + YftI*x[7]));
    c2 = (x[1] - (-YffI*x[4] - YftI*x[6] + YftR*x[7]));
    c3 = (x[2] - (YttR*x[5] + YtfR*x[6] - YtfI*x[7]));
    c4 = (x[3] - (-YttI*x[5] - YtfI*x[6] - YtfR*x[7]));
    c5 = (x[6]*x[6] + x[7]*x[7] - x[4]*x[5]);
    c6 = (x[8] - x[9] - atan2(x[7], x[6]));

    f += param[start + 18]*c1;
    f += param[start + 19]*c2;
    f += param[start + 20]*c3;
    f += param[start + 21]*c4;
    f += param[start + 22]*c5;
    f += param[start + 30]*c6;

    raug = param[start + 23];
    f += 0.5*raug*(c1*c1);
    f += 0.5*raug*(c2*c2);
    f += 0.5*raug*(c3*c3);
    f += 0.5*raug*(c4*c4);
    f += 0.5*raug*(c5*c5);
    f += 0.5*raug*(c6*c6);

    f *= scale;
    __syncthreads();

    return f;
}

__device__
void eval_grad_f_kernel(int n, double scale, double *x, double *g, double *param,
                        double YffR, double YffI,
                        double YftR, double YftI,
                        double YttR, double YttI,
                        double YtfR, double YtfI)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int I = blockIdx.x;

    double c1, c2, c3, c4, c5, c6;
    double g1, g2, g3, g4, g5, g6, g7, g8, g9, g10;
    double raug;

    c1 = (x[0] - (YffR*x[4] + YftR*x[6] + YftI*x[7]));
    c2 = (x[1] - (-YffI*x[4] - YftI*x[6] + YftR*x[7]));
    c3 = (x[2] - (YttR*x[5] + YtfR*x[6] - YtfI*x[7]));
    c4 = (x[3] - (-YttI*x[5] - YtfI*x[6] - YtfR*x[7]));
    c5 = (x[6]*x[6] + x[7]*x[7] - x[4]*x[5]);
    c6 = (x[8] - x[9] - atan2(x[7], x[6]));

    int start = 31*I;
    g1 = param[start] + param[start + 6]*(x[0] - param[start + 12]);
    g2 = param[start + 1] + param[start + 7]*(x[1] - param[start + 13]);
    g3 = param[start + 2] + param[start + 8]*(x[2] - param[start + 14]);
    g4 = param[start + 3] + param[start + 9]*(x[3] - param[start + 15]);
    g5 = param[start + 4] + param[start + 10]*(x[4] - param[start + 16]);
    g6 = param[start + 5] + param[start + 11]*(x[5] - param[start + 17]);

    g9 = param[start + 24] + param[start + 26]*(x[8] - param[start + 28]);
    g10 = param[start + 25] + param[start + 27]*(x[9] - param[start + 29]);

    raug = param[start + 23];
    g1 += param[start + 18] + raug*c1;
    g2 += param[start + 19] + raug*c2;
    g3 += param[start + 20] + raug*c3;
    g4 += param[start + 21] + raug*c4;

    g5 += param[start + 18]*(-YffR) + param[start + 19]*(YffI) + param[start + 22]*(-x[5]) +
          raug*(-YffR)*c1 + raug*(YffI)*c2 + raug*(-x[5])*c5;
    g6 += param[start + 20]*(-YttR) + param[start + 21]*(YttI) + param[start + 22]*(-x[4]) +
          raug*(-YttR)*c3 + raug*(YttI)*c4 + raug*(-x[4])*c5;
    g7 = param[start + 18]*(-YftR) + param[start + 19]*(YftI) + param[start + 20]*(-YtfR) +
         param[start + 21]*(YtfI) + param[start + 22]*(2*x[6]) +
         raug*(-YftR)*c1 + raug*(YftI)*c2 + raug*(-YtfR)*c3 +
         raug*(YtfI)*c4 + raug*(2*x[6])*c5;
    g8 = param[start + 18]*(-YftI) + param[start + 19]*(-YftR) + param[start + 20]*(YtfI) +
         param[start + 21]*(YtfR) + param[start + 22]*(2*x[7]) +
         raug*(-YftI)*c1 + raug*(-YftR)*c2 + raug*(YtfI)*c3 +
         raug*(YtfR)*c4 + raug*(2*x[7])*c5;

    g7 += (param[start + 30] + raug*c6)*(x[7] / (x[6]*x[6] + x[7]*x[7]));
    g8 += (-((param[start + 30] + raug*c6)*(x[6] / (x[6]*x[6] + x[7]*x[7]))));
    g9 += param[start + 30] + raug*c6;
    g10 += (-(param[start + 30] + raug*c6));

    if (tx == 0 && ty == 0) {
        g[0] = scale*g1;
        g[1] = scale*g2;
        g[2] = scale*g3;
        g[3] = scale*g4;
        g[4] = scale*g5;
        g[5] = scale*g6;
        g[6] = scale*g7;
        g[7] = scale*g8;
        g[8] = scale*g9;
        g[9] = scale*g10;
    }

    __syncthreads();

    return;
}

__device__
void eval_h_kernel(int n, double scale, double *x, double *A, double *param,
                   double YffR, double YffI,
                   double YftR, double YftI,
                   double YttR, double YttI,
                   double YtfR, double YtfI)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int start = 31*blockIdx.x;
    double alrect, altheta, raug, c5, c6, x6sq_plus_x7sq;

    alrect = param[start + 22];
    altheta = param[start + 30];
    raug = param[start + 23];
    c5 = (x[6]*x[6] + x[7]*x[7] - x[4]*x[5]);
    c6 = (x[8] - x[9] - atan2(x[7], x[6]));
    x6sq_plus_x7sq = x[6]*x[6] + x[7]*x[7];

    if (tx == 0 && ty == 0) {
        // 1st column
        A[0] = scale*(param[start + 6] + raug); // A[1,1]
        A[4] = scale*(raug*(-YffR)); // A[5,1]
        A[6] = scale*(raug*(-YftR)); // A[7,1]
        A[7] = scale*(raug*(-YftI)); // A[8,1]

        // 2nd columns
        A[n + 1] = scale*(param[start + 7] + raug); // A[2,2]
        A[n + 4] = scale*(raug*(YffI));  // A[5,2]
        A[n + 6] = scale*(raug*(YftI));  // A[7,2]
        A[n + 7] = scale*(raug*(-YftR)); // A[8,2]

        // 3rd column
        A[n*2 + 2] = scale*(param[start + 8] + raug); // A[3,3]
        A[n*2 + 5] = scale*(raug*(-YttR)); // A[6,3]
        A[n*2 + 6] = scale*(raug*(-YtfR)); // A[7,3]
        A[n*2 + 7] = scale*(raug*(YtfI));  // A[8,3]

        // 4th column
        A[n*3 + 3] = scale*(param[start + 9] + raug); // A[4,4]
        A[n*3 + 5] = scale*(raug*(YttI)); // A[6,4]
        A[n*3 + 6] = scale*(raug*(YtfI)); // A[7,4]
        A[n*3 + 7] = scale*(raug*(YtfR)); // A[8,4]

        // 5th column
        A[n*4 + 4] = scale*(param[start + 10] + raug*(YffR*YffR) + raug*(YffI*YffI) + raug*(x[5]*x[5])); // A[5,5]
        A[n*4 + 5] = scale*(-(alrect + raug*c5) + raug*(x[4]*x[5])); // A[6,5]
        A[n*4 + 6] = scale*(raug*(YffR*YftR) + raug*(YffI*YftI) + raug*((-x[5])*(2*x[6]))); // A[7,5]
        A[n*4 + 7] = scale*(raug*(YffR*YftI) + raug*(YffI*(-YftR)) + raug*((-x[5])*(2*x[7]))); // A[8,5]

        // 6th column
        A[n*5 + 5] = scale*(param[start + 11] + raug*(YttR*YttR) + raug*(YttI*YttI) + raug*(x[4]*x[4])); // A[6,6]
        A[n*5 + 6] = scale*(raug*(YttR*YtfR) + raug*(YttI*YtfI) + raug*((-x[4])*(2*x[6]))); // A[7,6]
        A[n*5 + 7] = scale*(raug*((-YttR)*YtfI) + raug*(YttI*YtfR) + raug*((-x[4])*(2*x[7]))); // A[8,6]

        // 7th column
        A[n*6 + 6] = scale*((alrect + raug*c5)*2 + raug*(YftR*YftR) + raug*(YftI*YftI) +
                    raug*(YtfR*YtfR) + raug*(YtfI*YtfI) + raug*((2*x[6])*(2*x[6]))); // A[7,7]
        A[n*6 + 6] += scale*((altheta + raug*c6)*((-2*x[6]*x[7]) / (x6sq_plus_x7sq * x6sq_plus_x7sq)));
        A[n*6 + 6] += scale*(raug*((x[7] / x6sq_plus_x7sq)*(x[7] / x6sq_plus_x7sq)));
        A[n*6 + 7] = scale*(raug*(YftR*YftI) + raug*(YftI*(-YftR)) + raug*((-YtfR)*YtfI) +
                    raug*(YtfI*YtfR) + raug*((2*x[6])*(2*x[7]))); // A[8,7]
        A[n*6 + 7] += scale*((altheta + raug*c6)*((x[6]*x[6] - x[7]*x[7])/(x6sq_plus_x7sq * x6sq_plus_x7sq)));
        A[n*6 + 7] += scale*(raug*((x[7]/x6sq_plus_x7sq)*(-x[6]/x6sq_plus_x7sq)));
        A[n*6 + 8] = scale*(raug*(x[7]/x6sq_plus_x7sq)); // A[9,7]
        A[n*6 + 9] = scale*(-raug*(x[7]/x6sq_plus_x7sq)); // A[10,7]

        // 8th column
        A[n*7 + 7] = scale*((alrect + raug*c5)*2 + raug*(YftI*YftI) + raug*(YftR*YftR) +
                    raug*(YtfI*YtfI) + raug*(YtfR*YtfR) + raug*((2*x[7])*(2*x[7]))); // A[8,8]
        A[n*7 + 7] += scale*((altheta + raug*c6)*((2*x[6]*x[7]) / (x6sq_plus_x7sq * x6sq_plus_x7sq)));
        A[n*7 + 7] += scale*(raug*((-x[6]/x6sq_plus_x7sq)*(-x[6]/x6sq_plus_x7sq)));
        A[n*7 + 8] = scale*(raug*(-x[6]/x6sq_plus_x7sq)); // A[9,8]
        A[n*7 + 9] = scale*(-raug*(-x[6]/x6sq_plus_x7sq)); // A[10,8]

        // 9th column
        A[n*8 + 8] = scale*(param[start + 26] + raug); // A[9,9]
        A[n*8 + 9] = scale*(-raug); // A[10,9]

        // 10th column
        A[n*9 + 9] = scale*(param[start + 27] + raug); // A[10,10]
    }

    __syncthreads();

    if (tx < n && ty == 0) {
        #pragma unroll
        for (int j = 0; j < n; j++) {
            if (tx > j) {
                A[n*tx + j] = A[n*j + tx];
            }
        }
    }
    /*
    if (tx > ty) {
        A[n*tx + ty] = A[n*ty + tx];
    }
    */

    __syncthreads();

    return;
}

__global__ void
//__launch_bounds__(32, 16)
auglag_kernel(int nbranches, int n, int major_iter, int max_auglag, int line_start,
              double scale, double mu_max,
              double *u_curr, double *v_curr, double *l_curr,
              double *rho, double *wRIij, double *param,
              double *_YffR, double *_YffI,
              double *_YftR, double *_YftI,
              double *_YttR, double *_YttI,
              double *_YtfR, double *_YtfI,
              double *frBound, double *toBound)
{
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int I = blockIdx.x;

    bool terminate;
    int max_feval, max_minor, it, status, minor_iter;
    double mu, eta, omega;
    double YffR, YffI, YftR, YftI, YttR, YttI, YtfR, YtfI;

    extern __shared__ double shmem[];
    double *x, *xl, *xu;

    x = shmem;
    xl = shmem + n;
    xu = shmem + 2*n;

    int pij_idx = line_start + 8*I;

    if (tx == 0 && ty == 0) {
        #pragma unroll
        for (int j = 0; j < n; j++) {
            xl[j] = -INF;
            xu[j] = INF;
        }

        xl[4] = frBound[2*I];
        xu[4] = frBound[2*I+1];
        xl[5] = toBound[2*I];
        xu[5] = toBound[2*I+1];
        xl[8] = -2*M_PI;
        xu[8] = 2*M_PI;
        xl[9] = -2*M_PI;
        xu[9] = 2*M_PI;

        x[0] = u_curr[pij_idx];
        x[1] = u_curr[pij_idx+1];
        x[2] = u_curr[pij_idx+2];
        x[3] = u_curr[pij_idx+3];
        x[4] = min(xu[4], max(xl[4], u_curr[pij_idx+4]));
        x[5] = min(xu[5], max(xl[5], u_curr[pij_idx+5]));
        x[6] = wRIij[2*I];
        x[7] = wRIij[2*I+1];
        x[8] = min(xu[8], max(xl[8], u_curr[pij_idx+6]));
        x[9] = min(xu[9], max(xl[9], u_curr[pij_idx+7]));
    }

    YffR = _YffR[I]; YffI = _YffI[I];
    YftR = _YftR[I]; YftI = _YftI[I];
    YttR = _YttR[I]; YttI = _YttI[I];
    YtfR = _YtfR[I]; YtfI = _YtfI[I];

    int start = 31*I;

    param[start] = l_curr[pij_idx];
    param[start + 1] = l_curr[pij_idx+1];
    param[start + 2] = l_curr[pij_idx+2];
    param[start + 3] = l_curr[pij_idx+3];
    param[start + 4] = l_curr[pij_idx+4];
    param[start + 5] = l_curr[pij_idx+5];
    param[start + 6] = rho[pij_idx];
    param[start + 7] = rho[pij_idx+1];
    param[start + 8] = rho[pij_idx+2];
    param[start + 9] = rho[pij_idx+3];
    param[start + 10] = rho[pij_idx+4];
    param[start + 11] = rho[pij_idx+5];
    param[start + 12] = v_curr[pij_idx];
    param[start + 13] = v_curr[pij_idx+1];
    param[start + 14] = v_curr[pij_idx+2];
    param[start + 15] = v_curr[pij_idx+3];
    param[start + 16] = v_curr[pij_idx+4];
    param[start + 17] = v_curr[pij_idx+5];

    param[start + 24] = l_curr[pij_idx+6];
    param[start + 25] = l_curr[pij_idx+7];
    param[start + 26] = rho[pij_idx+6];
    param[start + 27] = rho[pij_idx+7];
    param[start + 28] = v_curr[pij_idx+6];
    param[start + 29] = v_curr[pij_idx+7];

    if (major_iter == 1) {
        param[start + 23] = 10.0;
        mu = 10.0;
    } else {
        mu = param[start + 23];
    }

    __syncthreads();

    eta = 1 / pow(mu, 0.1);
    omega = 1 / mu;
    max_feval = 500;
    max_minor = 200;

    it = 0;
    terminate = false;

    double cviol1, cviol2, cviol3, cviol4, cviol5, cviol6, cnorm;

    while (!terminate) {
        it += 1;

        // Solve the branch problem.
        cdriver_auglag(n, max_feval, max_minor, &status, &minor_iter, scale,
                       x, xl, xu, param, YffR, YffI, YftR, YftI, YttR, YttI, YtfR, YtfI,
                       &eval_f_kernel, &eval_grad_f_kernel, &eval_h_kernel);

        // Check the termination condition.

        cviol1 = x[0] - (YffR*x[4] + YftR*x[6] + YftI*x[7]);
        cviol2 = x[1] - (-YffI*x[4] - YftI*x[6] + YftR*x[7]);
        cviol3 = x[2] - (YttR*x[5] + YtfR*x[6] - YtfI*x[7]);
        cviol4 = x[3] - (-YttI*x[5] - YtfI*x[6] - YtfR*x[7]);
        cviol5 = x[6]*x[6] + x[7]*x[7] - x[4]*x[5];
        cviol6 = x[8] - x[9] - atan2(x[7], x[6]);

        cnorm = max(abs(cviol1), max(abs(cviol2), max(abs(cviol3), max(abs(cviol4), max(abs(cviol5), abs(cviol6))))));

        if (cnorm <= eta) {
            if (cnorm <= 1e-6) {
                terminate = true;
            } else {
                if (tx == 0 && ty == 0) {
                    param[start + 18] += mu*cviol1;
                    param[start + 19] += mu*cviol2;
                    param[start + 20] += mu*cviol3;
                    param[start + 21] += mu*cviol4;
                    param[start + 22] += mu*cviol5;
                    param[start + 30] += mu*cviol6;
                }

                eta = eta / pow(mu, 0.9);
                omega  = omega / mu;
            }
        } else {
            mu = min(mu_max, mu*10);
            eta = 1 / pow(mu, 0.1);
            omega = 1 / mu;
            param[start + 23] = mu;
        }

        if (it >= max_auglag) {
            terminate = true;
        }

        __syncthreads();
    }

    u_curr[pij_idx] = x[0];
    u_curr[pij_idx+1] = x[1];
    u_curr[pij_idx+2] = x[2];
    u_curr[pij_idx+3] = x[3];
    u_curr[pij_idx+4] = x[4];
    u_curr[pij_idx+5] = x[5];
    u_curr[pij_idx+6] = x[8];
    u_curr[pij_idx+7] = x[9];
    wRIij[2*I] = x[6];
    wRIij[2*I+1] = x[7];
    param[start + 23] = mu;

    __syncthreads();

    return;
}

double eval_f(int I, int n, double scale, double *x, double *param,
              double YffR, double YffI,
              double YftR, double YftI,
              double YttR, double YttI,
              double YtfR, double YtfI)
{
    double f = 0, c1, c2, c3, c4, c5, c6, raug;

    int start = 31*I;
    for (int j = 0; j < 6; j++) {
        f += param[start + j]*x[j];
    }
    f += param[start + 24]*x[8];
    f += param[start + 25]*x[9];
    for (int j = 0; j < 6; j++) {
        f += 0.5*(param[start + 6+j]*((x[j] - param[start + 12+j])*(x[j] - param[start + 12+j])));
    }
    f += 0.5*(param[start + 26]*((x[8] - param[start + 28])*(x[8] - param[start + 28])));
    f += 0.5*(param[start + 27]*((x[9] - param[start + 29])*(x[9] - param[start + 29])));

    c1 = (x[0] - (YffR*x[4] + YftR*x[6] + YftI*x[7]));
    c2 = (x[1] - (-YffI*x[4] - YftI*x[6] + YftR*x[7]));
    c3 = (x[2] - (YttR*x[5] + YtfR*x[6] - YtfI*x[7]));
    c4 = (x[3] - (-YttI*x[5] - YtfI*x[6] - YtfR*x[7]));
    c5 = (x[6]*x[6] + x[7]*x[7] - x[4]*x[5]);
    c6 = (x[8] - x[9] - atan2(x[7], x[6]));

    f += param[start + 18]*c1;
    f += param[start + 19]*c2;
    f += param[start + 20]*c3;
    f += param[start + 21]*c4;
    f += param[start + 22]*c5;
    f += param[start + 30]*c6;

    raug = param[start + 23];
    f += 0.5*raug*(c1*c1);
    f += 0.5*raug*(c2*c2);
    f += 0.5*raug*(c3*c3);
    f += 0.5*raug*(c4*c4);
    f += 0.5*raug*(c5*c5);
    f += 0.5*raug*(c6*c6);

    f *= scale;

    return f;
}

void eval_grad_f(int I, int n, double scale, double *x, double *g, double *param,
                double YffR, double YffI,
                double YftR, double YftI,
                double YttR, double YttI,
                double YtfR, double YtfI)
{
    double c1, c2, c3, c4, c5, c6;
    double g1, g2, g3, g4, g5, g6, g7, g8, g9, g10;
    double raug;

    c1 = (x[0] - (YffR*x[4] + YftR*x[6] + YftI*x[7]));
    c2 = (x[1] - (-YffI*x[4] - YftI*x[6] + YftR*x[7]));
    c3 = (x[2] - (YttR*x[5] + YtfR*x[6] - YtfI*x[7]));
    c4 = (x[3] - (-YttI*x[5] - YtfI*x[6] - YtfR*x[7]));
    c5 = (x[6]*x[6] + x[7]*x[7] - x[4]*x[5]);
    c6 = (x[8] - x[9] - atan2(x[7], x[6]));

    int start = 31*I;
    g1 = param[start] + param[start + 6]*(x[0] - param[start + 12]);
    g2 = param[start + 1] + param[start + 7]*(x[1] - param[start + 13]);
    g3 = param[start + 2] + param[start + 8]*(x[2] - param[start + 14]);
    g4 = param[start + 3] + param[start + 9]*(x[3] - param[start + 15]);
    g5 = param[start + 4] + param[start + 10]*(x[4] - param[start + 16]);
    g6 = param[start + 5] + param[start + 11]*(x[5] - param[start + 17]);

    g9 = param[start + 24] + param[start + 26]*(x[8] - param[start + 28]);
    g10 = param[start + 25] + param[start + 27]*(x[9] - param[start + 29]);

    raug = param[start + 23];
    g1 += param[start + 18] + raug*c1;
    g2 += param[start + 19] + raug*c2;
    g3 += param[start + 20] + raug*c3;
    g4 += param[start + 21] + raug*c4;

    g5 += param[start + 18]*(-YffR) + param[start + 19]*(YffI) + param[start + 22]*(-x[5]) +
          raug*(-YffR)*c1 + raug*(YffI)*c2 + raug*(-x[5])*c5;
    g6 += param[start + 20]*(-YttR) + param[start + 21]*(YttI) + param[start + 22]*(-x[4]) +
          raug*(-YttR)*c3 + raug*(YttI)*c4 + raug*(-x[4])*c5;
    g7 = param[start + 18]*(-YftR) + param[start + 19]*(YftI) + param[start + 20]*(-YtfR) +
         param[start + 21]*(YtfI) + param[start + 22]*(2*x[6]) +
         raug*(-YftR)*c1 + raug*(YftI)*c2 + raug*(-YtfR)*c3 +
         raug*(YtfI)*c4 + raug*(2*x[6])*c5;
    g8 = param[start + 18]*(-YftI) + param[start + 19]*(-YftR) + param[start + 20]*(YtfI) +
         param[start + 21]*(YtfR) + param[start + 22]*(2*x[7]) +
         raug*(-YftI)*c1 + raug*(-YftR)*c2 + raug*(YtfI)*c3 +
         raug*(YtfR)*c4 + raug*(2*x[7])*c5;

    g7 += (param[start + 30] + raug*c6)*(x[7] / (x[6]*x[6] + x[7]*x[7]));
    g8 += (-((param[start + 30] + raug*c6)*(x[6] / (x[6]*x[6] + x[7]*x[7]))));
    g9 += param[start + 30] + raug*c6;
    g10 += (-(param[start + 30] + raug*c6));

    g[0] = scale*g1;
    g[1] = scale*g2;
    g[2] = scale*g3;
    g[3] = scale*g4;
    g[4] = scale*g5;
    g[5] = scale*g6;
    g[6] = scale*g7;
    g[7] = scale*g8;
    g[8] = scale*g9;
    g[9] = scale*g10;

    return;
}

void eval_h(int I, int n, double scale, double *x, double *A, double *param,
            double YffR, double YffI,
            double YftR, double YftI,
            double YttR, double YttI,
            double YtfR, double YtfI)
{
    double alrect, altheta, raug, c5, c6, x6sq_plus_x7sq;
    int start = 31*I;

    alrect = param[start + 22];
    altheta = param[start + 30];
    raug = param[start + 23];
    c5 = (x[6]*x[6] + x[7]*x[7] - x[4]*x[5]);
    c6 = (x[8] - x[9] - atan2(x[7], x[6]));
    x6sq_plus_x7sq = x[6]*x[6] + x[7]*x[7];

    // 1st column
    A[0] = scale*(param[start + 6] + raug); // A[1,1]
    A[4] = scale*(raug*(-YffR)); // A[5,1]
    A[6] = scale*(raug*(-YftR)); // A[7,1]
    A[7] = scale*(raug*(-YftI)); // A[8,1]

    // 2nd columns
    A[n + 1] = scale*(param[start + 7] + raug); // A[2,2]
    A[n + 4] = scale*(raug*(YffI));  // A[5,2]
    A[n + 6] = scale*(raug*(YftI));  // A[7,2]
    A[n + 7] = scale*(raug*(-YftR)); // A[8,2]

    // 3rd column
    A[n*2 + 2] = scale*(param[start + 8] + raug); // A[3,3]
    A[n*2 + 5] = scale*(raug*(-YttR)); // A[6,3]
    A[n*2 + 6] = scale*(raug*(-YtfR)); // A[7,3]
    A[n*2 + 7] = scale*(raug*(YtfI));  // A[8,3]

    // 4th column
    A[n*3 + 3] = scale*(param[start + 9] + raug); // A[4,4]
    A[n*3 + 5] = scale*(raug*(YttI)); // A[6,4]
    A[n*3 + 6] = scale*(raug*(YtfI)); // A[7,4]
    A[n*3 + 7] = scale*(raug*(YtfR)); // A[8,4]

    // 5th column
    A[n*4 + 4] = scale*(param[start + 10] + raug*(YffR*YffR) + raug*(YffI*YffI) + raug*(x[5]*x[5])); // A[5,5]
    A[n*4 + 5] = scale*(-(alrect + raug*c5) + raug*(x[4]*x[5])); // A[6,5]
    A[n*4 + 6] = scale*(raug*(YffR*YftR) + raug*(YffI*YftI) + raug*((-x[5])*(2*x[6]))); // A[7,5]
    A[n*4 + 7] = scale*(raug*(YffR*YftI) + raug*(YffI*(-YftR)) + raug*((-x[5])*(2*x[7]))); // A[8,5]

    // 6th column
    A[n*5 + 5] = scale*(param[start + 11] + raug*(YttR*YttR) + raug*(YttI*YttI) + raug*(x[4]*x[4])); // A[6,6]
    A[n*5 + 6] = scale*(raug*(YttR*YtfR) + raug*(YttI*YtfI) + raug*((-x[4])*(2*x[6]))); // A[7,6]
    A[n*5 + 7] = scale*(raug*((-YttR)*YtfI) + raug*(YttI*YtfR) + raug*((-x[4])*(2*x[7]))); // A[8,6]

    // 7th column
    A[n*6 + 6] = scale*((alrect + raug*c5)*2 + raug*(YftR*YftR) + raug*(YftI*YftI) +
                 raug*(YtfR*YtfR) + raug*(YtfI*YtfI) + raug*((2*x[6])*(2*x[6]))); // A[7,7]
    A[n*6 + 6] += scale*((altheta + raug*c6)*((-2*x[6]*x[7]) / (x6sq_plus_x7sq * x6sq_plus_x7sq)));
    A[n*6 + 6] += scale*(raug*((x[7] / x6sq_plus_x7sq)*(x[7] / x6sq_plus_x7sq)));
    A[n*6 + 7] = scale*(raug*(YftR*YftI) + raug*(YftI*(-YftR)) + raug*((-YtfR)*YtfI) +
                 raug*(YtfI*YtfR) + raug*((2*x[6])*(2*x[7]))); // A[8,7]
    A[n*6 + 7] += scale*((altheta + raug*c6)*((x[6]*x[6] - x[7]*x[7])/(x6sq_plus_x7sq * x6sq_plus_x7sq)));
    A[n*6 + 7] += scale*(raug*((x[7]/x6sq_plus_x7sq)*(-x[6]/x6sq_plus_x7sq)));
    A[n*6 + 8] = scale*(raug*(x[7]/x6sq_plus_x7sq)); // A[9,7]
    A[n*6 + 9] = scale*(-raug*(x[7]/x6sq_plus_x7sq)); // A[10,7]

    // 8th column
    A[n*7 + 7] = scale*((alrect + raug*c5)*2 + raug*(YftI*YftI) + raug*(YftR*YftR) +
                 raug*(YtfI*YtfI) + raug*(YtfR*YtfR) + raug*((2*x[7])*(2*x[7]))); // A[8,8]
    A[n*7 + 7] += scale*((altheta + raug*c6)*((2*x[6]*x[7]) / (x6sq_plus_x7sq * x6sq_plus_x7sq)));
    A[n*7 + 7] += scale*(raug*((-x[6]/x6sq_plus_x7sq)*(-x[6]/x6sq_plus_x7sq)));
    A[n*7 + 8] = scale*(raug*(-x[6]/x6sq_plus_x7sq)); // A[9,8]
    A[n*7 + 9] = scale*(-raug*(-x[6]/x6sq_plus_x7sq)); // A[10,8]

    // 9th column
    A[n*8 + 8] = scale*(param[start + 26] + raug); // A[9,9]
    A[n*8 + 9] = scale*(-raug); // A[10,9]

    // 10th column
    A[n*9 + 9] = scale*(param[start + 27] + raug); // A[10,10]

    for (int j = 0; j < n; j++) {
        for (int k = 0; k < j; k++) {
            A[n*j + k] = A[n*k + j];
        }
    }

    return;
}

void auglag(int nbranches, int n, int major_iter, int max_auglag,
            int line_start, double scale, double mu_max,
            double *u_curr, double *v_curr, double *l_curr,
            double *rho, double *wRIij, double *param,
            double *_YffR, double *_YffI,
            double *_YftR, double *_YftI,
            double *_YttR, double *_YttI,
            double *_YtfR, double *_YtfI,
            double *frBound, double *toBound)
{
    bool terminate;
    int max_feval, max_minor, it, status, minor_iter;
    double mu, eta, omega;
    double YffR, YffI, YftR, YftI, YttR, YttI, YtfR, YtfI;

    double *x, *xl, *xu;

    x = (double *)calloc(n, sizeof(double));
    xl = (double *)calloc(n, sizeof(double));
    xu = (double *)calloc(n, sizeof(double));

    for (int I = 0; I < nbranches; I++) {
        int pij_idx = line_start + 8*I;

        for (int j = 0; j < n; j++) {
            xl[j] = -INF;
            xu[j] = INF;
        }

        xl[4] = frBound[2*I];
        xu[4] = frBound[2*I+1];
        xl[5] = toBound[2*I];
        xu[5] = toBound[2*I+1];
        xl[8] = -2*M_PI;
        xu[8] = 2*M_PI;
        xl[9] = -2*M_PI;
        xu[9] = 2*M_PI;

        x[0] = u_curr[pij_idx];
        x[1] = u_curr[pij_idx+1];
        x[2] = u_curr[pij_idx+2];
        x[3] = u_curr[pij_idx+3];
        x[4] = min(xu[4], max(xl[4], u_curr[pij_idx+4]));
        x[5] = min(xu[5], max(xl[5], u_curr[pij_idx+5]));
        x[6] = wRIij[2*I];
        x[7] = wRIij[2*I+1];
        x[8] = min(xu[8], max(xl[8], u_curr[pij_idx+6]));
        x[9] = min(xu[9], max(xl[9], u_curr[pij_idx+7]));

        YffR = _YffR[I]; YffI = _YffI[I];
        YftR = _YftR[I]; YftI = _YftI[I];
        YttR = _YttR[I]; YttI = _YttI[I];
        YtfR = _YtfR[I]; YtfI = _YtfI[I];

        int start = 31*I;
        param[start] = l_curr[pij_idx];
        param[start + 1] = l_curr[pij_idx+1];
        param[start + 2] = l_curr[pij_idx+2];
        param[start + 3] = l_curr[pij_idx+3];
        param[start + 4] = l_curr[pij_idx+4];
        param[start + 5] = l_curr[pij_idx+5];
        param[start + 6] = rho[pij_idx];
        param[start + 7] = rho[pij_idx+1];
        param[start + 8] = rho[pij_idx+2];
        param[start + 9] = rho[pij_idx+3];
        param[start + 10] = rho[pij_idx+4];
        param[start + 11] = rho[pij_idx+5];
        param[start + 12] = v_curr[pij_idx];
        param[start + 13] = v_curr[pij_idx+1];
        param[start + 14] = v_curr[pij_idx+2];
        param[start + 15] = v_curr[pij_idx+3];
        param[start + 16] = v_curr[pij_idx+4];
        param[start + 17] = v_curr[pij_idx+5];

        param[start + 24] = l_curr[pij_idx+6];
        param[start + 25] = l_curr[pij_idx+7];
        param[start + 26] = rho[pij_idx+6];
        param[start + 27] = rho[pij_idx+7];
        param[start + 28] = v_curr[pij_idx+6];
        param[start + 29] = v_curr[pij_idx+7];

        if (major_iter == 1) {
            param[start + 23] = 10.0;
            mu = 10.0;
        } else {
            mu = param[start + 23];
        }

        eta = 1 / pow(mu, 0.1);
        omega = 1 / mu;
        max_feval = 500;
        max_minor = 200;

        it = 0;
        terminate = false;

        double cviol1, cviol2, cviol3, cviol4, cviol5, cviol6, cnorm;
        while (!terminate) {
            it += 1;

            // Solve the branch problem.

            driver_auglag(I, n, max_feval, max_minor, &status, &minor_iter, scale,
                          x, xl, xu, param, YffR, YffI, YftR, YftI, YttR, YttI, YtfR, YtfI,
                          &eval_f, &eval_grad_f, &eval_h);

            // Check the termination condition.
            cviol1 = x[0] - (YffR*x[4] + YftR*x[6] + YftI*x[7]);
            cviol2 = x[1] - (-YffI*x[4] - YftI*x[6] + YftR*x[7]);
            cviol3 = x[2] - (YttR*x[5] + YtfR*x[6] - YtfI*x[7]);
            cviol4 = x[3] - (-YttI*x[5] - YtfI*x[6] - YtfR*x[7]);
            cviol5 = x[6]*x[6] + x[7]*x[7] - x[4]*x[5];
            cviol6 = x[8] - x[9] - atan2(x[7], x[6]);

            cnorm = max(abs(cviol1), max(abs(cviol2), max(abs(cviol3), max(abs(cviol4), max(abs(cviol5), abs(cviol6))))));

            if (cnorm <= eta) {
                if (cnorm <= 1e-6) {
                    terminate = true;
                } else {
                    param[start + 18] += mu*cviol1;
                    param[start + 19] += mu*cviol2;
                    param[start + 20] += mu*cviol3;
                    param[start + 21] += mu*cviol4;
                    param[start + 22] += mu*cviol5;
                    param[start + 30] += mu*cviol6;

                    eta = eta / pow(mu, 0.9);
                    omega  = omega / mu;
                }
            } else {
                mu = min(mu_max, mu*10);
                eta = 1 / pow(mu, 0.1);
                omega = 1 / mu;
                param[start + 23] = mu;
            }

            if (it >= max_auglag) {
                terminate = true;
            }
        }

        u_curr[pij_idx] = x[0];
        u_curr[pij_idx+1] = x[1];
        u_curr[pij_idx+2] = x[2];
        u_curr[pij_idx+3] = x[3];
        u_curr[pij_idx+4] = x[4];
        u_curr[pij_idx+5] = x[5];
        u_curr[pij_idx+6] = x[8];
        u_curr[pij_idx+7] = x[9];
        wRIij[2*I] = x[6];
        wRIij[2*I+1] = x[7];
        param[start + 23] = mu;
    }

    free(x);
    free(xl);
    free(xu);

    return;
}
