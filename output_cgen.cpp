#include <iostream>
#include <stdlib.h>
#include <math.h>

// Necessary to declare constant for OPS.
double dx, dy, dt;
// Step for each dimension and time.
double nx, ny;

#define PREVIOUS 0
#define CURRENT 1
#define NEXT 2
#define M_PI 3.14159265358979323846
#define BORDER_SIZE 20
#define CPML 20
#define SOURCE_LOCATION_X 0
#define SOURCE_LOCATION_Y 0
#define RICKER_PEAK_FREQUENCY 5
#define MAX_VELOCITY 5000

#define OPS_2D
#include <ops_seq.h>
#include <ops_lib_cpp.h>
#include <sources.cpp>
#include <velocity-model.cpp>
#include <wave-propagation-ops.h>

using namespace std;

void wavePropagation(double *u_new, double const *u_current, double const *u_previous, double const *velocity, double const *wx, double const *wy, double const *zetax, double const *zetay, int const *idx)
{
    // Propagates the wave.
    u_new[OPS_ACC0(0, 0)] = (velocity[OPS_ACC3(0, 0)] * velocity[OPS_ACC3(0, 0)]) *
                                ((-205.0 / 72) * u_current[OPS_ACC1(0, 0)] +
                                 (8.0 / 5) * (u_current[OPS_ACC1(1, 0)] + u_current[OPS_ACC1(-1, 0)] + u_current[OPS_ACC1(0, 1)] + u_current[OPS_ACC1(0, -1)]) +
                                 (-0.2) * (u_current[OPS_ACC1(2, 0)] + u_current[OPS_ACC1(-2, 0)] + u_current[OPS_ACC1(0, 2)] + u_current[OPS_ACC1(0, -2)]) +
                                 (8.0 / 315) * (u_current[OPS_ACC1(3, 0)] + u_current[OPS_ACC1(-3, 0)] + u_current[OPS_ACC1(0, 3)] + u_current[OPS_ACC1(0, -3)]) +
                                 (-1.0 / 560.0) * (u_current[OPS_ACC1(4, 0)] + u_current[OPS_ACC1(-4, 0)] + u_current[OPS_ACC1(0, 4)] + u_current[OPS_ACC1(0, -4)])) +
                            2 * u_current[OPS_ACC1(0, 0)] - u_previous[OPS_ACC2(0, 0)];

    // Apply absorbent border conditions.
    u_new[OPS_ACC0(0, 0)] += velocity[OPS_ACC3(0, 0)] * velocity[OPS_ACC3(0, 0)] * dx * (wx[OPS_ACC4(1, 0)] - wx[OPS_ACC4(-1, 0)] + wy[OPS_ACC5(0, 1)] - wy[OPS_ACC5(0, -1)]) * 0.5 +
                             velocity[OPS_ACC3(0, 0)] * velocity[OPS_ACC3(0, 0)] * dx * dy * (zetax[OPS_ACC6(0, 0)] + zetay[OPS_ACC7(0, 0)]);

    // printf("p1=%lf | p2=%lf | p3=%lf | v=%lf | dx=%lf | wst=%lf | zetast=%lf | dx=%lf | dy=%lf | idx=%d %d\n", u_previous[OPS_ACC2(0, 0)], u_current[OPS_ACC1(0, 0)], u_new[OPS_ACC0(0, 0)], velocity[OPS_ACC3(0, 0)], dx, (wx[OPS_ACC4(1, 0)] - wx[OPS_ACC4(-1, 0)] + wy[OPS_ACC5(0, 1)] - wy[OPS_ACC5(0, -1)]), (zetax[OPS_ACC6(0, 0)] + zetay[OPS_ACC7(0, 0)]), dx, dy, idx[0], idx[1]);
}

/**
 * @brief  Injects the source into the wavefield.
 * @note
 * @param  *u_new: Field to be injected.
 * @param  *value: Value of the amplitude of the source to be injected.
 * @retval None
 */
void sourceInjection(double *u_new, double const *value)
{
    u_new[OPS_ACC0(0, 0)] += *value;
}

/**
 * @brief  Updates omega CPML parameter for X and Y dimension.
 * @note
 * @param  *wx: Omega X parameter.
 * @param  *wy: Omega Y parameter.
 * @param  *ax: A X auxiliar parameter.
 * @param  *ay: A Y auxiliar parameter.
 * @param  *bx: B X auxiliar parameter.
 * @param  *by: B Y auxiliar parameter.
 * @param  *u_next: Wave grid.
 * @retval None
 */
void update_omega(double *wx, double *wy, double const *ax, double const *ay, double const *bx, double const *by, double const *u_next, int const *idx)
{
    wx[OPS_ACC0(0, 0)] = bx[OPS_ACC4(0, 0)] * wx[OPS_ACC0(0, 0)] +
                         ax[OPS_ACC2(0, 0)] * (u_next[OPS_ACC6(1, 0)] - u_next[OPS_ACC6(-1, 0)]) / (2 * dx);
    wy[OPS_ACC1(0, 0)] = by[OPS_ACC5(0, 0)] * wy[OPS_ACC1(0, 0)] +
                         ay[OPS_ACC3(0, 0)] * (u_next[OPS_ACC6(0, 1)] - u_next[OPS_ACC6(0, -1)]) / (2 * dy);

    // printf("wy=%lf | by=%lf | ay=%lf | idx = %d %d\n", wy[OPS_ACC1(0, 0)], by[OPS_ACC5(0, 0)], ay[OPS_ACC3(0, 0)], idx[0], idx[1]);
}

/**
 * @brief  Updates Zeta CPML parameter for X and Y dimension.
 * @note
 * @param  *zetax: Zeta X parameter.
 * @param  *zetay: Zeta Y parameter.
 * @param  *wx: Omega X parameter.
 * @param  *wy: Omega Y parameter.
 * @param  *ax: A X auxiliar parameter.
 * @param  *ay: A Y auxiliar parameter.
 * @param  *bx: B X auxiliar parameter.
 * @param  *by: B Y auxiliar parameter.
 * @param  *u_next: Wave grid.
 * @retval None
 */
void update_zeta(double *zetax, double *zetay, double const *wx, double const *wy, double const *ax, double const *ay, double const *bx, double const *by, double const *u_next, int const *idx)
{
    zetax[OPS_ACC0(0, 0)] = bx[OPS_ACC6(0, 0)] * zetax[OPS_ACC0(0, 0)] +
                            ax[OPS_ACC4(0, 0)] * (u_next[OPS_ACC8(1, 0)] - (2 * u_next[OPS_ACC8(0, 0)]) + u_next[OPS_ACC8(-1, 0)]) / (dx * dx) +
                            ax[OPS_ACC4(0, 0)] * (wx[OPS_ACC2(1, 0)] - wx[OPS_ACC2(-1, 0)]) / (2 * dx);

    // printf("zetax=%lf | bx=%lf | ax=%lf | h=%lf | idx = %d %d\n", zetax[OPS_ACC0(0, 0)], bx[OPS_ACC6(0, 0)], ax[OPS_ACC4(0, 0)], dx, idx[0], idx[1]);

    zetay[OPS_ACC1(0, 0)] = by[OPS_ACC7(0, 0)] * zetay[OPS_ACC1(0, 0)] +
                            ay[OPS_ACC5(0, 0)] * (u_next[OPS_ACC8(0, 1)] - (2 * u_next[OPS_ACC8(0, 0)]) + u_next[OPS_ACC8(0, -1)]) / (dy * dy) +
                            ay[OPS_ACC5(0, 0)] * (wy[OPS_ACC3(0, 1)] - wy[OPS_ACC3(0, -1)]) / (2 * dy);

    // printf("zetay=%lf | by=%lf | ay=%lf | h=%lf | idx = %d %d\n", zetay[OPS_ACC1(0, 0)], by[OPS_ACC7(0, 0)], ay[OPS_ACC5(0, 0)], dy, idx[0], idx[1]);
}

/**
 * @brief  Initializes auxiliar CPML variables in X dimension.
 * @note
 * @param  *ax:
 * @param  *bx:
 * @param  *idx: Index being calculated.
 * @retval None
 */
void initializeAxBx(double *ax, double *bx, int const *idx)
{
    // Local variables
    double d_x, b_x, f_x, alfa_x;
    double Lx;

    Lx = CPML * dx;
    // This sould be outside the kernel.

    if (idx[0] < 0)
      f_x = Lx - (idx[0] + CPML) * dx;
    else
      if (idx[0] > nx - 2 * CPML)
        f_x = dx * (idx[0] - nx + 2 * CPML);
      else
        f_x = 0;

    alfa_x = M_PI * RICKER_PEAK_FREQUENCY * (Lx - f_x) / Lx;
    d_x = -3.0 * log(0.00001) / (2.0 * Lx) * MAX_VELOCITY * (f_x / Lx) * (f_x / Lx);
    b_x = exp(-dt * (d_x + alfa_x));

    ax[OPS_ACC0(0, 0)] = d_x * (b_x - 1) / (d_x + alfa_x);
    bx[OPS_ACC1(0, 0)] = exp(-dt * (d_x + alfa_x));
}

/**
 * @brief  Initializes auxiliar CPML variables in Y dimension.
 * @note
 * @param  *ay:
 * @param  *by:
 * @param  *idx:
 * @retval None
 */
void initializeAyBy(double *ay, double *by, int const *idx)
{
    // Local variables
    double d_y, b_y, f_y, alfa_y;
    double Ly;

    Ly = CPML * dy;
    // This sould be outside the kernel.

    if (idx[0] < 0)
    {
        f_y = Ly - (idx[0] + CPML) * dy;
    }
    else
      if (idx[0] > ny - 2 * CPML)
      {
          f_y = dy * (idx[0] - ny + 2 * CPML);
      }
      else
      {
          f_y = 0.0;
      }

    alfa_y = M_PI * RICKER_PEAK_FREQUENCY * (Ly - f_y) / Ly;
    d_y = -3.0 * log(0.00001) / (2.0 * Ly) * MAX_VELOCITY * (f_y / Ly) * (f_y / Ly);
    b_y = exp(-dt * (d_y + alfa_y));

    ay[OPS_ACC0(0, 0)] = d_y * (b_y - 1) / (d_y + alfa_y);
    by[OPS_ACC1(0, 0)] = exp(-dt * (d_y + alfa_y));
}

/**
 * @brief  Copy data from an array to another.
 * @note
 * @param  *to: Array will be copied here.
 * @param  *from: Array to be copied.
 * @retval None
 */
void makeCopy(double *to, double const *from)
{
    to[OPS_ACC0(0, 0)] = from[OPS_ACC1(0, 0)];
}

void save_array(char filename[], double *array, int size)
{
    FILE *fp;
    int i;

    fp = fopen(filename, "w");

    for (i = 0; i < size; i++)
    {
        fprintf(fp, "%lf ", array[i]);
    }
    fflush(fp);
    fclose(fp);
}

/**
 * @brief  Entry point for wave propagation example using OPS.
 * @note
 * @param  argc: Number of input arguments. Expected 7.
 * @param  *argv[]: xSize:
 *                  ySize:
 *                  xIntervals:
 *                  yIntervals:
 *                  tTotal:
 *                  tIntervals:
 * @retval Positive number for sucess. Negative for error.
 */
int propagationCPML(int xSize, int ySize, int xIntervals, int yIntervals, int tTotal, int tIntervals)
{
    // Local variables
    double **u, *source;
    // u represents the grid data, c the velocity for each location in space and source is a vector of the source magnitude in time.
    int i, j, t;
    // Counters.

    // OPS variables
    int size[2];
    // Size in each dimension.
    int base[] = {0, 0};
    // Base indices for each dimension of the block.
    int d_m[] = {-BORDER_SIZE, -BORDER_SIZE};
    // d_m - padding from the face in the negative direction for each dimension. (Halo)
    int d_p[] = {BORDER_SIZE, BORDER_SIZE};
    // d_p - padding from the face in the positive direction for each dimension. (Halo)
    ops_block wave_grid;
    // Block struct variable.
    ops_dat d_u_Previous, d_u_Current, d_u_Next, d_velocity;
    // Block data.

    // CPML variables
    double *wx, *wy, *zetax, *zetay;
    double *ax, *bx, *ay, *by;

    // OPS variables for CPML conditions.
    int s1d_0[] = {0};
    int s2d_00[] = {0, 0};
    int s2d_8th[] = {0, 0, 1, 0, -1, 0, 2, 0, -2, 0, 3, 0, -3, 0, 4, 0, -4, 0, 0, 1, 0, -1, 0, 2, 0, -2, 0, 3, 0, -3, 0, 4, 0, -4};
    int s2d_2px[] = {1, 0, -1, 0};
    int s2d_2py[] = {0, 1, 0, -1};
    int s2d_4pt[] = {1, 0, -1, 0, 0, 1, 0, -1};
    int s2d_5pt[] = {1, 0, 0, 0, -1, 0, 0, 1, 0, -1};
    int stride2d_x[] = {1, 0};
    // Stride only in X dimension.
    int stride2d_y[] = {0, 1};
    // Stride only in Y dimension.
    ops_dat d_wx, d_wy, d_zetax, d_zetay;
    ops_dat d_ax, d_ay, d_bx, d_by;
    ops_block wave_grid_1d;

    // Initialize OPS.
    ops_init(1, NULL, 1);

    nx = xIntervals + 2 * BORDER_SIZE;
    ny = yIntervals + 2 * BORDER_SIZE;

    // Set grid dimension.
    size[0] = xIntervals;
    size[1] = yIntervals;

    // Calculate size of each step for dimension and time.
    dx = (double)xSize / (double)xIntervals;
    dy = (double)ySize / (double)yIntervals;
    dt = (double)tTotal / (double)tIntervals;

    printf("dx=%lf | dy=%lf | dt=%lf\n", dx, dy, dt);
    printf("nx=%d | ny=%d | nT=%d | borders=%d\n", xIntervals, yIntervals, tIntervals, BORDER_SIZE);

    // Check CFL convergency conditions.
    if (dt / dx > 1 && dt / dy > 1)
    {
        cout << "Does not comply with CFL conditions." << endl;
        return -1;
    }

    ops_decl_const2("dx", 1, "double", &dx);
    ops_decl_const2("dy", 1, "double", &dy);
    ops_decl_const2("dt", 1, "double", &dt);
    ops_decl_const2("nx", 1, "double", &nx);
    ops_decl_const2("ny", 1, "double", &ny);

    int range_CPML[] = {1 - BORDER_SIZE, xIntervals + BORDER_SIZE - 1, 1 - BORDER_SIZE, yIntervals + BORDER_SIZE - 1};
    int range_CPML_aux_1D_X_init[] = {-BORDER_SIZE, xIntervals + BORDER_SIZE - 1, 0, 1};
    int range_CPML_aux_1D_Y_init[] = {0, 1, -BORDER_SIZE, yIntervals + BORDER_SIZE - 1};
    int range[] = {-BORDER_SIZE + 4, xIntervals + BORDER_SIZE - 4 - 1, -BORDER_SIZE + 4, yIntervals + BORDER_SIZE - 4 - 1};

    int whole_range[] = {-BORDER_SIZE, xIntervals + BORDER_SIZE, -BORDER_SIZE, yIntervals + BORDER_SIZE};

    // Allocate grids. 3 x 2D space represented in 1 dimension.
    u = (double **)malloc(3 * sizeof(double *));
    for (t = 0; t < 3; t++)
    {
        u[t] = (double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double));
    }

    // Allocate source time vector.
    source = rickerSource(dt, tIntervals, RICKER_PEAK_FREQUENCY);

    // save_array("output/source-ops.txt", source, tIntervals);

    // Allocates CPML Variables.
    wx = (double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double));
    wy = (double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double));
    zetax = (double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double));
    zetay = (double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double));
    ax = (double *)calloc((xIntervals + 2 * BORDER_SIZE), sizeof(double));
    bx = (double *)calloc((xIntervals + 2 * BORDER_SIZE), sizeof(double));
    ay = (double *)calloc((yIntervals + 2 * BORDER_SIZE), sizeof(double));
    by = (double *)calloc((yIntervals + 2 * BORDER_SIZE), sizeof(double));

    // Defines structured block.
    wave_grid = ops_decl_block(2, "wave_grid");
    wave_grid_1d = ops_decl_block(1, "wave_grid_1d");

    // Declares stencil
    ops_stencil S1D_0 = ops_decl_stencil(1, 1, s1d_0, "0");
    ops_stencil S2D_00 = ops_decl_stencil(2, 1, s2d_00, "0,0");
    ops_stencil S2D_8TH = ops_decl_stencil(2, 17, s2d_8th, "8th");
    ops_stencil S2D_2PX = ops_decl_stencil(2, 2, s2d_2px, "2px");
    ops_stencil S2D_2PY = ops_decl_stencil(2, 2, s2d_2py, "2py");
    ops_stencil S2D_4PT = ops_decl_stencil(2, 4, s2d_4pt, "4pt");
    ops_stencil S2D_5PT = ops_decl_stencil(2, 5, s2d_5pt, "5pt");
    ops_stencil S2D_00_STRIDE_X = ops_decl_strided_stencil(2, 1, s2d_00, stride2d_x, "s2D_00_stride2D_x");
    ops_stencil S2D_00_STRIDE_Y = ops_decl_strided_stencil(2, 1, s2d_00, stride2d_y, "s2D_00_stride2D_y");

    // Set dataset.
    d_m[0] = -BORDER_SIZE;
    d_m[1] = -BORDER_SIZE;
    d_p[0] = BORDER_SIZE;
    d_p[1] = BORDER_SIZE;
    d_u_Previous = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, u[PREVIOUS], "double", "u_previous");
    d_u_Current = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, u[CURRENT], "double", "u_current");
    d_u_Next = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, u[NEXT], "double", "u_next");

    // Declare CPML datasets.
    size[0] = xIntervals;
    size[1] = yIntervals;
    d_wx = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, wx, "double", "wx");
    d_wy = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, wy, "double", "wy");
    d_zetax = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, zetax, "double", "zetax");
    d_zetay = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, zetay, "double", "zetay");

    // Declare auxiliar CPML datasets
    d_m[0] = -BORDER_SIZE;
    d_m[1] = 0;
    d_p[0] = BORDER_SIZE;
    d_p[1] = 0;
    size[0] = xIntervals;
    size[1] = 1;
    d_ax = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, ax, "double", "d_ax");
    d_bx = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, bx, "double", "d_bx");
    d_m[0] = 0;
    d_m[1] = -BORDER_SIZE;
    d_p[0] = 0;
    d_p[1] = BORDER_SIZE;
    size[0] = 1;
    size[1] = yIntervals;
    d_ay = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, ay, "double", "d_ay");
    d_by = ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, by, "double", "d_by");

    int source_location_range[] = {SOURCE_LOCATION_X, SOURCE_LOCATION_X + 1, SOURCE_LOCATION_Y, SOURCE_LOCATION_Y + 1};
    double source_field_factor = getVelocityAtSourcePoint(FourHorizontalLayers, SOURCE_LOCATION_X, SOURCE_LOCATION_Y, xIntervals, yIntervals) * dt;
    source_field_factor = source_field_factor * source_field_factor / (dx * dy);
    double source_magnitude = 0;

    d_velocity = initializeFourHorizontalLayersModel(xIntervals, yIntervals, BORDER_SIZE, dt, dx);

    // ops_print_dat_to_txtfile(d_velocity, "output/velocity-ops.txt");

    // Initialize CPML auxiliar variables.
    ops_par_loop(initializeAxBx, "initializeAxBx", wave_grid, 2, range_CPML_aux_1D_X_init,
                 ops_arg_dat(d_ax, 1, S2D_00, "double", OPS_WRITE),
                 ops_arg_dat(d_bx, 1, S2D_00, "double", OPS_WRITE),
                 ops_arg_idx());

    ops_par_loop(initializeAyBy, "initializeAyBy", wave_grid, 2, range_CPML_aux_1D_X_init,
                 ops_arg_dat(d_ay, 1, S2D_00, "double", OPS_WRITE),
                 ops_arg_dat(d_by, 1, S2D_00, "double", OPS_WRITE),
                 ops_arg_idx());

    // ops_print_dat_to_txtfile(d_ax, "output/ax-ops.txt");
    // ops_print_dat_to_txtfile(d_bx, "output/bx-ops.txt");
    // ops_print_dat_to_txtfile(d_ay, "output/ay-ops.txt");
    // ops_print_dat_to_txtfile(d_by, "output/by-ops.txt");

    t = 0;
    do
    {
        // Propagates the wave.
        ops_par_loop(wavePropagation, "wave_propagation", wave_grid, 2, range,
                     ops_arg_dat(d_u_Next, 1, S2D_00, "double", OPS_RW),
                     ops_arg_dat(d_u_Current, 1, S2D_8TH, "double", OPS_READ),
                     ops_arg_dat(d_u_Previous, 1, S2D_00, "double", OPS_READ),
                     ops_arg_dat(d_velocity, 1, S2D_00, "double", OPS_READ),
                     ops_arg_dat(d_wx, 1, S2D_2PX, "double", OPS_READ),
                     ops_arg_dat(d_wy, 1, S2D_2PY, "double", OPS_READ),
                     ops_arg_dat(d_zetax, 1, S2D_00, "double", OPS_READ),
                     ops_arg_dat(d_zetay, 1, S2D_00, "double", OPS_READ),
                     ops_arg_idx());

        source_magnitude = source_field_factor * source[t];

        // Injects the source.
        ops_par_loop(sourceInjection, "source_injection", wave_grid, 2, source_location_range,
                     ops_arg_dat(d_u_Next, 1, S2D_00, "double", OPS_WRITE),
                     ops_arg_gbl(&source_magnitude, 1, "double", OPS_READ));

        // Updates Omega X and Y
        ops_par_loop(update_omega, "update_omega", wave_grid, 2, range_CPML,
                     ops_arg_dat(d_wx, 1, S2D_00, "double", OPS_RW),
                     ops_arg_dat(d_wy, 1, S2D_00, "double", OPS_RW),
                     ops_arg_dat(d_ax, 1, S2D_00_STRIDE_X, "double", OPS_READ),
                     ops_arg_dat(d_ay, 1, S2D_00_STRIDE_Y, "double", OPS_READ),
                     ops_arg_dat(d_bx, 1, S2D_00_STRIDE_X, "double", OPS_READ),
                     ops_arg_dat(d_by, 1, S2D_00_STRIDE_Y, "double", OPS_READ),
                     ops_arg_dat(d_u_Next, 1, S2D_4PT, "double", OPS_READ),
                     ops_arg_idx());

        // ops_print_dat_to_txtfile(d_wx, "output/omegax-ops.txt");

        // Updates Zeta X and Y
        ops_par_loop(update_zeta, "update_zeta", wave_grid, 2, range_CPML,
                     ops_arg_dat(d_zetax, 1, S2D_00, "double", OPS_WRITE),
                     ops_arg_dat(d_zetay, 1, S2D_00, "double", OPS_WRITE),
                     ops_arg_dat(d_wx, 1, S2D_2PX, "double", OPS_READ),
                     ops_arg_dat(d_wy, 1, S2D_2PY, "double", OPS_READ),
                     ops_arg_dat(d_ax, 1, S2D_00_STRIDE_X, "double", OPS_READ),
                     ops_arg_dat(d_ay, 1, S2D_00_STRIDE_Y, "double", OPS_READ),
                     ops_arg_dat(d_bx, 1, S2D_00_STRIDE_X, "double", OPS_READ),
                     ops_arg_dat(d_by, 1, S2D_00_STRIDE_Y, "double", OPS_READ),
                     ops_arg_dat(d_u_Next, 1, S2D_5PT, "double", OPS_READ),
                     ops_arg_idx());

        // Save  dat to file.
        // if (t == 0)
        // {
        //     ops_print_dat_to_txtfile(d_zetax, "output/zetax-ops.txt");
        // }

        // ops_print_dat_to_txtfile(d_u_Current, "output/ops-output.txt");

        // Transfer current to previous.
        ops_par_loop(makeCopy, "copy_current_to_previous", wave_grid, 2, whole_range,
                     ops_arg_dat(d_u_Previous, 1, S2D_00, "double", OPS_WRITE),
                     ops_arg_dat(d_u_Current, 1, S2D_00, "double", OPS_READ));

        // Transfer next to current.
        ops_par_loop(makeCopy, "copy_next_to_current", wave_grid, 2, whole_range,
                     ops_arg_dat(d_u_Current, 1, S2D_00, "double", OPS_WRITE),
                     ops_arg_dat(d_u_Next, 1, S2D_00, "double", OPS_READ));

        // Save  dat to file.
        // if (t % 1000 == 0)
        //     ops_print_dat_to_txtfile(d_u_Current, "output/u-ops.txt");

        // printf("%d\n", t);
        t++;
    }
    while (t < tIntervals);

    ops_print_dat_to_txtfile(d_u_Current, "output/u-ops.txt");

    ops_exit();

    // free(u[PREVIOUS]);
    // free(u[CURRENT]);
    // free(u[NEXT]);
    // free(source);
}
