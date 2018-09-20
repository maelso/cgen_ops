import cgen as c

code = c.Module([
    c.Include('iostream'),
    c.Include('stdlib.h'),
    c.Include('math.h'),
    c.Line(),

    c.LineComment('Necessary to declare constant for OPS.'),
    c.Value('double', 'dx, dy, dt'),
    c.LineComment('Step for each dimension and time.'),
    c.Value('double', 'nx, ny'),
    c.Line(),
    
    c.Define('PREVIOUS', '0'),
    c.Define('CURRENT', '1'),
    c.Define('NEXT', '2'),
    c.Define('M_PI', '3.14159265358979323846'),
    c.Define('BORDER_SIZE', '20'),
    c.Define('CPML', '20'),
    c.Define('SOURCE_LOCATION_X', '0'),
    c.Define('SOURCE_LOCATION_Y', '0'),
    c.Define('RICKER_PEAK_FREQUENCY', '5'),
    c.Define('MAX_VELOCITY', '5000'),
    c.Line(),

    c.Define('OPS_2D', ''),
    c.Include('ops_seq.h'),
    c.Include('ops_lib_cpp.h'),
    c.Include('sources.cpp'),
    c.Include('velocity-model.cpp'),
    c.Include('wave-propagation-ops.h'),    
    c.Line(),

    c.Statement('using namespace std'),
    c.Line()
    
    ])

code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'wavePropagation'), [c.Pointer(c.Value('double', 'u_new')), c.Const(c.Pointer(c.Value('double', 'u_current'))), c.Const(c.Pointer(c.Value('double', 'u_previous'))), c.Const(c.Pointer(c.Value('double', 'velocity'))), c.Const(c.Pointer(c.Value('double', 'wx'))), c.Const(c.Pointer(c.Value('double', 'wy'))), c.Const(c.Pointer(c.Value('double', 'zetax'))), c.Const(c.Pointer(c.Value('double', 'zetay'))), c.Const(c.Pointer(c.Value('int', 'idx')))]),
    c.Block([
        c.LineComment('Propagates the wave.'),
        c.Assign('u_new[OPS_ACC0(0, 0)]','(velocity[OPS_ACC3(0, 0)] * velocity[OPS_ACC3(0, 0)]) *\n\
                                ((-205.0 / 72) * u_current[OPS_ACC1(0, 0)] +\n\
                                 (8.0 / 5) * (u_current[OPS_ACC1(1, 0)] + u_current[OPS_ACC1(-1, 0)] + u_current[OPS_ACC1(0, 1)] + u_current[OPS_ACC1(0, -1)]) +\n\
                                 (-0.2) * (u_current[OPS_ACC1(2, 0)] + u_current[OPS_ACC1(-2, 0)] + u_current[OPS_ACC1(0, 2)] + u_current[OPS_ACC1(0, -2)]) +\n\
                                 (8.0 / 315) * (u_current[OPS_ACC1(3, 0)] + u_current[OPS_ACC1(-3, 0)] + u_current[OPS_ACC1(0, 3)] + u_current[OPS_ACC1(0, -3)]) +\n\
                                 (-1.0 / 560.0) * (u_current[OPS_ACC1(4, 0)] + u_current[OPS_ACC1(-4, 0)] + u_current[OPS_ACC1(0, 4)] + u_current[OPS_ACC1(0, -4)])) +\n\
                            2 * u_current[OPS_ACC1(0, 0)] - u_previous[OPS_ACC2(0, 0)]'),
        c.Line(),

        c.LineComment('Apply absorbent border conditions.'),
        c.Increment('u_new[OPS_ACC0(0, 0)]', 'velocity[OPS_ACC3(0, 0)] * velocity[OPS_ACC3(0, 0)] * dx * (wx[OPS_ACC4(1, 0)] - wx[OPS_ACC4(-1, 0)] + wy[OPS_ACC5(0, 1)] - wy[OPS_ACC5(0, -1)]) * 0.5 +\n\
                             velocity[OPS_ACC3(0, 0)] * velocity[OPS_ACC3(0, 0)] * dx * dy * (zetax[OPS_ACC6(0, 0)] + zetay[OPS_ACC7(0, 0)])'),
        c.Line(),
        c.LineComment('printf("p1=%lf | p2=%lf | p3=%lf | v=%lf | dx=%lf | wst=%lf | zetast=%lf | dx=%lf | dy=%lf | idx=%d %d\\n", u_previous[OPS_ACC2(0, 0)], u_current[OPS_ACC1(0, 0)], u_new[OPS_ACC0(0, 0)], velocity[OPS_ACC3(0, 0)], dx, (wx[OPS_ACC4(1, 0)] - wx[OPS_ACC4(-1, 0)] + wy[OPS_ACC5(0, 1)] - wy[OPS_ACC5(0, -1)]), (zetax[OPS_ACC6(0, 0)] + zetay[OPS_ACC7(0, 0)]), dx, dy, idx[0], idx[1]);')
    ])
))
code.append(c.Line())

code.append(c.MultilineComment('\
@brief  Injects the source into the wavefield.\n\
@note\n\
@param  *u_new: Field to be injected.\n\
@param  *value: Value of the amplitude of the source to be injected.\n\
@retval None'))
code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'sourceInjection'), [c.Pointer(c.Value('double', 'u_new')), c.Const(c.Pointer(c.Value('double', 'value')))]),
    c.Block([
    c.Increment('u_new[OPS_ACC0(0, 0)]', '*value')
    ])
))
code.append(c.Line())

code.append(c.MultilineComment('\
@brief  Updates omega CPML parameter for X and Y dimension.\n\
@note\n\
@param  *wx: Omega X parameter.\n\
@param  *wy: Omega Y parameter.\n\
@param  *ax: A X auxiliar parameter.\n\
@param  *ay: A Y auxiliar parameter.\n\
@param  *bx: B X auxiliar parameter.\n\
@param  *by: B Y auxiliar parameter.\n\
@param  *u_next: Wave grid.\n\
@retval None'))
code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'update_omega'), [c.Pointer(c.Value('double', 'wx')), c.Pointer(c.Value('double', 'wy')), c.Const(c.Pointer(c.Value('double', 'ax'))), c.Const(c.Pointer(c.Value('double', 'ay'))), c.Const(c.Pointer(c.Value('double', 'bx'))), c.Const(c.Pointer(c.Value('double', 'by'))), c.Const(c.Pointer(c.Value('double', 'u_next'))), c.Const(c.Pointer(c.Value('int', 'idx')))]),
    c.Block([
        c.Assign('wx[OPS_ACC0(0, 0)]', 'bx[OPS_ACC4(0, 0)] * wx[OPS_ACC0(0, 0)] +\n\
                         ax[OPS_ACC2(0, 0)] * (u_next[OPS_ACC6(1, 0)] - u_next[OPS_ACC6(-1, 0)]) / (2 * dx)'),
        c.Assign('wy[OPS_ACC1(0, 0)]', 'by[OPS_ACC5(0, 0)] * wy[OPS_ACC1(0, 0)] +\n\
                         ay[OPS_ACC3(0, 0)] * (u_next[OPS_ACC6(0, 1)] - u_next[OPS_ACC6(0, -1)]) / (2 * dy)'),
        c.Line(),

        c.LineComment('printf("wy=%lf | by=%lf | ay=%lf | idx = %d %d\\n", wy[OPS_ACC1(0, 0)], by[OPS_ACC5(0, 0)], ay[OPS_ACC3(0, 0)], idx[0], idx[1]);')
    ])
))
code.append(c.Line())

code.append(c.MultilineComment('\
@brief  Updates Zeta CPML parameter for X and Y dimension.\n\
@note\n\
@param  *zetax: Zeta X parameter.\n\
@param  *zetay: Zeta Y parameter.\n\
@param  *wx: Omega X parameter.\n\
@param  *wy: Omega Y parameter.\n\
@param  *ax: A X auxiliar parameter.\n\
@param  *ay: A Y auxiliar parameter.\n\
@param  *bx: B X auxiliar parameter.\n\
@param  *by: B Y auxiliar parameter.\n\
@param  *u_next: Wave grid.\n\
@retval None'))
code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'update_zeta'), [c.Pointer(c.Value('double', 'zetax')), c.Pointer(c.Value('double', 'zetay')), c.Const(c.Pointer(c.Value('double', 'wx'))), c.Const(c.Pointer(c.Value('double', 'wy'))), c.Const(c.Pointer(c.Value('double', 'ax'))), c.Const(c.Pointer(c.Value('double', 'ay'))), c.Const(c.Pointer(c.Value('double', 'bx'))), c.Const(c.Pointer(c.Value('double', 'by'))), c.Const(c.Pointer(c.Value('double', 'u_next'))), c.Const(c.Pointer(c.Value('int', 'idx')))]),
    c.Block([
        c.Assign('zetax[OPS_ACC0(0, 0)]', 'bx[OPS_ACC6(0, 0)] * zetax[OPS_ACC0(0, 0)] +\n\
                            ax[OPS_ACC4(0, 0)] * (u_next[OPS_ACC8(1, 0)] - (2 * u_next[OPS_ACC8(0, 0)]) + u_next[OPS_ACC8(-1, 0)]) / (dx * dx) +\n\
                            ax[OPS_ACC4(0, 0)] * (wx[OPS_ACC2(1, 0)] - wx[OPS_ACC2(-1, 0)]) / (2 * dx)'),
        c.Line(),

        c.LineComment('printf("zetax=%lf | bx=%lf | ax=%lf | h=%lf | idx = %d %d\\n", zetax[OPS_ACC0(0, 0)], bx[OPS_ACC6(0, 0)], ax[OPS_ACC4(0, 0)], dx, idx[0], idx[1]);'),
        c.Line(),

        c.Assign('zetay[OPS_ACC1(0, 0)]', 'by[OPS_ACC7(0, 0)] * zetay[OPS_ACC1(0, 0)] +\n\
                            ay[OPS_ACC5(0, 0)] * (u_next[OPS_ACC8(0, 1)] - (2 * u_next[OPS_ACC8(0, 0)]) + u_next[OPS_ACC8(0, -1)]) / (dy * dy) +\n\
                            ay[OPS_ACC5(0, 0)] * (wy[OPS_ACC3(0, 1)] - wy[OPS_ACC3(0, -1)]) / (2 * dy)'),
        c.Line(),

        c.LineComment('printf("zetay=%lf | by=%lf | ay=%lf | h=%lf | idx = %d %d\\n", zetay[OPS_ACC1(0, 0)], by[OPS_ACC7(0, 0)], ay[OPS_ACC5(0, 0)], dy, idx[0], idx[1]);')
    ])
))
code.append(c.Line())

code.append(c.MultilineComment('\
@brief  Initializes auxiliar CPML variables in X dimension.\n\
@note\n\
@param  *ax:\n\
@param  *bx:\n\
@param  *idx: Index being calculated.\n\
@retval None'))
code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'initializeAxBx'), [c.Pointer(c.Value('double', 'ax')), c.Pointer(c.Value('double', 'bx')), c.Const(c.Pointer(c.Value('int', 'idx')))]),
    c.Block([
        c.LineComment('Local variables'),
        c.Value('double', 'd_x, b_x, f_x, alfa_x'),
        c.Value('double', 'Lx'),
        c.Line(),

        c.Assign('Lx', 'CPML * dx'),
        c.LineComment('This sould be outside the kernel.'),
        c.Line(),

        c.make_multiple_ifs(
            [
                ['idx[0] < 0', c.Assign('f_x', 'Lx - (idx[0] + CPML) * dx')],
                ['idx[0] > nx - 2 * CPML', c.Assign('f_x', 'dx * (idx[0] - nx + 2 * CPML)')],
                ['', c.Assign('f_x', '0')]
            ], 'last'
        ),
        c.Line(),

        c.Assign('alfa_x', 'M_PI * RICKER_PEAK_FREQUENCY * (Lx - f_x) / Lx'),
        c.Assign('d_x', '-3.0 * log(0.00001) / (2.0 * Lx) * MAX_VELOCITY * (f_x / Lx) * (f_x / Lx)'),
        c.Assign('b_x', 'exp(-dt * (d_x + alfa_x))'),
        c.Line(),

        c.Assign('ax[OPS_ACC0(0, 0)]', 'd_x * (b_x - 1) / (d_x + alfa_x)'),
        c.Assign('bx[OPS_ACC1(0, 0)]', 'exp(-dt * (d_x + alfa_x))')
    ])
))
code.append(c.Line())

code.append(c.MultilineComment('\
@brief  Initializes auxiliar CPML variables in Y dimension.\n\
@note\n\
@param  *ay:\n\
@param  *by:\n\
@param  *idx:\n\
@retval None\
'))
code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'initializeAyBy'), [c.Pointer(c.Value('double', 'ay')), c.Pointer(c.Value('double', 'by')), c.Const(c.Pointer(c.Value('int', 'idx')))]),
    c.Block([
        c.LineComment('Local variables'),
        c.Value('double', 'd_y, b_y, f_y, alfa_y'),
        c.Value('double', 'Ly'),
        c.Line(),

        c.Assign('Ly', 'CPML * dy'),
        c.LineComment('This sould be outside the kernel.'),
        c.Line(),

        c.make_multiple_ifs(
            [
                ['idx[0] < 0', c.Block([
                    c.Assign('f_y', 'Ly - (idx[0] + CPML) * dy')
                    ])
                ],
                ['idx[0] > ny - 2 * CPML', c.Block([
                    c.Assign('f_y', 'dy * (idx[0] - ny + 2 * CPML)')
                    ])
                ],
                ['', c.Block([
                    c.Assign('f_y', '0.0')
                    ])
                ]
            ], 'last'
        ),
        c.Line(),

        c.Assign('alfa_y', 'M_PI * RICKER_PEAK_FREQUENCY * (Ly - f_y) / Ly'),
        c.Assign('d_y', '-3.0 * log(0.00001) / (2.0 * Ly) * MAX_VELOCITY * (f_y / Ly) * (f_y / Ly)'),
        c.Assign('b_y', 'exp(-dt * (d_y + alfa_y))'),
        c.Line(),

        c.Assign('ay[OPS_ACC0(0, 0)]', 'd_y * (b_y - 1) / (d_y + alfa_y)'),
        c.Assign('by[OPS_ACC1(0, 0)]', 'exp(-dt * (d_y + alfa_y))')
    ])
))
code.append(c.Line())

code.append(c.MultilineComment('\
@brief  Copy data from an array to another.\n\
@note\n\
@param  *to: Array will be copied here.\n\
@param  *from: Array to be copied.\n\
@retval None\
'))
code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'makeCopy'),[c.Pointer(c.Value('double', 'to')), c.Const(c.Pointer(c.Value('double', 'from')))]),
    c.Block([
        c.Assign('to[OPS_ACC0(0, 0)]', 'from[OPS_ACC1(0, 0)]')
    ])
))
code.append(c.Line())

code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('void', 'save_array'),[c.Value('char', 'filename[]'), c.Pointer(c.Value('double', 'array')), c.Value('int', 'size')]),
    c.Block([
        c.Pointer(c.Value('FILE', 'fp')),
        c.Value('int', 'i'),
        c.Line(),

        c.Assign('fp', 'fopen(filename, "w")'),
        c.Line(),

        c.For('i = 0', 'i < size', 'i++',
            c.Block([
                c.Statement('fprintf(fp, \"%lf \", array[i])')
            ])
        ),
        c.Statement('fflush(fp)'),
        c.Statement('fclose(fp)')
    ])
))
code.append(c.Line())

code.append(c.MultilineComment('\
@brief  Entry point for wave propagation example using OPS.\n\
@note\n\
@param  argc: Number of input arguments. Expected 7.\n\
@param  *argv[]: xSize:\n\
                 ySize:\n\
                 xIntervals:\n\
                 yIntervals:\n\
                 tTotal:\n\
                 tIntervals:\n\
@retval Positive number for sucess. Negative for error.\
'))
code.append(c.FunctionBody(
    c.FunctionDeclaration(c.Value('int', 'propagationCPML'), [c.Value('int', 'xSize'), c.Value('int', 'ySize'), c.Value('int', 'xIntervals'), c.Value('int', 'yIntervals'), c.Value('int', 'tTotal'), c.Value('int', 'tIntervals')]),
    c.Block([
        c.LineComment('Local variables'),
        # c.Pointer(c.Pointer(c.Value('double', 'u'))),
        # c.Pointer(c.Value('double', 'source')),
        c.Value('double', '**u, *source'),
        c.LineComment('u represents the grid data, c the velocity for each location in space and source is a vector of the source magnitude in time.'),
        c.Value('int', 'i, j, t'),
        c.LineComment('Counters.'),
        c.Line(),

        c.LineComment('OPS variables'),
        c.Value('int', 'size[2]'),
        c.LineComment('Size in each dimension.'),
        c.Initializer(c.Value('int', 'base[]'), '{0, 0}'),
        c.LineComment('Base indices for each dimension of the block.'),
        c.Initializer(c.Value('int', 'd_m[]'), '{-BORDER_SIZE, -BORDER_SIZE}'),
        c.LineComment('d_m - padding from the face in the negative direction for each dimension. (Halo)'),
        c.Initializer(c.Value('int', 'd_p[]'), '{BORDER_SIZE, BORDER_SIZE}'),
        c.LineComment('d_p - padding from the face in the positive direction for each dimension. (Halo)'),
        c.Value('ops_block', 'wave_grid'),
        c.LineComment('Block struct variable.'),
        c.Value('ops_dat', 'd_u_Previous, d_u_Current, d_u_Next, d_velocity'),
        c.LineComment('Block data.'),
        c.Line(),

        c.LineComment('CPML variables'),
        c.Value('double','*wx, *wy, *zetax, *zetay'),
        # c.Pointer(c.Value('double', 'wx')),
        # c.Pointer(c.Value('double', 'wy')),
        # c.Pointer(c.Value('double', 'zetax')),
        # c.Pointer(c.Value('double', 'zetay')),
        c.Value('double', '*ax, *bx, *ay, *by'),
        c.Line(),

        c.LineComment('OPS variables for CPML conditions.'),
        c.Initializer(c.Value('int', 's1d_0[]'), '{0}'),
        c.Initializer(c.Value('int', 's2d_00[]'), '{0, 0}'),
        c.Initializer(c.Value('int', 's2d_8th[]'), '{0, 0, 1, 0, -1, 0, 2, 0, -2, 0, 3, 0, -3, 0, 4, 0, -4, 0, 0, 1, 0, -1, 0, 2, 0, -2, 0, 3, 0, -3, 0, 4, 0, -4}'),
        c.Initializer(c.Value('int', 's2d_2px[]'), '{1, 0, -1, 0}'),
        c.Initializer(c.Value('int', 's2d_2py[]'), '{0, 1, 0, -1}'),
        c.Initializer(c.Value('int', 's2d_4pt[]'), '{1, 0, -1, 0, 0, 1, 0, -1}'),
        c.Initializer(c.Value('int', 's2d_5pt[]'), '{1, 0, 0, 0, -1, 0, 0, 1, 0, -1}'),
        c.Initializer(c.Value('int', 'stride2d_x[]'), '{1, 0}'),
        c.LineComment('Stride only in X dimension.'),
        c.Initializer(c.Value('int', 'stride2d_y[]'), '{0, 1}'),
        c.LineComment('Stride only in Y dimension.'),
        c.Value('ops_dat', 'd_wx, d_wy, d_zetax, d_zetay'),
        c.Value('ops_dat', 'd_ax, d_ay, d_bx, d_by'),
        c.Value('ops_block', 'wave_grid_1d'),
        c.Line(),

        c.LineComment('Initialize OPS.'),
        c.Statement('ops_init(1, NULL, 1)'),
        c.Line(),

        c.Assign('nx', 'xIntervals + 2 * BORDER_SIZE'),
        c.Assign('ny', 'yIntervals + 2 * BORDER_SIZE'),
        c.Line(),

        c.LineComment('Set grid dimension.'),
        c.Assign('size[0]', 'xIntervals'),
        c.Assign('size[1]', 'yIntervals'),
        c.Line(),

        c.LineComment('Calculate size of each step for dimension and time.'),
        c.Assign('dx', '(double)xSize / (double)xIntervals'),
        c.Assign('dy', '(double)ySize / (double)yIntervals'),
        c.Assign('dt', '(double)tTotal / (double)tIntervals'),
        c.Line(),

        c.Statement('printf("dx=%lf | dy=%lf | dt=%lf\\n", dx, dy, dt)'),
        c.Statement('printf("nx=%d | ny=%d | nT=%d | borders=%d\\n", xIntervals, yIntervals, tIntervals, BORDER_SIZE)'),
        c.Line(),

        c.LineComment('Check CFL convergency conditions.'),
        c.If('dt / dx > 1 && dt / dy > 1', c.block_if_necessary([
            c.Statement('cout << "Does not comply with CFL conditions." << endl'),
            c.Statement('return -1')])
        ),
        c.Line(),

        c.Statement('ops_decl_const2("dx", 1, "double", &dx)'),
        c.Statement('ops_decl_const2("dy", 1, "double", &dy)'),
        c.Statement('ops_decl_const2("dt", 1, "double", &dt)'),
        c.Statement('ops_decl_const2("nx", 1, "double", &nx)'),
        c.Statement('ops_decl_const2("ny", 1, "double", &ny)'),
        c.Line(),

        c.Initializer(c.Value('int', 'range_CPML[]'), '{1 - BORDER_SIZE, xIntervals + BORDER_SIZE - 1, 1 - BORDER_SIZE, yIntervals + BORDER_SIZE - 1}'),
        c.Initializer(c.Value('int', 'range_CPML_aux_1D_X_init[]'), '{-BORDER_SIZE, xIntervals + BORDER_SIZE - 1, 0, 1}'),
        c.Initializer(c.Value('int', 'range_CPML_aux_1D_Y_init[]'), '{0, 1, -BORDER_SIZE, yIntervals + BORDER_SIZE - 1}'),
        c.Initializer(c.Value('int', 'range[]'), '{-BORDER_SIZE + 4, xIntervals + BORDER_SIZE - 4 - 1, -BORDER_SIZE + 4, yIntervals + BORDER_SIZE - 4 - 1}'),
        c.Line(),

        c.Initializer(c.Value('int', 'whole_range[]'), '{-BORDER_SIZE, xIntervals + BORDER_SIZE, -BORDER_SIZE, yIntervals + BORDER_SIZE}'),
        c.Line(),

        c.LineComment('Allocate grids. 3 x 2D space represented in 1 dimension.'),
        c.Assign('u', '(double **)malloc(3 * sizeof(double *))'),
        c.For('t = 0', 't < 3', 't++',
        c.Block([
            c.Assign('u[t]', '(double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double))')
        ])
        ),
        c.Line(),

        c.LineComment('Allocate source time vector.'),
        c.Assign('source', 'rickerSource(dt, tIntervals, RICKER_PEAK_FREQUENCY)'),
        c.Line(),

        c.LineComment('save_array("output/source-ops.txt", source, tIntervals);'),
        c.Line(),
        
        c.LineComment('Allocates CPML Variables.'),
        c.Assign('wx', '(double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Assign('wy', '(double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Assign('zetax', '(double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Assign('zetay', '(double *)calloc((xIntervals + 2 * BORDER_SIZE) * (yIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Assign('ax', '(double *)calloc((xIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Assign('bx', '(double *)calloc((xIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Assign('ay', '(double *)calloc((yIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Assign('by', '(double *)calloc((yIntervals + 2 * BORDER_SIZE), sizeof(double))'),
        c.Line(),

        c.LineComment('Defines structured block.'),
        c.Assign('wave_grid', 'ops_decl_block(2, "wave_grid")'),
        c.Assign('wave_grid_1d', 'ops_decl_block(1, "wave_grid_1d")'),
        c.Line(),

        c.LineComment('Declares stencil'),
        c.Assign('ops_stencil S1D_0', 'ops_decl_stencil(1, 1, s1d_0, "0")'),
        c.Assign('ops_stencil S2D_00', 'ops_decl_stencil(2, 1, s2d_00, "0,0")'),
        c.Assign('ops_stencil S2D_8TH', 'ops_decl_stencil(2, 17, s2d_8th, "8th")'),
        c.Assign('ops_stencil S2D_2PX', 'ops_decl_stencil(2, 2, s2d_2px, "2px")'),
        c.Assign('ops_stencil S2D_2PY', 'ops_decl_stencil(2, 2, s2d_2py, "2py")'),
        c.Assign('ops_stencil S2D_4PT', 'ops_decl_stencil(2, 4, s2d_4pt, "4pt")'),
        c.Assign('ops_stencil S2D_5PT', 'ops_decl_stencil(2, 5, s2d_5pt, "5pt")'),
        c.Assign('ops_stencil S2D_00_STRIDE_X', 'ops_decl_strided_stencil(2, 1, s2d_00, stride2d_x, "s2D_00_stride2D_x")'),
        c.Assign('ops_stencil S2D_00_STRIDE_Y', 'ops_decl_strided_stencil(2, 1, s2d_00, stride2d_y, "s2D_00_stride2D_y")'),
        c.Line(),

        c.LineComment('Set dataset.'),
        c.Assign('d_m[0]', '-BORDER_SIZE'),
        c.Assign('d_m[1]', '-BORDER_SIZE'),
        c.Assign('d_p[0]', 'BORDER_SIZE'),
        c.Assign('d_p[1]', 'BORDER_SIZE'),
        c.Assign('d_u_Previous', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, u[PREVIOUS], "double", "u_previous")'),
        c.Assign('d_u_Current', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, u[CURRENT], "double", "u_current")'),
        c.Assign('d_u_Next', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, u[NEXT], "double", "u_next")'),
        c.Line(),

        c.LineComment('Declare CPML datasets.'),
        c.Assign('size[0]', 'xIntervals'),
        c.Assign('size[1]', 'yIntervals'),
        c.Assign('d_wx', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, wx, "double", "wx")'),
        c.Assign('d_wy', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, wy, "double", "wy")'),
        c.Assign('d_zetax', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, zetax, "double", "zetax")'),
        c.Assign('d_zetay', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, zetay, "double", "zetay")'),
        c.Line(),

        c.LineComment('Declare auxiliar CPML datasets'),
        c.Assign('d_m[0]', '-BORDER_SIZE'),
        c.Assign('d_m[1]', '0'),
        c.Assign('d_p[0]', 'BORDER_SIZE'),
        c.Assign('d_p[1]', '0'),
        c.Assign('size[0]', 'xIntervals'),
        c.Assign('size[1]', '1'),
        c.Assign('d_ax', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, ax, "double", "d_ax")'),
        c.Assign('d_bx', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, bx, "double", "d_bx")'),
        c.Assign('d_m[0]', '0'),
        c.Assign('d_m[1]', '-BORDER_SIZE'),
        c.Assign('d_p[0]', '0'),
        c.Assign('d_p[1]', 'BORDER_SIZE'),
        c.Assign('size[0]', '1'),
        c.Assign('size[1]', 'yIntervals'),
        c.Assign('d_ay', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, ay, "double", "d_ay")'),
        c.Assign('d_by', 'ops_decl_dat(wave_grid, 1, size, base, d_m, d_p, by, "double", "d_by")'),
        c.Line(),

        c.Initializer(c.Value('int', 'source_location_range[]'), '{SOURCE_LOCATION_X, SOURCE_LOCATION_X + 1, SOURCE_LOCATION_Y, SOURCE_LOCATION_Y + 1}'),
        c.Initializer(c.Value('double', 'source_field_factor'), 'getVelocityAtSourcePoint(FourHorizontalLayers, SOURCE_LOCATION_X, SOURCE_LOCATION_Y, xIntervals, yIntervals) * dt'),
        c.Assign('source_field_factor', 'source_field_factor * source_field_factor / (dx * dy)'),
        c.Initializer(c.Value('double', 'source_magnitude'), '0'),
        c.Line(),

        c.Assign('d_velocity', 'initializeFourHorizontalLayersModel(xIntervals, yIntervals, BORDER_SIZE, dt, dx)'),
        c.Line(),

        c.LineComment('ops_print_dat_to_txtfile(d_velocity, "output/velocity-ops.txt");'),
        c.Line(),

        c.LineComment('Initialize CPML auxiliar variables.'),

        c.Statement('ops_par_loop(initializeAxBx, "initializeAxBx", wave_grid, 2, range_CPML_aux_1D_X_init,\n\
                 ops_arg_dat(d_ax, 1, S2D_00, "double", OPS_WRITE),\n\
                 ops_arg_dat(d_bx, 1, S2D_00, "double", OPS_WRITE),\n\
                 ops_arg_idx())'),
        c.Line(),

        c.Statement('ops_par_loop(initializeAyBy, "initializeAyBy", wave_grid, 2, range_CPML_aux_1D_X_init,\n\
                 ops_arg_dat(d_ay, 1, S2D_00, "double", OPS_WRITE),\n\
                 ops_arg_dat(d_by, 1, S2D_00, "double", OPS_WRITE),\n\
                 ops_arg_idx())'),
        c.Line(),

        c.LineComment('ops_print_dat_to_txtfile(d_ax, "output/ax-ops.txt");'),
        c.LineComment('ops_print_dat_to_txtfile(d_bx, "output/bx-ops.txt");'),
        c.LineComment('ops_print_dat_to_txtfile(d_ay, "output/ay-ops.txt");'),
        c.LineComment('ops_print_dat_to_txtfile(d_by, "output/by-ops.txt");'),
        c.Line(),

        c.Assign('t', '0'),
        c.DoWhile('t < tIntervals', c.Block([
            c.LineComment('Propagates the wave.'),
            c.Statement('ops_par_loop(wavePropagation, "wave_propagation", wave_grid, 2, range,\n\
                     ops_arg_dat(d_u_Next, 1, S2D_00, "double", OPS_RW),\n\
                     ops_arg_dat(d_u_Current, 1, S2D_8TH, "double", OPS_READ),\n\
                     ops_arg_dat(d_u_Previous, 1, S2D_00, "double", OPS_READ),\n\
                     ops_arg_dat(d_velocity, 1, S2D_00, "double", OPS_READ),\n\
                     ops_arg_dat(d_wx, 1, S2D_2PX, "double", OPS_READ),\n\
                     ops_arg_dat(d_wy, 1, S2D_2PY, "double", OPS_READ),\n\
                     ops_arg_dat(d_zetax, 1, S2D_00, "double", OPS_READ),\n\
                     ops_arg_dat(d_zetay, 1, S2D_00, "double", OPS_READ),\n\
                     ops_arg_idx())'),
            c.Line(),

            c.Assign('source_magnitude', 'source_field_factor * source[t]'),
            c.Line(),

            c.LineComment('Injects the source.'),
            c.Statement('ops_par_loop(sourceInjection, "source_injection", wave_grid, 2, source_location_range,\n\
                     ops_arg_dat(d_u_Next, 1, S2D_00, "double", OPS_WRITE),\n\
                     ops_arg_gbl(&source_magnitude, 1, "double", OPS_READ))'),
            c.Line(),

            c.LineComment('Updates Omega X and Y'),
            c.Statement('ops_par_loop(update_omega, "update_omega", wave_grid, 2, range_CPML,\n\
                     ops_arg_dat(d_wx, 1, S2D_00, "double", OPS_RW),\n\
                     ops_arg_dat(d_wy, 1, S2D_00, "double", OPS_RW),\n\
                     ops_arg_dat(d_ax, 1, S2D_00_STRIDE_X, "double", OPS_READ),\n\
                     ops_arg_dat(d_ay, 1, S2D_00_STRIDE_Y, "double", OPS_READ),\n\
                     ops_arg_dat(d_bx, 1, S2D_00_STRIDE_X, "double", OPS_READ),\n\
                     ops_arg_dat(d_by, 1, S2D_00_STRIDE_Y, "double", OPS_READ),\n\
                     ops_arg_dat(d_u_Next, 1, S2D_4PT, "double", OPS_READ),\n\
                     ops_arg_idx())'),
            c.Line(),

            c.LineComment('ops_print_dat_to_txtfile(d_wx, "output/omegax-ops.txt");'),
            c.Line(),

            c.LineComment('Updates Zeta X and Y'),
            c.Statement('ops_par_loop(update_zeta, "update_zeta", wave_grid, 2, range_CPML,\n\
                     ops_arg_dat(d_zetax, 1, S2D_00, "double", OPS_WRITE),\n\
                     ops_arg_dat(d_zetay, 1, S2D_00, "double", OPS_WRITE),\n\
                     ops_arg_dat(d_wx, 1, S2D_2PX, "double", OPS_READ),\n\
                     ops_arg_dat(d_wy, 1, S2D_2PY, "double", OPS_READ),\n\
                     ops_arg_dat(d_ax, 1, S2D_00_STRIDE_X, "double", OPS_READ),\n\
                     ops_arg_dat(d_ay, 1, S2D_00_STRIDE_Y, "double", OPS_READ),\n\
                     ops_arg_dat(d_bx, 1, S2D_00_STRIDE_X, "double", OPS_READ),\n\
                     ops_arg_dat(d_by, 1, S2D_00_STRIDE_Y, "double", OPS_READ),\n\
                     ops_arg_dat(d_u_Next, 1, S2D_5PT, "double", OPS_READ),\n\
                     ops_arg_idx())'),
            c.Line(),

            c.LineComment('Save  dat to file.'),
            c.LineComment('if (t == 0)'),
            c.LineComment('{'),
            c.LineComment('    ops_print_dat_to_txtfile(d_zetax, "output/zetax-ops.txt");'),
            c.LineComment('}'),
            c.Line(),

            c.LineComment('ops_print_dat_to_txtfile(d_u_Current, "output/ops-output.txt");'),
            c.Line(),

            c.LineComment('Transfer current to previous.'),
            c.Statement('ops_par_loop(makeCopy, "copy_current_to_previous", wave_grid, 2, whole_range,\n\
                     ops_arg_dat(d_u_Previous, 1, S2D_00, "double", OPS_WRITE),\n\
                     ops_arg_dat(d_u_Current, 1, S2D_00, "double", OPS_READ))'),
            c.Line(),
            
            c.LineComment('Transfer next to current.'),
            c.Statement('ops_par_loop(makeCopy, "copy_next_to_current", wave_grid, 2, whole_range,\n\
                     ops_arg_dat(d_u_Current, 1, S2D_00, "double", OPS_WRITE),\n\
                     ops_arg_dat(d_u_Next, 1, S2D_00, "double", OPS_READ))'),
            c.Line(),

            c.LineComment('Save  dat to file.'),
            c.LineComment('if (t % 1000 == 0)'),
            c.LineComment('    ops_print_dat_to_txtfile(d_u_Current, "output/u-ops.txt");'),
            c.Line(),

            c.LineComment('printf("%d\\n", t);'),
            c.Statement('t++'),
        ])),
        c.Line(),

        c.Statement('ops_print_dat_to_txtfile(d_u_Current, "output/u-ops.txt")'),
        c.Line(),

        c.Statement('ops_exit()'),
        c.Line(),

        c.LineComment('free(u[PREVIOUS]);'),
        c.LineComment('free(u[CURRENT]);'),
        c.LineComment('free(u[NEXT]);'),
        c.LineComment('free(source);')
    ])
))

print(code)


# c.If("argc != 7", c.Statement("printf('Incorrect number of parameters.\n\t1)Dimensions X of the subsurface (eg. in Km or meters)\n\t2)Dimensions Y of the subsurface (eg. in Km or meters)\n\t3)Number of cells for X dimension\n\t4)Number of cells for Y dimension\n\t5)Real time propagation in seconds\n\t6)Number of time intervals.\n')"))
