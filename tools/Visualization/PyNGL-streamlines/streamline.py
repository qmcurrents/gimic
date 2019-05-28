#!/usr/bin/env python

from math import sqrt
import sys
import string
import Ngl
import numpy
from optparse import OptionParser, OptionGroup

#-------------------------------------------------------------------------------

def parse_input():

    # initialize parser

    parser = OptionParser()

    # define options

    group = OptionGroup(parser, 'Basic options')
    group.add_option('--data',
                     type='string',
                     action='store',
                     default=None,
                     help='vector file (au) to read [default: %default]')
    group.add_option('--molecule',
                     type='string',
                     action='store',
                     default=None,
                     help='xy file (bohr) name to read [not required, default: %default]')
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Fine tuning')
    group.add_option('--smax',
                     type='float',
                     action='store',
                     default=0.1,
                     help='this will affect the range of vector norms [default: %default]')
    group.add_option('--legend',
                     action='store_true',
                     default=False,
                     help='get a legend bar [default: %default]')
    group.add_option('--text_color',
                     type='string',
                     action='store',
                     default='Blue',
                     # default='White',
                     help='text color [default: %default]')
    group.add_option('--rgb_min',
                     type='string',
                     action='store',
                     default='0.9 0.9 0.9',
                     help='min color (RGB) [default: %default]')
    group.add_option('--rgb_max',
                     type='string',
                     action='store',
                     default='0.0 0.0 0.0',
                     help='max color (RGB) [default: %default]')
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if not options.data:
        sys.stderr.write('you need to provide --data to read')
        sys.stderr.write('use --help to see all options')
        sys.exit()

    return options

#-------------------------------------------------------------------------------

def read_molecule(options):
    molecule = []
    if options.molecule:
        for line in open(options.molecule).readlines():
            # atom label
            e = line.split()[0]
            # coordinates
            x = float(line.split()[1])
            y = float(line.split()[2])
            molecule.append([e, x, y])
    return molecule

#-------------------------------------------------------------------------------

def get_vector_res(options, x_min, x_max, y_min, y_max):

    vector_res = Ngl.Resources()

    vector_res.nglFrame                 = False
    vector_res.nglMaximize              = False

    vector_res.vpWidthF                 = 0.5
    vector_res.vpHeightF                = 0.5

    vector_res.vcGlyphStyle             = "CurlyVector"
    vector_res.vcMonoFillArrowFillColor = False
    vector_res.vcRefLengthF             = 0.2
    vector_res.vcMinDistanceF           = 0.01
    vector_res.vcLineArrowThicknessF    = 1.0
    vector_res.vcRefAnnoOn              = False

    vector_res.vfXCStartV               = x_min
    vector_res.vfXCEndV                 = x_max
    vector_res.vfYCStartV               = y_min
    vector_res.vfYCEndV                 = y_max

    line_thickness                      = 2.0
    vector_res.tmXBFormat               = "f"
    vector_res.tmYLFormat               = "f"
    vector_res.tmXBMajorThicknessF      = line_thickness
    vector_res.tmXBMinorThicknessF      = line_thickness/2.0
    vector_res.tmXTMajorThicknessF      = line_thickness
    vector_res.tmXTMinorThicknessF      = line_thickness/2.0
    vector_res.tmYLMajorThicknessF      = line_thickness
    vector_res.tmYLMinorThicknessF      = line_thickness/2.0
    vector_res.tmYRMajorThicknessF      = line_thickness
    vector_res.tmYRMinorThicknessF      = line_thickness/2.0
    vector_res.tmBorderThicknessF       = line_thickness

    if options.legend:
        vector_res.lbOrientation            = "Horizontal"
        vector_res.pmLabelBarOrthogonalPosF = 0.1
        vector_res.pmLabelBarHeightF        = 0.05
        vector_res.pmLabelBarWidthF         = 0.3
    else:
        vector_res.lbLabelBarOn = False

    return vector_res

#-------------------------------------------------------------------------------

def read_data(options):

    f = open(options.data, 'r')
    data = map(string.strip, f.readlines())
    f.close()

    x_l = []
    y_l = []
    u_l = []
    v_l = []
    s_l = []

    for line in data:
        # skip empty lines - gimic
        if not line.strip():
            continue
        else:
            x = float(line.split()[1])
            y = float(line.split()[0])
            u = float(line.split()[3])
            v = float(line.split()[2])

            if len(line.split()) == 5:
                s = float(line.split()[4])
            else:
                # norm
                s = (u**2.0 + v**2.0)**0.5

            # find smin and smax
            if s > options.smax:
                s = options.smax
            if s < -options.smax:
                s = -options.smax

            if (abs(u)) < 1.0e-20:
                u = 0.0
            if (abs(v)) < 1.0e-20:
                v = 0.0

            x_l.append(x)
            y_l.append(y)
            u_l.append(u)
            v_l.append(v)
            s_l.append(s)

    nr_points = int(sqrt(len(u_l)))

    # get min and max values in each direction
    x_min = min(x_l)
    x_max = max(x_l)
    y_min = min(y_l)
    y_max = max(y_l)

    u_array = []
    v_array = []
    s_array = []

    ipoint = 0
    # loop in x direction and collect
    for i in range(nr_points):   

        u_slice = []
        v_slice = []
        s_slice = []

        # loop in y direction and collect normalized vec
        for j in range(nr_points):

            norm = (u_l[ipoint]**2.0 + v_l[ipoint]**2.0)**0.5

            # normalize vector
            if abs(norm) > 1.0e-15:
                u_slice.append(u_l[ipoint]/norm)
                v_slice.append(v_l[ipoint]/norm)
            else:
                u_slice.append(u_l[ipoint])
                v_slice.append(v_l[ipoint])
            s_slice.append(s_l[ipoint])

            ipoint = ipoint + 1

        u_array.append(u_slice)
        v_array.append(v_slice)
        s_array.append(s_slice)

    u_array2 = numpy.array(u_array, numpy.float32)
    v_array2 = numpy.array(v_array, numpy.float32)
    s_array2 = numpy.array(s_array, numpy.float32)

    return u_array2, v_array2, s_array2, x_min, x_max, y_min, y_max

#-------------------------------------------------------------------------------

def main():

    options = parse_input()

    u_array2, v_array2, s_array2, x_min, x_max, y_min, y_max = read_data(options)

    # no jpg, jpeg, tiff, tif, bmp
    wks_type = 'ps'
    # wks_type = 'png'
    wks = Ngl.open_wks(wks_type, options.data)

    color_res  = Ngl.Resources()
    gs_res     = Ngl.Resources()
    text_res   = Ngl.Resources()

    rgb_min = []
    for s in options.rgb_min.split():
        rgb_min.append(float(s))

    rgb_max = []
    for s in options.rgb_max.split():
        rgb_max.append(float(s))

    l = []

    l.append([1.00, 1.00, 1.00])
    l.append([0.00, 0.00, 0.00])

    f_l = []
    n   = 100
    for f in range(n):
        f_l.append(float((f+1)*(1.0/n)))

    for f in f_l:
        r = rgb_min[0] + f*(rgb_max[0]-rgb_min[0])
        g = rgb_min[1] + f*(rgb_max[1]-rgb_min[1])
        b = rgb_min[2] + f*(rgb_max[2]-rgb_min[2])
        l.append([r, g, b])

    rgb_map = numpy.array(l, 'f')

    color_res.wkColorMap = rgb_map
    Ngl.set_values(wks, color_res)

    vector_res = get_vector_res(options, x_min, x_max, y_min, y_max)
    plot = Ngl.vector_scalar(wks, u_array2, v_array2, s_array2, vector_res)

    color_res.wkColorMap = "default"
    Ngl.set_values(wks, color_res)

    line_thickness            = 2.0

    gs_res.gsMarkerSizeF      = 10.0
    gs_res.gsMarkerThicknessF = 2.0
    gs_res.gsLineDashPattern  = 2
    gs_res.gsLineThicknessF   = line_thickness

    text_res.txFont        = 21
    text_res.txFontHeightF = 0.015
    text_res.txFontColor   = options.text_color

    for atom in read_molecule(options):

        gs_res.gsMarkerIndex = 16
        gs_res.gsMarkerColor = "White"
        Ngl.polymarker(wks, plot, [atom[2]], [atom[1]], gs_res)

        gs_res.gsMarkerIndex = 4
        gs_res.gsMarkerColor = options.text_color
        Ngl.polymarker(wks, plot, [atom[2]], [atom[1]], gs_res)

# no printing of labels for a better overview
#       Ngl.text(wks, plot, atom[0], [atom[2] + 1.0], [atom[1]], text_res)

    Ngl.frame(wks)
    del plot
    Ngl.end()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
