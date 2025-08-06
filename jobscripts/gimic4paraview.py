import paraview.simple as pv
import subprocess
import sys
import os
import re

BOHR2ANGST = 0.52917726

def molxyz2cml(workdir):
    mol_xyz_path = os.path.join(workdir, 'mol.xyz')

    if not os.path.isfile(mol_xyz_path):
        print(f"[Prepare .pvsm script] \'mol.xyz\' file not found: {mol_xyz_path}", file=sys.stderr)

    else:
        try:
            subprocess.run("module load openbabel",
                           shell=True, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run(f'obabel -ixyz {mol_xyz_path} -ocml > ' + os.path.join(workdir, 'mol.cml'),
                           shell=True, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception as e1:
            try:
                subprocess.run(f'obabel.exe -ixyz {mol_xyz_path} -ocml > ' + os.path.join(workdir, 'mol.cml'),
                               shell=True, check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except Exception as e2:
                print(f"[Prepare .pvsm script] obabel failed:\n  Linux attempt: {e1}\n  Windows attempt: {e2}",
                      file=sys.stderr)

        mol_cml_path = os.path.join(workdir, 'mol.cml')
        mol_bohr_cml_path = os.path.join(workdir, 'mol-bohr.cml')
        cml2bohr(mol_cml_path, mol_bohr_cml_path)

def cml2bohr(input_cml, output_cml):

    with open(input_cml, 'r') as infile, open(output_cml, 'w') as outfile:
        for line in infile:
            if '<atom id' in line:
                match = re.search(r'x3="([-.\d]+)"\s+y3="([-.\d]+)"\s+z3="([-.\d]+)"', line)

                if match:
                    x3, y3, z3 = map(float, match.groups())

                    line = re.sub(r'x3="[-.\d]+"', f'x3="{x3 / BOHR2ANGST:.6f}"', line)
                    line = re.sub(r'y3="[-.\d]+"', f'y3="{y3 / BOHR2ANGST:.6f}"', line)
                    line = re.sub(r'z3="[-.\d]+"', f'z3="{z3 / BOHR2ANGST:.6f}"', line)

            outfile.write(line)


def process_jmod_file(filename, name, renderView, isosurface_param,
                      display_flag=False):
    contours = []
    reader = pv.XMLImageDataReader(FileName=filename)
    pv.RenameSource(name, reader)

    for color, param, name in zip([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
                                  [isosurface_param, -isosurface_param],
                                  ['dia', 'para']):
        contour = pv.Contour(Input=reader)
        contour.ContourBy = ['scalars']
        contour.Isosurfaces = [param]
        contourDisplay = pv.Show(contour, renderView)
        pv.RenameSource(name, contour)

        contourDisplay.ColorArrayName = [None, '']
        contourDisplay.DiffuseColor = color
        contourDisplay.Visibility = display_flag
        contours.append(contour)

    return contours


def process_jvec_file(filename, name, renderView):
    reader_vec = pv.XMLImageDataReader(FileName=filename)
    pv.RenameSource(name + "_jvec", reader_vec)

    # Apply stream tracer
    streamTracer = pv.StreamTracer(Input=reader_vec,
                                Vectors=['POINTS', 'vectors'],
                                MaximumStreamlineLength=25.0)
    streamTracer.SeedType = 'Point Cloud'
    streamTracer.SeedType.Center = [0.0, 0.0, 0.0]
    streamTracer.SeedType.Radius = 1.0

    # Display stream tracer
    streamTracerDisplay = pv.Show(streamTracer, renderView)
    pv.RenameSource(name + "_streamTracer", streamTracer)
    streamTracerDisplay.Visibility = False

    # Apply tube filter
    tube = pv.Tube(Input=streamTracer)
    tube.Radius = 0.05

    # Display tube
    tubeDisplay = pv.Show(tube, renderView)
    tubeDisplay.ColorArrayName = ['POINTS', 'vectors']
    pv.RenameSource(name + "_tube", tube)

    # Set color preset
    pv.ColorBy(tubeDisplay, ('POINTS', 'vectors'))
    tubeDisplay.RescaleTransferFunctionToDataRange(False)

    # Get color transfer function for 'vectors'
    vectorsLUT = pv.GetColorTransferFunction('vectors')
    vectorsLUT.ApplyPreset('Black-Body Radiation', True)

    # Rescale to custom range 1e-6 - 0.1
    vectorsLUT.RescaleTransferFunction(1e-6, 0.1)

    # Hide tube and stream tracer
    tubeDisplay.Visibility = False

    return reader_vec, streamTracer, tube


def process_sigma_file(sigma_file, renderView, iso_param=0.2):

    name = os.path.basename(sigma_file).replace('.vtu', '')
    reader = pv.XMLUnstructuredGridReader(FileName=sigma_file)
    pv.RenameSource(name, reader)

    # Slice
    slice = pv.Slice(Input=reader)
    slice.SliceType = 'Plane'
    slice.SliceOffsetValues = [0.0]
    slice.SliceType.Origin = [0.0, 0.0, 0.0]
    slice.SliceType.Normal = [0.0, 1.0, 0.0]
    sliceDisplay = pv.Show(slice, renderView)
    pv.RenameSource(name + "_slice", slice)

    # Coloring
    sliceDisplay.ColorArrayName = ['POINTS', 'scalars']
    pv.ColorBy(sliceDisplay, ('POINTS', 'scalars'))

    # Get and configure LUT
    scalarsLUT = pv.GetColorTransferFunction('scalars')

    # Rescale to fixed custom range
    scalarsLUT.RescaleTransferFunction(-iso_param, iso_param)

    # Disable automatic rescaling
    scalarsLUT.AutomaticRescaleRangeMode = "Never"

    # Optional: invert
    scalarsLUT.InvertTransferFunction()
    # Lighting
    sliceDisplay.Ambient = 1.0
    sliceDisplay.Diffuse = 0.25

    pv.Hide3DWidgets(renderView)
    sliceDisplay.Visibility = True

    renderView.CameraPosition = [0, 100, 0]
    renderView.CameraFocalPoint = [0, 0, 0]
    renderView.CameraViewUp = [0, 0, 1]

    pv.Render()
    pv.SaveScreenshot(os.path.join(workdir, f"{name}_slice.png"), renderView, ImageResolution=[1920, 1080])

    sliceDisplay.Visibility = False



def create_pvsm(workdir, filename, mo_dir,
                cur_dens_isosurface_param=0.005, sigma_isosurface_param=0.2):
    print("[Prepare .pvsm script] Preparing .pvsm file for Paraview.", file=sys.stderr, flush=True)
    pv.ResetSession()

    try:
        molxyz2cml(workdir)
    except:
        print(f"[Prepare .pvsm script]Failed to convert {workdir}/mol.xyz to {workdir}/mol-bohr.cml", file=sys.stderr)

    renderView = pv.GetActiveViewOrCreate('RenderView')
    renderView.Background = [1.0, 1.0, 1.0]
    renderView.UseColorPaletteForBackground = False

    cml_path = os.path.abspath(os.path.join(workdir, "mol-bohr.cml"))
    cmlFile = pv.OpenDataFile(cml_path)
    cmlDisplay = pv.Show(cmlFile, renderView)

    contours = []
    jmod_path = os.path.abspath(os.path.join(workdir, 'jmod.vti'))
    jvec_path = os.path.abspath(os.path.join(workdir, 'jvec.vti'))

    if os.path.isfile(jmod_path):
        contours.append(process_jmod_file(jmod_path, "main", renderView, cur_dens_isosurface_param, True))
    else:
        print(f"[Prepare .pvsm script] Warning: {jmod_path} not found. Skipping jmod visualization.", file=sys.stderr)
    if os.path.isfile(jvec_path):
        reader_vec, streamTracer, tube = process_jvec_file(jvec_path, "main", renderView)
    else:
        print(f"[Prepare .pvsm script] Warning: {jvec_path} not found. Skipping jvec visualization.", file=sys.stderr)

    # Process each sigma_*.vtu file
    sigma_files = sorted([
        os.path.abspath(os.path.join(workdir, f))
        for f in os.listdir(workdir)
        if f.startswith("sigma_") and f.endswith(".vtu")
    ])
    for sigma_file in sigma_files:
        process_sigma_file(sigma_file, renderView, iso_param=sigma_isosurface_param)

    # Optional: orbitals
    try:
        if os.path.isdir(os.path.join(workdir, mo_dir)):
            try:
                folders = sorted([f for f in os.listdir(os.path.join(workdir, mo_dir)) if
                                  os.path.isdir(os.path.join(workdir, mo_dir, f))],
                                 key=lambda x: int(x.split('_')[0]))
            except:
                folders = os.listdir(os.path.join(workdir, mo_dir))

            for d in folders:
                contours.append(process_jmod_file(os.path.abspath(os.path.join(workdir, mo_dir, d, 'jmod.vti')),
                                                  d, renderView, cur_dens_isosurface_param))
                reader_vec, streamTracer, tube = process_jvec_file(os.path.abspath(os.path.join(workdir, mo_dir, d, 'jvec.vti')),
                                                                   d, renderView)
    except:
        print("[Prepare .pvsm script] Failed to process orbitals.", file=sys.stderr)

    pv.Render()
    pv.SaveScreenshot(os.path.join(os.path.abspath(workdir), "jmod.jpeg"), ImageResolution=[1920, 1080])
    pv.SaveState(os.path.abspath(workdir) + '/' + filename)


if __name__ == "__main__":

    default_workdir = '.' #os.getcwd()
    default_filename = 'cdens.pvsm'
    default_mo_dir = 'XDENSes'

    workdir = sys.argv[1] if len(sys.argv) > 1 else default_workdir
    filename = sys.argv[2] if len(sys.argv) > 2 else default_filename
    mo_dir = sys.argv[3] if len(sys.argv) > 3 else default_mo_dir

    create_pvsm(workdir, filename, mo_dir)
