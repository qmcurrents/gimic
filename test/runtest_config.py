def configure(options, input_files, extra_args):
    """
    This function is used by runtest to configure runtest
    at runtime for GIMIC specific launch command and file naming.
    """

    from os import path

    launcher = 'gimic'
    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    inp = input_files[0]

    full_command = '{0} {1}'.format(launcher_full_path, inp)

    output_prefix = None

    relative_reference_path = 'reference'

    return launcher, full_command, output_prefix, relative_reference_path
