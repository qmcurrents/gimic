# Apptainer Container

This is guide on [Apptainer](https://apptainer.org/) container use for Gimic.

## Building Container

Build container

```bash
cd  "gimic folder"
apptainer build --fakeroot gimic.sif container/gimic-apptainer.def
```

**Note**, you need to be in gimic folder or the build will faill. Apptainer does not support relative paths to definition file, yet, so this is important!

## Using Container

Container file is an executable by on its own and
is equivalent to `gimic` command.

```bash
gimic.sif gimic.inp > gimic.out
```

All commands can be called with `apptainer exec gimic.sif <command>` calls, e.g.

```bash
apptainer exec gimic.sif gimic gimic.inp > gimic.out
apptainer exec gimic.sif 3D-run.sh
```

This includes all executables that are located in

- `build/bin/`
- `jobscripts/`
- `tools/`

For mode detailed ways to call Gimic from the container refer
to [Apptainer documentation](https://apptainer.org/docs/user/main/quick_start.html).

### Note Certain Jobscripts Do not Work

Some jobscripts call external programs like `sbatch`.
These wont work with the container.
Because container by definition cannot access external programs.

These scripts are all located in `jobscripts` folder
and include `cluster` in their names like `current-profile-cluster.sh`.
But `local` versions of the commands work.

## Aliases for Gimic Commands

Consider also adding aliases, to make calling Gimic from the container easier, e.g.

```bash
alias gimic="path_to/gimic.sif"
alias 3D-run.sh="apptainer exec path_to/gimic.sif 3D-run.sh"
# etc.
```