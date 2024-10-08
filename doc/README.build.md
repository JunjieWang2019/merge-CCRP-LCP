Obtaining the software
======================

The authoritative location of the software is the following git
repository:
   <https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-ges-tm.git>

Each released version may be identified by a version control system tag in
the form `ges-tm-v${version}`.

An example:

```console
$ git clone https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-ges-tm.git
$ cd mpeg-pcc-ges-tm
$ git checkout ges-tm-v1.0-rc1
```

> It is strongly advised to obtain the software using the version control
> system rather than to download a zip (or other archive) of a particular
> release.  The build system uses the version control system to accurately
> identify the version being built.

Building
--------

The codec is supported on Linux, OSX and Windows platforms.  The build
configuration is managed using CMake.

> It is strongly advised to build the software in a separate build directory.

### Linux
```console
$ mkdir build
$ cd build
$ cmake ..
$ make
$ tmc3/tmc3 --help
```

### OSX
```console
$ mkdir build
$ cd build
$ cmake .. -G Xcode
$ xcodebuild
$ tmc3/tmc3 --help
```

As an alternative, the generated XCode project may be opened and built from
XCode itself.

### Windows
```console
> md build
> cd build
> cmake .. -G "Visual Studio 15 2017 Win64"
```

Open the generated visual studio solution to build it.
